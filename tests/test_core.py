from contextlib import redirect_stderr
from io import StringIO
import logging
import math
from pathlib import Path
import subprocess
import sys
from dataclasses import replace
from types import SimpleNamespace
from tempfile import TemporaryDirectory
from threading import Lock
from time import sleep
import unittest
from unittest.mock import patch

import pandas as pd

from bridge.rmlst.cli import parser as cli_parser
from bridge.rmlst.clustering import cluster_candidates
from bridge.rmlst.consensus import (
    MINIMAP2_SHORT_READ_PARAMETERS,
    _align_cluster_reads,
    alignment_min_mapping_quality,
    cliquesnv_heap_megabytes,
    cliquesnv_launcher,
    compatible_alleles,
    composite_consensus_possible,
    mask_weak_variants,
    pad_terminal_consensus,
    run_cluster_consensus,
)
from bridge.rmlst.haplotypes import (
    AlignmentEvidence,
    BaseEvidence,
    HaplotypeCandidate,
    alignment_evidence_from_sam,
    evaluate_haplotype,
    maximum_minor_allele_depth,
    read_cliquesnv_haplotypes,
)
from bridge.rmlst.kallisto import (
    derive_candidates,
    index_kmer_size,
    run_candidate_pseudobam,
    run_quant,
)
from bridge.rmlst.models import Config, ConsensusCall, FastaRecord
from bridge.rmlst.plotting import write_called_species_pi_pdf
from bridge.rmlst.pipeline import (
    INTERMEDIATE_OUTPUT_NAMES,
    _remove_intermediate_outputs,
)
from bridge.rmlst.reporting import write_direct_species_alleles
from bridge.rmlst.seqdef import add_rst_matches, allele_species
from bridge.rmlst.species import (
    COVERAGE_COLUMNS,
    demoted_profile_species,
    direct_profile_species,
    em_species,
    final_species_calls,
    merge_consensus_calls,
    species_coverage,
)
from bridge.tools import close_logging, configure_logging, run, stage


def call(cluster, locus, consensus, alleles, species, counts=10.0):
    return ConsensusCall(
        sample_id="sample",
        locus_id=locus,
        runtime_cluster_id=cluster,
        representative=f"{locus}_1",
        source_members=(f"{locus}_1",),
        est_counts=counts,
        assigned_fragments=10,
        exclusive_fragments=8,
        ambiguous_fragments=2,
        status="complete_exact_match" if "N" not in consensus else "unique_partial_match",
        consensus=consensus,
        callable_bases=sum(base in "ACGT" for base in consensus),
        allele_length=len(consensus),
        callable_breadth=sum(base in "ACGT" for base in consensus) / len(consensus),
        compatible_alleles=tuple(alleles),
        compatible_species=tuple(species),
        allele_search_scope="cluster_members_exact",
        consensus_path=None,
        blast_path=None,
    )


class CoreTests(unittest.TestCase):
    def test_intermediate_output_cleanup_preserves_final_and_unrelated_outputs(self):
        with TemporaryDirectory() as directory:
            output = Path(directory)
            for name in INTERMEDIATE_OUTPUT_NAMES:
                path = output / name
                if Path(name).suffix:
                    path.write_text("diagnostic\n", encoding="utf-8")
                else:
                    path.mkdir()
                    (path / "detail.txt").write_text(
                        "diagnostic\n", encoding="utf-8"
                    )

            retained = (
                "called_species_pi.pdf",
                "species_calls.csv",
                "species_allele_fastas",
                "rmlst_manifest.json",
                "rmlst.log",
                "unrelated-user-file.txt",
            )
            for name in retained:
                path = output / name
                if name == "species_allele_fastas":
                    path.mkdir()
                    (path / "Species_A.fasta").write_text(
                        ">a\nACGT\n", encoding="utf-8"
                    )
                else:
                    path.write_text("retained\n", encoding="utf-8")

            _remove_intermediate_outputs(output)

            self.assertTrue(all(
                not (output / name).exists()
                for name in INTERMEDIATE_OUTPUT_NAMES
            ))
            self.assertTrue(all((output / name).exists() for name in retained))

    def test_composite_consensus_requires_member_supplied_base_at_every_position(self):
        members = {
            "BACT000001_1": "ACGT",
            "BACT000001_2": "ATGA",
        }
        self.assertTrue(composite_consensus_possible("ATGT", members))
        self.assertFalse(composite_consensus_possible("ATGC", members))
        self.assertFalse(composite_consensus_possible("ACGT", members))
        self.assertFalse(composite_consensus_possible("ATNT", members))

    def test_cli_enables_composite_filter_by_default(self):
        required = [
            "--allele-target", "reference.fasta",
            "--rmlstseqdef", "seqdef.tsv",
            "--reads", "reads_1.fastq", "reads_2.fastq",
            "--output-folder", "output",
        ]
        self.assertTrue(
            cli_parser().parse_args(required).filter_composite_consensus
        )
        self.assertFalse(
            cli_parser().parse_args(
                [*required, "--no-filter-composite-consensus"]
            ).filter_composite_consensus
        )

    def test_cli_optionally_distinguishes_subspecies(self):
        required = [
            "--allele-target", "reference.fasta",
            "--rmlstseqdef", "seqdef.tsv",
            "--reads", "reads_1.fastq", "reads_2.fastq",
            "--output-folder", "output",
        ]
        self.assertFalse(
            cli_parser().parse_args(required).distinguish_subspecies
        )
        self.assertTrue(
            cli_parser().parse_args(
                [*required, "--distinguish-subspecies"]
            ).distinguish_subspecies
        )

    def test_allele_species_separates_named_and_blank_subspecies(self):
        with TemporaryDirectory() as directory:
            seqdef = Path(directory) / "seqdef.tsv"
            seqdef.write_text(
                "rST\tBACT000001\tspecies\tsubspecies\n"
                "1\t7\tMycobacterium tuberculosis\tbovis\n"
                "2\t7\tMycobacterium tuberculosis\t\n"
                "3\t8\tMycobacterium tuberculosis\tbovis\n",
                encoding="utf-8",
            )

            pooled = allele_species(seqdef)
            distinguished = allele_species(
                seqdef, distinguish_subspecies=True
            )

            self.assertEqual(
                pooled["BACT000001_7"], {"Mycobacterium tuberculosis"}
            )
            self.assertEqual(distinguished["BACT000001_7"], {
                "Mycobacterium tuberculosis subsp. bovis",
                "Mycobacterium tuberculosis [no subspecies annotation]",
            })
            self.assertEqual(distinguished["BACT000001_8"], {
                "Mycobacterium tuberculosis subsp. bovis"
            })

    def test_subspecies_mode_requires_subspecies_column(self):
        with TemporaryDirectory() as directory:
            seqdef = Path(directory) / "seqdef.tsv"
            seqdef.write_text(
                "rST\tBACT000001\tspecies\n1\t7\tSpecies A\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "subspecies column"):
                allele_species(seqdef, distinguish_subspecies=True)

    def test_called_species_plot_is_created_as_pdf(self):
        calls = pd.DataFrame([{
            "sample_id": "sample_1",
            "species": "Streptococcus pneumoniae",
            "reported": True,
            "pi": 0.75,
            "allele_evidence_unambiguous_loci": 18,
            "allele_evidence_ambiguous_loci": 35,
        }])
        with TemporaryDirectory() as directory:
            output = Path(directory) / "called_species_pi.pdf"
            returned = write_called_species_pi_pdf(calls, output, "sample_1")
            self.assertEqual(returned, output)
            self.assertTrue(output.read_bytes().startswith(b"%PDF"))

    def test_no_reported_species_still_creates_default_plot(self):
        calls = pd.DataFrame([{
            "sample_id": "sample_1",
            "species": "Streptococcus mitis",
            "reported": False,
            "pi": 0.25,
        }])
        with TemporaryDirectory() as directory:
            output = Path(directory) / "called_species_pi.pdf"
            write_called_species_pi_pdf(calls, output, "sample_1")
            self.assertTrue(output.read_bytes().startswith(b"%PDF"))

    def test_cli_accepts_native_cliquesnv_time_limit(self):
        args = cli_parser().parse_args([
            "--allele-target", "reference.fasta",
            "--rmlstseqdef", "seqdef.tsv",
            "--reads", "reads_1.fastq", "reads_2.fastq",
            "--output-folder", "output",
            "--cliquesnv-time-limit-seconds", "37",
        ])
        self.assertEqual(args.cliquesnv_time_limit_seconds, 37)
        self.assertEqual(args.minor_allele_min_depth, 5)
        self.assertEqual(args.short_read_aligner, "minimap2")
        self.assertFalse(args.verbose)

    def test_cli_accepts_verbose(self):
        args = cli_parser().parse_args([
            "--allele-target", "reference.fasta",
            "--rmlstseqdef", "seqdef.tsv",
            "--reads", "reads_1.fastq", "reads_2.fastq",
            "--output-folder", "output",
            "--verbose",
        ])
        self.assertTrue(args.verbose)

    def test_logging_keeps_full_file_and_filters_default_console(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            for verbose in (False, True):
                with self.subTest(verbose=verbose):
                    stream = StringIO()
                    log_path = root / f"verbose_{verbose}.log"
                    with redirect_stderr(stream):
                        handlers = configure_logging(log_path, verbose)
                        try:
                            run([
                                sys.executable, "-c",
                                "print('processed 100000 reads')",
                            ], "kallisto")
                            with stage("kallisto"):
                                logging.getLogger("bridge.rmlst").info(
                                    "bootstrap_parse_status=running"
                                )
                        finally:
                            close_logging(handlers)

                    console = stream.getvalue()
                    full_log = log_path.read_text(encoding="utf-8")
                    self.assertIn("stage=kallisto status=started", console)
                    self.assertIn("stage=kallisto status=completed", console)
                    self.assertEqual(
                        "processed 100000 reads" in console,
                        verbose,
                    )
                    self.assertEqual(
                        "bootstrap_parse_status=running" in console,
                        verbose,
                    )
                    self.assertIn("command=", full_log)
                    self.assertIn("processed 100000 reads", full_log)
                    self.assertIn("bootstrap_parse_status=running", full_log)

    def test_failed_captured_tool_logs_stdout_and_stderr(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            log_path = root / "rmlst.log"
            handlers = configure_logging(log_path)
            try:
                completed = SimpleNamespace(
                    returncode=2,
                    stdout="inspect stdout detail\n",
                    stderr="inspect stderr detail\n",
                )
                with patch("bridge.rmlst.kallisto.subprocess.run", return_value=completed):
                    with self.assertRaises(subprocess.CalledProcessError) as raised:
                        index_kmer_size(root / "reference.idx", "kallisto")
            finally:
                close_logging(handlers)

            self.assertEqual(raised.exception.returncode, 2)
            full_log = log_path.read_text(encoding="utf-8")
            self.assertIn("stdout: inspect stdout detail", full_log)
            self.assertIn("stderr: inspect stderr detail", full_log)

    def test_failed_streamed_tool_keeps_stdout_and_stderr_in_log(self):
        with TemporaryDirectory() as directory:
            log_path = Path(directory) / "rmlst.log"
            handlers = configure_logging(log_path)
            try:
                with self.assertRaises(subprocess.CalledProcessError) as raised:
                    run([
                        sys.executable,
                        "-c",
                        (
                            "import sys; "
                            "print('mash stdout detail'); "
                            "print('mash stderr detail', file=sys.stderr); "
                            "raise SystemExit(255)"
                        ),
                    ], "mash")
            finally:
                close_logging(handlers)

            self.assertEqual(raised.exception.returncode, 255)
            full_log = log_path.read_text(encoding="utf-8")
            self.assertIn("mash stdout detail", full_log)
            self.assertIn("mash stderr detail", full_log)

    def test_mash_sketch_retries_singularity_mount_collision(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            config = Config(
                allele_target=root / "reference.fasta",
                seqdef=root / "seqdef.tsv",
                reads=(root / "reads_1.fastq", root / "reads_2.fastq"),
                sample_id="sample",
                output=root / "output",
                tempdir=root / "work",
            )
            names = ("BACT000004_1", "BACT000004_2")
            calls = pd.DataFrame({
                "locus_id": ["BACT000004", "BACT000004"],
                "Template": list(names),
            })
            records = {
                name: FastaRecord(name, name, "ACGT") for name in names
            }
            sketch = config.tempdir / "mash" / "BACT000004_candidates.msh"
            attempts = 0

            def fake_run(command, tool, *, stdout_path=None):
                nonlocal attempts
                attempts += 1
                if attempts == 1:
                    sketch.write_text("partial", encoding="utf-8")
                    raise subprocess.CalledProcessError(
                        255,
                        command,
                        output=(
                            "FATAL: container creation failed: failed to mount "
                            "squashfs filesystem: file exists\n"
                        ),
                    )
                self.assertFalse(sketch.exists())

            completed = SimpleNamespace(returncode=0, stdout="", stderr="")
            with (
                patch("bridge.rmlst.clustering.executable", return_value="mash"),
                patch("bridge.rmlst.clustering.run", side_effect=fake_run),
                patch(
                    "bridge.rmlst.clustering.subprocess.run",
                    return_value=completed,
                ),
            ):
                result = cluster_candidates(calls, records, config)

            self.assertEqual(attempts, 2)
            self.assertEqual(len(result), 2)

    def test_mash_triangle_retries_silent_exit_255(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            config = Config(
                allele_target=root / "reference.fasta",
                seqdef=root / "seqdef.tsv",
                reads=(root / "reads_1.fastq", root / "reads_2.fastq"),
                sample_id="sample",
                output=root / "output",
                tempdir=root / "work",
            )
            names = ("BACT000005_1", "BACT000005_2")
            calls = pd.DataFrame({
                "locus_id": ["BACT000005", "BACT000005"],
                "Template": list(names),
            })
            records = {
                name: FastaRecord(name, name, "ACGT") for name in names
            }
            failed = SimpleNamespace(returncode=255, stdout="", stderr="")
            succeeded = SimpleNamespace(returncode=0, stdout="", stderr="")

            with (
                patch("bridge.rmlst.clustering.executable", return_value="mash"),
                patch("bridge.rmlst.clustering.run"),
                patch(
                    "bridge.rmlst.clustering.subprocess.run",
                    side_effect=(failed, succeeded),
                ) as triangle,
                self.assertLogs("mash", level="INFO") as captured,
            ):
                result = cluster_candidates(calls, records, config)

            self.assertEqual(triangle.call_count, 2)
            self.assertEqual(len(result), 2)
            self.assertTrue(any(
                "retrying triangle once" in message
                for message in captured.output
            ))

    def test_cli_selects_minimap2_backend(self):
        args = cli_parser().parse_args([
            "--allele-target", "reference.fasta",
            "--rmlstseqdef", "seqdef.tsv",
            "--reads", "reads_1.fastq", "reads_2.fastq",
            "--output-folder", "output",
            "--short-read-aligner", "minimap2",
            "--minimap2", "/tools/minimap2",
        ])
        self.assertEqual(args.short_read_aligner, "minimap2")
        self.assertEqual(args.minimap2, "/tools/minimap2")

    def test_cli_accepts_supplied_kallisto_index(self):
        args = cli_parser().parse_args([
            "--allele-target", "reference.fasta",
            "--rmlstseqdef", "seqdef.tsv",
            "--reads", "reads_1.fastq", "reads_2.fastq",
            "--output-folder", "output",
            "--kallisto-index", "reference.k41.idx",
        ])
        self.assertEqual(args.kallisto_index, Path("reference.k41.idx"))

    def test_supplied_kallisto_index_bypasses_index_construction(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            supplied = root / "reference.k41.idx"
            supplied.write_bytes(b"test-index")
            config = Config(
                allele_target=root / "reference.fasta",
                seqdef=root / "seqdef.tsv",
                reads=(root / "r1.fastq", root / "r2.fastq"),
                sample_id="sample", output=root / "output",
                tempdir=root / "work", kallisto_index=supplied,
                bootstraps=2,
            )
            with (
                patch("bridge.rmlst.kallisto.executable", return_value="/tools/kallisto"),
                patch("bridge.rmlst.kallisto.ensure_index") as ensure,
                patch("bridge.rmlst.kallisto.index_kmer_size", return_value=41) as inspect,
                patch("bridge.rmlst.kallisto.run") as runner,
            ):
                quant, returned_index, kmer_size = run_quant(
                    config, config.allele_target
                )

            ensure.assert_not_called()
            inspect.assert_called_once_with(supplied, "/tools/kallisto")
            self.assertEqual(returned_index, supplied)
            self.assertEqual(kmer_size, 41)
            self.assertIn("quant_index_", quant.name)
            self.assertEqual(
                runner.call_args.args[0][0:4],
                ["/tools/kallisto", "quant", "-i", str(supplied)],
            )

    def test_supplied_index_kmer_size_is_used_for_candidate_index(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            config = Config(
                allele_target=root / "reference.fasta",
                seqdef=root / "seqdef.tsv",
                reads=(root / "r1.fastq", root / "r2.fastq"),
                sample_id="sample", output=root / "output",
                tempdir=root / "work",
            )
            candidates = pd.DataFrame({"Template": ["BACT000001_1"]})
            records = {
                "BACT000001_1": FastaRecord(
                    "BACT000001_1", "BACT000001_1", "A" * 100
                )
            }
            with (
                patch("bridge.rmlst.kallisto.executable", return_value="/tools/kallisto"),
                patch("bridge.rmlst.kallisto.run") as runner,
            ):
                run_candidate_pseudobam(
                    config, candidates, records, index_kmer_size=41
                )

            index_command = runner.call_args_list[0].args[0]
            self.assertIn("candidate_alleles.k41.kallisto.idx", index_command[3])
            self.assertEqual(index_command[4:6], ["-k", "41"])

    def test_minimap2_backend_uses_validated_short_read_parameters(self):
        config = Config(
            allele_target=Path("reference.fasta"), seqdef=Path("seqdef.tsv"),
            reads=(Path("r1.fastq"), Path("r2.fastq")), sample_id="sample",
            output=Path("output"), tempdir=Path("work"),
            short_read_aligner="minimap2",
        )
        with patch("bridge.rmlst.consensus.run") as runner:
            _align_cluster_reads(
                config, "/tools/minimap2", Path("representative.fasta"),
                Path("r1.fastq"), Path("r2.fastq"), Path("alignment.sam"), 2,
            )

        runner.assert_called_once_with([
            "/tools/minimap2", *MINIMAP2_SHORT_READ_PARAMETERS,
            "-t", "2", "representative.fasta", "r1.fastq", "r2.fastq",
        ], "minimap2", stdout_path=Path("alignment.sam"))
        self.assertEqual(alignment_min_mapping_quality(config), 10)
        self.assertEqual(
            alignment_min_mapping_quality(
                replace(config, short_read_aligner="bwa-mem2")
            ),
            20,
        )

    def test_cliquesnv_heap_uses_override_and_worker_count(self):
        with patch(
            "bridge.rmlst.consensus.detected_memory_limit_bytes",
            return_value=(64 * 1024 ** 3, "physical_memory"),
        ) as detect:
            heap_mb, total_gb, source = cliquesnv_heap_megabytes(2, 16)

        self.assertEqual(heap_mb, 6553)
        self.assertEqual(total_gb, 16)
        self.assertEqual(source, "configured")
        detect.assert_not_called()

    def test_cliquesnv_heap_uses_detected_memory_when_unset(self):
        with patch(
            "bridge.rmlst.consensus.detected_memory_limit_bytes",
            return_value=(10 * 1024 ** 3, "cgroup:test"),
        ):
            heap_mb, total_gb, source = cliquesnv_heap_megabytes(5)

        self.assertEqual(heap_mb, 1638)
        self.assertEqual(total_gb, 10)
        self.assertEqual(source, "cgroup:test")

    def test_cliquesnv_launcher_places_heap_options_correctly(self):
        self.assertEqual(
            cliquesnv_launcher("/usr/bin/cliquesnv", 2048),
            ("/usr/bin/cliquesnv", "-Xms512m", "-Xmx2048m"),
        )
        with patch("bridge.rmlst.consensus.executable", return_value="/usr/bin/java"):
            self.assertEqual(
                cliquesnv_launcher("/apps/clique-snv.jar", 2048),
                (
                    "/usr/bin/java", "-Xms512m", "-Xmx2048m", "-jar",
                    "/apps/clique-snv.jar",
                ),
            )

    def test_species_fastas_contain_consensuses_for_reported_direct_calls(self):
        calls = pd.DataFrame([
            {
                "sample_id": "sample", "species": "Species A", "reported": True,
                "call_method": "exact_allele_profile",
            },
            {
                "sample_id": "sample", "species": "Species B", "reported": True,
                "call_method": "em_fallback",
            },
        ])
        groups = pd.DataFrame([
            {
                "locus_id": "BACT000001",
                "consensus_group_id": "BACT000001_CG0001",
                "consensus_status": "complete_exact_match",
                "consensus": "ACGTGT",
                "compatible_species": "Species A;Species B",
                "compatible_alleles": "BACT000001_7",
            },
            {
                "locus_id": "BACT000002",
                "consensus_group_id": "BACT000002_CG0001",
                "consensus_status": "unique_partial_match",
                "consensus": "ACNNGT",
                "compatible_species": "Species A",
                "compatible_alleles": "BACT000002_8",
            },
            {
                "locus_id": "BACT000003",
                "consensus_group_id": "BACT000003_CG0001",
                "consensus_status": "ambiguous_partial_match",
                "consensus": "NNGTGT",
                "compatible_species": "Species A",
                "compatible_alleles": "BACT000003_9;BACT000003_10;BACT000003_11",
            },
        ])
        with TemporaryDirectory() as directory:
            updated, manifest = write_direct_species_alleles(
                calls, groups, Path(directory)
            )
            direct = updated.loc[updated["species"].eq("Species A")].iloc[0]
            em = updated.loc[updated["species"].eq("Species B")].iloc[0]
            fasta = Path(direct["allele_fasta"])
            self.assertTrue(fasta.is_file())
            contents = fasta.read_text()
            self.assertIn(">BACT000001_CG0001 locus=BACT000001 species=Species A", contents)
            self.assertIn("database_match=exact allele=BACT000001_7", contents)
            self.assertIn("database_match=partial allele=BACT000002_8", contents)
            self.assertIn("database_match=partial allele_count=3", contents)
            self.assertIn("ACNNGT", contents)
            self.assertEqual(int(direct["allele_fasta_allele_count"]), 3)
            self.assertEqual(em["allele_fasta"], "")
            self.assertEqual(manifest["species"].tolist(), ["Species A"])
            self.assertEqual(
                manifest["alleles"].tolist(),
                ["BACT000001_CG0001;BACT000002_CG0001;BACT000003_CG0001"],
            )

    def test_consensus_workers_share_the_global_thread_budget(self):
        clusters = pd.DataFrame([
            {
                "sample_id": "sample", "locus_id": "BACT000001",
                "runtime_cluster_id": f"BACT000001_RC{i:04d}",
                "representative": f"BACT000001_{i}", "est_counts": 1.0,
            }
            for i in range(1, 5)
        ])
        assignments = pd.DataFrame([
            {
                "runtime_cluster_id": row.runtime_cluster_id,
                "Member_Template": row.representative,
            }
            for row in clusters.itertuples(index=False)
        ])
        records = {
            row.representative: FastaRecord(row.representative, row.representative, "ACGT")
            for row in clusters.itertuples(index=False)
        }
        active = 0
        maximum = 0
        lock = Lock()

        def fake_generate(row, *args):
            nonlocal active, maximum
            with lock:
                active += 1
                maximum = max(maximum, active)
            sleep(0.03)
            with lock:
                active -= 1
            cluster = str(row.runtime_cluster_id)
            representative = str(row.representative)
            return [{
                "sample_id": "sample", "locus_id": "BACT000001",
                "runtime_cluster_id": cluster, "representative": representative,
                "source_members": (representative,), "est_counts": 1.0,
                "assigned_fragments": 0, "exclusive_fragments": 0,
                "ambiguous_fragments": 0, "status": "insufficient_reads",
                "consensus": "", "callable_bases": 0, "allele_length": 4,
                "callable_breadth": 0.0, "compatible_alleles": (),
                "allele_search_scope": "not_searched", "consensus_path": None,
                "blast_path": None, "directory": Path("/tmp") / cluster,
            }]

        config = Config(
            allele_target=Path("reference.fasta"), seqdef=Path("seqdef.tsv"),
            reads=(Path("r1.fastq"), Path("r2.fastq")), sample_id="sample",
            output=Path("output"), tempdir=Path("work"), threads=4,
            memory_limit_gb=4,
            bwa_mem2="/usr/bin/true", samtools="/usr/bin/true",
            cliquesnv="/usr/bin/true", viral_consensus="/usr/bin/true",
        )
        with patch("bridge.rmlst.consensus._generate_cluster_consensus", side_effect=fake_generate):
            calls = run_cluster_consensus(
                clusters, assignments, {}, pd.DataFrame(), records, {},
                Path("reference.fasta"), config,
            )
        self.assertEqual(len(calls), 4)
        self.assertEqual(maximum, 2)

    def test_subconsensus_proportion_is_strictly_over_one_percent(self):
        evidence = AlignmentEvidence(
            depth=(4000, 4000, 4000, 4000),
            bases=tuple({} for _ in range(4)),
        )
        at_floor = evaluate_haplotype(
            HaplotypeCandidate("h1", "ACGT", 0.01, Path("calls.fasta")),
            "ACGT", evidence, 2000,
            min_proportion=0.01, min_supporting_fragments=20,
            min_depth=10, min_breadth=0.8,
        )
        above_floor = evaluate_haplotype(
            HaplotypeCandidate("h2", "ACGT", 0.0101, Path("calls.fasta")),
            "ACGT", evidence, 2000,
            min_proportion=0.01, min_supporting_fragments=20,
            min_depth=10, min_breadth=0.8,
        )
        self.assertFalse(at_floor.retained)
        self.assertIn("proportion_not_over_threshold", at_floor.reasons)
        self.assertTrue(above_floor.retained)

    def test_above_one_percent_is_still_rejected_when_depth_is_insufficient(self):
        evidence = AlignmentEvidence(
            depth=(400, 400, 400, 400),
            bases=tuple({} for _ in range(4)),
        )
        evaluated = evaluate_haplotype(
            HaplotypeCandidate("h1", "ACGT", 0.02, Path("calls.fasta")),
            "ACGT", evidence, 2000,
            min_proportion=0.01, min_supporting_fragments=20,
            min_depth=10, min_breadth=0.8,
        )
        self.assertFalse(evaluated.retained)
        self.assertNotIn("proportion_not_over_threshold", evaluated.reasons)
        self.assertIn("insufficient_projected_depth", evaluated.reasons)
        self.assertIn("insufficient_projected_breadth", evaluated.reasons)

    def test_tenfold_projected_depth_meets_the_default_floor(self):
        evidence = AlignmentEvidence(
            depth=(500, 500, 500, 500),
            bases=tuple({} for _ in range(4)),
        )
        evaluated = evaluate_haplotype(
            HaplotypeCandidate("h1", "ACGT", 0.02, Path("calls.fasta")),
            "ACGT", evidence, 2000,
            min_proportion=0.01, min_supporting_fragments=20,
            min_depth=10, min_breadth=0.8,
        )
        self.assertTrue(evaluated.retained)
        self.assertEqual(evaluated.projected_median_depth, 10)

    def test_subconsensus_requires_absolute_fragment_support_independently(self):
        evidence = AlignmentEvidence(
            depth=(4000, 4000, 4000, 4000),
            bases=tuple({} for _ in range(4)),
        )
        evaluated = evaluate_haplotype(
            HaplotypeCandidate("h1", "ACGT", 0.02, Path("calls.fasta")),
            "ACGT", evidence, 500,
            min_proportion=0.01, min_supporting_fragments=20,
            min_depth=10, min_breadth=0.8,
        )
        self.assertFalse(evaluated.retained)
        self.assertEqual(evaluated.supporting_fragments, 10)
        self.assertIn("insufficient_supporting_fragments", evaluated.reasons)

    def test_variant_requires_both_strands_and_interior_read_support(self):
        weak_bases = tuple({} for _ in range(4))
        weak_bases[1]["C"] = BaseEvidence(total=40, forward=40, reverse=0, interior=40)
        strong_bases = tuple({} for _ in range(4))
        strong_bases[1]["C"] = BaseEvidence(total=40, forward=20, reverse=20, interior=30)
        candidate = HaplotypeCandidate("h1", "ACGT", 0.02, Path("calls.fasta"))
        arguments = dict(
            min_proportion=0.01, min_supporting_fragments=20,
            min_depth=10, min_breadth=0.8,
        )
        weak = evaluate_haplotype(
            candidate, "AAGT", AlignmentEvidence((4000,) * 4, weak_bases),
            2000, **arguments,
        )
        strong = evaluate_haplotype(
            candidate, "AAGT", AlignmentEvidence((4000,) * 4, strong_bases),
            2000, **arguments,
        )
        self.assertFalse(weak.retained)
        self.assertIn("weak_variant_evidence_at_2", weak.reasons)
        self.assertTrue(strong.retained)

    def test_weak_variants_are_masked_and_terminal_truncations_are_padded(self):
        self.assertEqual(
            mask_weak_variants(
                "ACGT", ("weak_variant_evidence_at_2", "unrelated_reason")
            ),
            "ANGT",
        )
        self.assertEqual(
            pad_terminal_consensus("NCG", {"BACT000001_1": "ACGT"}),
            "NCGN",
        )
        self.assertEqual(
            pad_terminal_consensus("AC", {"BACT000001_1": "ACGT"}, max_padding=1),
            "AC",
        )

    def test_cliquesnv_frequency_and_sam_evidence_are_parsed(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            fasta = root / "cs01_region.fasta"
            fasta.write_text(
                ">gapped_0_0.750000000\nAC-G\n"
                ">sample_0_0.010100000\nACGT\n"
                ">reconstructed_1 Freq=0.250000\nACGA\n",
                encoding="utf-8",
            )
            with self.assertLogs("cliquesnv", level="WARNING") as captured:
                candidates = read_cliquesnv_haplotypes([fasta])
            self.assertEqual(len(candidates), 2)
            self.assertAlmostEqual(candidates[0].proportion, 0.0101)
            self.assertAlmostEqual(candidates[1].proportion, 0.25)
            self.assertIn("reason=invalid_sequence", "\n".join(captured.output))

            sam = root / "alignment.sam"
            sam.write_text(
                "@SQ\tSN:ref\tLN:4\n"
                "forward\t0\tref\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
                "reverse\t16\tref\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n",
                encoding="utf-8",
            )
            evidence = alignment_evidence_from_sam(
                sam, 4, min_end_distance=0
            )
            self.assertEqual(evidence.depth, (2, 2, 2, 2))
            self.assertEqual(evidence.bases[1]["C"].forward, 1)
            self.assertEqual(evidence.bases[1]["C"].reverse, 1)

    def test_minor_allele_gate_uses_depth_only(self):
        bases = tuple({} for _ in range(2))
        bases[0]["A"] = BaseEvidence(
            total=9997, forward=9997, reverse=0, interior=0
        )
        bases[0]["C"] = BaseEvidence(
            total=3, forward=3, reverse=0, interior=0
        )
        bases[1]["G"] = BaseEvidence(total=100, forward=50, reverse=50, interior=90)
        evidence = AlignmentEvidence((10000, 100), bases)

        self.assertEqual(maximum_minor_allele_depth(evidence), 3)

    def test_no_minor_allele_routes_to_viralconsensus(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            cluster_id = "BACT000001_RC0001"
            clusters = pd.DataFrame([{
                "sample_id": "sample", "locus_id": "BACT000001",
                "runtime_cluster_id": cluster_id,
                "representative": "BACT000001_1", "est_counts": 100.0,
            }])
            assignments = pd.DataFrame([{
                "runtime_cluster_id": cluster_id,
                "Member_Template": "BACT000001_1",
            }])
            records = {
                "BACT000001_1": FastaRecord(
                    "BACT000001_1", "BACT000001_1", "ACGT"
                ),
            }
            read_stats = pd.DataFrame([{
                "runtime_cluster_id": cluster_id, "assigned_fragments": 100,
                "exclusive_fragments": 100, "ambiguous_fragments": 0,
            }])
            config = Config(
                allele_target=root / "reference.fasta", seqdef=root / "seqdef.tsv",
                reads=(root / "r1.fastq", root / "r2.fastq"), sample_id="sample",
                output=root / "output", tempdir=root / "work", threads=2,
                memory_limit_gb=4,
                bwa_mem2="/usr/bin/true", samtools="/usr/bin/true",
                cliquesnv="/usr/bin/true", viral_consensus="/usr/bin/true",
                blastn="/usr/bin/true", makeblastdb="/usr/bin/true",
            )
            bases = tuple({} for _ in range(4))
            for index, base in enumerate("ACGT"):
                bases[index][base] = BaseEvidence(total=100)
            evidence = AlignmentEvidence((100,) * 4, bases)
            called_tools = []

            def fake_run(command, tool, *, stdout_path=None):
                called_tools.append(tool)
                if tool == "viralconsensus":
                    output = Path(command[command.index("-o") + 1])
                    output.write_text(">consensus\nACGT\n", encoding="utf-8")

            with (
                patch("bridge.rmlst.consensus.run", side_effect=fake_run),
                patch(
                    "bridge.rmlst.consensus.alignment_evidence_from_sam",
                    return_value=evidence,
                ),
            ):
                calls = run_cluster_consensus(
                    clusters, assignments,
                    {cluster_id: (root / "cluster_1.fastq", root / "cluster_2.fastq")},
                    read_stats, records, {"BACT000001_1": {"Species A"}},
                    config.allele_target, config,
                )

            self.assertEqual(len(calls), 1)
            self.assertEqual(calls[0].consensus_method, "viralconsensus")
            self.assertEqual(calls[0].status, "complete_exact_match")
            self.assertIn("viralconsensus", called_tools)
            self.assertNotIn("cliquesnv", called_tools)
            self.assertNotIn("samtools", called_tools)

            mixed_bases = tuple({} for _ in range(4))
            mixed_bases[0]["A"] = BaseEvidence(total=95)
            mixed_bases[0]["C"] = BaseEvidence(total=5)
            forced_evidence = AlignmentEvidence((100,) * 4, mixed_bases)
            forced_config = replace(
                config,
                tempdir=root / "work_forced",
                cliquesnv="/not/installed/cliquesnv",
                viralconsensus_only=True,
            )
            called_tools.clear()
            with (
                patch("bridge.rmlst.consensus.run", side_effect=fake_run),
                patch(
                    "bridge.rmlst.consensus.alignment_evidence_from_sam",
                    return_value=forced_evidence,
                ),
            ):
                forced_calls = run_cluster_consensus(
                    clusters, assignments,
                    {cluster_id: (
                        root / "cluster_1.fastq", root / "cluster_2.fastq"
                    )},
                    read_stats, records, {"BACT000001_1": {"Species A"}},
                    forced_config.allele_target, forced_config,
                )

            self.assertEqual(forced_calls[0].consensus_method, "viralconsensus")
            self.assertIn("viralconsensus", called_tools)
            self.assertNotIn("cliquesnv", called_tools)

    def test_cliquesnv_non_subconsensus_results_fall_back_to_viralconsensus(self):
        scenarios = {
            "single_full_proportion": (
                ">cluster_0_1.000000000\nACGT\n",
                "reason=cliquesnv_single_full_proportion_haplotype",
            ),
            "no_retained_haplotype": (
                ">cluster_0_0.990000000\nAACGT\n"
                ">cluster_1_0.010000000\nAAGT\n",
                "reason=no_cliquesnv_haplotype_passed_filters",
            ),
        }
        for name, (clique_fasta, expected_log) in scenarios.items():
            with self.subTest(name=name), TemporaryDirectory() as directory:
                root = Path(directory)
                cluster_id = "BACT000001_RC0001"
                clusters = pd.DataFrame([{
                    "sample_id": "sample", "locus_id": "BACT000001",
                    "runtime_cluster_id": cluster_id,
                    "representative": "BACT000001_1", "est_counts": 100.0,
                }])
                assignments = pd.DataFrame([{
                    "runtime_cluster_id": cluster_id,
                    "Member_Template": "BACT000001_1",
                }])
                records = {
                    "BACT000001_1": FastaRecord(
                        "BACT000001_1", "BACT000001_1", "AAGT"
                    ),
                }
                read_stats = pd.DataFrame([{
                    "runtime_cluster_id": cluster_id,
                    "assigned_fragments": 100,
                    "exclusive_fragments": 100,
                    "ambiguous_fragments": 0,
                }])
                config = Config(
                    allele_target=root / "reference.fasta",
                    seqdef=root / "seqdef.tsv",
                    reads=(root / "r1.fastq", root / "r2.fastq"),
                    sample_id="sample", output=root / "output",
                    tempdir=root / "work", threads=2, memory_limit_gb=4,
                    bwa_mem2="/usr/bin/true", samtools="/usr/bin/true",
                    cliquesnv="/usr/bin/true",
                    viral_consensus="/usr/bin/true",
                    blastn="/usr/bin/true", makeblastdb="/usr/bin/true",
                )
                bases = tuple({} for _ in range(4))
                bases[1]["A"] = BaseEvidence(total=95)
                bases[1]["C"] = BaseEvidence(total=5)
                evidence = AlignmentEvidence((100,) * 4, bases)
                called_tools = []

                def fake_run(command, tool, *, stdout_path=None):
                    called_tools.append(tool)
                    if tool == "cliquesnv":
                        outdir = Path(command[command.index("-outDir") + 1])
                        outdir.mkdir(parents=True, exist_ok=True)
                        (outdir / "cs01_region.fasta").write_text(
                            clique_fasta, encoding="utf-8"
                        )
                    elif tool == "viralconsensus":
                        output = Path(command[command.index("-o") + 1])
                        output.write_text(
                            ">consensus\nAAGT\n", encoding="utf-8"
                        )

                with (
                    patch("bridge.rmlst.consensus.run", side_effect=fake_run),
                    patch(
                        "bridge.rmlst.consensus.alignment_evidence_from_sam",
                        return_value=evidence,
                    ),
                    self.assertLogs("bridge.rmlst", level="INFO") as captured,
                ):
                    calls = run_cluster_consensus(
                        clusters, assignments,
                        {cluster_id: (
                            root / "cluster_1.fastq",
                            root / "cluster_2.fastq",
                        )},
                        read_stats, records,
                        {"BACT000001_1": {"Species A"}},
                        config.allele_target, config,
                    )

                self.assertEqual(len(calls), 1)
                self.assertEqual(calls[0].consensus, "AAGT")
                self.assertEqual(calls[0].consensus_method, "viralconsensus")
                self.assertEqual(calls[0].status, "complete_exact_match")
                self.assertEqual(called_tools.count("cliquesnv"), 1)
                self.assertEqual(called_tools.count("viralconsensus"), 1)
                self.assertTrue(any(
                    expected_log in message for message in captured.output
                ))

    def test_retained_haplotypes_become_separate_weighted_consensus_calls(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            cluster_id = "BACT000001_RC0001"
            clusters = pd.DataFrame([{
                "sample_id": "sample", "locus_id": "BACT000001",
                "runtime_cluster_id": cluster_id,
                "representative": "BACT000001_1", "est_counts": 100.0,
            }])
            assignments = pd.DataFrame([
                {"runtime_cluster_id": cluster_id, "Member_Template": "BACT000001_1"},
                {"runtime_cluster_id": cluster_id, "Member_Template": "BACT000001_2"},
            ])
            records = {
                "BACT000001_1": FastaRecord("BACT000001_1", "BACT000001_1", "AAGT"),
                "BACT000001_2": FastaRecord("BACT000001_2", "BACT000001_2", "ACGT"),
            }
            read_stats = pd.DataFrame([{
                "runtime_cluster_id": cluster_id, "assigned_fragments": 2000,
                "exclusive_fragments": 1600, "ambiguous_fragments": 400,
            }])
            config = Config(
                allele_target=root / "reference.fasta", seqdef=root / "seqdef.tsv",
                reads=(root / "r1.fastq", root / "r2.fastq"), sample_id="sample",
                output=root / "output", tempdir=root / "work", threads=2,
                memory_limit_gb=4,
                bwa_mem2="/usr/bin/true", samtools="/usr/bin/true",
                cliquesnv="/usr/bin/true", viral_consensus="/usr/bin/true",
                cliquesnv_time_limit_seconds=37,
                blastn="/usr/bin/true",
                makeblastdb="/usr/bin/true",
            )
            strong_bases = tuple({} for _ in range(4))
            strong_bases[1]["A"] = BaseEvidence(
                total=3920, forward=1960, reverse=1960, interior=3500
            )
            strong_bases[1]["C"] = BaseEvidence(
                total=80, forward=40, reverse=40, interior=70
            )
            evidence = AlignmentEvidence((4000,) * 4, strong_bases)

            cliquesnv_commands = []

            def fake_run(command, tool, *, stdout_path=None):
                if tool == "cliquesnv":
                    cliquesnv_commands.append(command)
                    outdir = Path(command[command.index("-outDir") + 1])
                    outdir.mkdir(parents=True, exist_ok=True)
                    (outdir / "cs01_region.fasta").write_text(
                        ">cluster_0_0.980000000\nAAGT\n"
                        ">cluster_1_0.020000000\nACGT\n",
                        encoding="utf-8",
                    )

            with (
                patch("bridge.rmlst.consensus.run", side_effect=fake_run),
                patch(
                    "bridge.rmlst.consensus.alignment_evidence_from_sam",
                    return_value=evidence,
                ),
            ):
                calls = run_cluster_consensus(
                    clusters, assignments,
                    {cluster_id: (root / "cluster_1.fastq", root / "cluster_2.fastq")},
                    read_stats, records,
                    {"BACT000001_1": {"Species A"}, "BACT000001_2": {"Species B"}},
                    config.allele_target, config,
                )

            self.assertEqual(len(calls), 2)
            self.assertEqual(
                [call.subconsensus_id for call in calls],
                [f"{cluster_id}_SC0001", f"{cluster_id}_SC0002"],
            )
            self.assertEqual([call.status for call in calls], ["complete_exact_match"] * 2)
            self.assertAlmostEqual(calls[0].est_counts, 98.0)
            self.assertAlmostEqual(calls[1].est_counts, 2.0)
            self.assertAlmostEqual(calls[1].supporting_fragments, 40.0)
            self.assertTrue(calls[1].variant_evidence_passed)
            self.assertEqual(len(cliquesnv_commands), 1)
            self.assertEqual(
                cliquesnv_commands[0][cliquesnv_commands[0].index("-tl") + 1],
                "37",
            )

    def test_exact_cluster_member_bypasses_medoid_length_validation(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            cluster_id = "BACT000001_RC0001"
            clusters = pd.DataFrame([{
                "sample_id": "sample", "locus_id": "BACT000001",
                "runtime_cluster_id": cluster_id,
                "representative": "BACT000001_1", "est_counts": 100.0,
            }])
            assignments = pd.DataFrame([
                {"runtime_cluster_id": cluster_id, "Member_Template": "BACT000001_1"},
                {"runtime_cluster_id": cluster_id, "Member_Template": "BACT000001_2"},
            ])
            records = {
                "BACT000001_1": FastaRecord("BACT000001_1", "BACT000001_1", "AAGT"),
                "BACT000001_2": FastaRecord("BACT000001_2", "BACT000001_2", "AACGT"),
            }
            read_stats = pd.DataFrame([{
                "runtime_cluster_id": cluster_id, "assigned_fragments": 100,
                "exclusive_fragments": 100, "ambiguous_fragments": 0,
            }])
            config = Config(
                allele_target=root / "reference.fasta", seqdef=root / "seqdef.tsv",
                reads=(root / "r1.fastq", root / "r2.fastq"), sample_id="sample",
                output=root / "output", tempdir=root / "work", threads=2,
                memory_limit_gb=4,
                bwa_mem2="/usr/bin/true", samtools="/usr/bin/true",
                cliquesnv="/usr/bin/true", viral_consensus="/usr/bin/true",
                blastn="/usr/bin/true",
                makeblastdb="/usr/bin/true",
            )
            mixed_bases = tuple({} for _ in range(4))
            mixed_bases[0]["A"] = BaseEvidence(total=95)
            mixed_bases[0]["C"] = BaseEvidence(total=5)
            evidence = AlignmentEvidence((100,) * 4, mixed_bases)

            def fake_run(command, tool, *, stdout_path=None):
                if tool == "cliquesnv":
                    outdir = Path(command[command.index("-outDir") + 1])
                    outdir.mkdir(parents=True, exist_ok=True)
                    (outdir / "cs01_region.fasta").write_text(
                        ">cluster_0_0.980000000\nAACGT\n"
                        ">cluster_1_0.020000000\nAAGT\n",
                        encoding="utf-8",
                    )

            with (
                patch("bridge.rmlst.consensus.run", side_effect=fake_run),
                patch(
                    "bridge.rmlst.consensus.alignment_evidence_from_sam",
                    return_value=evidence,
                ),
            ):
                calls = run_cluster_consensus(
                    clusters, assignments,
                    {cluster_id: (root / "cluster_1.fastq", root / "cluster_2.fastq")},
                    read_stats, records,
                    {"BACT000001_1": {"Species A"}, "BACT000001_2": {"Species A"}},
                    config.allele_target, config,
                )

            self.assertEqual(len(calls), 1)
            self.assertEqual(calls[0].status, "complete_exact_match")
            self.assertEqual(calls[0].compatible_alleles, ("BACT000001_2",))
            self.assertEqual(calls[0].allele_length, 5)
            self.assertTrue(calls[0].variant_evidence_passed)

    def test_exact_match_bypass_still_requires_depth_and_support(self):
        evaluated = evaluate_haplotype(
            HaplotypeCandidate("h1", "AACGT", 0.5, Path("calls.fasta")),
            "AAGT", AlignmentEvidence((5,) * 4, tuple({} for _ in range(4))),
            20, min_proportion=0.01, min_supporting_fragments=20,
            min_depth=10, min_breadth=0.8, validate_variants=False,
        )
        self.assertFalse(evaluated.retained)
        self.assertIn("insufficient_supporting_fragments", evaluated.reasons)
        self.assertIn("insufficient_projected_depth", evaluated.reasons)
        self.assertIn("insufficient_projected_breadth", evaluated.reasons)
        self.assertNotIn("unphaseable_length_change", evaluated.reasons)

    def test_weak_cluster_variant_is_masked_before_allele_search(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            cluster_id = "BACT000001_RC0001"
            clusters = pd.DataFrame([{
                "sample_id": "sample", "locus_id": "BACT000001",
                "runtime_cluster_id": cluster_id,
                "representative": "BACT000001_1", "est_counts": 100.0,
            }])
            assignments = pd.DataFrame([
                {"runtime_cluster_id": cluster_id, "Member_Template": "BACT000001_1"},
                {"runtime_cluster_id": cluster_id, "Member_Template": "BACT000001_2"},
            ])
            records = {
                "BACT000001_1": FastaRecord("BACT000001_1", "BACT000001_1", "AAGT"),
                "BACT000001_2": FastaRecord("BACT000001_2", "BACT000001_2", "ATGT"),
            }
            read_stats = pd.DataFrame([{
                "runtime_cluster_id": cluster_id, "assigned_fragments": 100,
                "exclusive_fragments": 100, "ambiguous_fragments": 0,
            }])
            config = Config(
                allele_target=root / "reference.fasta", seqdef=root / "seqdef.tsv",
                reads=(root / "r1.fastq", root / "r2.fastq"), sample_id="sample",
                output=root / "output", tempdir=root / "work", threads=2,
                memory_limit_gb=4,
                bwa_mem2="/usr/bin/true", samtools="/usr/bin/true",
                cliquesnv="/usr/bin/true", viral_consensus="/usr/bin/true",
                blastn="/usr/bin/true",
                makeblastdb="/usr/bin/true",
            )
            weak_bases = tuple({} for _ in range(4))
            weak_bases[1]["A"] = BaseEvidence(
                total=60, forward=30, reverse=30, interior=50
            )
            weak_bases[1]["C"] = BaseEvidence(
                total=40, forward=40, reverse=0, interior=40
            )
            evidence = AlignmentEvidence((100,) * 4, weak_bases)

            def fake_run(command, tool, *, stdout_path=None):
                if tool == "cliquesnv":
                    outdir = Path(command[command.index("-outDir") + 1])
                    outdir.mkdir(parents=True, exist_ok=True)
                    (outdir / "cs01_region.fasta").write_text(
                        ">cluster_0_0.980000000\nACGT\n"
                        ">cluster_1_0.020000000\nAAGT\n",
                        encoding="utf-8",
                    )

            def fake_search(pending, *args):
                self.assertEqual(len(pending), 1)
                self.assertEqual(pending[0]["consensus"], "ANGT")
                pending[0]["compatible_alleles"] = ("BACT000001_2",)
                pending[0]["status"] = "unique_partial_match"

            with (
                patch("bridge.rmlst.consensus.run", side_effect=fake_run),
                patch(
                    "bridge.rmlst.consensus.alignment_evidence_from_sam",
                    return_value=evidence,
                ),
                patch(
                    "bridge.rmlst.consensus._batch_full_reference_search",
                    side_effect=fake_search,
                ),
            ):
                calls = run_cluster_consensus(
                    clusters, assignments,
                    {cluster_id: (root / "cluster_1.fastq", root / "cluster_2.fastq")},
                    read_stats, records,
                    {"BACT000001_1": {"Species A"}, "BACT000001_2": {"Species A"}},
                    config.allele_target, config,
                )

            self.assertEqual(len(calls), 1)
            self.assertEqual(calls[0].status, "unique_partial_match")
            self.assertEqual(calls[0].consensus, "ANGT")
            self.assertTrue(calls[0].variant_evidence_passed)

    def test_vectorized_bootstrap_candidate_parsing(self):
        with TemporaryDirectory() as directory:
            quant = Path(directory)
            header = "target_id\tlength\teff_length\test_counts\ttpm\n"
            (quant / "abundance.tsv").write_text(
                header
                + "BACT000001_1\t400\t300\t8\t80\n"
                + "BACT000001_2\t400\t300\t2\t20\n",
                encoding="utf-8",
            )
            (quant / "bs_abundance_0.tsv").write_text(
                header
                + "BACT000001_1\t400\t300\t10\t100\n"
                + "BACT000001_2\t400\t300\t0\t0\n",
                encoding="utf-8",
            )
            (quant / "bs_abundance_1.tsv").write_text(
                header
                + "BACT000001_1\t400\t300\t6\t75\n"
                + "BACT000001_2\t400\t300\t2\t25\n",
                encoding="utf-8",
            )
            config = SimpleNamespace(
                sample_id="sample", bootstraps=2, support_threshold=2,
                use_bootstrap_mean=True,
            )
            calls = derive_candidates(config, quant)
            self.assertEqual(calls["Template"].tolist(), ["BACT000001_1"])
            self.assertEqual(float(calls.loc[0, "est_counts"]), 8.0)
            self.assertEqual(int(calls.loc[0, "bootstrap_support"]), 2)
            self.assertEqual(int(calls.loc[0, "effective_bootstrap_threshold"]), 2)
            self.assertFalse(bool(calls.loc[0, "bootstrap_threshold_relaxed"]))

    def test_empty_loci_relax_to_their_highest_bootstrap_support(self):
        with TemporaryDirectory() as directory:
            quant = Path(directory)
            header = "target_id\tlength\teff_length\test_counts\ttpm\n"
            targets = [
                ("BACT000001_1", 10, 50), ("BACT000001_2", 2, 10),
                ("BACT000002_1", 5, 25), ("BACT000002_2", 4, 20),
                ("BACT000003_1", 1, 5), ("BACT000003_2", 0, 0),
            ]
            (quant / "abundance.tsv").write_text(
                header + "".join(
                    f"{target}\t400\t300\t{count}\t{tpm}\n"
                    for target, count, tpm in targets
                ),
                encoding="utf-8",
            )
            bootstrap_counts = [
                (10, 1, 5, 4, 0, 0),
                (8, 0, 4, 3, 0, 0),
                (9, 0, 0, 0, 0, 0),
            ]
            for replicate, counts in enumerate(bootstrap_counts):
                (quant / f"bs_abundance_{replicate}.tsv").write_text(
                    header + "".join(
                        f"{target}\t400\t300\t{count}\t0\n"
                        for (target, _, _), count in zip(targets, counts)
                    ),
                    encoding="utf-8",
                )
            config = SimpleNamespace(
                sample_id="sample", bootstraps=3, support_threshold=3,
                use_bootstrap_mean=True,
            )
            calls = derive_candidates(config, quant).set_index("Template")

            self.assertEqual(
                set(calls.index),
                {"BACT000001_1", "BACT000002_1", "BACT000002_2", "BACT000003_1"},
            )
            self.assertEqual(
                int(calls.loc["BACT000001_1", "effective_bootstrap_threshold"]), 3
            )
            self.assertFalse(
                bool(calls.loc["BACT000001_1", "bootstrap_threshold_relaxed"])
            )
            self.assertEqual(
                int(calls.loc["BACT000002_1", "effective_bootstrap_threshold"]), 2
            )
            self.assertTrue(
                bool(calls.loc["BACT000002_1", "bootstrap_threshold_relaxed"])
            )
            self.assertEqual(
                int(calls.loc["BACT000003_1", "effective_bootstrap_threshold"]), 0
            )

    def test_rst_matches_are_limited_but_total_count_is_retained(self):
        with TemporaryDirectory() as directory:
            seqdef = Path(directory) / "seqdef.tsv"
            rows = ["rST\tBACT000001\tBACT000002\tspecies"]
            rows.extend(
                f"{number}\t7\t9\tSpecies A" for number in range(1, 13)
            )
            rows.append("99\t8\t9\tSpecies A")
            seqdef.write_text("\n".join(rows) + "\n", encoding="utf-8")
            calls = pd.DataFrame([{"species": "Species A", "reported": True}])
            groups = pd.DataFrame([
                {
                    "locus_id": "BACT000001",
                    "consensus_status": "unique_partial_match",
                    "compatible_alleles": "BACT000001_7",
                    "compatible_species": "Species A",
                },
                {
                    "locus_id": "BACT000002",
                    "consensus_status": "complete_exact_match",
                    "compatible_alleles": "BACT000002_9",
                    "compatible_species": "Species A",
                },
            ])
            matched = add_rst_matches(calls, groups, seqdef, limit=10)
            self.assertEqual(matched.loc[0, "matching_rST_count"], 12)
            self.assertEqual(
                matched.loc[0, "matching_rSTs"], "1;2;3;4;5;6;7;8;9;10"
            )

    def test_rst_matching_uses_species_assigned_options_at_each_locus(self):
        with TemporaryDirectory() as directory:
            seqdef = Path(directory) / "seqdef.tsv"
            seqdef.write_text(
                "rST\tBACT000001\tBACT000002\tspecies\n"
                "1\t1\t3\tSpecies A\n"
                "2\t2\t4\tSpecies A\n"
                "3\t1\t4\tSpecies A\n"
                "4\t9\t8\tSpecies B\n",
                encoding="utf-8",
            )
            calls = pd.DataFrame([{"species": "Species A", "reported": True}])
            groups = pd.DataFrame([
                {
                    "locus_id": "BACT000001",
                    "consensus_status": "ambiguous_partial_match",
                    "compatible_alleles": "BACT000001_1;BACT000001_2;BACT000001_9",
                    "compatible_species": "Species A;Species B",
                },
                {
                    "locus_id": "BACT000002",
                    "consensus_status": "ambiguous_partial_match",
                    "compatible_alleles": "BACT000002_3;BACT000002_4;BACT000002_8",
                    "compatible_species": "Species A;Species B",
                },
            ])
            matched = add_rst_matches(calls, groups, seqdef)
            self.assertEqual(matched.loc[0, "matching_rSTs"], "1;2;3")
            self.assertEqual(int(matched.loc[0, "matching_rST_count"]), 3)

    def test_rst_matching_respects_subspecies_and_blank_category(self):
        with TemporaryDirectory() as directory:
            seqdef = Path(directory) / "seqdef.tsv"
            seqdef.write_text(
                "rST\tBACT000001\tspecies\tsubspecies\n"
                "1\t7\tMycobacterium tuberculosis\tbovis\n"
                "2\t7\tMycobacterium tuberculosis\t\n",
                encoding="utf-8",
            )
            calls = pd.DataFrame([
                {
                    "species": "Mycobacterium tuberculosis subsp. bovis",
                    "reported": True,
                },
                {
                    "species": (
                        "Mycobacterium tuberculosis "
                        "[no subspecies annotation]"
                    ),
                    "reported": True,
                },
            ])
            groups = pd.DataFrame([
                {
                    "locus_id": "BACT000001",
                    "consensus_status": "complete_exact_match",
                    "compatible_alleles": "BACT000001_7",
                    "compatible_species": (
                        "Mycobacterium tuberculosis subsp. bovis;"
                        "Mycobacterium tuberculosis "
                        "[no subspecies annotation]"
                    ),
                },
            ])

            matched = add_rst_matches(
                calls, groups, seqdef, distinguish_subspecies=True
            )

            self.assertEqual(matched["matching_rSTs"].tolist(), ["1", "2"])

    def test_partial_consensus_requires_all_callable_bases(self):
        found = compatible_alleles(
            "ACNNGT",
            {"L_1": "ACCCGT", "L_2": "ACTTGT", "L_3": "ATCCGT"},
        )
        self.assertEqual(found, ("L_1", "L_2"))

    def test_partial_consensus_accepts_bounded_terminal_allele_extensions(self):
        consensus = "NACGTN"
        found = compatible_alleles(
            consensus,
            {
                "L_exact": "AACGTA",
                "L_left_extension": "TTTAACGTA",
                "L_right_extension_15": "AACGTA" + "T" * 15,
                "L_too_long": "AACGTA" + "T" * 16,
                "L_internal_mismatch": "AATGTA",
            },
        )
        self.assertEqual(
            found,
            ("L_exact", "L_left_extension", "L_right_extension_15"),
        )

    def test_shorter_allele_can_only_trim_non_callable_consensus_ends(self):
        alleles = {"L_short": "AACGT"}
        self.assertEqual(
            compatible_alleles("NACGTNNN", alleles),
            ("L_short",),
        )
        self.assertEqual(
            compatible_alleles("NACGTAAA", alleles),
            (),
        )
        self.assertEqual(
            compatible_alleles(
                "NACGTNNN", alleles, max_terminal_difference=2
            ),
            (),
        )

    def test_terminal_compatibility_does_not_allow_internal_indels(self):
        self.assertEqual(
            compatible_alleles("NACGTN", {"L_internal_insertion": "AACCGTA"}),
            (),
        )

    def test_identical_consensus_uses_winner_assignment_and_sums_evidence(self):
        groups = merge_consensus_calls([
            call("L_RC1", "L", "ACGT", ["L_1"], ["Species A"], 4),
            call("L_RC2", "L", "ACGT", ["L_2"], ["Species B"], 6),
        ])
        self.assertEqual(len(groups), 1)
        self.assertEqual(groups.iloc[0]["merged_cluster_count"], 2)
        self.assertEqual(groups.iloc[0]["combined_est_counts"], 10)
        self.assertEqual(groups.iloc[0]["compatible_alleles"], "L_2")
        self.assertEqual(groups.iloc[0]["compatible_species"], "Species B")

    def test_distinct_consensuses_are_preserved(self):
        groups = merge_consensus_calls([
            call("L_RC1", "L", "ACGT", ["L_1"], ["Species A"]),
            call("L_RC2", "L", "ACGA", ["L_2"], ["Species B"]),
        ])
        self.assertEqual(len(groups), 2)

    def test_two_species_can_be_directly_supported_together(self):
        rows = []
        for locus in ("L1", "L2", "L3"):
            rows.extend([
                {
                    "locus_id": locus,
                    "consensus_status": "complete_exact_match",
                    "compatible_species": "Species A",
                    "combined_est_counts": 10.0,
                },
                {
                    "locus_id": locus,
                    "consensus_status": "complete_exact_match",
                    "compatible_species": "Species B",
                    "combined_est_counts": 5.0,
                },
            ])
        groups = pd.DataFrame(rows)
        self.assertEqual(direct_profile_species(groups, 3), {"Species A": 3, "Species B": 3})
        em = em_species(groups, 3).set_index("species")
        self.assertTrue(em.loc["Species A", "detected"])
        self.assertTrue(em.loc["Species B", "detected"])
        self.assertAlmostEqual(em.loc["Species A", "pi"], 2 / 3)

    def test_ambiguous_profile_uses_em_instead_of_direct_call(self):
        groups = pd.DataFrame([
            {
                "locus_id": locus,
                "consensus_status": "ambiguous_partial_match",
                "compatible_species": "Species A;Species B",
                "combined_est_counts": 10.0,
            }
            for locus in ("L1", "L2", "L3")
        ])
        self.assertEqual(direct_profile_species(groups, 3), {})
        em = em_species(groups, 3)
        self.assertEqual(set(em["species"]), {"Species A", "Species B"})
        self.assertFalse(em["detected"].any())

    def test_composite_consensus_is_counted_as_ambiguous_evidence(self):
        groups = pd.DataFrame([
            {
                "locus_id": "L1",
                "consensus_status": "complete_exact_match",
                "compatible_species": "Species A",
                "combined_est_counts": 10.0,
                "composite_from_members": True,
            },
            {
                "locus_id": "L2",
                "consensus_status": "complete_exact_match",
                "compatible_species": "Species A",
                "combined_est_counts": 5.0,
                "composite_from_members": False,
            },
        ])
        em = em_species(groups, 1).set_index("species").loc["Species A"]
        self.assertEqual(int(em["allele_evidence_unambiguous_loci"]), 1)
        self.assertEqual(int(em["allele_evidence_ambiguous_loci"]), 1)
        self.assertEqual(float(em["unambiguous_est_counts"]), 5.0)
        self.assertEqual(float(em["ambiguous_est_counts"]), 10.0)

    def test_em_output_includes_every_locus_used_by_em(self):
        groups = pd.DataFrame([
            {
                "locus_id": "L1", "consensus_status": "complete_exact_match",
                "compatible_species": "Species A", "combined_est_counts": 10.0,
            },
            {
                "locus_id": "L2", "consensus_status": "ambiguous_partial_match",
                "compatible_species": "Species A;Species B",
                "combined_est_counts": 10.0,
            },
            {
                "locus_id": "L3", "consensus_status": "no_compatible_match",
                "compatible_species": "Species A", "combined_est_counts": 10.0,
            },
            {
                "locus_id": "L4", "consensus_status": "alignment_failed",
                "compatible_species": "Species A;Species B",
                "combined_est_counts": 10.0,
            },
        ])

        em = em_species(groups, 1).set_index("species")
        species_a = em.loc["Species A"]
        self.assertEqual(int(species_a["allele_evidence_unambiguous_loci"]), 1)
        self.assertEqual(int(species_a["allele_evidence_ambiguous_loci"]), 1)
        self.assertEqual(int(species_a["em_evidence_unambiguous_loci"]), 2)
        self.assertEqual(int(species_a["em_evidence_ambiguous_loci"]), 2)

    def test_species_requires_an_unambiguous_locus_to_be_reported(self):
        groups = pd.DataFrame([
            {
                "locus_id": "L1",
                "consensus_status": "ambiguous_partial_match",
                "compatible_species": "Species A;Species B",
                "combined_est_counts": 10.0,
            },
            {
                "locus_id": "L2",
                "consensus_status": "ambiguous_partial_match",
                "compatible_species": "Species A;Species C",
                "combined_est_counts": 10.0,
            },
            {
                "locus_id": "L3",
                "consensus_status": "ambiguous_partial_match",
                "compatible_species": "Species A;Species D",
                "combined_est_counts": 10.0,
            },
        ])
        self.assertEqual(direct_profile_species(groups, 3), {})
        final, _ = final_species_calls("sample", groups, 3)
        species_a = final.loc[final["species"].eq("Species A")].iloc[0]
        self.assertEqual(
            int(species_a["allele_evidence_unambiguous_loci"]), 0
        )
        self.assertFalse(bool(species_a["reported"]))

    def test_empty_candidate_set_is_insufficient_material_not_a_crash(self):
        from bridge.rmlst.errors import InsufficientMaterial

        with self.assertRaises(InsufficientMaterial) as caught:
            raise InsufficientMaterial("no candidates", stage="kallisto")
        self.assertEqual(caught.exception.stage, "kallisto")
        # Existing handlers catch ValueError; that must keep working.
        self.assertIsInstance(caught.exception, ValueError)

    def test_run_records_per_tool_runtime(self):
        from bridge import tools

        tools.reset_runtime_accounting()
        self.assertEqual(tools.runtime_accounting()["tools"], {})
        tools.run(
            [sys.executable, "-c", "pass"], "probe_tool",
        )
        accounting = tools.runtime_accounting()
        self.assertEqual(accounting["tools"]["probe_tool"]["calls"], 1)
        self.assertGreaterEqual(accounting["tools"]["probe_tool"]["seconds"], 0.0)
        tools.reset_runtime_accounting()

    def test_failed_tool_still_records_its_runtime(self):
        """A tool that dies still consumed wall clock; the breakdown must show it."""
        from bridge import tools

        tools.reset_runtime_accounting()
        with self.assertRaises(subprocess.CalledProcessError):
            tools.run(
                [sys.executable, "-c", "import sys; sys.exit(3)"], "failing_tool",
            )
        self.assertEqual(
            tools.runtime_accounting()["tools"]["failing_tool"]["calls"], 1
        )
        tools.reset_runtime_accounting()

    def test_padding_is_skipped_when_nothing_is_longer_than_the_consensus(self):
        """The guard must not change the answer, only skip work that cannot help."""
        alleles = {"a1": "ACGTACGT", "a2": "ACGTACGA"}
        full_length = "ACGTACGT"
        self.assertEqual(
            pad_terminal_consensus(full_length, alleles), full_length
        )
        # A genuinely short consensus still pads.
        padded = pad_terminal_consensus("CGTACGT", alleles)
        self.assertEqual(len(padded), 8)
        self.assertIn("N", padded)

    def _mitis_group_geometry(self, bystander=False):
        """One dominant organism plus a near neighbour on a single exclusive locus.

        This is the shape that over-called Streptococcus pneumoniae on a
        synthetic mixture panel: the real organism carried many exclusive loci,
        while the neighbour had exactly one exclusive locus and rode the rest of
        the way on evidence that is merely *compatible* with it. The profile-wide
        intersection collapses to the neighbour because
        ``{X} & {X, Y} & {X, Y} == {X}``.

        With ``bystander=True`` a third mitis-group name rides along on the
        shared loci with **no** exclusive locus of its own - the shape of the
        EM-fallback background, which outnumbered the real demotions 94 to 1 on
        the panel.
        """
        shared = "Streptococcus mitis;Streptococcus pneumoniae"
        if bystander:
            shared = shared + ";Streptococcus pseudopneumoniae"
        rows = [
            # The organism actually present: four loci that name it alone.
            {
                "locus_id": f"L{index}", "consensus_status": "complete_exact_match",
                "compatible_species": "Streptococcus mitis",
                "combined_est_counts": 500.0,
            }
            for index in range(1, 5)
        ]
        rows.append({
            "locus_id": "L5", "consensus_status": "complete_exact_match",
            "compatible_species": "Streptococcus pneumoniae",
            "combined_est_counts": 1.0,
        })
        rows.extend([
            {
                "locus_id": f"L{index}", "consensus_status": "ambiguous_partial_match",
                "compatible_species": shared, "combined_est_counts": 5.0,
            }
            for index in range(6, 9)
        ])
        return pd.DataFrame(rows)

    def test_one_exclusive_locus_does_not_report_a_near_neighbour(self):
        groups = self._mitis_group_geometry()

        # The old bar - a single unambiguous locus - admits the neighbour.
        lenient = direct_profile_species(groups, 3, min_unambiguous_loci=1)
        self.assertIn("Streptococcus pneumoniae", lenient)

        # The default bar does not, and the real organism is untouched.
        direct = direct_profile_species(groups, 3)
        self.assertIn("Streptococcus mitis", direct)
        self.assertNotIn("Streptococcus pneumoniae", direct)

        final, _ = final_species_calls("sample", groups, 3)
        by_species = final.set_index("species")
        self.assertTrue(bool(by_species.loc["Streptococcus mitis", "reported"]))
        self.assertFalse(bool(by_species.loc["Streptococcus pneumoniae", "reported"]))

    def test_profile_loci_do_not_overstate_a_call(self):
        """exact_profile_loci counts compatibility; the unambiguous count is the evidence."""
        groups = self._mitis_group_geometry()
        final, _ = final_species_calls("sample", groups, 3)
        mitis = final.set_index("species").loc["Streptococcus mitis"]
        self.assertEqual(int(mitis["exact_profile_loci"]), 7)
        self.assertEqual(int(mitis["exact_profile_unambiguous_loci"]), 4)

    def test_demoted_near_neighbour_is_marked_as_an_alternative_name(self):
        """A second name for one organism must not read as a second organism."""
        groups = self._mitis_group_geometry()
        final, _ = final_species_calls("sample", groups, 3)
        neighbour = final.set_index("species").loc["Streptococcus pneumoniae"]
        self.assertFalse(bool(neighbour["reported"]))
        self.assertEqual(neighbour["alternative_to"], "Streptococcus mitis")
        reported = final.loc[final["reported"].astype(bool)]
        self.assertEqual(reported["alternative_to"].tolist(), [""] * len(reported))

    def test_compatible_only_species_is_not_an_alternative_name(self):
        """No exclusive locus means no second name - it is background, not a synonym.

        The panel separates the two classes completely: every real demotion
        carried exactly one exclusive locus, every EM-fallback row carried none.
        Labelling on a shared genus alone erased that distinction and would feed
        the reconciler 94 spurious GENUS_ONLY signals for each real one.
        """
        groups = self._mitis_group_geometry(bystander=True)

        # It never names itself alone, so it is not a profile candidate at all.
        self.assertNotIn(
            "Streptococcus pseudopneumoniae",
            demoted_profile_species(groups, 3),
        )
        self.assertNotIn(
            "Streptococcus pseudopneumoniae",
            direct_profile_species(groups, 3),
        )

        final, _ = final_species_calls("sample", groups, 3)
        by_species = final.set_index("species")

        # The bystander shares the genus with a reported S. mitis and must still
        # come back unlabelled...
        bystander = by_species.loc["Streptococcus pseudopneumoniae"]
        self.assertFalse(bool(bystander["reported"]))
        self.assertEqual(bystander["alternative_to"], "")

        # ...while the one exclusive locus still earns S. pneumoniae its label.
        neighbour = by_species.loc["Streptococcus pneumoniae"]
        self.assertFalse(bool(neighbour["reported"]))
        self.assertEqual(neighbour["alternative_to"], "Streptococcus mitis")

    def test_demoted_species_alone_in_its_genus_is_not_an_alternative_name(self):
        """An alternative name needs something to be an alternative *to*."""
        groups = pd.DataFrame([
            {
                "locus_id": "L1", "consensus_status": "complete_exact_match",
                "compatible_species": "Solo species", "combined_est_counts": 5.0,
            },
        ] + [
            {
                "locus_id": f"L{index}", "consensus_status": "ambiguous_partial_match",
                "compatible_species": "Solo species;Solo other",
                "combined_est_counts": 5.0,
            }
            for index in range(2, 4)
        ])
        # One exclusive locus: clears the profile geometry (3 compatible loci,
        # intersection collapses to this species) but falls short of the
        # reporting bar, which is what makes it a demotion rather than a call.
        self.assertIn("Solo species", demoted_profile_species(groups, 3))

        final, _ = final_species_calls("sample", groups, 3)
        row = final.set_index("species").loc["Solo species"]
        self.assertFalse(bool(row["reported"]))
        self.assertEqual(row["alternative_to"], "")

    def test_a_call_whose_every_locus_carried_no_reads_is_not_made(self):
        """Loci with zero estimated counts are not evidence.

        ``_profile_candidates`` applies the same ``combined_est_counts > 0``
        filter as ``em_species``; otherwise a negative control can report a
        species on loci carrying no reads. The trade-off is that a species whose
        every locus has zero estimated counts is not called, even if another
        arm corroborates it.
        """
        groups = pd.DataFrame([
            {
                "locus_id": f"L{index}", "consensus_status": "complete_exact_match",
                "compatible_species": "Species A", "combined_est_counts": 0.0,
            }
            for index in range(1, 4)
        ])
        final, em = final_species_calls("sample", groups, 3)
        self.assertNotIn("Species A", list(final["species"]))
        self.assertTrue(em.empty)

    def test_one_locus_with_reads_is_enough_to_keep_the_others_as_evidence(self):
        """The filter is per locus, so a real call is not lost to a quiet locus.

        Only the loci that carried nothing are discarded. A species with reads
        somewhere keeps every locus that had them, and the bar is applied to
        that count - which is what separates this from an all-or-nothing rule.
        """
        groups = pd.DataFrame(
            [
                {
                    "locus_id": f"L{index}",
                    "consensus_status": "complete_exact_match",
                    "compatible_species": "Species A",
                    "combined_est_counts": 12.0,
                }
                for index in range(1, 4)
            ]
            + [{
                "locus_id": "L4", "consensus_status": "complete_exact_match",
                "compatible_species": "Species A", "combined_est_counts": 0.0,
            }]
        )
        final, _ = final_species_calls("sample", groups, 3)
        species_a = final.set_index("species").loc["Species A"]
        self.assertTrue(bool(species_a["reported"]))
        self.assertEqual(int(species_a["exact_profile_unambiguous_loci"]), 3)

    def test_generic_species_uses_definitive_genus_candidates_as_ambiguous_evidence(self):
        groups = pd.DataFrame([
            {
                "locus_id": "L1", "consensus_status": "complete_exact_match",
                "compatible_species": "Neisseria meningitidis",
                "combined_est_counts": 1.0,
            },
            {
                "locus_id": "L2", "consensus_status": "complete_exact_match",
                "compatible_species": "Neisseria gonorrhoeae",
                "combined_est_counts": 1.0,
            },
            {
                "locus_id": "L3", "consensus_status": "complete_exact_match",
                "compatible_species": "Neisseria sp.",
                "combined_est_counts": 10.0,
            },
        ])
        em = em_species(groups, 1).set_index("species")
        self.assertNotIn("Neisseria sp.", em.index)
        self.assertEqual(
            int(em.loc[
                "Neisseria meningitidis", "allele_evidence_ambiguous_loci"
            ]),
            1,
        )
        self.assertEqual(
            int(em.loc[
                "Neisseria gonorrhoeae", "allele_evidence_ambiguous_loci"
            ]),
            1,
        )

    def test_generic_species_is_retained_without_definitive_genus_evidence(self):
        groups = pd.DataFrame([{
            "locus_id": "L1", "consensus_status": "complete_exact_match",
            "compatible_species": "Neisseria sp.", "combined_est_counts": 10.0,
        }])
        em = em_species(groups, 1)
        self.assertEqual(em["species"].tolist(), ["Neisseria sp."])

    def test_no_species_evidence_produces_empty_stable_outputs(self):
        groups = pd.DataFrame([{
            "locus_id": "L1", "consensus_status": "no_compatible_match",
            "compatible_species": "", "combined_est_counts": 10.0,
        }])
        final, em = final_species_calls("sample", groups, 3)
        self.assertTrue(final.empty)
        self.assertTrue(em.empty)
        self.assertIn("reported", final.columns)

    def test_zero_abundance_loci_are_not_evidence(self):
        """A locus the EM gave no reads cannot support a species call.

        This is what let a negative control report Escherichia coli: four loci
        named it, none carried a single estimated read. Consensus can be
        assembled for a locus that is then assigned no abundance - that says a
        sequence was seen, not that an organism was present.
        """
        groups = pd.DataFrame([
            {
                "locus_id": f"L{index}", "consensus_status": "complete_exact_match",
                "compatible_species": "Ghost species", "combined_est_counts": 0.0,
            }
            for index in range(1, 6)
        ])
        self.assertEqual(direct_profile_species(groups, 3), {})
        final, _ = final_species_calls("sample", groups, 3)
        if not final.empty:
            self.assertFalse(bool(final.set_index("species")
                                  .loc["Ghost species", "reported"]))

    def test_a_species_with_reads_is_unaffected_by_the_zero_filter(self):
        groups = pd.DataFrame([
            {
                "locus_id": f"L{index}", "consensus_status": "complete_exact_match",
                "compatible_species": "Real species", "combined_est_counts": 40.0,
            }
            for index in range(1, 6)
        ])
        self.assertIn("Real species", direct_profile_species(groups, 3))


class SpeciesCoverageTests(unittest.TestCase):
    """How well assembled the evidence behind a call was.

    These pin the contract the reconciler now leans on: with the confidence
    tiers gone, coverage is one of the few remaining quality signals, so it has
    to describe exactly the loci the call rests on and must never invent a zero.
    """

    A = "Staphylococcus aureus"
    B = "Staphylococcus epidermidis"

    def group(self, locus, species, *, status="complete_exact_match",
              counts=100.0, breadth=0.9, depth=50.0, bases=900, length=1000,
              fragments=200, ordinal=1):
        return {
            "sample_id": "sample", "locus_id": locus,
            "consensus_group_id": f"{locus}_CG{ordinal:04d}",
            "combined_est_counts": counts, "consensus_status": status,
            "compatible_species": ";".join(species),
            "compatible_species_count": len(species),
            "callable_breadth": breadth, "projected_median_depth": depth,
            "callable_bases": bases, "allele_length": length,
            "supporting_fragments": fragments, "composite_from_members": False,
        }

    def test_only_loci_that_name_the_species_alone_are_counted(self):
        """Shared evidence is excluded, or a well-covered relative inflates it."""
        groups = pd.DataFrame([
            self.group("L1", [self.A], breadth=0.9, depth=50.0, bases=900),
            self.group("L2", [self.A], breadth=0.7, depth=30.0, bases=700),
            self.group("L3", [self.A], breadth=1.0, depth=80.0, bases=1000),
            # Shared with a near neighbour: must not contribute.
            self.group("L4", [self.A, self.B], breadth=0.1, depth=1.0, bases=100),
        ])
        coverage = species_coverage(groups)[self.A]
        self.assertAlmostEqual(coverage["mean_callable_breadth"], (0.9 + 0.7 + 1.0) / 3)
        self.assertAlmostEqual(coverage["min_callable_breadth"], 0.7)
        self.assertAlmostEqual(coverage["median_projected_depth"], 50.0)
        self.assertAlmostEqual(coverage["min_projected_depth"], 30.0)
        self.assertEqual(coverage["total_callable_bases"], 2600)
        self.assertEqual(coverage["total_allele_length"], 3000)
        self.assertEqual(coverage["supporting_fragments"], 600)

    def test_the_loci_measured_are_the_loci_the_call_reports(self):
        """Coverage and allele_evidence_unambiguous_loci describe one set.

        The whole point of computing both from _scored_evidence. If they ever
        diverge the report says a call rests on N loci while quoting the
        breadth of some other collection of them.
        """
        groups = pd.DataFrame([
            self.group("L1", [self.A], breadth=0.9),
            self.group("L2", [self.A], breadth=0.8),
            self.group("L3", [self.A], breadth=0.7),
            self.group("L4", [self.A, self.B], breadth=0.6),
            self.group("L5", [self.B], breadth=0.5),
            self.group("L6", [self.B], breadth=0.5),
            self.group("L7", [self.B], breadth=0.5),
        ])
        final = final_species_calls("sample", groups, 3)[0].set_index("species")
        self.assertEqual(final.loc[self.A, "allele_evidence_unambiguous_loci"], 3)
        self.assertAlmostEqual(final.loc[self.A, "mean_callable_breadth"],
                               (0.9 + 0.8 + 0.7) / 3)

    def test_one_value_per_locus_not_per_consensus_group(self):
        """A within-species mixture must not double-count its locus.

        Two haplotypes at one locus are two consensus groups. The
        best-assembled one speaks for the locus; summing both would report more
        callable bases than the locus contains.
        """
        groups = pd.DataFrame([
            self.group("L1", [self.A], breadth=0.9, bases=900, fragments=200, ordinal=1),
            self.group("L1", [self.A], breadth=0.4, bases=400, fragments=20, ordinal=2),
            self.group("L2", [self.A], breadth=0.8, bases=800, fragments=100),
        ])
        coverage = species_coverage(groups)[self.A]
        self.assertAlmostEqual(coverage["mean_callable_breadth"], (0.9 + 0.8) / 2)
        self.assertEqual(coverage["total_callable_bases"], 1700)
        self.assertEqual(coverage["supporting_fragments"], 300)

    def test_an_unmeasured_depth_is_skipped_rather_than_poisoning_the_aggregate(self):
        """projected_median_depth is NaN for some consensus methods."""
        groups = pd.DataFrame([
            self.group("L1", [self.A], depth=float("nan")),
            self.group("L2", [self.A], depth=60.0),
            self.group("L3", [self.A], depth=40.0),
        ])
        coverage = species_coverage(groups)[self.A]
        self.assertAlmostEqual(coverage["median_projected_depth"], 50.0)
        self.assertAlmostEqual(coverage["min_projected_depth"], 40.0)

    def test_all_depths_unmeasured_stays_blank(self):
        groups = pd.DataFrame([
            self.group(f"L{index}", [self.A], depth=float("nan"))
            for index in range(1, 4)
        ])
        coverage = species_coverage(groups)[self.A]
        self.assertTrue(math.isnan(coverage["median_projected_depth"]))
        self.assertTrue(math.isnan(coverage["min_projected_depth"]))
        # Breadth was measured, so it is still reported.
        self.assertFalse(math.isnan(coverage["mean_callable_breadth"]))

    def test_a_zero_abundance_locus_is_not_evidence(self):
        """Same rule that stopped a negative control reporting E. coli."""
        groups = pd.DataFrame([
            self.group("L1", [self.A], counts=0.0, breadth=0.1),
            self.group("L2", [self.A], counts=50.0, breadth=0.9),
        ])
        self.assertAlmostEqual(
            species_coverage(groups)[self.A]["mean_callable_breadth"], 0.9)

    def test_no_unambiguous_evidence_is_blank_not_zero(self):
        """Blank means "not measured"; zero would mean "measured, and none"."""
        groups = pd.DataFrame([
            self.group(f"L{index}", [self.A, self.B]) for index in range(1, 4)
        ])
        self.assertEqual(species_coverage(groups), {})
        final = final_species_calls("sample", groups, 3)[0]
        for column in COVERAGE_COLUMNS:
            self.assertTrue(final[column].isna().all(), column)

    def test_the_coverage_columns_are_in_the_header_even_with_no_calls(self):
        """A sample that called nothing still gets the full header.

        Downstream concatenation skips a file whose header differs, so a column
        that appears only when there is something to report would silently drop
        the sample from the run table.
        """
        empty = pd.DataFrame(columns=[
            "locus_id", "consensus_status", "compatible_species",
            "combined_est_counts", "composite_from_members",
        ])
        final = final_species_calls("sample", empty, 3)[0]
        self.assertTrue(final.empty)
        for column in COVERAGE_COLUMNS:
            self.assertIn(column, list(final.columns))


class BlastDatabaseCacheTests(unittest.TestCase):
    """Where the shared, unlocked BLAST cache is built.

    It defaults to the parent of the output folder, which is right when -o is
    <run_root>/<SAMPLE>. The dual-marker driver writes to
    <run_root>/<SAMPLE>/rmlst, one level deeper, so without an explicit setting
    every sample would silently build its own copy and the up-front makeblastdb
    that keeps concurrent jobs from corrupting it would protect nothing.
    """

    def config(self, output, **kwargs):
        return Config(
            allele_target=Path("ref.fasta"), seqdef=Path("seqdef.txt"),
            reads=(Path("r1.fastq"), Path("r2.fastq")), sample_id="sample",
            output=output, tempdir=output / "work", **kwargs,
        )

    def test_the_default_is_the_parent_of_the_output_folder(self):
        config = self.config(Path("/runs/RUN1/SAMPLE_A"))
        self.assertIsNone(config.blastdb_cache_dir)
        cache_root = config.blastdb_cache_dir or config.output.parent
        self.assertEqual(cache_root, Path("/runs/RUN1"))

    def test_an_explicit_cache_dir_wins_over_the_parent(self):
        config = self.config(Path("/runs/RUN1/SAMPLE_A/rmlst"),
                             blastdb_cache_dir=Path("/runs/RUN1"))
        cache_root = config.blastdb_cache_dir or config.output.parent
        self.assertEqual(cache_root, Path("/runs/RUN1"))
        self.assertNotEqual(cache_root, config.output.parent)

    def test_the_cli_accepts_and_resolves_the_flag(self):
        args = cli_parser().parse_args([
            "--allele-target", "ref.fasta", "--rmlstseqdef", "seqdef.txt",
            "--sample-id", "sample", "-o", "out",
            "--reads", "r1.fastq", "r2.fastq",
            "--blastdb-cache-dir", "cache",
        ])
        self.assertEqual(args.blastdb_cache_dir, Path("cache"))

    def test_omitting_the_flag_leaves_it_unset(self):
        args = cli_parser().parse_args([
            "--allele-target", "ref.fasta", "--rmlstseqdef", "seqdef.txt",
            "--sample-id", "sample", "-o", "out",
            "--reads", "r1.fastq", "r2.fastq",
        ])
        self.assertIsNone(args.blastdb_cache_dir)


class RetainedOutputTests(unittest.TestCase):
    def test_consensus_groups_survives_a_normal_run(self):
        """It is the per-locus detail behind every coverage aggregate.

        Deleting it made the figures in species_calls.csv unauditable without
        re-running the sample.
        """
        self.assertNotIn("consensus_groups.csv", INTERMEDIATE_OUTPUT_NAMES)

    def test_the_other_diagnostics_are_still_removed(self):
        for name in ("kallisto_candidates.csv", "runtime_clusters.csv",
                     "em_fallback.csv", "cluster_consensus"):
            self.assertIn(name, INTERMEDIATE_OUTPUT_NAMES)


if __name__ == "__main__":
    unittest.main()
