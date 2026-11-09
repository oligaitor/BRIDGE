"""The 16S identity/coverage stage (bridge.marker16s.alignment_metrics).

Every binary is faked: the tests check the arithmetic, the file plumbing and
the failure handling, not minimap2.
"""

from __future__ import annotations

import csv
import json
import tempfile
import unittest
from collections import Counter
from pathlib import Path
from unittest.mock import patch

from bridge.rmlst.consensus import MINIMAP2_SHORT_READ_PARAMETERS as RMLST_PARAMETERS
from bridge.marker16s import alignment_metrics as am
from bridge.marker16s.calls import CALLS_COLUMNS, Marker16SCall, write_marker16s_calls
from bridge.marker16s.quant import read_run_info
from bridge.marker16s.reference import load_reference


def header(genus: str, species: str, accession: str) -> str:
    return (f"{accession};K__Bacteria;P__Bacillota;C__Bacilli;O__Lactobacillales;"
            f"F__Streptococcaceae;G__{genus};S__{genus}_{species}")


REFS = [
    (header("Streptococcus", "pyogenes", "NR_1"), "ACGT" * 25),        # 100 bp
    (header("Streptococcus", "pyogenes", "NR_2"), "ACGT" * 30),        # 120 bp
    (header("Streptococcus", "mitis", "NR_3"), "ACGA" * 25),
    (header("Klebsiella", "pneumoniae", "NR_4"), "TTGA" * 25),
    (header("Bacillus", "cereus", "NR_5"), "GGCC" * 25),
]


def call(genus, species, strength):
    return Marker16SCall(
        sample_id="s", taxon=f"{genus} {species}" if species else genus, genus=genus,
        species=species, rank="species" if species else "genus", strength=strength,
        abundance_frac=0.5, genus_frac=0.5, reads_assigned=100.0, raw_reads=400.0,
        copy_number=4.0, within_genus_share=0.9, margin=0.8, n_references=2,
        competing_species=[], reason="test")


class PureHelperTests(unittest.TestCase):

    def test_parameters_match_the_rmlst_arm(self):
        self.assertEqual(am.MINIMAP2_SHORT_READ_PARAMETERS, RMLST_PARAMETERS)

    def test_read_identity_counts_alignment_columns(self):
        self.assertEqual(am.read_identity("50M2I48M", 4), (100, 0.96))
        self.assertEqual(am.read_identity("10S90M", 0), (90, 1.0))
        self.assertEqual(am.read_identity("*", 0), (0, 0.0))
        columns, identity = am.read_identity("30=2X68=", 2)
        self.assertEqual(columns, 100)
        self.assertAlmostEqual(identity, 0.98)

    def test_primary_alignments_are_summarised_per_genus(self):
        genus_of = {"NR_1": "Streptococcus", "NR_2": "Streptococcus", "NR_4": "Klebsiella"}
        sam = [
            "@SQ\tSN:NR_1\tLN:100",
            "r1\t0\tNR_1\t1\t0\t100M\t*\t0\t0\tACGT\tIIII\tNM:i:2",
            "r2\t0\tNR_1\t1\t0\t50M\t*\t0\t0\tACGT\tIIII\tNM:i:0",
            "r3\t0\tNR_2\t1\t0\t100M\t*\t0\t0\tACGT\tIIII\tNM:i:10",
            "r4\t0\tNR_4\t1\t0\t100M\t*\t0\t0\tACGT\tIIII\tNM:i:0",
            "r5\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII",
            "r6\t0\tNR_9\t1\t0\t100M\t*\t0\t0\tACGT\tIIII\tNM:i:0",   # not a called genus
        ]
        stats = am.parse_primary_alignments(sam, genus_of)
        n, identity, refs = stats["Streptococcus"]
        self.assertEqual(n, 3)
        self.assertAlmostEqual(identity, (0.98 + 1.0 + 0.90) / 3)
        self.assertEqual(refs.most_common(1)[0], ("NR_1", 2))
        self.assertEqual(stats["Klebsiella"][0], 1)
        self.assertNotIn("Bacillus", stats)

    def test_depth_parser(self):
        breadth, depth = am.parse_depth(["r\t1\t3", "r\t2\t0", "r\t3\t5", "r\t4\t0"], 4)
        self.assertEqual(breadth, 0.5)
        self.assertEqual(depth, 2.0)
        self.assertEqual(am.parse_depth([], 0), (0.0, 0.0))

    def test_blast_top_hit_by_bitscore_with_query_coverage(self):
        rows = [
            "c\tNR_2\t97.50\t120\t3\t0\t1\t120\t1\t120\t150\t120\t180.0",
            "c\tNR_1\t99.00\t100\t1\t0\t21\t120\t1\t100\t150\t100\t190.5",
        ]
        self.assertEqual(am.parse_blast_top_hit(rows), ("NR_1", 99.0, (120 - 21 + 1) / 150 * 100))
        self.assertIsNone(am.parse_blast_top_hit(["garbage"]))

    def test_genus_subsets_are_written_from_one_pass(self):
        tmp = Path(tempfile.mkdtemp())
        fasta = tmp / "mimt.fna"
        fasta.write_text("".join(f">{h}\n{s[:50]}\n{s[50:]}\n" for h, s in REFS))
        reference = load_reference(fasta)
        combined, per_genus, lengths = am.write_genus_subsets(
            fasta, reference, ["Streptococcus", "Klebsiella"], tmp / "aln")
        self.assertEqual(sorted(per_genus), ["Klebsiella", "Streptococcus"])
        self.assertEqual(len(lengths), 4)
        self.assertEqual(lengths[REFS[1][0]], 120)
        self.assertEqual(combined.read_text().count(">"), 4)
        self.assertEqual(per_genus["Klebsiella"].read_text().count(">"), 1)
        # multi-line sequences are re-joined
        self.assertIn("ACGT" * 25 + "\n", combined.read_text())
        single = am.write_single_record(per_genus["Streptococcus"], REFS[1][0], tmp / "one.fa")
        self.assertEqual(single.read_text().count(">"), 1)

    def test_missing_tool_skips_the_stage_with_a_warning(self):
        with patch.object(am, "executable", side_effect=FileNotFoundError("Executable not found: minimap2")), \
             self.assertLogs("marker16s", level="WARNING") as logs:
            self.assertIsNone(am.resolve_tools())
        self.assertIn("missing_tool", logs.output[0])

    def test_metric_cells_are_blank_where_nothing_was_measured(self):
        blank = am.metric_cells(None)
        self.assertEqual(set(blank), set(am.METRIC_COLUMNS))
        self.assertTrue(all(value == "" for value in blank.values()))
        cells = am.metric_cells(am.GenusAlignmentMetrics(
            "Streptococcus", 3, 250, 0.98765, "NR_1", 0.9712, 41.26, 1490, 99.531, 98.25, "NR_1"))
        self.assertEqual(cells["mean_read_identity"], "0.9877")
        self.assertEqual(cells["consensus_pident"], "99.53")
        self.assertEqual(cells["consensus_qcov"], "98.2")
        self.assertEqual(cells["alignment_status"], "ok")
        no_reads = am.metric_cells(am.GenusAlignmentMetrics("X", 2, status="no_reads"))
        self.assertEqual(no_reads["n_aligned_reads"], "")
        self.assertEqual(no_reads["alignment_status"], "no_reads")


class StageWithFakedToolsTests(unittest.TestCase):
    """Drive alignment_metrics_for_calls end to end with `run` faked."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.fasta = self.tmp / "mimt.fna"
        self.fasta.write_text("".join(f">{h}\n{s}\n" for h, s in REFS))
        self.reference = load_reference(self.fasta)
        self.reads = [self.tmp / "R1.fq", self.tmp / "R2.fq"]
        for path in self.reads:
            path.write_text("@r\nACGT\n+\nIIII\n")
        self.tools = am.AlignmentTools("minimap2", "samtools", "viral_consensus", "blastn", "makeblastdb")
        self.called: list[tuple[str, list[str]]] = []

    def fake_run(self, command, tool, *, stdout_path=None):
        self.called.append((tool, list(map(str, command))))
        command = list(map(str, command))
        if tool == "minimap2":
            Path(stdout_path).write_text("@HD\tVN:1.6\n")
        elif tool == "samtools" and command[1] == "sort":
            Path(command[command.index("-o") + 1]).write_bytes(b"BAM")
        elif tool == "samtools" and command[1] == "index":
            Path(command[2] + ".bai").write_bytes(b"BAI")
        elif tool == "samtools" and command[1] == "view" and stdout_path is not None:
            strep1, strep2 = REFS[0][0], REFS[1][0]
            Path(stdout_path).write_text(
                f"r1\t0\t{strep1}\t1\t0\t100M\t*\t0\t0\tA\tI\tNM:i:1\n"
                f"r2\t0\t{strep1}\t1\t0\t100M\t*\t0\t0\tA\tI\tNM:i:3\n"
                f"r3\t0\t{strep2}\t1\t0\t100M\t*\t0\t0\tA\tI\tNM:i:0\n")
        elif tool == "samtools" and command[1] == "view":            # -h ... -o top.sam
            Path(command[command.index("-o") + 1]).write_text("@HD\n")
        elif tool == "samtools" and command[1] == "depth":
            Path(stdout_path).write_text("".join(f"x\t{i}\t{2 if i <= 90 else 0}\n" for i in range(1, 101)))
        elif tool == "viralconsensus":
            Path(command[command.index("-o") + 1]).write_text(">consensus\n" + "ACGT" * 24 + "\n")
        elif tool == "makeblastdb":
            pass
        elif tool == "blastn":
            out = Path(command[command.index("-out") + 1])
            out.write_text(f"consensus\t{REFS[0][0]}\t99.00\t96\t1\t0\t1\t96\t1\t96\t96\t100\t170.0\n")
        else:
            raise AssertionError(f"unexpected tool {tool}: {command}")

    def test_metrics_for_a_called_genus(self):
        calls = [call("Streptococcus", "pyogenes", "species_strong"),
                 call("Klebsiella", None, "genus_only"),
                 call("Bacillus", None, "silent_weak")]
        with patch.object(am, "run", side_effect=self.fake_run):
            metrics = am.alignment_metrics_for_calls(
                calls, self.reference, self.fasta, self.reads, self.tmp / "alignment",
                tools=self.tools, threads=2)
        self.assertEqual(sorted(metrics), ["Klebsiella", "Streptococcus"])   # silent genus not aligned
        strep = metrics["Streptococcus"]
        self.assertEqual(strep.status, "ok")
        self.assertEqual(strep.n_references, 3)
        self.assertEqual(strep.n_aligned_reads, 3)
        self.assertAlmostEqual(strep.mean_read_identity, (0.99 + 0.97 + 1.0) / 3)
        self.assertEqual(strep.top_reference, REFS[0][0])
        self.assertAlmostEqual(strep.ref_breadth, 0.9)
        self.assertAlmostEqual(strep.mean_depth, 1.8)
        self.assertEqual(strep.consensus_length, 96)
        self.assertEqual(strep.consensus_pident, 99.0)
        self.assertEqual(strep.consensus_qcov, 100.0)
        self.assertEqual(strep.consensus_top_hit, REFS[0][0])
        self.assertEqual(metrics["Klebsiella"].status, "no_reads")
        self.assertEqual(metrics["Klebsiella"].n_references, 1)

        tools_used = [tool for tool, _ in self.called]
        self.assertEqual(tools_used[:3], ["minimap2", "samtools", "samtools"])   # align, sort, index
        minimap = self.called[0][1]
        self.assertEqual(minimap[1:1 + len(RMLST_PARAMETERS)], list(RMLST_PARAMETERS))
        self.assertIn("-F", self.called[3][1])                                   # primary only
        self.assertIn("0x904", self.called[3][1])
        self.assertTrue(any(tool == "blastn" and "-max_target_seqs" in cmd for tool, cmd in self.called))
        self.assertTrue((self.tmp / "alignment" / "Streptococcus" / "consensus.fasta").is_file())

        # ...and the calls file carries them, blank for the silent genus
        path = write_marker16s_calls(calls, self.tmp / "marker16s_calls.csv", metrics)
        with open(path, newline="") as handle:
            rows = {row["genus"]: row for row in csv.DictReader(handle)}
        self.assertEqual(list(rows["Streptococcus"]), CALLS_COLUMNS)
        self.assertEqual(rows["Streptococcus"]["consensus_pident"], "99.00")
        self.assertEqual(rows["Streptococcus"]["ref_breadth"], "0.9000")
        self.assertEqual(rows["Streptococcus"]["n_aligned_reads"], "3")
        self.assertEqual(rows["Klebsiella"]["alignment_status"], "no_reads")
        self.assertEqual(rows["Klebsiella"]["consensus_pident"], "")
        self.assertEqual(rows["Bacillus"]["alignment_status"], "")

    def test_an_alignment_failure_marks_every_genus_and_raises_nothing(self):
        import subprocess

        def failing_run(command, tool, *, stdout_path=None):
            raise subprocess.CalledProcessError(1, command)

        calls = [call("Streptococcus", "pyogenes", "species_strong")]
        with patch.object(am, "run", side_effect=failing_run):
            metrics = am.alignment_metrics_for_calls(
                calls, self.reference, self.fasta, self.reads, self.tmp / "alignment",
                tools=self.tools)
        self.assertEqual(metrics["Streptococcus"].status, "failed:align")
        self.assertEqual(am.metric_cells(metrics["Streptococcus"])["alignment_status"], "failed:align")

    def test_nothing_called_means_nothing_aligned(self):
        with patch.object(am, "run", side_effect=self.fake_run):
            metrics = am.alignment_metrics_for_calls(
                [call("Bacillus", None, "silent_weak")], self.reference, self.fasta,
                self.reads, self.tmp / "alignment", tools=self.tools)
        self.assertEqual(metrics, {})
        self.assertEqual(self.called, [])


class RunInfoTests(unittest.TestCase):

    def test_read_counts_come_from_the_quant_folder(self):
        tmp = Path(tempfile.mkdtemp())
        (tmp / "quant").mkdir()
        (tmp / "quant" / "run_info.json").write_text(json.dumps({
            "n_targets": 29090, "n_processed": 22539, "n_pseudoaligned": 2337,
            "n_unique": 282, "p_pseudoaligned": 10.4, "p_unique": 1.3,
            "kallisto_version": "0.52.0"}))
        self.assertEqual(read_run_info(tmp), {
            "n_processed": 22539, "n_pseudoaligned": 2337, "n_unique": 282,
            "p_pseudoaligned": 10.4, "p_unique": 1.3})
        self.assertEqual(read_run_info(tmp / "missing"), {})


class SingleReferenceSamTests(unittest.TestCase):
    """ViralConsensus needs a SAM whose header declares exactly one reference."""

    TARGET = REFS[0][0]
    OTHER = REFS[2][0]

    def write(self, tmp: Path, text: str) -> Path:
        path = tmp / "all.sam"
        path.write_text(text, encoding="utf-8")
        return path

    def test_keeps_only_the_target_sq_line(self):
        sam = (f"@HD\tVN:1.6\tSO:coordinate\n"
               f"@SQ\tSN:{self.TARGET}\tLN:100\n"
               f"@SQ\tSN:{self.OTHER}\tLN:100\n"
               f"@SQ\tSN:{self.TARGET}X\tLN:100\n"      # prefix match must not survive
               f"@PG\tID:minimap2\tPN:minimap2\n"
               f"r1\t0\t{self.TARGET}\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n")
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            out = tmp / "one.sam"
            n = am.single_reference_sam(self.write(tmp, sam), self.TARGET, out)
            lines = out.read_text(encoding="utf-8").splitlines()
        sq = [line for line in lines if line.startswith("@SQ")]
        self.assertEqual(sq, [f"@SQ\tSN:{self.TARGET}\tLN:100"])
        self.assertIn("@HD\tVN:1.6\tSO:coordinate", lines)
        self.assertIn("@PG\tID:minimap2\tPN:minimap2", lines)
        self.assertEqual(n, 1)

    def test_mate_on_another_reference_is_cleared(self):
        sam = (f"@SQ\tSN:{self.TARGET}\tLN:100\n@SQ\tSN:{self.OTHER}\tLN:100\n"
               f"r1\t65\t{self.TARGET}\t1\t60\t4M\t{self.OTHER}\t40\t0\tACGT\tIIII\n"
               f"r2\t99\t{self.TARGET}\t5\t60\t4M\t=\t20\t19\tACGT\tIIII\n")
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            out = tmp / "one.sam"
            am.single_reference_sam(self.write(tmp, sam), self.TARGET, out)
            records = [line.split("\t") for line in out.read_text(encoding="utf-8").splitlines()
                       if not line.startswith("@")]
        self.assertEqual(records[0][6:8], ["*", "0"])
        self.assertEqual(records[1][6:8], ["=", "20"])      # same-reference mate untouched
        self.assertEqual(records[0][9], "ACGT")               # sequence untouched


if __name__ == "__main__":
    unittest.main()
