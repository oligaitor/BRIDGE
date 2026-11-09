"""Intermediate files (bridge.clean, and the tidy-up in bridge.sample).

The point of the module is that it removes only what nothing downstream
reads, so most of these tests pin the files that must SURVIVE: a tree that has
been cleaned must still finalise.
"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from bridge import clean, pbs
from bridge import sample as sample_mod
from bridge.cli import main as cli_main
from bridge.reads import ReadPair
from bridge.sample import ArmFailed, run_sample
from bridge.settings import ARMS, load_settings

# What finalise / amr-report read, plus the small files worth keeping.
KEPT = (
    "rmlst/species_calls.csv",
    "rmlst/rmlst_manifest.json",
    "16s/marker16s_calls.csv",
    "16s/marker16s_manifest.json",
    "16s/quant/run_info.json",
    "16s/quant/abundance.tsv",
    "16s/alignment/Escherichia/consensus.fasta",
    "16s/alignment/Escherichia/blast.tsv",
    "16s/alignment/Escherichia/depth.tsv",
    "16s/alignment/Escherichia/top_reference.fasta",
    "amr/amr_hits.csv",
    "amr/amr_manifest.json",
    "amr/kma.res",
    "amr/kma.mapstat",
    "amr/kma.fsa",
    "bacterial_ID/identification_summary.csv",
)
REMOVED = (
    "16s/alignment/alignment.bam",
    "16s/alignment/alignment.bam.bai",
    "16s/alignment/primary.sam",
    "16s/alignment/references.fasta",
    "16s/alignment/Escherichia/references.fasta",
    "16s/alignment/Escherichia/top_reference.sam",
    "16s/alignment/Escherichia/blastdb/references.nsq",
    "amr/kma.aln",
    "amr/kma.mat.gz",
)


def populate(sample_dir: Path, names=KEPT + REMOVED, size: int = 100) -> None:
    for name in names:
        path = sample_dir / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"x" * size)


class CleanRunTests(unittest.TestCase):

    def setUp(self):
        self.out = Path(tempfile.mkdtemp()) / "results"
        for name in ("S1", "NEG_1"):
            populate(self.out / name)
        cache = self.out / clean.BLAST_CACHE
        cache.mkdir(parents=True)
        (cache / "reference_abc.nsq").write_bytes(b"x" * 1000)
        (self.out / "identification_summary_all.csv").write_text("a\n")
        (self.out / "oe").mkdir()
        (self.out / "oe" / "S1_RUN.o123").write_text("log\n")

    def test_a_dry_run_measures_and_deletes_nothing(self):
        report = clean.clean_run(self.out)
        self.assertEqual(report.samples, 2)
        self.assertEqual(sum(report.sizes.values()), 2 * len(REMOVED) * 100)
        for sample in ("S1", "NEG_1"):
            for name in KEPT + REMOVED:
                self.assertTrue((self.out / sample / name).is_file(), name)
        text = clean.describe(report, self.out)
        self.assertIn("Would remove", text)
        self.assertIn("Nothing was deleted", text)
        self.assertIn("16s/alignment/<genus>/blastdb", text)

    def test_apply_removes_the_intermediates_and_keeps_what_finalise_reads(self):
        clean.clean_run(self.out, apply=True)
        for sample in ("S1", "NEG_1"):
            for name in KEPT:
                self.assertTrue((self.out / sample / name).is_file(), f"lost {name}")
            for name in REMOVED:
                self.assertFalse((self.out / sample / name).exists(), f"kept {name}")
            self.assertFalse((self.out / sample / "16s/alignment/Escherichia/blastdb").exists())
        self.assertTrue((self.out / "identification_summary_all.csv").is_file())
        self.assertTrue((self.out / "oe" / "S1_RUN.o123").is_file())

    def test_a_second_pass_finds_nothing(self):
        clean.clean_run(self.out, apply=True)
        report = clean.clean_run(self.out)
        self.assertEqual(report.samples, 0)
        self.assertEqual(report.sizes, {})

    def test_the_blast_cache_goes_only_when_asked(self):
        clean.clean_run(self.out, apply=True)
        self.assertTrue((self.out / clean.BLAST_CACHE).is_dir())
        clean.clean_run(self.out, apply=True, blast_cache=True)
        self.assertFalse((self.out / clean.BLAST_CACHE).exists())

    def test_a_sample_whose_job_is_still_running_is_left_alone(self):
        jobids = self.out / "oe" / ".jobids"
        jobids.mkdir(parents=True)
        (jobids / "S1").write_text("123.pbs\n")
        (jobids / "NEG_1").write_text("124.pbs\n")
        with patch.object(pbs, "active_jobs", return_value={"123.pbs"}):
            report = clean.clean_run(self.out, apply=True, blast_cache=True)
        self.assertEqual(report.skipped, ["S1"])
        self.assertTrue((self.out / "S1" / "amr" / "kma.mat.gz").is_file())
        self.assertFalse((self.out / "NEG_1" / "amr" / "kma.mat.gz").exists())
        # the shared database is still in use by that job
        self.assertTrue((self.out / clean.BLAST_CACHE).is_dir())

    def test_the_command_reports_and_only_deletes_with_apply(self):
        with patch("sys.stdout"):
            self.assertEqual(cli_main(["clean", "--out", str(self.out)]), 0)
        self.assertTrue((self.out / "S1" / "16s/alignment/alignment.bam").is_file())
        with patch("sys.stdout"):
            self.assertEqual(cli_main(["clean", "--out", str(self.out), "--apply"]), 0)
        self.assertFalse((self.out / "S1" / "16s/alignment/alignment.bam").exists())

    def test_a_missing_tree_is_an_error_not_a_crash(self):
        with patch("sys.stderr"):
            self.assertEqual(cli_main(["clean", "--out", str(self.out / "nope")]), 1)


class SampleTidyTests(unittest.TestCase):
    """New runs tidy each arm as it finishes, as bridge.rmlst does for rMLST."""

    def setUp(self):
        self.settings = load_settings()
        self.tmp = Path(tempfile.mkdtemp())
        r1, r2 = self.tmp / "S_R1.fq", self.tmp / "S_R2.fq"
        r1.write_text("@r\nACGT\n+\nIIII\n")
        r2.write_text("@r\nACGT\n+\nIIII\n")
        self.pair = ReadPair("S", r1, r2)
        self.out = self.tmp / "out"
        self.sample_dir = self.out / "S"

    def _arms(self, *, amr_error=None):
        """Patch the three arms; the 16S and AMR ones leave a full set of files."""
        def fake_16s(argv):
            populate(self.sample_dir, [n for n in KEPT + REMOVED if n.startswith("16s/")])
            return 0

        def fake_amr(*args, **kwargs):
            populate(self.sample_dir, [n for n in KEPT + REMOVED if n.startswith("amr/")])
            if amr_error:
                raise amr_error
            return self.sample_dir / "amr" / "amr_hits.csv"

        return (patch.object(sample_mod, "run_pipeline"),
                patch.object(sample_mod, "marker16s_main", side_effect=fake_16s),
                patch("bridge.amr.kma.run_arm", side_effect=fake_amr))

    def test_intermediates_are_removed_by_default(self):
        a, b, c = self._arms()
        with a, b, c, patch("sys.stdout"):
            run_sample(self.settings, self.pair, self.out, arms=ARMS)
        for name in REMOVED:
            self.assertFalse((self.sample_dir / name).exists(), name)
        for name in KEPT:
            if name.startswith(("16s/", "amr/")):
                self.assertTrue((self.sample_dir / name).is_file(), name)

    def test_keep_temp_keeps_them(self):
        a, b, c = self._arms()
        with a, b, c, patch("sys.stdout"):
            run_sample(self.settings, self.pair, self.out, arms=ARMS, keep_temp=True)
        for name in REMOVED:
            self.assertTrue((self.sample_dir / name).exists(), name)

    def test_a_failed_arm_keeps_its_files_for_diagnosis(self):
        a, b, c = self._arms(amr_error=RuntimeError("kma died"))
        with a, b, c, patch("sys.stdout"):
            with self.assertRaises(ArmFailed):
                run_sample(self.settings, self.pair, self.out, arms=ARMS)
        self.assertTrue((self.sample_dir / "amr" / "kma.mat.gz").is_file())
        self.assertFalse((self.sample_dir / "16s/alignment/alignment.bam").exists())


class PbsKeepTempTests(unittest.TestCase):

    def test_keep_temp_reaches_the_job_script_only_when_asked(self):
        settings = load_settings()
        tmp = Path(tempfile.mkdtemp())
        pair = ReadPair("S1", tmp / "S1_R1.fq", tmp / "S1_R2.fq")
        tier = settings.hpc.tiers[0]
        plain = pbs.render_sample_script(settings, pair, tmp, "RUN1", tier,
                                         arms=ARMS, size_bytes=1)
        kept = pbs.render_sample_script(settings, pair, tmp, "RUN1", tier,
                                        arms=ARMS, size_bytes=1, keep_temp=True)
        self.assertNotIn("--keep-temp", plain)
        self.assertIn("--keep-temp", kept)
        # still one continued shell command
        self.assertIn('bridge_S1" \\\n    --keep-temp', kept)


if __name__ == "__main__":
    unittest.main()
