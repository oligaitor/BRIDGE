"""One outcome per sample (bridge.outcome) - the rule, finalise step 5, both CLI modes."""

from __future__ import annotations

import csv
import io
import json
import shutil
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

from bridge.finalise import describe, finalise_run
from bridge.outcome import (
    OUTCOMES_NAME, OutcomeSettings, build_outcomes, outcome_settings, read_sample_list,
    record_sample_list,
)
from bridge.reads import ReadPair
from bridge.settings import load_settings

import test_finalise as fin
from test_finalise import rmlst_row, s16_row

SETTINGS = OutcomeSettings()


class TreeCase(unittest.TestCase):
    """A synthetic results tree, as in test_finalise (its tests not re-run here)."""

    setUp = fin.FinaliseTests.setUp
    sample = fin.FinaliseTests.sample


def id_row(sample, organism, status, rank="species"):
    return {"sample": sample, "organism": organism, "status": status,
            "reported_rank": rank, "basis": "rmlst"}


def read(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return {row["sample"]: row for row in csv.DictReader(handle)}


class RuleTests(unittest.TestCase):
    """build_outcomes on its own: the order the outcomes are tried in."""

    def outcome(self, sample="S", rows=(), reads=1_000_000, **kwargs):
        rows_out = build_outcomes("RUN1", [sample], id_rows=list(rows), amr_rows=[],
                                  depth={sample: reads}, settings=SETTINGS, **kwargs)
        return rows_out[0]

    def test_a_reported_organism_is_detected_and_named(self):
        row = self.outcome(rows=[id_row("S", "Proteus mirabilis", "report"),
                                 id_row("S", "Streptococcus", "report", rank="genus"),
                                 id_row("S", "Escherichia", "review", rank="genus")])
        self.assertEqual(row["outcome"], "detected")
        self.assertEqual(row["statement"],
                         "Bacteria detected: Proteus mirabilis, Streptococcus sp.")
        self.assertEqual((row["n_report"], row["n_review"]), ("2", "1"))

    def test_only_review_override_calls_are_held_for_review(self):
        row = self.outcome(rows=[id_row("S", "Massilia", "review_override", "genus"),
                                 id_row("S", "Escherichia", "review", "genus")])
        self.assertEqual(row["outcome"], "held_for_review")
        self.assertEqual(row["statement"], "Result pending scientist review")

    def test_only_background_calls_are_a_negative(self):
        row = self.outcome(rows=[id_row("S", "Escherichia", "review", "genus")])
        self.assertEqual(row["outcome"], "not_detected")
        self.assertEqual(row["statement"], "No bacteria present above reporting thresholds")
        self.assertIn("held_calls_not_reported", row["flags"])

    def test_no_call_is_a_negative(self):
        row = self.outcome(rows=[{"sample": "S", "organism": "", "status": "no_call"}])
        self.assertEqual(row["outcome"], "not_detected")
        self.assertEqual(row["flags"], "")

    def test_a_shallow_negative_is_insufficient_but_a_shallow_positive_is_not(self):
        shallow = self.outcome(rows=[{"sample": "S", "status": "no_call"}], reads=9_999)
        self.assertEqual(shallow["outcome"], "insufficient_reads")
        self.assertEqual(self.outcome(rows=[{"sample": "S", "status": "no_call"}],
                                      reads=10_000)["outcome"], "not_detected")
        positive = self.outcome(rows=[id_row("S", "Proteus mirabilis", "report")], reads=500)
        self.assertEqual(positive["outcome"], "detected")
        self.assertIn("low_depth", positive["flags"])

    def test_zero_reads_beats_a_failed_arm(self):
        row = self.outcome(reads=0, failures={"S": "no rmlst/species_calls.csv"})
        self.assertEqual(row["outcome"], "insufficient_reads")

    def test_a_failure_ignores_rows_left_by_an_earlier_finalise(self):
        row = self.outcome(rows=[id_row("S", "Proteus mirabilis", "report")],
                           failures={"S": "exit 1"})
        self.assertEqual((row["outcome"], row["reason"]), ("analysis_failed", "exit 1"))

    def test_no_rows_and_no_depth_is_a_failure_that_says_so(self):
        row = self.outcome(reads=None)
        self.assertEqual(row["outcome"], "analysis_failed")
        self.assertEqual(row["reason"], "no identification rows; depth unknown")
        self.assertIn("depth_unknown", row["flags"])

    def test_controls_and_samples_without_a_folder(self):
        self.assertEqual(self.outcome("NEGMATRIX_1", control_prefixes=("NEG",))["outcome"],
                         "negative_control")
        self.assertEqual(self.outcome("Water", controls={"Water"})["outcome"],
                         "negative_control")
        missing = self.outcome("Gone", present={"Other"})
        self.assertEqual((missing["outcome"], missing["statement"]), ("not_run", "Not analysed"))

    def test_samples_off_the_list_are_flagged(self):
        row = self.outcome(rows=[{"sample": "S", "status": "no_call"}], listed={"T"})
        self.assertIn("not_on_sample_list", row["flags"])

    def test_amr_outcome(self):
        rows = build_outcomes("RUN1", ["A", "B", "C", "D"], id_rows=[], depth={},
                              settings=SETTINGS,
                              amr_rows=[{"sample": "A", "status": "report"},
                                        {"sample": "A", "status": "review"},
                                        {"sample": "B", "status": "review"},
                                        {"sample": "C", "status": "no_call"}])
        self.assertEqual([row["amr_outcome"] for row in rows],
                         ["genes_reported", "genes_held_for_review",
                          "no_gene_above_floor", "not_run"])

    def test_statements_come_from_the_config(self):
        settings = outcome_settings(load_settings().cfg)
        self.assertEqual(settings.min_read_pairs, 10_000)
        self.assertEqual(settings.statements["not_detected"],
                         "No bacteria present above reporting thresholds")
        edited = OutcomeSettings(statements={**SETTINGS.statements,
                                             "detected": "Found {organisms} {oops}"})
        row = build_outcomes("RUN1", ["S"], id_rows=[id_row("S", "Proteus mirabilis", "report")],
                             amr_rows=[], depth={}, settings=edited)[0]
        # a typo in the wording must not stop a run: the template is used as written
        self.assertEqual(row["statement"], "Found {organisms} {oops}")


class SampleListTests(unittest.TestCase):
    def test_read_sample_list_takes_the_r1_parent_and_never_touches_the_reads(self):
        tmp = Path(tempfile.mkdtemp())
        path = tmp / "sample.list"
        path.write_text("# a comment\n/scratch/runs/R/reads/P1/P1_S1_R1_001.fastq.gz\r\n\n"
                        "/scratch/runs/R/reads/P2/P2_S2_R1_001.fastq.gz\n"
                        "/scratch/runs/R/reads/P1/P1_S1_R1_001.fastq.gz\n")
        self.assertEqual(read_sample_list(path), ["P1", "P2"])

    def test_record_sample_list_only_ever_adds(self):
        tmp = Path(tempfile.mkdtemp())
        first = [ReadPair("P1", Path("/r/P1/P1_R1.fq.gz"), Path("/r/P1/P1_R2.fq.gz")),
                 ReadPair("P2", Path("/r/P2/P2_R1.fq.gz"), Path("/r/P2/P2_R2.fq.gz"))]
        path, added = record_sample_list(tmp, first)
        self.assertEqual(added, 2)
        # a re-run of a subset plus one new sample
        _, added = record_sample_list(tmp, [first[1], ReadPair(
            "P3", Path("/r/P3/P3_R1.fq.gz"), Path("/r/P3/P3_R2.fq.gz"))])
        self.assertEqual(added, 1)
        self.assertEqual(read_sample_list(path), ["P1", "P2", "P3"])

    def test_a_colleagues_sample_list_is_read_but_never_written(self):
        from bridge.outcome import sample_list_path
        tmp = Path(tempfile.mkdtemp())
        theirs = tmp / "sample.list"
        theirs.write_text("/r/X1/X1_R1.fq.gz\n")
        self.assertEqual(sample_list_path(tmp), theirs)          # read when ours is absent
        path, _ = record_sample_list(tmp, [ReadPair(
            "P1", Path("/r/P1/P1_R1.fq.gz"), Path("/r/P1/P1_R2.fq.gz"))])
        self.assertEqual(path.name, "sample_bacteria.list")
        self.assertEqual(theirs.read_text(), "/r/X1/X1_R1.fq.gz\n")  # untouched
        self.assertEqual(sample_list_path(tmp), path)            # ours wins once it exists


class FinaliseOutcomeTests(TreeCase):
    """Step 5 end to end, on the same kind of synthetic tree as test_finalise."""

    def build_run(self):
        empty_rmlst = lambda name: [rmlst_row(name, "Escherichia coli", reported=False, est=20.0)]
        self.sample("A", [rmlst_row("A", "Proteus mirabilis")],
                    [s16_row("A", "Proteus mirabilis", "Proteus", "mirabilis", "species_strong")])
        self.sample("N", empty_rmlst("N"),
                    [s16_row("N", "Paracoccus yeei", "Paracoccus", "yeei", "species_strong"),
                     s16_row("N", "Moraxella tetraodonis", "Moraxella", "tetraodonis",
                             "species_strong"),
                     s16_row("N", "Sphingomonas", "Sphingomonas", "", "genus_only")])
        # one fallback call, and it is the blank's Escherichia: background only
        self.sample("M", empty_rmlst("M"),
                    [s16_row("M", "Escherichia", "Escherichia", "", "genus_only")])
        self.sample("Q", [], [])                                  # arms ran, called nothing
        self.sample("L", [], [], total=5_000)                     # too shallow to be negative
        self.sample("Neg_1", empty_rmlst("Neg_1"),
                    [s16_row("Neg_1", "Escherichia", "Escherichia", "", "genus_only", raw=70.0)],
                    total=22_539)
        # zero reads: kallisto failed, so no rMLST output; the AMR manifest says 0
        (self.root / "Z" / "amr").mkdir(parents=True)
        (self.root / "Z" / "amr" / "amr_manifest.json").write_text(
            json.dumps({"total_read_pairs": 0}))
        # rMLST crashed on a deep sample: 16S output only
        (self.root / "X" / "16s" / "quant").mkdir(parents=True)
        (self.root / "X" / "16s" / "marker16s_calls.csv").write_text(
            "sample_id,taxon,genus,species,rank,strength,abundance_frac,genus_frac,"
            "reads_assigned,raw_reads,copy_number,within_genus_share,margin,"
            "n_references,competing_species,reason\n")
        (self.root / "X" / "16s" / "quant" / "run_info.json").write_text(
            json.dumps({"n_processed": 400_000}))
        (self.root / "sample.list").write_text("".join(
            f"/scratch/r/reads/{name}/{name}_S1_L001_R1_001.fastq.gz\n"
            for name in ("A", "N", "M", "Q", "L", "Z", "X", "Neg_1", "Missing")))

    def test_every_sample_on_the_list_gets_exactly_one_outcome(self):
        self.build_run()
        result = finalise_run(self.settings, self.root, "RUN1")
        self.assertIsNone(result.outcome_error)
        rows = read(self.root / OUTCOMES_NAME)
        self.assertEqual({name: row["outcome"] for name, row in rows.items()}, {
            "A": "detected", "N": "held_for_review", "M": "not_detected",
            "Q": "not_detected", "L": "insufficient_reads", "Z": "insufficient_reads",
            "X": "analysis_failed", "Neg_1": "negative_control", "Missing": "not_run"})
        self.assertEqual(rows["A"]["statement"], "Bacteria detected: Proteus mirabilis")
        self.assertIn("held_calls_not_reported", rows["M"]["flags"])
        self.assertEqual(rows["Z"]["total_read_pairs"], "0")
        self.assertIn("arm incomplete", rows["X"]["reason"])
        self.assertEqual(result.outcome.counts["insufficient_reads"], 2)
        clinical = read(self.root / "clinical_summary_all.csv")
        self.assertEqual(clinical["A"]["organism"], "Proteus mirabilis")
        self.assertNotIn("Missing", clinical)
        self.assertNotIn("Neg_1", clinical)
        text = describe(result, self.root)
        self.assertIn(OUTCOMES_NAME, text)
        self.assertIn("held_for_review", text)

    def test_without_a_sample_list_the_folders_are_the_run(self):
        self.sample("A", [rmlst_row("A", "Proteus mirabilis")],
                    [s16_row("A", "Proteus mirabilis", "Proteus", "mirabilis", "species_strong")])
        result = finalise_run(self.settings, self.root, "RUN1")
        self.assertIsNone(result.outcome.sample_list)
        self.assertEqual(list(read(self.root / OUTCOMES_NAME)), ["A"])
        self.assertIn("no sample list", describe(result, self.root))

    def test_backup_folders_in_the_tree_are_not_samples(self):
        self.build_run()
        for name in ("pre_pooled_20261012", "pooled_261012"):
            (self.root / name).mkdir()
            (self.root / name / "identification_summary_all.csv").write_text("sample\n")
        finalise_run(self.settings, self.root, "RUN1")
        rows = read(self.root / OUTCOMES_NAME)
        self.assertNotIn("pre_pooled_20261012", rows)
        self.assertNotIn("pooled_261012", rows)
        self.assertEqual(rows["X"]["outcome"], "analysis_failed")   # a real sample still is
        # ... and without a list, a folder with arm output is still a sample
        (self.root / "sample.list").unlink()
        finalise_run(self.settings, self.root, "RUN1")
        rows = read(self.root / OUTCOMES_NAME)
        self.assertNotIn("pooled_261012", rows)
        self.assertEqual(rows["Z"]["outcome"], "insufficient_reads")

    def test_a_stale_summary_does_not_rescue_a_sample_that_now_fails(self):
        self.sample("A", [rmlst_row("A", "Proteus mirabilis")],
                    [s16_row("A", "Proteus mirabilis", "Proteus", "mirabilis", "species_strong")])
        finalise_run(self.settings, self.root, "RUN1")
        (self.root / "A" / "rmlst" / "species_calls.csv").unlink()   # arm lost on a re-run
        finalise_run(self.settings, self.root, "RUN1")
        self.assertEqual(read(self.root / OUTCOMES_NAME)["A"]["outcome"], "analysis_failed")

    def test_an_outcome_failure_leaves_the_other_tables(self):
        from unittest.mock import patch
        self.sample("A", [rmlst_row("A", "Proteus mirabilis")],
                    [s16_row("A", "Proteus mirabilis", "Proteus", "mirabilis", "species_strong")])
        with patch("bridge.outcome.outcomes_from_tree", side_effect=RuntimeError("boom")):
            result = finalise_run(self.settings, self.root, "RUN1")
        self.assertIn("boom", result.outcome_error)
        self.assertTrue((self.root / "identification_summary_all.csv").is_file())
        from bridge.cli import _finalise_failed
        self.assertTrue(_finalise_failed(result))


class CliTests(TreeCase):
    def run_cli(self, *argv):
        from bridge.cli import main
        buffer = io.StringIO()
        with redirect_stdout(buffer):
            code = main(list(argv))
        return code, buffer.getvalue()

    def test_outcomes_out_rebuilds_from_the_tree(self):
        FinaliseOutcomeTests.build_run(self)
        finalise_run(self.settings, self.root, "RUN1")
        (self.root / OUTCOMES_NAME).unlink()
        code, text = self.run_cli("outcomes", "--out", str(self.root), "--run-id", "RUN1")
        self.assertEqual(code, 0, text)
        rows = read(self.root / OUTCOMES_NAME)
        self.assertEqual(rows["N"]["outcome"], "held_for_review")
        self.assertEqual(rows["X"]["outcome"], "analysis_failed")
        self.assertEqual(rows["Missing"]["outcome"], "not_run")

    def test_outcomes_tables_works_from_the_run_tables_alone(self):
        FinaliseOutcomeTests.build_run(self)
        finalise_run(self.settings, self.root, "RUN1")
        tables = Path(tempfile.mkdtemp()) / "RUN1_tables"
        tables.mkdir()
        for name in ("identification_summary_all.csv", "sample.list"):
            shutil.copy(self.root / name, tables / name)
        dest = Path(tempfile.mkdtemp())
        depth = dest / "depth.csv"
        depth.write_text("sample,total_read_pairs\nZ,0\nL,5000\n")
        code, text = self.run_cli("outcomes", "--tables", str(tables), "--dest", str(dest),
                                  "--depth-csv", str(depth))
        self.assertEqual(code, 0, text)
        self.assertFalse((tables / OUTCOMES_NAME).exists())   # --dest respected
        rows = read(dest / OUTCOMES_NAME)
        self.assertEqual(rows["A"]["run_id"], "RUN1")           # from the table, not the folder
        self.assertEqual(rows["Z"]["outcome"], "insufficient_reads")
        # no folders to look at: never "not_run", and without depth a gap is a failure
        self.assertEqual(rows["Missing"]["outcome"], "analysis_failed")
        self.assertEqual(rows["Missing"]["reason"], "no identification rows; depth unknown")
        self.assertEqual(rows["N"]["outcome"], "held_for_review")

    def test_out_needs_a_run_id(self):
        from contextlib import redirect_stderr
        with redirect_stderr(io.StringIO()):
            code, _ = self.run_cli("outcomes", "--out", str(self.root))
        self.assertEqual(code, 1)


if __name__ == "__main__":
    unittest.main()
