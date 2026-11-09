"""Splitting held_for_review (bridge.held_split), and its place in the outcome rule."""

from __future__ import annotations

import csv
import io
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

from bridge.held_split import (
    BACKGROUND, HELD, UNCONFIRMED, RunContext, SplitSettings, build_recurrence,
    recurrent_elsewhere, split_held,
)
from bridge.outcome import OutcomeSettings, build_outcomes, outcome_settings
from bridge.settings import load_settings


def call(sample, genus, reads, identity="0.995", status="review_override"):
    return {"sample": sample, "organism": genus, "genus": genus, "status": status,
            "basis": "16s_fallback", "16s_raw_reads": str(reads),
            "16s_read_identity": identity}


def recurrence(*rows):
    return tuple({"run_id": run, "genus": genus, "prevalence": str(p)} for run, genus, p in rows)


def run_of(*samples):
    """Other samples on the run: (sample, [genera])."""
    return [call(name, genus, 10, status="review") for name, genera in samples for genus in genera]


class SplitTests(unittest.TestCase):
    def split(self, calls, others=(), settings=SplitSettings(), run_id="RUN1"):
        context = RunContext.build(run_id, list(calls) + run_of(*others), settings)
        return split_held("S", calls, context, settings)

    def test_one_dominant_high_identity_call_is_unconfirmed_dna(self):
        outcome, genus, reason = self.split([call("S", "Stenotrophomonas", 1000),
                                             call("S", "Moraxella", 50)])
        self.assertEqual((outcome, genus), (UNCONFIRMED, "Stenotrophomonas"))
        self.assertIn("95%", reason)

    def test_a_dominant_call_needs_identity_reads_and_share(self):
        self.assertEqual(self.split([call("S", "Steno", 1000, "0.98"), call("S", "X", 10)])[0], HELD)
        self.assertEqual(self.split([call("S", "Steno", 40), call("S", "X", 1)])[0], HELD)
        self.assertEqual(self.split([call("S", "Steno", 700), call("S", "X", 300)])[0], HELD)
        self.assertEqual(self.split([call("S", "Steno", 1000, ""), call("S", "X", 10)])[0], HELD)

    def test_a_dominant_genus_common_on_its_run_is_background(self):
        # A genus in most of the run's samples is a batch contaminant, however dominant
        others = [(f"O{i}", ["Blastomonas"]) for i in range(7)] + [(f"P{i}", ["Other"]) for i in range(3)]
        outcome, _, reason = self.split([call("S", "Blastomonas", 16000), call("S", "Aquabacterium", 900)],
                                        others=others + [("Q", ["Aquabacterium"])] * 0)
        self.assertEqual(outcome, HELD)       # Aquabacterium is not common here
        outcome, _, reason = self.split([call("S", "Blastomonas", 16000), call("S", "Aquabacterium", 900)],
                                        others=[(f"O{i}", ["Blastomonas", "Aquabacterium"]) for i in range(7)]
                                        + [(f"P{i}", ["Other"]) for i in range(3)])
        self.assertEqual(outcome, BACKGROUND)
        self.assertIn("Blastomonas: common on run", reason)

    def test_recurrence_on_other_runs_or_a_weak_call_makes_background(self):
        settings = SplitSettings(recurrence=recurrence(
            ("R2", "Moraxella", 0.5), ("R3", "Moraxella", 0.2), ("R4", "Moraxella", 0.1),
            ("R2", "Paracoccus", 0.9), ("R3", "Paracoccus", 0.9), ("R4", "Paracoccus", 0.12)))
        outcome, _, reason = self.split([call("S", "Moraxella", 69, "0.967"),
                                         call("S", "Paracoccus", 46, "0.970"),
                                         call("S", "Rivihabitans", 42, "0.973")], settings=settings)
        self.assertEqual(outcome, BACKGROUND)
        self.assertIn("Rivihabitans: weak (read identity 0.973)", reason)

    def test_one_strong_call_that_recurs_nowhere_keeps_it_held(self):
        settings = SplitSettings(recurrence=recurrence(*[(f"R{i}", "Clostridium", 0.3) for i in range(5)]))
        outcome, _, reason = self.split([call("S", "Massilia", 119428, "0.993"),
                                         call("S", "Microbacterium", 94326, "0.997"),
                                         call("S", "Clostridium", 14043, "0.943")], settings=settings)
        self.assertEqual(outcome, HELD)
        self.assertIn("Massilia is not background-compatible", reason)

    def test_a_dominant_call_escapes_cross_run_recurrence(self):
        settings = SplitSettings(recurrence=recurrence(*[(f"R{i}", "Pseudomonas", 0.2) for i in range(6)]))
        self.assertEqual(self.split([call("S", "Pseudomonas", 2926), call("S", "Moraxella", 20)],
                                    settings=settings)[0], UNCONFIRMED)

    def test_the_run_is_never_compared_with_itself(self):
        settings = SplitSettings(recurrence=recurrence(
            ("RUN1", "Moraxella", 0.9), ("R2", "Moraxella", 0.9), ("R3", "Moraxella", 0.9)))
        self.assertEqual(recurrent_elsewhere("RUN1", settings), set())     # only 2 other runs
        self.assertEqual(recurrent_elsewhere("R9", settings), {"Moraxella"})


class OutcomeWiringTests(unittest.TestCase):
    def outcome(self, rows, reads=1_000_000, settings=None):
        settings = settings or OutcomeSettings(split=SplitSettings())
        return build_outcomes("RUN1", ["S"], id_rows=rows, amr_rows=[], depth={"S": reads},
                              settings=settings)[0]

    def test_unconfirmed_dna_names_the_genus_in_the_statement(self):
        row = self.outcome([call("S", "Stenotrophomonas", 1000), call("S", "Moraxella", 20)])
        self.assertEqual(row["outcome"], "dna_detected_unconfirmed")
        self.assertTrue(row["statement"].startswith(
            "Bacterial DNA detected, most consistent with Stenotrophomonas; not identified"))
        self.assertEqual(row["reported_organisms"], "")

    def test_background_profile_reads_as_a_negative_and_is_flagged(self):
        row = self.outcome([call("S", "Moraxella", 60, "0.97"), call("S", "Paracoccus", 50, "0.97")])
        self.assertEqual(row["outcome"], "background_profile")
        self.assertEqual(row["statement"], "No bacteria present above reporting thresholds")
        self.assertIn("background_profile", row["flags"])

    def test_a_shallow_background_profile_is_insufficient(self):
        row = self.outcome([call("S", "Moraxella", 60, "0.97"), call("S", "Paracoccus", 50, "0.97")],
                           reads=5_000)
        self.assertEqual(row["outcome"], "insufficient_reads")
        self.assertIn("background_profile", row["flags"])

    def test_the_split_can_be_turned_off(self):
        row = self.outcome([call("S", "Stenotrophomonas", 1000), call("S", "Moraxella", 20)],
                           settings=OutcomeSettings())       # split disabled
        self.assertEqual(row["outcome"], "held_for_review")

    def test_settings_come_from_the_config(self):
        split = outcome_settings(load_settings().cfg).split
        self.assertTrue(split.enabled)
        self.assertEqual((split.run_prevalence_min, split.cross_runs_min, split.min_read_identity),
                         (0.30, 3, 0.99))
        self.assertEqual(split.recurrence, ())            # no table configured


class RecurrenceCliTests(unittest.TestCase):
    def test_recurrence_table_from_run_folders(self):
        from bridge.cli import main
        tmp = Path(tempfile.mkdtemp())
        for run, samples in (("R1", {"A": ["Moraxella"], "B": ["Moraxella", "Proteus"],
                                     "NEG_1": ["Moraxella"]}),
                             ("R2", {"C": ["Paracoccus"]})):
            (tmp / run).mkdir()
            with open(tmp / run / "identification_summary_all.csv", "w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["sample", "run_id", "genus", "status"])
                writer.writeheader()
                for sample, genera in samples.items():
                    for genus in genera:
                        writer.writerow({"sample": sample, "run_id": run, "genus": genus,
                                         "status": "review"})
        buffer = io.StringIO()
        with redirect_stdout(buffer):
            code = main(["recurrence", str(tmp / "R1"), str(tmp / "R2"), "--out", str(tmp / "rec.csv")])
        self.assertEqual(code, 0)
        rows = {(r["run_id"], r["genus"]): r for r in csv.DictReader(open(tmp / "rec.csv"))}
        # the blank is not a sample: Moraxella in 2 of R1's 2 samples
        self.assertEqual(rows[("R1", "Moraxella")]["prevalence"], "1.0000")
        self.assertEqual(rows[("R1", "Proteus")]["n_samples"], "2")
        self.assertEqual(rows[("R2", "Paracoccus")]["n_samples_with"], "1")


if __name__ == "__main__":
    unittest.main()
