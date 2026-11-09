"""bridge finalise on a synthetic results tree (bridge.finalise)."""

from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path

from bridge.finalise import (
    detect_controls, finalise_run, find_arm_outputs, sample_dirs,
)
from bridge.settings import load_settings

RMLST_HEADER = (
    "sample_id,species,reported,call_method,exact_profile_loci,"
    "exact_profile_unambiguous_loci,pi,allele_evidence_unambiguous_loci,"
    "allele_evidence_ambiguous_loci,em_evidence_unambiguous_loci,"
    "em_evidence_ambiguous_loci,unambiguous_est_counts,ambiguous_est_counts,"
    "mean_callable_breadth,min_callable_breadth,median_projected_depth,"
    "min_projected_depth,total_callable_bases,total_allele_length,"
    "supporting_fragments,alternative_to,matching_rSTs,matching_rST_count,"
    "allele_fasta,allele_fasta_allele_count\n"
)
S16_HEADER = (
    "sample_id,taxon,genus,species,rank,strength,abundance_frac,genus_frac,"
    "reads_assigned,raw_reads,copy_number,within_genus_share,margin,"
    "n_references,competing_species,reason\n"
)


def rmlst_row(sample, species, reported=True, est=5000.0):
    method = "exact_allele_profile" if reported else "em_fallback"
    return (f"{sample},{species},{reported},{method},30,25,0.9,25,3,25,3,{est},100.0,"
            f"0.95,0.8,100.0,20.0,20000.0,21000.0,4000.0,,1,1,,0\n")


def s16_row(sample, taxon, genus, species, strength, raw=7000.0):
    rank = "species" if species else "genus"
    return (f"{sample},{taxon},{genus},{species},{rank},{strength},0.9,0.95,"
            f"{raw / 7:.3f},{raw},7.0,0.95,0.9,2,{taxon}; other,test\n")


class FinaliseTests(unittest.TestCase):
    def setUp(self):
        self.settings = load_settings()
        self.root = Path(tempfile.mkdtemp())
        (self.root / "oe").mkdir()

    def sample(self, name, rmlst_rows, s16_rows, total=1_000_000, flat=False):
        folder = self.root / name
        rmlst_dir = folder if flat else folder / "rmlst"
        s16_dir = folder / ("marker16s" if flat else "16s")
        (s16_dir / "quant").mkdir(parents=True)
        rmlst_dir.mkdir(parents=True, exist_ok=True)
        (rmlst_dir / "species_calls.csv").write_text(RMLST_HEADER + "".join(rmlst_rows))
        (s16_dir / "marker16s_calls.csv").write_text(S16_HEADER + "".join(s16_rows))
        (s16_dir / "quant" / "run_info.json").write_text(json.dumps({"n_processed": total}))

    def test_a_full_run_reconciles_everything_but_the_control(self):
        self.sample("A", [rmlst_row("A", "Escherichia coli")],
                    [s16_row("A", "Escherichia coli", "Escherichia", "coli", "species_strong")])
        self.sample("B", [rmlst_row("B", "Klebsiella pneumoniae")],
                    [s16_row("B", "Klebsiella pneumoniae", "Klebsiella", "pneumoniae", "species_strong")])
        self.sample("C_flat", [rmlst_row("C_flat", "Staphylococcus aureus")],
                    [s16_row("C_flat", "Staphylococcus", "Staphylococcus", "", "genus_only")], flat=True)
        self.sample("Neg_1", [rmlst_row("Neg_1", "Escherichia coli", reported=False, est=20.0)],
                    [s16_row("Neg_1", "Escherichia", "Escherichia", "", "genus_only", raw=70.0)],
                    total=22_539)
        (self.root / "no_arms").mkdir()

        self.assertEqual(detect_controls(self.root, self.settings.hpc.controls_prefixes), ["Neg_1"])
        self.assertEqual([p.name for p in sample_dirs(self.root)], ["A", "B", "C_flat", "Neg_1", "no_arms"])
        rmlst, s16 = find_arm_outputs(self.root / "C_flat")
        self.assertEqual(rmlst, self.root / "C_flat" / "species_calls.csv")
        self.assertEqual(s16, self.root / "C_flat" / "marker16s" / "marker16s_calls.csv")

        result = finalise_run(self.settings, self.root, "RUN1")

        self.assertEqual(result.controls, ["Neg_1"])
        self.assertEqual(result.baseline_path, self.root / "rpm_baselines.csv")
        self.assertGreater(result.baseline_rows, 0)
        self.assertEqual(result.reconciled, ["A", "B", "C_flat"])
        self.assertIn(("Neg_1", "negative control"), result.skipped)
        self.assertTrue(any(name == "no_arms" for name, _ in result.skipped))
        self.assertEqual(result.failed, [])
        self.assertEqual(result.combined_samples, 3)
        for name in ("A", "B", "C_flat"):
            self.assertTrue((self.root / name / "bacterial_ID" / "identification_summary.csv").is_file())
            self.assertTrue((self.root / name / "bacterial_ID" / "reconciled_calls.json").is_file())

        with open(result.combined_path, newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        by_sample = {row["sample"]: row for row in rows}
        self.assertEqual(set(by_sample), {"A", "B", "C_flat"})
        self.assertEqual(by_sample["A"]["run_id"], "RUN1")
        # A's E. coli sits against the control's Escherichia; B's Klebsiella has no background.
        self.assertEqual(by_sample["A"]["rpm_band"], "background")
        self.assertEqual(by_sample["B"]["rpm_band"], "no_control_background")
        self.assertEqual(by_sample["B"]["status"], "report")

        # what the JSON record carried
        record = json.loads((self.root / "A" / "bacterial_ID" / "reconciled_calls.json").read_text())
        self.assertEqual(record["total_raw_reads"], 1_000_000)

    def test_without_a_control_the_run_is_reconciled_ungated(self):
        self.sample("A", [rmlst_row("A", "Escherichia coli")],
                    [s16_row("A", "Escherichia coli", "Escherichia", "coli", "species_strong")])
        result = finalise_run(self.settings, self.root, "RUN1")
        self.assertEqual(result.controls, [])
        self.assertIsNone(result.baseline_path)
        self.assertEqual(result.reconciled, ["A"])
        self.assertEqual(result.combined_samples, 1)

    def test_several_16s_fallback_calls_are_held_as_review_override(self):
        """End to end: config default -> reconciler flag -> summary status."""
        # rMLST resolved nothing (a sub-bar row only); 16S names three genera.
        empty_rmlst = [rmlst_row("N", "Escherichia coli", reported=False, est=20.0)]
        self.sample("N", empty_rmlst,
                    [s16_row("N", "Paracoccus yeei", "Paracoccus", "yeei", "species_strong"),
                     s16_row("N", "Moraxella tetraodonis", "Moraxella", "tetraodonis",
                             "species_strong"),
                     s16_row("N", "Sphingomonas", "Sphingomonas", "", "genus_only")])
        self.sample("S", empty_rmlst,
                    [s16_row("S", "Streptococcus pyogenes", "Streptococcus", "pyogenes",
                             "species_strong")])
        result = finalise_run(self.settings, self.root, "RUN1")
        with open(result.combined_path, newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        held = [row for row in rows if row["sample"] == "N"]
        self.assertEqual(len(held), 3)
        for row in held:
            self.assertEqual((row["basis"], row["status"]), ("16s_fallback", "review_override"))
            self.assertIn("multiple_16s_fallback", row["flags"])
        single = [row for row in rows if row["sample"] == "S"]
        self.assertEqual([row["status"] for row in single], ["report"])

    def test_a_background_fallback_call_does_not_count_towards_the_rule(self):
        """Through the real gate: Escherichia is in the blank, Moraxella is not."""
        empty_rmlst = [rmlst_row("M", "Escherichia coli", reported=False, est=20.0)]
        self.sample("M", empty_rmlst,
                    [s16_row("M", "Escherichia", "Escherichia", "", "genus_only"),
                     s16_row("M", "Moraxella tetraodonis", "Moraxella", "tetraodonis",
                             "species_strong")])
        self.sample("Neg_1", [rmlst_row("Neg_1", "Escherichia coli", reported=False, est=20.0)],
                    [s16_row("Neg_1", "Escherichia", "Escherichia", "", "genus_only", raw=70.0)],
                    total=22_539)
        result = finalise_run(self.settings, self.root, "RUN1")
        with open(result.combined_path, newline="", encoding="utf-8") as handle:
            status = {row["organism"]: row["status"] for row in csv.DictReader(handle)}
        # Two fallback calls, but only one eligible for report: the rule stays quiet.
        self.assertEqual(status, {"Escherichia": "review", "Moraxella tetraodonis": "report"})

    def test_explicit_control_name_wins_over_prefixes(self):
        self.sample("A", [rmlst_row("A", "Escherichia coli")],
                    [s16_row("A", "Escherichia coli", "Escherichia", "coli", "species_strong")])
        self.sample("Water", [rmlst_row("Water", "Escherichia coli", reported=False, est=20.0)],
                    [s16_row("Water", "Escherichia", "Escherichia", "", "genus_only", raw=70.0)],
                    total=22_539)
        result = finalise_run(self.settings, self.root, "RUN1", controls=["Water"])
        self.assertEqual(result.controls, ["Water"])
        self.assertEqual(result.reconciled, ["A"])
        self.assertIn(("Water", "negative control"), result.skipped)


if __name__ == "__main__":
    unittest.main()
