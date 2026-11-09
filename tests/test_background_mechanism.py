"""The mechanism-attributed contamination gate.

Reagent background is a fraction of a library, so RPM-r is its comparator.
Cross-contamination from a run positive is an absolute count, so the blank's
read count is. The numbers here model a 22,539-read blank against samples of
~1-3 M reads.
"""

from __future__ import annotations

import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from bridge.reconcile import baseline as baseline_mod
from bridge.reconcile.config import load_config
from bridge.reconcile.rpm import (
    BASELINE_COLUMNS, MECHANISM_CROSS, MECHANISM_REAGENT, Baseline, CohortEvidence,
    as_baseline, baselines_from_rows, build_baselines, evaluate, load_baselines,
    rpm_detail, ungated, write_baselines,
)
from bridge.reconcile.summary_csv import SUMMARY_COLUMNS, _rank, summary_rows

CONFIG = load_config(Path(__file__).resolve().parents[1] / "config.yaml")
NTC_TOTAL = 22_539


class MechanismRoutingTests(unittest.TestCase):

    def test_cross_contamination_is_gated_on_absolute_reads(self):
        """077X S. aureus: 725 est_counts vs 29.4 in the blank.

        RPM-r says 0.29 (the blank is 86x shallower); abs-r says 24.6. The
        blank's S. aureus is carry-over from this very organism elsewhere on
        the run, so the count is the right comparator.
        """
        blank = Baseline(rpm_ntc=1305.29, reads_ntc=29.42, total_reads_ntc=NTC_TOTAL,
                         mechanism=MECHANISM_CROSS, source_sample="SAMPLE_077X",
                         source_reads=724.9)
        result = evaluate("Staphylococcus aureus", "rmlst", 724.9, 1_943_956, blank)
        self.assertEqual(result.rule, "abs_r")
        self.assertTrue(result.passed)
        self.assertAlmostEqual(result.abs_r, 24.64, places=1)
        self.assertAlmostEqual(result.rpm_r, 0.286, places=2)     # still reported
        self.assertEqual(result.band, "moderate")
        self.assertIn("abs-r 24.6 >= 10", result.reason)
        self.assertIn("cross_contamination", result.reason)
        self.assertEqual(result.source_sample, "SAMPLE_077X")

    def test_reagent_background_is_unchanged(self):
        """082 Bacillus tuaregi: 46,036 reads at 24 M depth vs 15 in the blank.

        Kitome in a deep, low-biomass sample: 3,069x the blank's count, but
        RPM-r 2.9. Reagent stays on RPM-r and this stays background.
        """
        blank = Baseline(rpm_ntc=665.513, reads_ntc=15.0, total_reads_ntc=NTC_TOTAL)
        result = evaluate("Bacillus tuaregi", "marker16s", 46_036.0, 24_014_606, blank)
        self.assertEqual(result.rule, "rpm_r")
        self.assertFalse(result.passed)
        self.assertEqual(result.band, "background")
        self.assertAlmostEqual(result.abs_r, 3069.07, places=1)  # computed, not used
        self.assertIn("RPM-r", result.reason)
        self.assertIn("[reagent]", result.reason)

    def test_cross_without_a_read_count_falls_back_to_rpm_r_and_says_so(self):
        blank = Baseline(rpm_ntc=200.0, mechanism=MECHANISM_CROSS)
        result = evaluate("X", "rmlst", 5000, 1_000_000, blank)
        self.assertEqual(result.rule, "rpm_r")
        self.assertIsNone(result.abs_r)
        self.assertIn("reads_ntc unknown", result.reason)

    def test_the_sample_read_floor_applies_under_abs_r(self):
        blank = Baseline(rpm_ntc=44.0, reads_ntc=0.5, mechanism=MECHANISM_CROSS)
        result = evaluate("X", "marker16s", 8.0, 1_000_000, blank, min_sample_reads=10)
        self.assertEqual(result.rule, "abs_r")
        self.assertAlmostEqual(result.abs_r, 8.0)     # reads / max(reads_ntc, 1)
        self.assertFalse(result.passed)
        self.assertIn("floor 10", result.reason)
        enough = evaluate("X", "marker16s", 12.0, 1_000_000, blank, min_sample_reads=10)
        self.assertTrue(enough.passed)

    def test_a_bare_float_is_a_reagent_baseline(self):
        """Callers may hand in the control RPM alone."""
        self.assertEqual(as_baseline(200.0), Baseline(rpm_ntc=200.0))
        self.assertIsNone(as_baseline(None))
        result = evaluate("X", "rmlst", 5000, 1_000_000, 200.0)
        self.assertEqual(result.rule, "rpm_r")
        self.assertEqual(result.mechanism, MECHANISM_REAGENT)
        self.assertIsNone(result.abs_r)

    def test_absent_from_controls_is_the_floor_rule(self):
        result = evaluate("X", "rmlst", 100, 1_000_000, None)
        self.assertEqual(result.rule, "floor")
        self.assertEqual(result.band, "no_control_background")

    def test_ungated_computes_rpm_and_tests_nothing(self):
        result = ungated("X", "rmlst", 5000, 1_000_000)
        self.assertEqual(result.rule, "not_applied")
        self.assertAlmostEqual(result.rpm, 5000.0)
        self.assertIsNone(result.rpm_r)
        self.assertTrue(result.passed)
        self.assertEqual(result.band, "ungated")
        detail = rpm_detail(result)
        self.assertEqual(detail["rule"], "not_applied")
        self.assertEqual(detail["band"], "ungated")

    def test_the_detail_carries_both_comparators(self):
        blank = Baseline(rpm_ntc=3484.54, reads_ntc=78.54, total_reads_ntc=NTC_TOTAL,
                         mechanism=MECHANISM_CROSS, source_sample="SAMPLE_001")
        detail = rpm_detail(evaluate("Enterobacter hormaechei", "marker16s",
                                     15_488.6, 715_078, blank))
        for key in ("rpm_ntc", "reads_ntc", "total_reads_ntc", "abs_r", "mechanism",
                    "rule", "source_sample"):
            self.assertIn(key, detail)
        self.assertAlmostEqual(detail["abs_r"], 197.2, places=0)
        self.assertEqual(detail["rule"], "abs_r")
        self.assertNotIn("other_arm_background", detail)


class BaselineFileTests(unittest.TestCase):

    def test_a_legacy_six_column_file_loads_as_reagent(self):
        path = Path(tempfile.mkdtemp()) / "rpm_baselines.csv"
        path.write_text(
            "run_id,arm,organism,rpm_ntc,n_controls,source\n"
            "RUN1,rmlst,Salmonella enterica,1305.29,1,max_across_controls\n"
            "RUN1,marker16s,Acidovorax,803.008,1,max_across_controls\n",
            encoding="utf-8")
        baselines = load_baselines(path)
        entry = baselines[("RUN1", "rmlst", "salmonella enterica")]
        self.assertEqual(entry.mechanism, MECHANISM_REAGENT)
        self.assertIsNone(entry.reads_ntc)
        self.assertAlmostEqual(entry.rpm_ntc, 1305.29)
        self.assertEqual(entry.n_controls, 1)

    def test_the_new_schema_round_trips(self):
        rows = build_baselines(
            {"Neg_A": {"Salmonella enterica": (516.0, NTC_TOTAL)}},
            run_id="RUN1", arm="marker16s",
            cohort=CohortEvidence({"Salmonella enterica": ("SAMPLE_003", 97_155.4)}))
        self.assertEqual(list(rows[0]), BASELINE_COLUMNS)
        row = rows[0]
        self.assertEqual(row["mechanism"], MECHANISM_CROSS)
        self.assertEqual(row["source_sample"], "SAMPLE_003")
        self.assertEqual(row["reads_ntc"], "516")
        self.assertEqual(row["total_reads_ntc"], NTC_TOTAL)
        self.assertEqual(row["control_ids"], "Neg_A")
        self.assertAlmostEqual(float(row["ntc_over_source_frac"]), 0.00531, places=4)
        self.assertEqual(row["note"], "")

        path = Path(tempfile.mkdtemp()) / "rpm_baselines.csv"
        write_baselines(rows, path)
        loaded = load_baselines(path)[("RUN1", "marker16s", "salmonella enterica")]
        self.assertEqual(loaded, baselines_from_rows(rows)[("RUN1", "marker16s", "salmonella enterica")])
        self.assertEqual(loaded.mechanism, MECHANISM_CROSS)
        self.assertAlmostEqual(loaded.reads_ntc, 516.0)
        self.assertEqual(loaded.total_reads_ntc, NTC_TOTAL)
        self.assertEqual(loaded.source_sample, "SAMPLE_003")

    def test_a_high_ntc_over_source_fraction_is_noted(self):
        rows = build_baselines(
            {"Neg_A": {"Escherichia": (1000.0, NTC_TOTAL)}}, run_id="R", arm="marker16s",
            cohort=CohortEvidence(genera={"Escherichia": ("SAMPLE_070", 12_535.5)}))
        self.assertEqual(rows[0]["mechanism"], MECHANISM_CROSS)
        self.assertIn("high for cross-contamination", rows[0]["note"])

    def test_without_cohort_evidence_everything_is_reagent(self):
        rows = build_baselines({"Neg_A": {"Salmonella enterica": (516.0, NTC_TOTAL)}},
                               run_id="R", arm="rmlst")
        self.assertEqual(rows[0]["mechanism"], MECHANISM_REAGENT)
        self.assertEqual(rows[0]["source_sample"], "")

    def test_the_control_that_set_the_maximum_supplies_the_counts(self):
        rows = build_baselines({
            "NEG1": {"Escherichia coli": (10.0, 1_000_000)},
            "NEG2": {"Escherichia coli": (500.0, 1_000_000)},
        }, run_id="RUN1", arm="rmlst")
        self.assertEqual(rows[0]["reads_ntc"], "500")
        self.assertEqual(rows[0]["n_controls"], 2)
        self.assertEqual(rows[0]["control_ids"], "NEG1; NEG2")

    def test_cohort_lookup_falls_from_species_to_genus(self):
        evidence = CohortEvidence(species={"Enterobacter hormaechei": ("001", 15_488.6)},
                                  genera={"Enterobacter": ("001", 15_488.6)})
        self.assertEqual(evidence.source_for("Enterobacter hormaechei"), ("001", 15_488.6))
        self.assertEqual(evidence.source_for("Enterobacter quasihormaechei"), ("001", 15_488.6))
        self.assertEqual(evidence.source_for("Enterobacter"), ("001", 15_488.6))
        self.assertIsNone(evidence.source_for("Bacillus tuaregi"))


RMLST_HEADER = ("sample_id,species,reported,call_method,exact_profile_loci,"
                "exact_profile_unambiguous_loci,pi,allele_evidence_unambiguous_loci,"
                "allele_evidence_ambiguous_loci,em_evidence_unambiguous_loci,"
                "em_evidence_ambiguous_loci,unambiguous_est_counts,ambiguous_est_counts,"
                "alternative_to\n")
S16_HEADER = ("sample_id,taxon,genus,species,rank,strength,abundance_frac,genus_frac,"
              "reads_assigned,raw_reads,copy_number,within_genus_share,margin,"
              "n_references,competing_species,reason\n")


class CohortScanTests(unittest.TestCase):
    """The credible-source rule on a synthetic run: a blank, a strong positive,
    a marginal call, a junk folder."""

    def setUp(self):
        self.root = Path(tempfile.mkdtemp())

    def sample(self, name, rmlst_rows, s16_rows, depth):
        folder = self.root / name
        (folder / "rmlst").mkdir(parents=True)
        (folder / "16s" / "quant").mkdir(parents=True)
        (folder / "rmlst" / "species_calls.csv").write_text(RMLST_HEADER + "".join(rmlst_rows))
        (folder / "16s" / "marker16s_calls.csv").write_text(S16_HEADER + "".join(s16_rows))
        (folder / "16s" / "quant" / "run_info.json").write_text(json.dumps({"n_processed": depth}))

    @staticmethod
    def rmlst(sample, species, reported, est):
        return f"{sample},{species},{reported},x,30,25,0.9,25,3,25,3,{est},0.0,\n"

    @staticmethod
    def s16(sample, taxon, genus, raw):
        return f"{sample},{taxon},{genus},,genus,genus_only,0.5,0.5,{raw/7:.3f},{raw},7.0,0.9,0.9,1,{taxon},t\n"

    def test_a_strong_positive_makes_the_blanks_entry_cross_contamination(self):
        # blank: Salmonella 516 reads (16S) / 29.4 est (rMLST) on 22,539 reads
        self.sample("Neg_1",
                    [self.rmlst("Neg_1", "Salmonella enterica", False, 29.42)],
                    [self.s16("Neg_1", "Salmonella", "Salmonella", 515.7)], NTC_TOTAL)
        # 003: rMLST-reported Salmonella, 102,190 est / 97,155 16S reads on 1.83 M
        self.sample("SAMPLE_003",
                    [self.rmlst("SAMPLE_003", "Salmonella enterica", True, 102_190.0)],
                    [self.s16("SAMPLE_003", "Salmonella", "Salmonella", 97_155.4)], 1_831_048)
        # a marginal Corynebacterium call that fails RPM-r on BOTH arms must not flip its genus
        self.sample("Neg_2",
                    [self.rmlst("Neg_2", "Corynebacterium resistens", False, 40.0)],
                    [self.s16("Neg_2", "Corynebacterium", "Corynebacterium", 60.0)], NTC_TOTAL)
        self.sample("SAMPLE_marginal",
                    [self.rmlst("SAMPLE_marginal", "Corynebacterium resistens", True, 50.0)],
                    [self.s16("SAMPLE_marginal", "Corynebacterium", "Corynebacterium", 80.0)],
                    2_000_000)
        (self.root / "oe").mkdir()
        (self.root / "NTC_ignored").mkdir()

        rows = baseline_mod.build_run_baselines(self.root, ["Neg_1", "Neg_2"], "R")
        by_key = {(r["arm"], r["organism"]): r for r in rows}

        salmonella_16s = by_key[("marker16s", "Salmonella")]
        self.assertEqual(salmonella_16s["mechanism"], MECHANISM_CROSS)
        self.assertEqual(salmonella_16s["source_sample"], "SAMPLE_003")
        self.assertAlmostEqual(float(salmonella_16s["source_reads"]), 97_155.4)
        salmonella_rmlst = by_key[("rmlst", "Salmonella enterica")]
        self.assertEqual(salmonella_rmlst["mechanism"], MECHANISM_CROSS)
        self.assertAlmostEqual(float(salmonella_rmlst["source_reads"]), 102_190.0)

        self.assertEqual(by_key[("rmlst", "Corynebacterium resistens")]["mechanism"], MECHANISM_REAGENT)
        self.assertEqual(by_key[("marker16s", "Corynebacterium")]["mechanism"], MECHANISM_REAGENT)

        # and the gate then behaves accordingly
        baselines = baselines_from_rows(rows)
        verdict = evaluate("Salmonella enterica", "marker16s", 97_155.4, 1_831_048,
                           baselines[("R", "marker16s", "salmonella")])
        self.assertEqual(verdict.rule, "abs_r")
        self.assertTrue(verdict.passed)
        old = evaluate("Salmonella enterica", "marker16s", 97_155.4, 1_831_048,
                       baselines[("R", "marker16s", "salmonella")].rpm_ntc)
        self.assertFalse(old.passed)

    def test_no_cohort_reproduces_the_reagent_only_file(self):
        self.sample("Neg_1", [self.rmlst("Neg_1", "Salmonella enterica", False, 29.42)],
                    [self.s16("Neg_1", "Salmonella", "Salmonella", 515.7)], NTC_TOTAL)
        self.sample("S", [self.rmlst("S", "Salmonella enterica", True, 102_190.0)],
                    [self.s16("S", "Salmonella", "Salmonella", 97_155.4)], 1_831_048)
        rows = baseline_mod.build_run_baselines(self.root, ["Neg_1"], "R", cohort=False)
        self.assertTrue(all(r["mechanism"] == MECHANISM_REAGENT for r in rows))
        self.assertEqual(baseline_mod.sample_dirs(self.root, ["Neg_1"]), [self.root / "S"])


class SummaryTableTests(unittest.TestCase):

    CALL = {"organism": "Salmonella enterica", "genus": "Salmonella",
            "reported_rank": "species", "basis": "rmlst",
            "rmlst_taxon": "Salmonella enterica", "16s_taxon": "Salmonella enterica",
            "16s_margin": 1.0, "flags": [], "notes": ""}

    def record(self, rpm, flags=None):
        return {"sample": "SAMPLE_003", "run_id": "RUN1", "outcome": "reported",
                "flags": flags or [], "calls": [self.CALL], "rpm": rpm}

    def test_the_new_columns_are_written(self):
        for column in ("abs_r_rmlst", "abs_r_16s", "bg_mechanism"):
            self.assertIn(column, SUMMARY_COLUMNS)
        rmlst = rpm_detail(evaluate("Salmonella enterica", "rmlst", 102_190.0, 1_831_048,
                                    Baseline(1305.29, reads_ntc=29.42, mechanism=MECHANISM_CROSS)))
        s16 = rpm_detail(evaluate("Salmonella enterica", "marker16s", 97_155.4, 1_831_048,
                                  Baseline(22_881.8, reads_ntc=515.7, mechanism=MECHANISM_CROSS)))
        row = summary_rows(self.record({"rmlst": {"Salmonella enterica": rmlst},
                                        "marker16s": {"Salmonella enterica": s16}}))[0]
        self.assertEqual(row["bg_mechanism"], MECHANISM_CROSS)
        self.assertEqual(float(row["abs_r_16s"]), 188.4)      # 97,155.4 / 515.7
        self.assertEqual(float(row["abs_r_rmlst"]), 3473.0)    # 102,190 / 29.42
        self.assertEqual(row["rpm_band"], "high")
        self.assertEqual(row["status"], "report")

    def test_a_record_without_rule_still_ranks_by_rpm_r(self):
        self.assertEqual(_rank({"band": "moderate", "rpm_r": 12.0}), (1, 3))
        self.assertEqual(_rank({"band": "moderate", "rule": "abs_r", "rpm_r": None}), (1, 3))
        self.assertEqual(_rank({"band": "ungated", "rule": "not_applied"}), (2, 2))

    def test_sample_level_flags_reach_every_organism_row(self):
        detail = rpm_detail(ungated("Salmonella enterica", "rmlst", 102_190.0, 1_831_048))
        row = summary_rows(self.record({"rmlst": {"Salmonella enterica": detail}, "marker16s": {}},
                                       flags=["rpm_gate_not_applied"]))[0]
        self.assertIn("rpm_gate_not_applied", row["flags"])
        self.assertEqual(row["rpm_band"], "ungated")
        self.assertEqual(float(row["rpm_rmlst"]), 55810.0)
        self.assertEqual(row["rpm_r_rmlst"], "")
        self.assertEqual(row["status"], "report")


class GateWiringWithoutBaselineTests(unittest.TestCase):

    def test_no_baseline_computes_rpm_and_flags_the_sample(self):
        from bridge.reconcile.adapters.rmlst import RmlstRow, rmlst_results
        from bridge.reconcile.reconcile_cli import _gate_results

        row = RmlstRow(
            sample_id="s", species="Klebsiella pneumoniae", reported=True,
            call_method="exact_allele_profile", pi=0.5, exact_profile_loci=25,
            exact_profile_unambiguous_loci=20, unambiguous_loci=20, ambiguous_loci=3,
            em_evidence_unambiguous_loci=20, em_evidence_ambiguous_loci=5,
            unambiguous_est_counts=1000.0, ambiguous_est_counts=100.0,
            mean_callable_breadth=0.94, min_callable_breadth=0.61,
            median_projected_depth=212.0, min_projected_depth=18.0,
            total_callable_bases=21847.0, total_allele_length=23100.0,
            supporting_fragments=4192.0, alternative_to="")
        results = rmlst_results([row], CONFIG)
        gated, detail, flags = _gate_results(
            results, "rmlst", {"klebsiella pneumoniae": 50_000.0}, None, "R", 1_000_000, CONFIG)
        self.assertEqual(flags, ["rpm_gate_not_applied"])
        self.assertAlmostEqual(detail["Klebsiella pneumoniae"]["rpm"], 50_000.0)
        self.assertEqual(detail["Klebsiella pneumoniae"]["rule"], "not_applied")
        self.assertNotIn("rpm_below_background", gated[0].flags)

    def test_no_read_total_still_short_circuits(self):
        from bridge.reconcile.reconcile_cli import _gate_results
        gated, detail, flags = _gate_results([], "rmlst", {}, None, "R", None, CONFIG)
        self.assertEqual((gated, detail, flags), ([], {}, ["rpm_gate_not_applied"]))


if __name__ == "__main__":
    unittest.main()
