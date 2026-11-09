"""Tests for the adapters, the RPM-r gate and the rMLST-primary reconciler.

Most of these pin a rule to the failure it exists to prevent.
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from bridge.reconcile.adapters.rmlst import (
    RmlstRow, load_rmlst_calls, rmlst_results,
)
from bridge.reconcile.config import load_config
from bridge.reconcile.rank_marker16s import Marker16SRow, rank_marker16s
from bridge.reconcile.reconcile import (
    CALL_BASES,
    MULTI_FALLBACK_FLAG,
    flag_multi_fallback,
    reconcile_sample,
    sample_outcome,
)
from bridge.reconcile.rpm import (
    MIN_NTC_RPM, build_baselines, compute_rpm, evaluate, lookup_baseline,
)
from bridge.reconcile.strength import (
    BlastStrength, BlastStrengthResult, RmlstStrength, RmlstStrengthResult,
)

CONFIG = load_config(Path(__file__).resolve().parents[1] / "config.yaml")


def rmlst_row(species, reported=True, alternative_to="", unamb=20, pi=0.5,
              unamb_counts=1000.0, amb_counts=100.0, em_unamb=20, em_amb=5,
              breadth=0.94, depth=212.0):
    return RmlstRow(
        sample_id="s", species=species, reported=reported,
        call_method="exact_allele_profile" if reported else "em_fallback",
        pi=pi, exact_profile_loci=unamb + 5, exact_profile_unambiguous_loci=unamb,
        unambiguous_loci=unamb, ambiguous_loci=3,
        em_evidence_unambiguous_loci=em_unamb, em_evidence_ambiguous_loci=em_amb,
        unambiguous_est_counts=unamb_counts, ambiguous_est_counts=amb_counts,
        mean_callable_breadth=breadth, min_callable_breadth=0.61,
        median_projected_depth=depth, min_projected_depth=18.0,
        total_callable_bases=21847.0, total_allele_length=23100.0,
        supporting_fragments=4192.0,
        alternative_to=alternative_to,
    )


def s16_row(taxon, genus, species, strength, margin=0.9, reads=5000.0):
    return Marker16SRow(
        sample_id="s", taxon=taxon, genus=genus, species=species,
        rank="species" if species else "genus", strength=strength,
        abundance_frac=0.5, genus_frac=0.5, reads_assigned=reads,
        raw_reads=reads * 4, copy_number=4.0, within_genus_share=0.9,
        margin=margin, n_references=5, competing_species=f"{taxon}; other sp",
        reason="test",
    )


class RmlstAdapterTests(unittest.TestCase):

    def test_reported_species_becomes_species_supported(self):
        results = rmlst_results([rmlst_row("Klebsiella pneumoniae")], CONFIG)
        self.assertEqual(len(results), 1)
        self.assertIs(results[0].strength, RmlstStrength.SPECIES_SUPPORTED)
        self.assertEqual(results[0].taxon.species_full, "Klebsiella pneumoniae")

    def test_demoted_alternative_becomes_genus_only_not_silent(self):
        """A second name for an organism already called still corroborates the genus.

        Dropping it would discard a real observation; promoting it would
        reinstate the over-call the reporting bar was introduced to remove.
        """
        results = rmlst_results([
            rmlst_row("Streptococcus mitis"),
            rmlst_row("Streptococcus oralis", reported=False,
                      alternative_to="Streptococcus mitis", unamb=1),
        ], CONFIG)
        by_strength = {r.strength: r for r in results}
        self.assertIn(RmlstStrength.GENUS_ONLY, by_strength)
        demoted = by_strength[RmlstStrength.GENUS_ONLY]
        self.assertEqual(demoted.taxon.genus, "Streptococcus")
        self.assertIn("rmlst_alternative_name", demoted.flags)

    def test_a_demoted_name_keeps_its_species(self):
        """The name IS the information; withholding it emptied the row.

        `basis=rmlst_alternative` already says the row is not a reported call,
        and a reader needs "Streptococcus oralis, an alternative
        name for Streptococcus mitis" - not "Streptococcus".
        """
        results = rmlst_results([
            rmlst_row("Streptococcus mitis"),
            rmlst_row("Streptococcus oralis", reported=False,
                      alternative_to="Streptococcus mitis", unamb=1),
        ], CONFIG)
        demoted = {r.strength: r for r in results}[RmlstStrength.GENUS_ONLY]
        self.assertEqual(demoted.taxon.species_full, "Streptococcus oralis")

    def test_the_em_evidence_loci_are_parsed(self):
        row = rmlst_row("X y", em_unamb=17, em_amb=4)
        self.assertEqual(row.em_evidence_unambiguous_loci, 17)
        self.assertEqual(row.em_evidence_ambiguous_loci, 4)

    def test_coverage_is_carried_on_the_row(self):
        row = rmlst_row("X y", breadth=0.83, depth=97.0)
        self.assertAlmostEqual(row.mean_callable_breadth, 0.83)
        self.assertAlmostEqual(row.median_projected_depth, 97.0)

    def test_em_background_is_dropped_entirely(self):
        results = rmlst_results([
            rmlst_row("Klebsiella pneumoniae"),
            rmlst_row("Citrobacter freundii", reported=False, unamb=0),
        ], CONFIG)
        self.assertEqual([r.taxon.species_full for r in results],
                         ["Klebsiella pneumoniae"])

    def test_several_species_per_sample_are_all_emitted(self):
        """The rMLST arm reports per species, so a sample can carry several."""
        results = rmlst_results([
            rmlst_row("Klebsiella pneumoniae"),
            rmlst_row("Staphylococcus aureus"),
        ], CONFIG)
        self.assertEqual(len(results), 2)

    def test_total_est_counts_includes_ambiguous_evidence(self):
        row = rmlst_row("X y", unamb_counts=900.0, amb_counts=100.0)
        self.assertEqual(row.total_est_counts, 1000.0)


class Marker16SRankingTests(unittest.TestCase):

    def test_species_strong_carries_the_species(self):
        result = rank_marker16s(
            s16_row("Klebsiella pneumoniae", "Klebsiella", "pneumoniae",
                    "species_strong"), CONFIG)
        self.assertIs(result.strength, BlastStrength.SPECIES_STRONG)
        self.assertEqual(result.top_taxon.species_full, "Klebsiella pneumoniae")

    def test_genus_only_withholds_the_species(self):
        result = rank_marker16s(
            s16_row("Streptococcus", "Streptococcus", "", "genus_only",
                    margin=0.03), CONFIG)
        self.assertIs(result.strength, BlastStrength.GENUS_ONLY)
        self.assertEqual(result.top_taxon.genus, "Streptococcus")
        self.assertIsNone(result.top_taxon.species)
        self.assertIn("marker16s_within_genus_ambiguous", result.flags)

    def test_identity_fields_are_none_when_absent_and_filled_when_present(self):
        """A calls file without the alignment stage carries no identity;
        one with it carries the consensus's BLAST figures and the top
        reference's breadth. Never a fabricated 0."""
        from dataclasses import replace
        row = s16_row("Klebsiella pneumoniae", "Klebsiella", "pneumoniae", "species_strong")
        result = rank_marker16s(row, CONFIG)
        self.assertIsNone(result.top_pident)
        self.assertIsNone(result.top_qcovs)
        self.assertIsNone(result.coverage)
        self.assertEqual(result.coverage_matched_by, "")
        # The margin IS meaningful - it is the within-genus abundance gap.
        self.assertIsNotNone(result.margin)

        measured = rank_marker16s(replace(row, consensus_pident=99.53, consensus_qcov=98.2,
                                          ref_breadth=0.971, mean_depth=41.2), CONFIG)
        self.assertEqual(measured.top_pident, 99.53)
        self.assertEqual(measured.top_qcovs, 98.2)
        self.assertEqual(measured.coverage, 0.971)
        self.assertEqual(measured.coverage_matched_by, "16s_alignment")
        self.assertIs(measured.strength, BlastStrength.SPECIES_STRONG)   # evidence, not a verdict

    def test_unknown_strength_is_an_error_not_a_silent_default(self):
        with self.assertRaises(ValueError):
            rank_marker16s(
                s16_row("X", "X", "y", "probably_fine"), CONFIG)

    def test_the_row_index_round_trips_through_the_adapter(self):
        """It is how the report gets back to this organism's own read counts.

        Matching on the name instead silently misses whenever the reconciler
        reports an organism under a different label than the arm used.
        """
        from bridge.reconcile.adapters.marker16s import marker16s_results

        rows = [
            s16_row("Klebsiella pneumoniae", "Klebsiella", "pneumoniae",
                    "species_strong"),
            s16_row("Trace", "Trace", "", "silent_weak"),   # dropped
            s16_row("Enterococcus", "Enterococcus", "", "genus_only"),
        ]
        results = marker16s_results(rows, CONFIG)
        self.assertEqual([r.row_index for r in results], [0, 2])
        # The index must address the ORIGINAL list, not the filtered one.
        self.assertEqual(rows[results[1].row_index].taxon, "Enterococcus")


class RpmTests(unittest.TestCase):

    def test_rpm_is_per_million_of_total_reads(self):
        self.assertAlmostEqual(compute_rpm(50, 1_000_000), 50.0)

    def test_organism_above_its_background_passes(self):
        result = evaluate("Escherichia coli", "rmlst", reads=5000,
                          total_reads=1_000_000, rpm_ntc=200.0)
        self.assertTrue(result.passed)
        self.assertAlmostEqual(result.rpm, 5000.0)

    def test_organism_at_its_background_fails(self):
        result = evaluate("Cutibacterium acnes", "rmlst", reads=300,
                          total_reads=1_000_000, rpm_ntc=200.0)
        self.assertFalse(result.passed)
        self.assertLess(result.rpm_r, 10.0)

    def test_very_clean_control_does_not_let_everything_through(self):
        """A near-zero denominator would make every ratio enormous.

        This is what the floor guards: without it, an organism at a trivial
        1 RPM against a control level of 0.01 would score RPM-r = 100 and be
        reported.
        """
        floored = evaluate("Rare organism", "rmlst", reads=1,
                           total_reads=1_000_000, rpm_ntc=0.01,
                           min_ntc_rpm=MIN_NTC_RPM)
        self.assertFalse(floored.passed)
        self.assertIn("floored", floored.reason)

        unfloored = evaluate("Rare organism", "rmlst", reads=1,
                             total_reads=1_000_000, rpm_ntc=0.01,
                             min_ntc_rpm=0.0)
        self.assertTrue(unfloored.passed)  # the failure mode being prevented

    def test_shallow_control_is_excluded_from_the_baseline(self):
        """A 79-read control makes one stray read worth 12,658 RPM.

        A floor cannot help - that value is far ABOVE it. Such a denominator
        would block genuine organisms, so the control is excluded instead and
        the pooled baseline is used.
        """
        rows = build_baselines({
            "NEG_shallow": {"Escherichia coli": (1.0, 79)},
            "NEG_good": {"Escherichia coli": (10.0, 1_000_000)},
        }, run_id="RUN1", arm="rmlst")
        self.assertEqual(len(rows), 1)
        self.assertAlmostEqual(float(rows[0]["rpm_ntc"]), 10.0)
        self.assertEqual(rows[0]["n_controls"], 1)
        # The control's own figures travel with the RPM.
        self.assertEqual(rows[0]["reads_ntc"], "10")
        self.assertEqual(rows[0]["total_reads_ntc"], 1_000_000)
        self.assertEqual(rows[0]["control_ids"], "NEG_good")
        self.assertEqual(rows[0]["mechanism"], "reagent")

    def test_all_controls_too_shallow_yields_no_baseline(self):
        rows = build_baselines({"NEG": {"E coli": (1.0, 79)}},
                               run_id="RUN1", arm="rmlst")
        self.assertEqual(rows, [])

    def test_organism_absent_from_controls_uses_the_absolute_floor(self):
        result = evaluate("Neisseria meningitidis", "rmlst", reads=100,
                          total_reads=1_000_000, rpm_ntc=None)
        self.assertTrue(result.passed)
        self.assertIsNone(result.rpm_r)
        self.assertEqual(result.band, "no_control_background")

    def test_baseline_takes_the_worst_control_not_the_average(self):
        """One heavily contaminated control is the case an average would hide."""
        rows = build_baselines({
            "NEG1": {"Escherichia coli": (10.0, 1_000_000)},
            "NEG2": {"Escherichia coli": (500.0, 1_000_000)},
        }, run_id="RUN1", arm="rmlst")
        self.assertAlmostEqual(float(rows[0]["rpm_ntc"]), 500.0)

    def test_run_specific_baseline_wins_over_the_pool(self):
        """Reagent lots change between runs, so the run's own control governs."""
        baselines = {("RUN1", "rmlst", "escherichia coli"): 900.0,
                     ("*", "rmlst", "escherichia coli"): 10.0}
        self.assertEqual(
            lookup_baseline(baselines, "RUN1", "rmlst", "Escherichia coli"), 900.0)
        self.assertEqual(
            lookup_baseline(baselines, "RUN9", "rmlst", "Escherichia coli"), 10.0)

    def test_banding_is_coarse_enough_to_absorb_run_to_run_noise(self):
        """est_counts move up to ~27% between runs; a band must not flip on that.

        Checked away from the band edges, which is the point of banding: a
        value sitting exactly on a boundary can always be pushed across it, so
        the guarantee is that ordinary noise does not move a *typical* result.
        """
        for reads in (400, 4_000, 40_000):
            with self.subTest(reads=reads):
                low = evaluate("X", "rmlst", reads, 1_000_000, 200.0)
                high = evaluate("X", "rmlst", reads * 1.27, 1_000_000, 200.0)
                self.assertEqual(low.band, high.band)


class GateWiringTests(unittest.TestCase):
    """The two result types name their taxon field differently.

    ``RmlstStrengthResult.taxon`` vs ``BlastStrengthResult.top_taxon``. Code
    that gates both must read both, or the 16S arm crashes the moment a
    baseline is supplied - which is invisible until someone actually turns the
    contamination gate on.
    """

    def test_both_result_types_expose_a_readable_taxon(self):
        from bridge.reconcile.reconcile_cli import _gate_results

        rmlst = rmlst_results([rmlst_row("Klebsiella pneumoniae")], CONFIG)
        s16 = [rank_marker16s(
            s16_row("Klebsiella pneumoniae", "Klebsiella", "pneumoniae",
                    "species_strong"), CONFIG)]
        baselines = {("R", "rmlst", "klebsiella pneumoniae"): 10.0,
                     ("R", "marker16s", "klebsiella pneumoniae"): 10.0}
        reads = {"klebsiella pneumoniae": 50_000.0}

        for arm, results in (("rmlst", rmlst), ("marker16s", s16)):
            with self.subTest(arm=arm):
                gated, detail, flags = _gate_results(
                    results, arm, reads, baselines, "R", 1_000_000, CONFIG)
                self.assertEqual(len(gated), 1)
                self.assertIn("Klebsiella pneumoniae", detail)
                # detail carries the JSON-shaped verdict, not the dataclass.
                self.assertTrue(detail["Klebsiella pneumoniae"]["passed"])

    def test_failing_the_gate_annotates_rather_than_silences(self):
        """A below-background organism STAYS in the call set.

        The RPM band IS the verdict a reader acts on, so an organism judged to be kitome is more useful
        visible and marked than absent - otherwise "not in the table" means both
        "nothing there" and "found, and rejected".
        """
        from bridge.reconcile.reconcile_cli import _gate_results

        results = rmlst_results([rmlst_row("Cutibacterium acnes")], CONFIG)
        gated, detail, _ = _gate_results(
            results, "rmlst", {"cutibacterium acnes": 100.0},
            {("R", "rmlst", "cutibacterium acnes"): 5000.0},
            "R", 1_000_000, CONFIG)
        self.assertEqual(len(gated), 1)
        self.assertIs(gated[0].strength, RmlstStrength.SPECIES_SUPPORTED)
        self.assertIn("rpm_below_background", gated[0].flags)
        # The measurement travels to the report so a reviewer sees why.
        self.assertIn("RPM-r", gated[0].reason)
        self.assertFalse(detail["Cutibacterium acnes"]["passed"])
        self.assertEqual(detail["Cutibacterium acnes"]["band"], "background")

    def test_a_demoted_name_is_gated_on_its_own_reads(self):
        """The genus fallback in _gate_results must be able to hit.

        Were `rmlst_reads` keyed by species binomial alone, a lookup by genus
        could never hit, and the row would reach the report flagged
        `rpm_gate_no_read_count`, never checked against background at all.
        """
        from bridge.reconcile.reconcile_cli import _gate_results, _rmlst_reads_by_organism

        rows = [rmlst_row("Streptococcus mitis", unamb_counts=6000.0,
                          amb_counts=0.0),
                rmlst_row("Streptococcus oralis", reported=False,
                          alternative_to="Streptococcus mitis", unamb=1,
                          unamb_counts=400.0, amb_counts=0.0)]
        reads = _rmlst_reads_by_organism(rows)
        # Species keys win; the genus total is only a fallback.
        self.assertEqual(reads["streptococcus mitis"], 6000.0)
        self.assertEqual(reads["streptococcus oralis"], 400.0)
        self.assertEqual(reads["streptococcus"], 6400.0)

        gated, detail, _ = _gate_results(
            rmlst_results(rows, CONFIG), "rmlst", reads,
            {("R", "rmlst", "streptococcus mitis"): 1.0}, "R", 1_000_000, CONFIG)
        for result in gated:
            self.assertNotIn("rpm_gate_no_read_count", result.flags)
        self.assertIn("Streptococcus oralis", detail)

    def test_a_genus_only_name_falls_back_to_the_genus_total(self):
        """A name with no species key of its own still gets measured.

        A genus total is >= any species within it, so as a numerator it makes
        the gate stricter, never looser.
        """
        from bridge.reconcile.reconcile_cli import _rmlst_reads_by_organism

        reads = _rmlst_reads_by_organism([
            rmlst_row("Staphylococcus aureus", unamb_counts=900.0, amb_counts=100.0),
            rmlst_row("Staphylococcus epidermidis", unamb_counts=400.0, amb_counts=0.0),
        ])
        self.assertEqual(reads["staphylococcus"], 1400.0)

    def test_an_unestimated_abundance_contributes_no_genus_total(self):
        """None means "not estimated"; it must not become a 0 in a sum."""
        from bridge.reconcile.reconcile_cli import _rmlst_reads_by_organism

        reads = _rmlst_reads_by_organism([
            rmlst_row("Staphylococcus aureus", unamb_counts=None, amb_counts=None),
        ])
        self.assertIsNone(reads["staphylococcus aureus"])
        self.assertNotIn("staphylococcus", reads)


class SampleOutcomeTests(unittest.TestCase):
    """Only reported calls decide the outcome.

    An `rmlst_alternative` row is a second NAME for an organism already in the
    list, and a `16s_only` row is an observation the governing arm did not make.
    Letting either reach the count would report one mitis-group isolate as a
    mixed infection.
    """

    def _rmlst(self, name, loci=30):
        from bridge.reconcile.taxonomy import parse_taxon
        return RmlstStrengthResult(
            strength=RmlstStrength.SPECIES_SUPPORTED,
            taxon=parse_taxon(name, CONFIG), support=90.0, n_loci=loci)

    def _s16(self, genus, species="", strength=BlastStrength.SPECIES_STRONG):
        label = f"{genus} {species}" if species else genus
        return rank_marker16s(
            s16_row(label, genus, species, strength.value), CONFIG)

    def test_nothing_at_all_is_a_no_call(self):
        calls = reconcile_sample([], [], CONFIG)
        self.assertEqual(calls, [])
        self.assertEqual(sample_outcome(calls), "no_call")

    def test_one_organism_is_reported(self):
        calls = reconcile_sample([self._rmlst("Klebsiella pneumoniae")], [], CONFIG)
        self.assertEqual(sample_outcome(calls), "reported")

    def test_two_organisms_are_polymicrobial(self):
        calls = reconcile_sample(
            [self._rmlst("Klebsiella pneumoniae"),
             self._rmlst("Staphylococcus aureus")], [], CONFIG)
        self.assertEqual(sample_outcome(calls), "polymicrobial")
        for call in calls:
            self.assertIn("polymicrobial", call.flags)

    def test_a_16s_only_row_does_not_make_a_sample_polymicrobial(self):
        """One organism plus reagent Halomonas is not a mixed infection."""
        calls = reconcile_sample(
            [self._rmlst("Legionella gratiana")],
            [self._s16("Legionella", strength=BlastStrength.GENUS_ONLY),
             self._s16("Halomonas", "elongata")], CONFIG)
        self.assertEqual(sample_outcome(calls), "reported")
        by_basis = {c.basis: c for c in calls}
        self.assertEqual(sorted(by_basis), ["16s_only", "rmlst"])
        self.assertEqual(by_basis["16s_only"].organism, "Halomonas elongata")
        self.assertNotIn("polymicrobial", by_basis["rmlst"].flags)

    def test_a_demoted_second_name_does_not_make_a_sample_polymicrobial(self):
        calls = reconcile_sample(rmlst_results([
            rmlst_row("Streptococcus mitis"),
            rmlst_row("Streptococcus oralis", reported=False,
                      alternative_to="Streptococcus mitis", unamb=1),
        ], CONFIG), [], CONFIG)
        self.assertEqual(sample_outcome(calls), "reported")
        self.assertEqual(sorted(c.basis for c in calls),
                         ["rmlst", "rmlst_alternative"])

    def test_reported_calls_sort_before_informational_rows(self):
        calls = reconcile_sample(
            [self._rmlst("Legionella gratiana")],
            [self._s16("Aeromonas", "hydrophila")], CONFIG)
        self.assertTrue(calls[0].is_call)
        self.assertFalse(calls[-1].is_call)
        self.assertEqual(CALL_BASES, {"rmlst", "16s_fallback"})


class MultiFallbackTests(unittest.TestCase):
    """Several report-eligible 16S fallback calls in one sample are held, not removed.

    E.g. rMLST resolved nothing, and 16S named Moraxella, Paracoccus and
    Rivihabitans, none of which the control held.
    """

    def _s16(self, genus, species=""):
        label = f"{genus} {species}" if species else genus
        return rank_marker16s(
            s16_row(label, genus, species, BlastStrength.SPECIES_STRONG.value), CONFIG)

    @staticmethod
    def _background(calls, *organisms):
        """What the RPM gate leaves on a call it failed, before reconciliation."""
        for call in calls:
            if call.organism in organisms:
                call.flags.append("16s:rpm_below_background")

    def test_several_eligible_fallback_calls_are_all_flagged(self):
        calls = reconcile_sample([], [self._s16("Moraxella", "tetraodonis"),
                                      self._s16("Paracoccus", "yeei"),
                                      self._s16("Rivihabitans", "pingtungensis")], CONFIG)
        self.assertEqual(flag_multi_fallback(calls, 2), 3)
        self.assertEqual(len(calls), 3)
        for call in calls:
            self.assertEqual(call.basis, "16s_fallback")
            self.assertIn(MULTI_FALLBACK_FLAG, call.flags)
            self.assertIn("one of 3 report-eligible 16S fallback calls", call.notes)

    def test_background_calls_are_neither_counted_nor_flagged(self):
        """One eligible call beside two background ones stays report."""
        calls = reconcile_sample([], [self._s16("Streptococcus", "sanguinis"),
                                      self._s16("Shigella", "sonnei"),
                                      self._s16("Stutzerimonas", "stutzeri")], CONFIG)
        self._background(calls, "Shigella sonnei", "Stutzerimonas stutzeri")
        self.assertEqual(flag_multi_fallback(calls, 2), 1)
        self.assertFalse(any(MULTI_FALLBACK_FLAG in call.flags for call in calls))

    def test_only_the_eligible_calls_are_flagged(self):
        calls = reconcile_sample([], [self._s16("Moraxella", "tetraodonis"),
                                      self._s16("Paracoccus", "yeei"),
                                      self._s16("Shigella", "sonnei")], CONFIG)
        self._background(calls, "Shigella sonnei")
        self.assertEqual(flag_multi_fallback(calls, 2), 2)
        flagged = {call.organism for call in calls if MULTI_FALLBACK_FLAG in call.flags}
        self.assertEqual(flagged, {"Moraxella tetraodonis", "Paracoccus yeei"})

    def test_all_background_is_left_to_the_gate(self):
        """Three fallback calls, all background - already review."""
        calls = reconcile_sample([], [self._s16("Escherichia", "ruysiae"),
                                      self._s16("Salmonella", "enterica"),
                                      self._s16("Shigella", "sonnei")], CONFIG)
        self._background(calls, *(call.organism for call in calls))
        self.assertEqual(flag_multi_fallback(calls, 2), 0)
        self.assertFalse(any(MULTI_FALLBACK_FLAG in call.flags for call in calls))

    def test_a_single_fallback_call_is_left_alone(self):
        calls = reconcile_sample([], [self._s16("Moraxella", "tetraodonis")], CONFIG)
        flag_multi_fallback(calls, 2)
        self.assertNotIn(MULTI_FALLBACK_FLAG, calls[0].flags)

    def test_rmlst_governed_samples_are_never_flagged(self):
        """16s_only rows are already informational; only fallback calls count."""
        from bridge.reconcile.taxonomy import parse_taxon
        rmlst = RmlstStrengthResult(
            strength=RmlstStrength.SPECIES_SUPPORTED,
            taxon=parse_taxon("Legionella gratiana", CONFIG), support=90.0, n_loci=30)
        calls = reconcile_sample([rmlst], [self._s16("Halomonas", "elongata"),
                                           self._s16("Paracoccus", "yeei")], CONFIG)
        self.assertEqual(flag_multi_fallback(calls, 2), 0)
        self.assertFalse(any(MULTI_FALLBACK_FLAG in call.flags for call in calls))

    def test_threshold_zero_disables_the_rule(self):
        calls = reconcile_sample([], [self._s16("Moraxella", "tetraodonis"),
                                      self._s16("Paracoccus", "yeei")], CONFIG)
        flag_multi_fallback(calls, 0)
        self.assertFalse(any(MULTI_FALLBACK_FLAG in call.flags for call in calls))


class MiniYamlTests(unittest.TestCase):
    """The built-in YAML reader used where PyYAML is absent.

    Its whole justification is that it must agree with PyYAML on this project's
    config, and refuse rather than guess on anything else - these are clinical
    thresholds.
    """

    def test_matches_pyyaml_on_the_real_config(self):
        try:
            import yaml
        except ModuleNotFoundError:
            self.skipTest("PyYAML not installed; nothing to compare against")
        from bridge.reconcile.config import _mini_yaml
        text = (Path(__file__).resolve().parents[1] / "config.yaml").read_text(
            encoding="utf-8")
        self.assertEqual(_mini_yaml(text), yaml.safe_load(text))

    def test_scalar_types_are_preserved(self):
        from bridge.reconcile.config import _mini_yaml
        parsed = _mini_yaml(
            "a:\n"
            "  i: 25\n"
            "  f: 0.05\n"
            "  t: true\n"
            "  f2: false\n"
            "  n: null\n"
            "  s: /data/x.fna\n"
            '  q: "~"\n'
            "  empty: {}\n"
        )
        self.assertEqual(parsed["a"], {
            "i": 25, "f": 0.05, "t": True, "f2": False, "n": None,
            "s": "/data/x.fna", "q": "~", "empty": {},
        })
        self.assertIsInstance(parsed["a"]["i"], int)
        self.assertIsInstance(parsed["a"]["f"], float)

    def test_comments_and_blank_lines_ignored(self):
        from bridge.reconcile.config import _mini_yaml
        self.assertEqual(
            _mini_yaml("# lead\n\nk:\n  v: 1  # trailing\n"), {"k": {"v": 1}})

    def test_hash_inside_a_quoted_value_is_not_a_comment(self):
        from bridge.reconcile.config import _mini_yaml
        self.assertEqual(_mini_yaml('k: "a # b"\n'), {"k": "a # b"})

    def test_lists_are_refused_not_guessed(self):
        from bridge.reconcile.config import ConfigError, _mini_yaml
        with self.assertRaises(ConfigError) as caught:
            _mini_yaml("k:\n  - one\n  - two\n")
        self.assertIn("PyYAML", str(caught.exception))

    def test_unsupported_constructs_are_refused(self):
        from bridge.reconcile.config import ConfigError, _mini_yaml
        for text in ("k: [1, 2]\n", "k: &anchor v\n", "k: |\n  block\n"):
            with self.subTest(text=text):
                with self.assertRaises(ConfigError):
                    _mini_yaml(text)

    def test_config_loads_without_pyyaml(self):
        """Simulate an environment without PyYAML."""
        import bridge.reconcile.config as config_module
        original = config_module.yaml
        config_module.yaml = None
        try:
            cfg = config_module.load_config(
                Path(__file__).resolve().parents[1] / "config.yaml")
            self.assertEqual(cfg.require("rpm.rpm_r_min"), 10.0)
            self.assertEqual(cfg.require("marker16s.species_dominance"), 0.70)
        finally:
            config_module.yaml = original


class RmlstPrimaryTests(unittest.TestCase):
    """rMLST governs; 16S is the fallback.

    This replaces the 3x3 confidence matrix. The matrix treated the two arms as
    co-equal witnesses, which is what sent clean rMLST species calls to review
    whenever a single 16S gene had not independently named the same organism.
    """

    def _rmlst(self, name, loci=30):
        from bridge.reconcile.taxonomy import parse_taxon
        return RmlstStrengthResult(
            strength=RmlstStrength.SPECIES_SUPPORTED,
            taxon=parse_taxon(name, CONFIG), support=90.0, n_loci=loci)

    def _s16(self, taxon, genus, species="",
             strength=BlastStrength.SPECIES_STRONG, row_index=0):
        return rank_marker16s(
            s16_row(taxon, genus, species, strength.value), CONFIG,
            row_index=row_index)

    def test_every_rmlst_species_becomes_a_reported_call(self):
        calls = reconcile_sample(
            [self._rmlst("Staphylococcus aureus"),
             self._rmlst("Legionella gratiana")], [], CONFIG)
        self.assertEqual(sorted(c.organism for c in calls),
                         ["Legionella gratiana", "Staphylococcus aureus"])
        for call in calls:
            self.assertEqual(call.basis, "rmlst")
            self.assertEqual(call.reported_rank.value, "species")

    def test_a_16s_row_in_a_called_genus_attaches_rather_than_adding_a_row(self):
        calls = reconcile_sample(
            [self._rmlst("Legionella gratiana")],
            [self._s16("Legionella", "Legionella",
                       strength=BlastStrength.GENUS_ONLY)], CONFIG)
        self.assertEqual(len(calls), 1)
        self.assertEqual(calls[0].basis, "rmlst")
        self.assertEqual(calls[0].s16_top_hit, "Legionella")
        self.assertIsNotNone(calls[0].s16_row_index)

    def test_one_16s_genus_row_corroborates_every_species_in_that_genus(self):
        """16S saw the genus and could not split it.

        That supports all of the species rMLST resolved there equally - not one
        picked arbitrarily, and not none of them.
        """
        calls = reconcile_sample(
            [self._rmlst("Staphylococcus aureus"),
             self._rmlst("Staphylococcus epidermidis")],
            [self._s16("Staphylococcus", "Staphylococcus",
                       strength=BlastStrength.GENUS_ONLY)], CONFIG)
        self.assertEqual(len(calls), 2)
        for call in calls:
            with self.subTest(organism=call.organism):
                self.assertEqual(call.s16_top_hit, "Staphylococcus")
                self.assertNotIn("no_16s_corroboration", call.flags)

    def test_a_16s_genus_rmlst_did_not_call_becomes_an_informational_row(self):
        """Nothing vanishes: PANDAv7 covers a restricted set of target genera,
        so 16S can legitimately see an organism rMLST structurally cannot."""
        calls = reconcile_sample(
            [self._rmlst("Legionella gratiana")],
            [self._s16("Halomonas elongata", "Halomonas", "elongata")], CONFIG)
        by_basis = {c.basis: c for c in calls}
        self.assertEqual(sorted(by_basis), ["16s_only", "rmlst"])
        self.assertFalse(by_basis["16s_only"].is_call)
        self.assertIn("target-genus", by_basis["16s_only"].notes)

    def test_a_silent_rmlst_arm_hands_the_call_set_to_16s(self):
        """rMLST resolved nothing at all."""
        calls = reconcile_sample(
            [], [self._s16("Streptococcus dysgalactiae", "Streptococcus",
                           "dysgalactiae")], CONFIG)
        self.assertEqual(len(calls), 1)
        self.assertEqual(calls[0].basis, "16s_fallback")
        self.assertTrue(calls[0].is_call)
        self.assertEqual(calls[0].organism, "Streptococcus dysgalactiae")
        self.assertEqual(calls[0].reported_rank.value, "species")

    def test_a_16s_only_sample_still_reports_at_genus_where_16s_declined(self):
        calls = reconcile_sample(
            [], [self._s16("Streptococcus", "Streptococcus",
                           strength=BlastStrength.GENUS_ONLY)], CONFIG)
        self.assertEqual(calls[0].reported_rank.value, "genus")
        self.assertEqual(calls[0].organism, "Streptococcus")

    def test_both_arms_silent_is_no_call(self):
        self.assertEqual(reconcile_sample([], [], CONFIG), [])

    def test_a_compound_16s_genus_joins_the_plain_rmlst_genus(self):
        """SILVA-style 'Escherichia-Shigella' must not stand up a second row
        for an organism rMLST already reported."""
        calls = reconcile_sample(
            [self._rmlst("Escherichia coli")],
            [self._s16("Escherichia-Shigella", "Escherichia-Shigella",
                       strength=BlastStrength.GENUS_ONLY)], CONFIG)
        self.assertEqual(len(calls), 1)
        self.assertEqual(calls[0].basis, "rmlst")
        self.assertEqual(calls[0].organism, "Escherichia coli")

    def test_an_rmlst_call_with_no_16s_partner_is_flagged_not_demoted(self):
        calls = reconcile_sample([self._rmlst("Klebsiella pneumoniae")], [],
                                 CONFIG)
        self.assertIn("no_16s_corroboration", calls[0].flags)
        self.assertTrue(calls[0].is_call)
        self.assertEqual(calls[0].basis, "rmlst")

    def test_a_16s_species_disagreeing_with_rmlst_no_longer_overrules_it(self):
        """Closes the defect where one gene overruled ~53 loci.

        The matrix demoted a species call to genus level and flagged
        `species_disagreement` whenever the arms named different species in one
        genus. rMLST resolves those; 16S by construction cannot.
        """
        calls = reconcile_sample(
            [self._rmlst("Streptococcus mitis")],
            [self._s16("Streptococcus oralis", "Streptococcus", "oralis")],
            CONFIG)
        self.assertEqual(len(calls), 1)
        self.assertEqual(calls[0].organism, "Streptococcus mitis")
        self.assertEqual(calls[0].reported_rank.value, "species")
        self.assertNotIn("species_disagreement", calls[0].flags)


class UnestimatedAbundanceTests(unittest.TestCase):
    """A blank abundance means "not estimated". It must never read as zero reads.

    The rMLST arm leaves ``pi`` and the ``est_counts`` cells blank for a species
    called from its allele profile whose consensus groups never reached EM. Read
    as 0.0, 0 reads would fail the RPM-r floor and silence a real organism.

    The header here deliberately has no coverage columns. Such a file must
    still parse, with those fields coming back blank rather than as zeros.
    """

    HEADER = ("sample_id,species,reported,call_method,exact_profile_loci,"
              "exact_profile_unambiguous_loci,pi,allele_evidence_unambiguous_loci,"
              "allele_evidence_ambiguous_loci,em_evidence_unambiguous_loci,"
              "em_evidence_ambiguous_loci,unambiguous_est_counts,"
              "ambiguous_est_counts,alternative_to\n")

    def _load(self, row):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "species_calls.csv"
            path.write_text(self.HEADER + row, encoding="utf-8")
            return load_rmlst_calls(path)

    def test_blank_abundance_cells_are_none_not_zero(self):
        rows = self._load(
            "S1,Streptococcus dysgalactiae,True,exact_allele_profile,"
            "31,13,,13,0,0,0,,,\n")
        self.assertEqual(len(rows), 1)
        self.assertIsNone(rows[0].pi)
        self.assertIsNone(rows[0].total_est_counts)

    def test_a_measured_zero_is_still_zero(self):
        """The distinction only means anything if a real 0 survives as 0."""
        rows = self._load(
            "s,Cutibacterium acnes,True,exact_allele_profile,8,8,0.0,8,0,0,0,0.0,0.0,\n")
        self.assertEqual(rows[0].pi, 0.0)
        self.assertEqual(rows[0].total_est_counts, 0.0)

    def test_a_pre_3_0_0_file_parses_with_blank_coverage(self):
        """An older species_calls.csv has no coverage columns at all.

        They must come back as None, not 0.0 - the same rule, for the same
        reason: 0.0 asserts a measurement that was never made.
        """
        rows = self._load(
            "s,Klebsiella pneumoniae,True,exact_allele_profile,40,40,0.5,40,3,40,5,"
            "900.0,100.0,\n")
        row = rows[0]
        for field in ("mean_callable_breadth", "min_callable_breadth",
                      "median_projected_depth", "min_projected_depth",
                      "total_callable_bases", "total_allele_length",
                      "supporting_fragments"):
            with self.subTest(field=field):
                self.assertIsNone(getattr(row, field))
        # ...while the columns that WERE present still parse.
        self.assertEqual(row.em_evidence_unambiguous_loci, 40)
        self.assertEqual(row.total_est_counts, 1000.0)

    def test_unknown_abundance_is_not_silenced_by_the_rpm_gate(self):
        """Ungated and flagged, rather than destroyed."""
        from bridge.reconcile.reconcile_cli import _gate_results
        from bridge.reconcile.taxonomy import parse_taxon
        result = RmlstStrengthResult(
            strength=RmlstStrength.SPECIES_SUPPORTED,
            taxon=parse_taxon("Streptococcus dysgalactiae", CONFIG),
            support=None, n_loci=13, reason="13 unambiguous loci")
        gated, detail, flags = _gate_results(
            [result], "rmlst", {"streptococcus dysgalactiae": None},
            {}, "RUN", 1_000_000, CONFIG)
        self.assertIs(gated[0].strength, RmlstStrength.SPECIES_SUPPORTED)
        self.assertIn("rpm_gate_no_read_count", gated[0].flags)
        self.assertNotIn("rpm_below_background", gated[0].flags)

    def test_a_measured_zero_is_still_marked_below_background(self):
        """A real 0 still fails the gate - it is now flagged, not deleted."""
        from bridge.reconcile.reconcile_cli import _gate_results
        from bridge.reconcile.taxonomy import parse_taxon
        result = RmlstStrengthResult(
            strength=RmlstStrength.SPECIES_SUPPORTED,
            taxon=parse_taxon("Cutibacterium acnes", CONFIG),
            support=0.0, n_loci=8, reason="8 unambiguous loci")
        gated, detail, flags = _gate_results(
            [result], "rmlst", {"cutibacterium acnes": 0.0},
            {}, "RUN", 1_000_000, CONFIG)
        self.assertIn("rpm_below_background", gated[0].flags)
        self.assertFalse(detail["Cutibacterium acnes"]["passed"])
        # ...and is distinguishable from the unestimated case above, which is
        # flagged rpm_gate_no_read_count and never reaches a verdict at all.
        self.assertNotIn("rpm_gate_no_read_count", gated[0].flags)


class PerOrganismEvidenceTests(unittest.TestCase):
    """Each row carries ITS organism's numbers, not the sample's.

    In a polymicrobial sample the read count is precisely the figure that has
    to differ between organisms for the result to be worth anything. Evidence is fetched by the arm ROW INDEX rather than by name:
    matching on the name silently misses whenever the reconciler labels an
    organism differently from the arm.
    """

    def _evidence(self, rmlst_rows, s16_rows=()):
        from bridge.reconcile.adapters.marker16s import marker16s_results
        from bridge.reconcile.reconcile_cli import _evidence_for

        s16_rows = list(s16_rows)
        calls = reconcile_sample(rmlst_results(rmlst_rows, CONFIG),
                                 marker16s_results(s16_rows, CONFIG), CONFIG)
        return {c.organism: (c, _evidence_for(c, rmlst_rows, s16_rows))
                for c in calls}

    def test_two_organisms_of_one_sample_carry_different_numbers(self):
        found = self._evidence([
            rmlst_row("Staphylococcus aureus", unamb=40, pi=0.8,
                      unamb_counts=8000.0, breadth=0.96, depth=310.0),
            rmlst_row("Staphylococcus epidermidis", unamb=22, pi=0.2,
                      unamb_counts=2500.0, breadth=0.81, depth=95.0),
        ])
        aureus = found["Staphylococcus aureus"][1]
        epi = found["Staphylococcus epidermidis"][1]
        self.assertEqual(aureus["unambiguous_est_counts"], 8000.0)
        self.assertEqual(epi["unambiguous_est_counts"], 2500.0)
        self.assertEqual(aureus["allele_evidence_unambiguous_loci"], 40)
        self.assertEqual(epi["allele_evidence_unambiguous_loci"], 22)
        self.assertAlmostEqual(aureus["mean_callable_breadth"], 0.96)
        self.assertAlmostEqual(epi["mean_callable_breadth"], 0.81)
        self.assertAlmostEqual(aureus["pi"], 0.8)
        self.assertAlmostEqual(epi["pi"], 0.2)

    def test_a_demoted_name_finds_its_own_row_despite_the_label(self):
        """The case a name lookup gets wrong.

        C. resistens arrives labelled `Corynebacterium` under the old
        species-stripping behaviour; either way the row index is exact.
        """
        rows = [rmlst_row("Streptococcus mitis", unamb=40, unamb_counts=6000.0),
                rmlst_row("Streptococcus oralis", reported=False,
                          alternative_to="Streptococcus mitis", unamb=1,
                          unamb_counts=400.0)]
        found = self._evidence(rows)
        oralis = [e for c, e in found.values() if c.basis == "rmlst_alternative"]
        self.assertEqual(len(oralis), 1)
        self.assertEqual(oralis[0]["unambiguous_est_counts"], 400.0)
        self.assertEqual(oralis[0]["allele_evidence_unambiguous_loci"], 1)
        self.assertEqual(oralis[0]["alternative_to"], "Streptococcus mitis")

    def test_a_losing_sibling_cannot_contaminate_the_reported_organism(self):
        """C. auriscanis had no abundance; C. resistens holds 95.9%.

        Merging across a genus would let the sibling's missing data demote the
        organism that was measured perfectly well. Per-organism evidence is
        the structural fix - there is no shared bucket to pollute.
        """
        found = self._evidence([
            rmlst_row("Corynebacterium resistens", unamb=49,
                      unamb_counts=40000.0, amb_counts=0.0),
            rmlst_row("Corynebacterium auriscanis", unamb=7,
                      unamb_counts=None, amb_counts=None),
        ])
        self.assertEqual(sorted(found), ["Corynebacterium auriscanis",
                                         "Corynebacterium resistens"])
        resistens = found["Corynebacterium resistens"][1]
        auriscanis = found["Corynebacterium auriscanis"][1]
        self.assertEqual(resistens["unambiguous_est_counts"], 40000.0)
        self.assertEqual(resistens["allele_evidence_unambiguous_loci"], 49)
        # The sibling's blank stays the sibling's, and stays blank.
        self.assertIsNone(auriscanis["unambiguous_est_counts"])

    def test_the_16s_numbers_reach_an_rmlst_governed_row(self):
        """The arm still runs on every sample; it just no longer votes."""
        s16 = [s16_row("Staphylococcus", "Staphylococcus", "", "genus_only",
                       reads=1600.0)]
        found = self._evidence([rmlst_row("Staphylococcus aureus")], s16)
        evidence = found["Staphylococcus aureus"][1]
        self.assertEqual(evidence["16s_reads_assigned"], 1600.0)
        self.assertEqual(evidence["16s_raw_reads"], 6400.0)

    def test_an_organism_with_no_rmlst_row_gets_no_rmlst_numbers(self):
        """A 16s_only row must not borrow another organism's evidence."""
        s16 = [s16_row("Halomonas elongata", "Halomonas", "elongata",
                       "species_strong", reads=78.0)]
        found = self._evidence([rmlst_row("Legionella gratiana")], s16)
        halomonas = found["Halomonas elongata"][1]
        self.assertNotIn("unambiguous_est_counts", halomonas)
        self.assertNotIn("mean_callable_breadth", halomonas)
        self.assertEqual(halomonas["16s_reads_assigned"], 78.0)


class SummaryCsvTests(unittest.TestCase):
    """The reported result as a flat table."""

    CALL = {
        "organism": "Legionella gratiana", "genus": "Legionella",
        "reported_rank": "species", "basis": "rmlst",
        "rmlst_taxon": "Legionella gratiana",
        "call_method": "exact_allele_profile", "pi": 0.941,
        "exact_profile_loci": 45, "exact_profile_unambiguous_loci": 40,
        "allele_evidence_unambiguous_loci": 40,
        "allele_evidence_ambiguous_loci": 3,
        "em_evidence_unambiguous_loci": 40, "em_evidence_ambiguous_loci": 5,
        "unambiguous_est_counts": 5000.0, "ambiguous_est_counts": 100.0,
        "alternative_to": "",
        "mean_callable_breadth": 0.94, "min_callable_breadth": 0.61,
        "median_projected_depth": 212.0, "min_projected_depth": 18.0,
        "total_callable_bases": 21847.0, "total_allele_length": 23100.0,
        "supporting_fragments": 4192.0,
        "16s_taxon": "Legionella", "16s_margin": 0.21,
        "16s_abundance_frac": 0.42, "16s_reads_assigned": 300.0,
        "16s_raw_reads": 1200.0,
        "flags": ["x"], "notes": "n",
    }

    RECORD = {
        "sample": "S1", "run_id": "R1", "outcome": "reported",
        "rmlst_reads": 5100.0, "reads_16s": 1200.0,
        "total_raw_reads": 2_000_000,
        "calls": [CALL],
        "flags": [],
        "rpm": {"rmlst": {"Legionella gratiana":
                          {"rpm": 1000.0, "rpm_r": 12.5, "band": "moderate",
                           "passed": True, "reason": "ok",
                           "rule": "rpm_r", "mechanism": "reagent", "abs_r": 30.0}},
                "marker16s": {}},
    }

    def test_the_gate_columns_carry_both_comparators(self):
        from bridge.reconcile.summary_csv import summary_rows
        row = summary_rows(self.RECORD)[0]
        self.assertEqual(float(row["rpm_r_rmlst"]), 12.5)
        self.assertEqual(float(row["abs_r_rmlst"]), 30.0)
        self.assertEqual(row["abs_r_16s"], "")
        self.assertEqual(row["bg_mechanism"], "reagent")

    def test_one_row_per_organism(self):
        from bridge.reconcile.summary_csv import summary_rows
        rows = summary_rows(self.RECORD)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["organism"], "Legionella gratiana")
        self.assertEqual(rows[0]["basis"], "rmlst")

    def test_the_removed_columns_are_really_gone(self):
        from bridge.reconcile.summary_csv import SUMMARY_COLUMNS, summary_rows
        row = summary_rows(self.RECORD)[0]
        for column in ("confidence", "gate_passed", "rmlst_reads",
                       "reads_16s", "total_raw_reads", "rmlst_abundance_pct",
                       "rmlst_loci", "sample_outcome"):
            with self.subTest(column=column):
                self.assertNotIn(column, SUMMARY_COLUMNS)
                self.assertNotIn(column, row)

    def test_the_16s_columns_are_spelled_16s_not_s16(self):
        """The arm is called 16S everywhere else in the output, incl. rpm_16s."""
        from bridge.reconcile.summary_csv import SUMMARY_COLUMNS
        self.assertFalse([c for c in SUMMARY_COLUMNS if c.startswith("s16")])
        for column in ("16s_taxon", "16s_margin", "16s_abundance_frac",
                       "16s_reads_assigned", "16s_raw_reads"):
            self.assertIn(column, SUMMARY_COLUMNS)

    def test_per_organism_evidence_reaches_the_table(self):
        from bridge.reconcile.summary_csv import summary_rows
        row = summary_rows(self.RECORD)[0]
        self.assertEqual(row["rmlst_taxon"], "Legionella gratiana")
        self.assertEqual(row["call_method"], "exact_allele_profile")
        self.assertEqual(float(row["pi"]), 0.941)
        self.assertEqual(row["allele_evidence_unambiguous_loci"], "40")
        self.assertEqual(row["em_evidence_ambiguous_loci"], "5")
        self.assertEqual(float(row["unambiguous_est_counts"]), 5000.0)
        self.assertEqual(float(row["mean_callable_breadth"]), 0.94)
        self.assertEqual(float(row["median_projected_depth"]), 212.0)
        self.assertEqual(float(row["supporting_fragments"]), 4192.0)
        self.assertEqual(float(row["16s_raw_reads"]), 1200.0)
        self.assertEqual(row["16s_margin"], "0.21")
        self.assertEqual(row["rpm_band"], "moderate")
        self.assertEqual(row["rpm_r_rmlst"], "12.5")
        # No single "rpm" column: an organism has one per arm, and they
        # are not the same measurement.
        self.assertNotIn("rpm", row)

    def test_two_organisms_of_one_sample_get_different_numbers(self):
        """Sample-level totals are not repeated against every organism.
        """
        from bridge.reconcile.summary_csv import summary_rows
        rows = summary_rows({**self.RECORD, "outcome": "polymicrobial", "calls": [
            self.CALL,
            {**self.CALL, "organism": "Staphylococcus aureus",
             "genus": "Staphylococcus", "rmlst_taxon": "Staphylococcus aureus",
             "pi": 0.059, "unambiguous_est_counts": 210.0,
             "allele_evidence_unambiguous_loci": 12,
             "mean_callable_breadth": 0.55},
        ]})
        self.assertEqual({r["unambiguous_est_counts"] for r in rows},
                         {"5000", "210"})
        self.assertEqual({r["allele_evidence_unambiguous_loci"] for r in rows},
                         {"40", "12"})
        self.assertEqual(len({r["mean_callable_breadth"] for r in rows}), 2)

    def test_a_blank_est_count_stays_blank_not_zero(self):
        """0.0 asserts "measured, and there was none"; blank says "not measured"."""
        from bridge.reconcile.summary_csv import summary_rows
        row = summary_rows({**self.RECORD, "calls": [{
            **self.CALL, "pi": None, "unambiguous_est_counts": None,
            "ambiguous_est_counts": None, "mean_callable_breadth": None,
            "median_projected_depth": None,
        }]})[0]
        for column in ("pi", "unambiguous_est_counts", "ambiguous_est_counts",
                       "mean_callable_breadth", "median_projected_depth"):
            with self.subTest(column=column):
                self.assertEqual(row[column], "")
        # ...while the locus counts beside them survive.
        self.assertEqual(row["allele_evidence_unambiguous_loci"], "40")

    def test_a_sample_with_no_calls_still_gets_a_row(self):
        """'Not in the table' must not be ambiguous between nothing-found and never-ran."""
        from bridge.reconcile.summary_csv import summary_rows
        rows = summary_rows({**self.RECORD, "calls": [], "outcome": "no_call"})
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["status"], "no_call")
        self.assertEqual(rows[0]["reported_rank"], "none")
        self.assertEqual(rows[0]["sample"], "S1")

    def test_rpm_falls_back_from_species_to_genus(self):
        from bridge.reconcile.summary_csv import summary_rows
        record = {**self.RECORD, "rpm": {
            "rmlst": {}, "marker16s": {"Legionella": {
                "rpm": 50.0, "rpm_r": None, "band": "no_control_background",
                "passed": True, "reason": "absent"}}}}
        self.assertEqual(summary_rows(record)[0]["rpm_band"],
                         "no_control_background")

    def test_a_trace_abundance_keeps_its_precision(self):
        """A pi of 9.58e-05 must not round to zero in the table."""
        from bridge.reconcile.summary_csv import summary_rows
        row = summary_rows({**self.RECORD,
                            "calls": [{**self.CALL, "pi": 0.0000958}]})[0]
        self.assertEqual(float(row["pi"]), 0.0000958)
        self.assertNotEqual(float(row["pi"]), 0.0)

    def test_status_is_report_for_a_call_and_review_for_an_observation(self):
        """A one-word reading of `basis`, for filtering a run to its findings."""
        from bridge.reconcile.summary_csv import summary_rows
        rows = summary_rows({**self.RECORD, "calls": [
            self.CALL,
            {**self.CALL, "organism": "Klebsiella pneumoniae",
             "genus": "Klebsiella", "basis": "16s_fallback"},
            {**self.CALL, "organism": "Streptococcus oralis",
             "genus": "Streptococcus", "basis": "rmlst_alternative"},
            {**self.CALL, "organism": "Halomonas elongata",
             "genus": "Halomonas", "basis": "16s_only"},
        ]})
        self.assertEqual({r["basis"]: r["status"] for r in rows}, {
            "rmlst": "report", "16s_fallback": "report",
            "rmlst_alternative": "review", "16s_only": "review",
        })

    def test_a_below_background_call_is_review(self):
        """The RPM half of the rule. It can only fire where a control exists."""
        from bridge.reconcile.summary_csv import summary_rows
        row = summary_rows({**self.RECORD, "calls": [
            {**self.CALL, "flags": ["rmlst:rpm_below_background"]}]})[0]
        self.assertEqual(row["basis"], "rmlst")
        self.assertEqual(row["status"], "review")

    def test_several_16s_fallback_calls_are_review_override(self):
        """Held separately, not dropped, and not confused with background."""
        from bridge.reconcile.summary_csv import summary_rows
        flagged = {**self.CALL, "basis": "16s_fallback",
                   "flags": ["multiple_16s_fallback"]}
        rows = summary_rows({**self.RECORD, "calls": [
            {**flagged, "organism": "Paracoccus", "genus": "Paracoccus"},
            {**flagged, "organism": "Moraxella tetraodonis", "genus": "Moraxella"},
            {**flagged, "organism": "Shigella sonnei", "genus": "Shigella",
             "flags": ["16s:rpm_below_background"]},
        ]})
        self.assertEqual([r["status"] for r in rows],
                         ["review_override", "review_override", "review"])
        self.assertEqual(["multiple_16s_fallback" in r["flags"] for r in rows],
                         [True, True, False])

    def test_status_never_contradicts_basis(self):
        """It is derived, not decided - so it cannot drift from the columns."""
        from bridge.reconcile.summary_csv import summary_rows
        for row in summary_rows(self.RECORD):
            self.assertEqual(row["status"], "report")
            self.assertIn(row["basis"], {"rmlst", "16s_fallback"})

    def test_an_informational_row_is_labelled_by_its_basis(self):
        from bridge.reconcile.summary_csv import summary_rows
        rows = summary_rows({**self.RECORD, "calls": [
            self.CALL,
            {"organism": "Halomonas elongata", "genus": "Halomonas",
             "reported_rank": "species", "basis": "16s_only",
             "16s_taxon": "Halomonas elongata", "16s_raw_reads": 310.0,
             "flags": [], "notes": "seen by the 16S arm only"},
        ]})
        informational = [r for r in rows if r["basis"] == "16s_only"][0]
        self.assertEqual(informational["organism"], "Halomonas elongata")
        # It borrows nothing from the rMLST row beside it.
        self.assertEqual(informational["pi"], "")
        self.assertEqual(informational["allele_evidence_unambiguous_loci"], "")
        self.assertEqual(float(informational["16s_raw_reads"]), 310.0)


class GovernorAndKeyingTests(unittest.TestCase):
    """The summary table must not hide the stricter of two disagreeing arms.

    Both cases are ones where the CSV and the JSON could tell different
    stories.
    """

    BASE = {
        "sample": "S1", "run_id": "RUN1",
        "outcome": "reported",
        "rmlst_reads": 43150.0, "reads_16s": 19650.0,
        "total_raw_reads": 3355879, "flags": [],
    }

    def test_the_stricter_arm_governs_the_headline_columns(self):
        """C. resistens: rMLST said 'no background', 16S said RPM-r 10.8."""
        from bridge.reconcile.summary_csv import summary_rows
        row = summary_rows({
            **self.BASE,
            "calls": [{"organism": "Corynebacterium resistens",
                       "genus": "Corynebacterium", "reported_rank": "species",
                       "basis": "rmlst",
                       "rmlst_taxon": "Corynebacterium resistens",
                       "allele_evidence_unambiguous_loci": 17, "pi": 0.959,
                       "16s_taxon": "Corynebacterium resistens",
                       "16s_margin": 0.961, "flags": [], "notes": ""}],
            "rpm": {
                "rmlst": {"Corynebacterium resistens": {
                    "rpm": 12103.6, "rpm_r": None,
                    "band": "no_control_background", "passed": True,
                    "reason": "absent from controls"}},
                "marker16s": {"Corynebacterium resistens": {
                    "rpm": 2782.95, "rpm_r": 10.824, "band": "moderate",
                    "passed": True, "reason": "RPM-r 10.8 >= 10.0"}},
            },
        })[0]
        # The marginal measured verdict governs, not the permissive unmeasured one.
        self.assertEqual(row["rpm_band"], "moderate")
        self.assertEqual(row["rpm_band_arm"], "marker16s")
        # ...and BOTH arms stay visible, so nothing is traded away for it.
        self.assertEqual(float(row["rpm_rmlst"]), 12100.0)
        self.assertEqual(row["rpm_r_rmlst"], "")
        self.assertEqual(float(row["rpm_r_16s"]), 10.82)

    def test_a_genus_level_call_finds_its_arms_species_keyed_rpm(self):
        """Clostridium: the gate filed it under 'Clostridium nigeriense'.

        The gate keys on the arm's finest resolution; the reconciler reports the
        rank it can defend. A name-only lookup shows a blank RPM for a figure
        that had in fact been computed.
        """
        from bridge.reconcile.summary_csv import summary_rows
        row = summary_rows({
            **self.BASE, "sample": "S2",
            "calls": [{"organism": "Clostridium", "genus": "Clostridium",
                       "reported_rank": "genus", "basis": "16s_only",
                       "rmlst_taxon": None,
                       "16s_taxon": "Clostridium nigeriense",
                       "16s_margin": 1.0, "flags": [], "notes": ""}],
            "rpm": {"rmlst": {}, "marker16s": {"Clostridium nigeriense": {
                "rpm": 52.4088, "rpm_r": None,
                "band": "no_control_background", "passed": True,
                "reason": "absent from controls"}}},
        })[0]
        self.assertEqual(float(row["rpm_16s"]), 52.41)
        self.assertEqual(row["rpm_band"], "no_control_background")

    def test_the_arm_that_measured_a_ratio_governs(self):
        """An absent measurement is not a clean result; it says nothing at all."""
        from bridge.reconcile.summary_csv import _governing_band
        arm, _ = _governing_band({
            "rmlst": {"band": "high", "rpm_r": 500.0},
            "marker16s": {"band": "no_control_background", "rpm_r": None},
        })
        self.assertEqual(arm, "rmlst")

    def test_a_failed_gate_governs_over_everything(self):
        from bridge.reconcile.summary_csv import _governing_band
        arm, _ = _governing_band({
            "rmlst": {"band": "high", "rpm_r": 500.0},
            "marker16s": {"band": "background", "rpm_r": 2.0}})
        self.assertEqual(arm, "marker16s")

    def test_two_measured_verdicts_resolve_to_the_worse_one(self):
        from bridge.reconcile.summary_csv import _governing_band
        arm, _ = _governing_band({
            "rmlst": {"band": "high", "rpm_r": 500.0},
            "marker16s": {"band": "moderate", "rpm_r": 10.8}})
        self.assertEqual(arm, "marker16s")


class CrossArmBackgroundTests(unittest.TestCase):
    """'No control background' is an arm-local fact stated as a general one."""

    def test_the_other_arms_control_evidence_is_named(self):
        from bridge.reconcile.rpm import arms_with_background
        baselines = {("R", "marker16s", "corynebacterium"): 257.11}
        self.assertEqual(
            arms_with_background(baselines, "R", "Corynebacterium resistens",
                                 "Corynebacterium"),
            ["marker16s"])

    def test_a_genus_call_matches_a_species_the_control_saw_within_it(self):
        """Reported as Staphylococcus, blank held S. aureus at 1,305 RPM.

        lookup_baseline will not make this match, and should not - a species
        level is not a safe denominator for a genus call. A flag is not a
        denominator, and here the fact is simply that the blank contained a
        staphylococcus.
        """
        from bridge.reconcile.rpm import arms_with_background
        self.assertEqual(
            arms_with_background({("R", "rmlst", "staphylococcus aureus"): 1305.29},
                                 "R", "Staphylococcus", "Staphylococcus"),
            ["rmlst"])

    def test_an_organism_in_no_control_reports_no_arms(self):
        from bridge.reconcile.rpm import arms_with_background
        self.assertEqual(
            arms_with_background({("R", "marker16s", "escherichia"): 10.0},
                                 "R", "Legionella gratiana", "Legionella"),
            [])

    def test_a_different_genus_is_not_matched_by_prefix(self):
        """Staphylococcus must not match Staphylococcaceae-adjacent names."""
        from bridge.reconcile.rpm import arms_with_background
        self.assertEqual(
            arms_with_background({("R", "rmlst", "staphylococcus aureus"): 100.0},
                                 "R", "Staphylothermus marinus", "Staphylothermus"),
            [])

    def test_the_flag_reaches_the_gated_result_and_the_row(self):
        from bridge.reconcile.reconcile_cli import _gate_results
        from bridge.reconcile.summary_csv import summary_rows

        results = rmlst_results([rmlst_row("Corynebacterium resistens")], CONFIG)
        gated, detail, _ = _gate_results(
            results, "rmlst", {"corynebacterium resistens": 40000.0},
            {("R", "marker16s", "corynebacterium"): 257.11},
            "R", 3_355_879, CONFIG)
        self.assertIn("rpm_control_background_in_marker16s", gated[0].flags)
        self.assertEqual(
            detail["Corynebacterium resistens"]["other_arm_background"],
            ["marker16s"])

        row = summary_rows({
            "sample": "s", "run_id": "R", "outcome": "reported", "flags": [],
            "calls": [{"organism": "Corynebacterium resistens",
                       "genus": "Corynebacterium", "reported_rank": "species",
                       "basis": "rmlst", "rmlst_taxon": "Corynebacterium resistens",
                       "allele_evidence_unambiguous_loci": 17, "pi": 0.959,
                       "16s_taxon": None, "16s_margin": None,
                       "flags": [], "notes": ""}],
            "rpm": {"rmlst": detail, "marker16s": {}},
        })[0]
        self.assertIn("rpm_control_background_in_marker16s", row["flags"])

    def test_a_normal_verdict_carries_no_cross_arm_key(self):
        """The record must not grow a field where there is nothing to say."""
        from bridge.reconcile.rpm import evaluate as _evaluate, rpm_detail
        detail = rpm_detail(_evaluate("X", "rmlst", 5000.0, 1_000_000, 200.0), [])
        self.assertNotIn("other_arm_background", detail)


class GateIsImmuneToCopyNumberTests(unittest.TestCase):
    """RPM-r compares one organism to itself, so copy number must cancel.

    It does cancel while both sides resolve to the same rank. It stops
    cancelling the moment they do not - and they routinely do not, because a
    negative control is shallow and resolves to genus where the sample resolves
    to species. That is the genus fallback working as designed, and it must not
    become a back door for the copy-number correction.
    """

    def _rows(self, path, rows):
        header = ("sample_id,taxon,genus,species,rank,strength,abundance_frac,"
                  "genus_frac,reads_assigned,raw_reads,copy_number,"
                  "within_genus_share,margin,n_references,competing_species,reason\n")
        path.write_text(header + "".join(rows), encoding="utf-8")
        return path

    def test_the_adapter_hands_the_gate_raw_reads(self):
        import tempfile
        from bridge.reconcile.adapters.marker16s import load_marker16s_calls

        path = self._rows(
            Path(tempfile.mkdtemp()) / "marker16s_calls.csv",
            ["S,Corynebacterium resistens,Corynebacterium,resistens,species,"
             "species_strong,0.88,0.88,13074.953,39224.860,3.000,1,1,2,x,y\n"])
        row = load_marker16s_calls(path)[0]
        self.assertAlmostEqual(row.raw_reads, 39224.860)
        self.assertAlmostEqual(row.reads_assigned, 13074.953)

    def test_a_species_sample_against_a_genus_control_is_not_inflated(self):
        """Sample resolves C. resistens (3.00 copies); the blank only resolves
        Corynebacterium (4.08). On corrected reads the ratio picks up a spurious
        4.08/3.00 = 1.36x and reads 14.61; on raw reads it is 10.82.
        """
        from bridge.reconcile.rpm import compute_rpm, evaluate

        sample_total, control_total = 3_355_879, 22_539
        raw_sample, raw_control = 39_224.860, 24.341
        corrected_sample, corrected_control = 13_074.953, 6.010

        on_raw = evaluate(
            "Corynebacterium resistens", "marker16s", raw_sample, sample_total,
            compute_rpm(raw_control, control_total))
        on_corrected = evaluate(
            "Corynebacterium resistens", "marker16s", corrected_sample,
            sample_total, compute_rpm(corrected_control, control_total))

        self.assertAlmostEqual(on_raw.rpm_r, 10.82, places=1)
        self.assertAlmostEqual(on_corrected.rpm_r, 14.61, places=1)
        # The gap is exactly the ratio of the two divisors - no read differs.
        # Derived from these same figures rather than from the copy_number
        # column, which is reported to 3 decimals: the control's 6.010 implies
        # 4.050, not the 4.080 printed beside it.
        divisor_sample = raw_sample / corrected_sample
        divisor_control = raw_control / corrected_control
        self.assertAlmostEqual(on_corrected.rpm_r / on_raw.rpm_r,
                               divisor_control / divisor_sample, places=3)

    def test_matching_ranks_cancel_either_way(self):
        """Where both sides resolve alike the choice is invisible, which is why
        this went unnoticed until a control resolved more coarsely than a
        sample."""
        from bridge.reconcile.rpm import compute_rpm, evaluate

        for divisor in (1.0, 4.2, 15.0):
            with self.subTest(divisor=divisor):
                result = evaluate(
                    "X", "marker16s", 40_000.0 / divisor, 3_000_000,
                    compute_rpm(400.0 / divisor, 20_000))
                self.assertAlmostEqual(result.rpm_r, 0.6666, places=3)

    def test_the_baseline_builder_reads_the_same_column(self):
        """If the two ever disagree every RPM-r in the run is wrong by the
        ratio of two copy numbers, and nothing would say so."""
        import tempfile
        from bridge.reconcile.baseline import marker16s_profile

        path = self._rows(
            Path(tempfile.mkdtemp()) / "marker16s_calls.csv",
            ["NEG,Corynebacterium,Corynebacterium,,none,silent_weak,0.01,0.01,"
             "6.010,24.341,4.080,0.49,0.35,1,x,y\n"])
        profile = marker16s_profile(path)
        self.assertAlmostEqual(profile["Corynebacterium"], 24.341)


# At the BOTTOM, deliberately: anywhere higher, `python tests/test_adapters.py`
# would exit before the remaining classes were defined.
if __name__ == "__main__":
    unittest.main()
