"""Tests for the 16S marker arm.

The header cases here are real shapes taken from
``MIMt-16S_M2c_25_04_tax_bacteria_tax_in_id.fna``, not invented ones - the
parser's whole job is to survive that file, so the tests are only worth
anything if they use what is actually in it.
"""

from __future__ import annotations

import unittest
from pathlib import Path

import pandas as pd

from bridge.marker16s.aggregate import (
    DEFAULT_COPY_NUMBER, aggregate_taxa, copy_number_coverage, genus_totals,
)
from bridge.marker16s.calls import marker16s_calls
from bridge.marker16s.reference import parse_lineage_header


def header(genus: str, species_field: str, accession: str = "NR_000001.1") -> str:
    return (f"{accession};K__Bacteria;P__Pseudomonadota;C__Gammaproteobacteria;"
            f"O__Enterobacterales;F__Enterobacteriaceae;G__{genus};"
            f"S__{species_field}")


class ReferenceParsingTests(unittest.TestCase):

    def test_plain_lineage(self):
        taxon = parse_lineage_header(header("Klebsiella", "Klebsiella_pneumoniae"))
        self.assertEqual(taxon.genus, "Klebsiella")
        self.assertEqual(taxon.species, "pneumoniae")
        self.assertEqual(taxon.species_full, "Klebsiella pneumoniae")
        self.assertEqual(taxon.accession, "NR_000001.1")

    def test_strain_suffix_is_not_part_of_the_epithet(self):
        """22,438 of 29,090 MIMt entries carry a _strain_ suffix."""
        taxon = parse_lineage_header(
            header("Puzihella", "Puzihella_rosea_strain_MI-28"))
        self.assertEqual(taxon.species, "rosea")

    def test_isolate_and_bare_strain_codes(self):
        for field, expected in [
            ("Stigonema_mamillosum_isolate_N01", "mamillosum"),
            ("Geodermatophilus_obscurus_DSM_43160", "obscurus"),
            ("Borreliella_bavariensis_PBi", "bavariensis"),
            ("Picosynechococcus_fontinalis_FG-1", "fontinalis"),
        ]:
            with self.subTest(field=field):
                self.assertEqual(parse_lineage_header(header("G", field)).species,
                                 expected)

    def test_curated_genus_wins_over_a_stale_species_string(self):
        """G__ and S__ disagree for reclassified organisms.

        ``G__Agrilutibacter;S__Lysobacter_terrestris_strain_II4`` is real. The
        G__ field is curated and current; the species string still carries the
        old genus. Taking the epithet positionally means the stale genus there
        costs nothing.
        """
        taxon = parse_lineage_header(
            header("Agrilutibacter", "Lysobacter_terrestris_strain_II4"))
        self.assertEqual(taxon.genus, "Agrilutibacter")
        self.assertEqual(taxon.species, "terrestris")

    def test_nomenclatural_authority_is_stripped_from_the_genus(self):
        taxon = parse_lineage_header(header(
            "Methylacidiphilum_(ex_Ratnadevi_et_al._2023)",
            "Methylacidiphilum_kamchatkense_Kam1"))
        self.assertEqual(taxon.genus, "Methylacidiphilum")

    def test_unresolved_species_yields_genus_only_not_a_guess(self):
        for field in ["Escherichia_sp._K12", "Escherichia_uncultured_clone1"]:
            with self.subTest(field=field):
                taxon = parse_lineage_header(header("Escherichia", field))
                self.assertEqual(taxon.genus, "Escherichia")
                self.assertIsNone(taxon.species)
                self.assertEqual(taxon.label, "Escherichia")

    def test_malformed_header_does_not_raise(self):
        """A reference the arm cannot read is reported, not crashed on."""
        taxon = parse_lineage_header(">not;a;lineage")
        self.assertIsNone(taxon.genus)
        self.assertIsNone(taxon.species)


class AggregationTests(unittest.TestCase):

    def _abundance(self, rows):
        refs, records = {}, []
        for accession, genus, species, counts in rows:
            head = header(genus, f"{genus}_{species}", accession)
            refs[head] = parse_lineage_header(head)
            records.append(dict(target_id=head, length=1500, eff_length=1350,
                                est_counts=counts, tpm=0.0))
        return pd.DataFrame(records), refs

    def test_references_of_one_taxon_are_summed(self):
        frame, refs = self._abundance([
            ("A", "Klebsiella", "pneumoniae", 100.0),
            ("B", "Klebsiella", "pneumoniae", 300.0),
        ])
        taxa = aggregate_taxa(frame, refs, min_reads=1.0)
        self.assertEqual(len(taxa), 1)
        self.assertEqual(taxa[0].raw_reads, 400.0)
        self.assertEqual(taxa[0].n_references, 2)

    def test_copy_number_correction_divides_multi_copy_organisms_down(self):
        """The rrnDB correction: 16S copy number spans 1-15 per genome.

        Two organisms with equal raw reads but different copy number must not
        report equal abundance - that is the bias the correction exists for.
        """
        frame, refs = self._abundance([
            ("A", "Klebsiella", "pneumoniae", 800.0),   # 8 copies
            ("B", "Mycobacterium", "tuberculosis", 800.0),  # 1 copy
        ])
        table = {"klebsiella pneumoniae": 8.0, "mycobacterium tuberculosis": 1.0}
        taxa = {t.taxon: t for t in aggregate_taxa(frame, refs, copy_numbers=table)}
        self.assertAlmostEqual(taxa["Klebsiella pneumoniae"].reads, 100.0)
        self.assertAlmostEqual(taxa["Mycobacterium tuberculosis"].reads, 800.0)
        # The raw figure is preserved so the correction stays auditable.
        self.assertEqual(taxa["Klebsiella pneumoniae"].raw_reads, 800.0)

    def test_missing_copy_number_falls_back_to_the_database_mean(self):
        frame, refs = self._abundance([("A", "Nocardia", "farcinica", 420.0)])
        taxa = aggregate_taxa(frame, refs, copy_numbers={"escherichia coli": 7.0})
        self.assertAlmostEqual(taxa[0].copy_number, DEFAULT_COPY_NUMBER)

    def test_unparsable_reference_is_reported_not_dropped(self):
        frame = pd.DataFrame([dict(target_id="mystery", length=1, eff_length=1,
                                   est_counts=50.0, tpm=0.0)])
        taxa = aggregate_taxa(frame, {}, min_reads=1.0)
        self.assertEqual([t.taxon for t in taxa], ["unassigned"])


class CallTests(unittest.TestCase):

    def _taxa(self, rows):
        refs, records = {}, []
        for genus, species, counts in rows:
            head = header(genus, f"{genus}_{species}", f"{genus}_{species}")
            refs[head] = parse_lineage_header(head)
            records.append(dict(target_id=head, length=1500, eff_length=1350,
                                est_counts=counts, tpm=0.0))
        return aggregate_taxa(pd.DataFrame(records), refs, min_reads=0.5)

    def test_dominant_species_is_named(self):
        calls = marker16s_calls("s", self._taxa([("Klebsiella", "pneumoniae", 5000.0)]))
        self.assertEqual(calls[0].strength, "species_strong")
        self.assertEqual(calls[0].species, "pneumoniae")
        self.assertEqual(calls[0].rank, "species")

    def test_mitis_group_split_reports_genus_only(self):
        """The arm must decline where 16S genuinely cannot resolve.

        Naming S. pneumoniae from evidence that cannot exclude S. mitis or
        S. oralis is the classic 16S over-call. The margin rule catches it with
        no curated exception list - the
        mass genuinely splits, so the call stops at the genus.
        """
        calls = marker16s_calls("s", self._taxa([
            ("Streptococcus", "mitis", 400.0),
            ("Streptococcus", "oralis", 380.0),
            ("Streptococcus", "pneumoniae", 350.0),
        ]))
        self.assertEqual(len(calls), 1)
        self.assertEqual(calls[0].strength, "genus_only")
        self.assertIsNone(calls[0].species)
        self.assertEqual(calls[0].taxon, "Streptococcus")
        self.assertLess(calls[0].margin, 0.20)

    def test_escherichia_shigella_style_tie_reports_genus_only(self):
        calls = marker16s_calls("s", self._taxa([
            ("Escherichia", "coli", 1000.0),
            ("Escherichia", "fergusonii", 900.0),
        ]))
        self.assertEqual(calls[0].strength, "genus_only")

    def test_trace_genus_is_silent_not_called(self):
        """A genus can hold 100% of a near-empty profile; that is not a finding."""
        calls = {c.genus: c for c in marker16s_calls("s", self._taxa([
            ("Klebsiella", "pneumoniae", 9000.0),
            ("Cutibacterium", "acnes", 2.0),
        ]))}
        self.assertEqual(calls["Cutibacterium"].strength, "silent_weak")
        self.assertEqual(calls["Klebsiella"].strength, "species_strong")

    def test_absolute_read_floor_applies_even_at_high_proportion(self):
        calls = marker16s_calls("s", self._taxa([("Gardnerella", "vaginalis", 4.0)]))
        self.assertEqual(calls[0].strength, "silent_weak")
        self.assertIn("reads", calls[0].reason)

    def test_the_silent_reason_names_only_the_clause_that_failed(self):
        """A genus silenced on reads alone must not report a passing
        fraction as the reason. Each clause is stated only when it fired."""
        # Gardnerella alone: 100% of the profile, 0.95 corrected reads -> reads only.
        reads_only = marker16s_calls("s", self._taxa([("Gardnerella", "vaginalis", 4.0)]))[0]
        self.assertIn("corrected reads", reads_only.reason)
        self.assertNotIn("frac", reads_only.reason)
        # Cutibacterium at 50 raw (11.9 corrected) beside 9000 Klebsiella: 0.55%
        # of the profile but comfortably above the read floor -> frac only.
        calls = {c.genus: c for c in marker16s_calls("s", self._taxa([
            ("Klebsiella", "pneumoniae", 9000.0),
            ("Cutibacterium", "acnes", 50.0),
        ]))}
        frac_only = calls["Cutibacterium"]
        self.assertEqual(frac_only.strength, "silent_weak")
        self.assertIn("frac", frac_only.reason)
        self.assertNotIn("reads", frac_only.reason)

    def test_one_call_per_genus(self):
        calls = marker16s_calls("s", self._taxa([
            ("Streptococcus", "mitis", 500.0),
            ("Streptococcus", "oralis", 400.0),
            ("Klebsiella", "pneumoniae", 3000.0),
        ]))
        self.assertEqual(sorted(c.genus for c in calls),
                         ["Klebsiella", "Streptococcus"])

    def test_genus_totals_sum_members(self):
        taxa = self._taxa([("Streptococcus", "mitis", 300.0),
                           ("Streptococcus", "oralis", 100.0)])
        self.assertAlmostEqual(genus_totals(taxa)["Streptococcus"],
                               sum(t.reads for t in taxa))

class CopyNumberProvenanceTests(unittest.TestCase):
    """A table that matched nothing must not look like a table that worked.

    Without a real table every taxon is divided by DEFAULT_COPY_NUMBER, and a
    constant divisor cancels out of every ratio the arm computes - the output
    still says "copy_number_corrected" while correcting nothing.
    """

    def _abundance(self, rows):
        refs, records = {}, []
        for accession, genus, species, counts in rows:
            head = header(genus, f"{genus}_{species}", accession)
            refs[head] = parse_lineage_header(head)
            records.append(dict(target_id=head, length=1500, eff_length=1350,
                                est_counts=counts, tpm=0.0))
        return pd.DataFrame(records), refs

    def test_the_resolution_level_is_recorded_per_taxon(self):
        frame, refs = self._abundance([
            ("A", "Escherichia", "coli", 100.0),
            ("B", "Staphylococcus", "caprae", 100.0),
            ("C", "Nocardia", "farcinica", 100.0),
        ])
        table = {"escherichia coli": 7.0,        # species hit
                 "staphylococcus": 5.9}          # genus hit only
        taxa = {t.taxon: t for t in aggregate_taxa(frame, refs, copy_numbers=table)}
        self.assertEqual(taxa["Escherichia coli"].copy_number_source, "species")
        self.assertEqual(taxa["Staphylococcus caprae"].copy_number_source, "genus")
        self.assertEqual(taxa["Nocardia farcinica"].copy_number_source, "default")

    def test_no_table_marks_every_taxon_as_default(self):
        frame, refs = self._abundance([("A", "Escherichia", "coli", 100.0)])
        taxa = aggregate_taxa(frame, refs)
        self.assertEqual(taxa[0].copy_number_source, "default")
        self.assertAlmostEqual(taxa[0].copy_number, DEFAULT_COPY_NUMBER)

    def test_coverage_is_weighted_by_reads_not_by_taxon(self):
        """One dominant organism resolved beats fifty trace ones missed."""
        frame, refs = self._abundance([
            ("A", "Escherichia", "coli", 9900.0),
            ("B", "Nocardia", "farcinica", 50.0),
            ("C", "Gordonia", "bronchialis", 50.0),
        ])
        taxa = aggregate_taxa(frame, refs, copy_numbers={"escherichia coli": 7.0})
        coverage = copy_number_coverage(taxa)
        self.assertEqual(coverage["species"], 1)
        self.assertEqual(coverage["default"], 2)
        self.assertAlmostEqual(coverage["reads_resolved_frac"], 0.99)

    def test_a_uniform_divisor_changes_no_within_sample_ratio(self):
        """Why the missing table mattered, stated as an assertion.

        With one divisor for everything the corrected fractions are identical
        to the raw ones, so nothing downstream - margin, within-genus share,
        RPM-r - sees any effect at all.
        """
        frame, refs = self._abundance([
            ("A", "Escherichia", "coli", 700.0),
            ("B", "Mycobacterium", "tuberculosis", 300.0),
        ])
        uncorrected = {t.taxon: t.abundance_frac
                       for t in aggregate_taxa(frame, refs)}
        self.assertAlmostEqual(uncorrected["Escherichia coli"], 0.7)
        self.assertAlmostEqual(uncorrected["Mycobacterium tuberculosis"], 0.3)

        # A real table moves them, and by a lot: E. coli carries ~7 rrn copies
        # against M. tuberculosis's 1, so equal 16S mass is a 7x cell-count gap.
        corrected = {t.taxon: t.abundance_frac for t in aggregate_taxa(
            frame, refs, copy_numbers={"escherichia coli": 7.0,
                                       "mycobacterium tuberculosis": 1.0})}
        self.assertAlmostEqual(corrected["Escherichia coli"], 0.25)
        self.assertAlmostEqual(corrected["Mycobacterium tuberculosis"], 0.75)


class CopyNumberTableResolutionTests(unittest.TestCase):
    """Which table the 16S CLI uses, and when it refuses."""

    class _Args:
        def __init__(self, explicit=None, disabled=False):
            self.copy_number_table = explicit
            self.no_copy_number_table = disabled

    def test_an_explicit_table_that_does_not_exist_is_an_error(self):
        """Not a silent fallback: passing the flag means you meant it."""
        from bridge.marker16s.cli import resolve_copy_number_table
        with self.assertRaises(FileNotFoundError):
            resolve_copy_number_table(self._Args(Path("no/such/table.tsv")))

    def test_opting_out_beats_an_existing_default(self):
        from bridge.marker16s.cli import resolve_copy_number_table
        self.assertIsNone(resolve_copy_number_table(self._Args(disabled=True)))

    def test_an_explicit_table_is_used_when_it_exists(self):
        import tempfile
        from bridge.marker16s.cli import (
            load_copy_numbers, resolve_copy_number_table,
        )
        path = Path(tempfile.mkdtemp()) / "cn.tsv"
        path.write_text("# comment\nescherichia coli\t7.0\nbroken line\n"
                        "bad\tvalue\n", encoding="utf-8")
        self.assertEqual(resolve_copy_number_table(self._Args(path)), path)
        # And the loader skips the junk rather than failing on it.
        self.assertEqual(load_copy_numbers(path), {"escherichia coli": 7.0})

    def test_an_empty_table_is_an_error_not_a_silent_no_op(self):
        import tempfile
        from bridge.marker16s.cli import load_copy_numbers
        path = Path(tempfile.mkdtemp()) / "empty.tsv"
        path.write_text("# nothing but comments\n", encoding="utf-8")
        with self.assertRaises(ValueError):
            load_copy_numbers(path)
class RawGenusTotalsTests(unittest.TestCase):
    """Corrected totals answer composition; raw totals answer depth.

    The margin rule wants corrected reads - two species at equal RAW reads but
    1 vs 8 copies are not equally abundant. The RPM-r gate and any detection
    floor want raw ones. Both are available, and confusing them is how a
    copy-number table ends up changing a threshold it has no business touching.
    """

    def _taxa(self, rows):
        refs, records = {}, []
        for accession, genus, species, counts in rows:
            head = header(genus, f"{genus}_{species}", accession)
            refs[head] = parse_lineage_header(head)
            records.append(dict(target_id=head, length=1500, eff_length=1350,
                                est_counts=counts, tpm=0.0))
        return pd.DataFrame(records), refs

    def test_raw_totals_ignore_the_correction(self):
        frame, refs = self._taxa([("A", "Sphingorhabdus", "flavimaris", 27.9),
                                  ("B", "Escherichia", "coli", 1000.0)])
        taxa = aggregate_taxa(frame, refs, min_reads=0.5,
                              copy_numbers={"sphingorhabdus": 1.8,
                                            "escherichia coli": 7.0})
        raw = genus_totals(taxa, raw=True)
        corrected = genus_totals(taxa)
        self.assertAlmostEqual(raw["Sphingorhabdus"], 27.9)
        self.assertAlmostEqual(corrected["Sphingorhabdus"], 27.9 / 1.8)
        self.assertAlmostEqual(raw["Escherichia"], 1000.0)
        self.assertAlmostEqual(corrected["Escherichia"], 1000.0 / 7.0)

    def test_corrected_never_exceeds_raw(self):
        """Copy numbers are all >= 1, so moving a threshold from corrected to
        raw reads can only loosen it."""
        frame, refs = self._taxa([("A", "Bacillus", "cereus", 500.0)])
        taxa = aggregate_taxa(frame, refs, copy_numbers={"bacillus": 9.0})
        self.assertLessEqual(genus_totals(taxa)["Bacillus"],
                             genus_totals(taxa, raw=True)["Bacillus"])

    def test_the_margin_rule_still_uses_corrected_reads(self):
        """Composition is exactly where the correction belongs."""
        frame, refs = self._taxa([("A", "Bacillus", "cereus", 800.0),
                                  ("B", "Bacillus", "subtilis", 800.0)])
        calls = marker16s_calls("s", aggregate_taxa(
            frame, refs, copy_numbers={"bacillus cereus": 1.0,
                                       "bacillus subtilis": 8.0}))
        # 800/1 vs 800/8 -> 8:1, so the genus is NOT evenly split.
        self.assertEqual(calls[0].strength, "species_strong")
        self.assertEqual(calls[0].species, "cereus")
        self.assertGreater(calls[0].within_genus_share, 0.85)


if __name__ == "__main__":
    unittest.main()
