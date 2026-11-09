"""Tiers, allele levels and phenotype certainty (bridge.amr.calls, bridge.amr.phenotype)."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from bridge.amr.calls import allele_level, classify, collapse_genes, est_differences, tier
from bridge.amr.database import load_family_overrides, load_phenotypes
from bridge.amr.phenotype import annotate, is_uncertain
from bridge.amr.settings import DATA_DIR, Thresholds

from amr_fixtures import write_phenotypes

TH = Thresholds()


def hit(identity=100.0, coverage=100.0, depth=20.0, fragments=100, query_identity=None,
        family="blaTEM", template="blaTEM-1B_1_AY458016", length=861):
    return {"template": template, "gene": template.rsplit("_", 2)[0],
            "family": family, "template_identity": str(identity),
            "template_coverage": str(coverage), "depth": str(depth),
            "fragments": str(fragments), "template_length": str(length),
            "query_identity": str(identity if query_identity is None else query_identity),
            "flags": ""}


class TierTests(unittest.TestCase):
    def test_tiers(self):
        self.assertEqual(tier(hit(), TH), "confident")
        self.assertEqual(tier(hit(coverage=75.0), TH), "partial")
        self.assertEqual(tier(hit(fragments=3), TH), "partial")       # too little material
        self.assertEqual(tier(hit(depth=1.0), TH), "partial")
        self.assertEqual(tier(hit(identity=85.0), TH), "below")
        self.assertEqual(tier(hit(coverage=40.0), TH), "below")

    def test_identity_is_judged_where_the_gene_is_covered(self):
        # Ureaplasma tet(M): 99.77% where covered, over 89% of it
        ureaplasma = hit(identity=88.75, coverage=88.96, query_identity=99.77, depth=77, fragments=1514)
        self.assertEqual(tier(ureaplasma, TH), "confident")
        self.assertEqual(allele_level(ureaplasma, TH), "nearest")
        # Serratia's oqxB homologue: genuinely divergent where covered
        self.assertEqual(tier(hit(identity=56.42, coverage=65.91, query_identity=85.61), TH), "below")

    def test_allele_levels(self):
        self.assertEqual(allele_level(hit(), TH), "exact")
        self.assertEqual(allele_level(hit(identity=99.88), TH), "nearest")
        self.assertEqual(allele_level(hit(coverage=98.0), TH), "nearest")
        # KMA's two decimals: 99.995 would print as 100.00
        self.assertEqual(allele_level(hit(identity=100.00), TH), "exact")

    def test_est_differences_over_the_covered_part(self):
        self.assertEqual(est_differences(hit(query_identity=99.77, identity=99.77)), 2)
        self.assertEqual(est_differences(hit()), 0)

    def test_two_alleles_of_one_family_are_flagged(self):
        hits = [hit(), hit(template="blaTEM-1A_1_HM749966", identity=99.9),
                hit(family="tet(M)", template="tet(M)_1_X92947"),
                hit(family="blaTEM", template="blaTEM-2_1_X", identity=50.0)]
        classify(hits, TH)
        self.assertIn("family_multiple_alleles", hits[0]["flags"])
        self.assertIn("family_multiple_alleles", hits[1]["flags"])
        self.assertEqual(hits[2]["flags"], "")
        self.assertEqual(hits[3]["tier"], "below")
        self.assertEqual(hits[3]["allele_level"], "")
        # re-classifying does not duplicate the flag
        classify(hits, TH)
        self.assertEqual(hits[0]["flags"].count("family_multiple_alleles"), 1)

    def test_two_submissions_of_one_allele_are_not_multiple_alleles(self):
        hits = classify([hit(), hit(template="blaTEM-1B_2_XX000002")], TH)
        self.assertFalse(any("family_multiple_alleles" in h["flags"] for h in hits))

    def test_collapse_keeps_the_best_template_and_sums_reads(self):
        raw = [hit(template="tet(M)_2_B", family="tet(M)", identity=99.0, fragments=40, depth=4),
               hit(template="tet(M)_1_A", family="tet(M)", fragments=100, depth=10),
               hit(template="tet(M)_3_C", family="tet(M)", identity=50.0, fragments=5)]
        for h in raw:
            h["gene"] = "tet(M)"
        out = collapse_genes(raw, TH)
        kept = [h for h in out if h["tier"] != "below"]
        self.assertEqual(len(kept), 1)
        self.assertEqual(kept[0]["template"], "tet(M)_1_A")      # exact beats nearest
        self.assertEqual(kept[0]["fragments"], "140")            # below-floor row not added
        self.assertEqual(kept[0]["depth"], "14.00")
        self.assertIn("templates_merged:2", kept[0]["flags"])
        self.assertEqual(len(out), 2)                            # the below row is kept as is
        self.assertEqual(raw[0]["fragments"], "40")               # input untouched


class PhenotypeCertaintyTests(unittest.TestCase):
    def setUp(self):
        overrides = load_family_overrides(DATA_DIR / "family_overrides.tsv")
        self.db = load_phenotypes(write_phenotypes(Path(tempfile.mkdtemp())), overrides)

    def test_exact_shv12_is_certain_and_names_its_esbl_spectrum(self):
        h = classify([hit(family="blaSHV", template="blaSHV-12_1_AF462395")], TH)[0]
        annotate(h, self.db)
        self.assertIn("Cefotaxime", h["antibiotics"])
        self.assertEqual(h["family_spectrum_varies"], "True")
        self.assertFalse(is_uncertain(h))

    def test_nearest_shv1_is_uncertain_and_lists_the_esbl_drugs_as_possible(self):
        h = classify([hit(family="blaSHV", template="blaSHV-1_1_AF148850",
                          identity=99.7)], TH)[0]
        annotate(h, self.db)
        self.assertNotIn("Cefotaxime", h["antibiotics"])
        self.assertIn("cefotaxime", h["family_possible_antibiotics"])
        self.assertTrue(is_uncertain(h))

    def test_nearest_in_a_family_that_agrees_with_itself_is_certain(self):
        h = classify([hit(identity=99.7)], TH)[0]         # blaTEM-1B vs 1A: same spectrum
        annotate(h, self.db)
        self.assertFalse(is_uncertain(h))

    def test_template_without_a_phenotype_entry_is_flagged(self):
        h = hit(template="blaNOPE-1_1_X", family="blaNOPE")
        annotate(h, self.db)
        self.assertIn("no_phenotype_entry", h["flags"])
        self.assertEqual(h["antibiotics"], "")


if __name__ == "__main__":
    unittest.main()
