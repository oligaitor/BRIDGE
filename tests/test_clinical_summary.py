"""The clinician view (bridge.clinical): organisms, resistance by class, statements."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from amr_fixtures import amr_sample, read_csv, settings_with_db
from bridge.clinical import build_clinical, by_class, gene_label
from bridge.outcome import OutcomeSettings, build_outcomes


def gene(sample, name, drug_class, attribution, host="", **extra):
    return {"sample": sample, "gene": name, "drug_class": drug_class, "status": "report",
            "attribution": attribution, "host_organism": host, **extra}


class ClinicalTests(unittest.TestCase):
    def build(self, id_rows, amr_rows, phenotype_rows=(), depth=None, samples=None):
        samples = samples or sorted({row["sample"] for row in id_rows})
        outcomes = build_outcomes("RUN1", samples, id_rows=id_rows, amr_rows=amr_rows,
                                  depth=depth or {}, settings=OutcomeSettings(),
                                  control_prefixes=("NEG",))
        return build_clinical(outcomes, id_rows, amr_rows, list(phenotype_rows))

    def test_one_row_per_reported_organism_with_its_own_genes(self):
        ids = [{"sample": "S", "organism": "Enterococcus faecium", "status": "report",
                "reported_rank": "species"},
               {"sample": "S", "organism": "Streptococcus", "status": "report",
                "reported_rank": "genus"},
               {"sample": "S", "organism": "Escherichia", "status": "review",
                "reported_rank": "genus"}]
        amr = [gene("S", "msr(C)", "Macrolide, Streptogramin B", "intrinsic_host",
                    "Enterococcus faecium", host_verdict="intrinsic"),
               gene("S", "vanA", "Glycopeptide", "unassigned"),
               gene("S", "aac(6')-Ii", "Aminoglycoside", "unassigned",
                    phenotype_uncertain="True"),
               {**gene("S", "blaTEM-1A", "Beta-lactam", "unassigned"), "status": "review"}]
        pheno = [{"sample": "S", "organism": "Enterococcus faecium", "agent": "cefotaxime",
                  "eucast_expected": "True"},
                 {"sample": "S", "organism": "Enterococcus faecium", "agent": "aztreonam",
                  "eucast_expected": "True"},
                 {"sample": "S", "organism": "Enterococcus faecium", "agent": "erythromycin",
                  "eucast_expected": ""}]
        rows = self.build(ids, amr, pheno)
        self.assertEqual([row["organism"] for row in rows],
                         ["Enterococcus faecium", "Streptococcus sp."])
        faecium, strep = rows
        # msr(C) is species-typical for E. faecium: kept apart from resistance
        self.assertEqual(faecium["resistance_genes"], "No acquired resistance genes detected")
        self.assertEqual(faecium["species_typical_genes"], "Macrolide, Streptogramin B: msr(C)")
        self.assertEqual(faecium["expected_resistance"], "aztreonam, cefotaxime")
        self.assertEqual(strep["resistance_genes"], "No acquired resistance genes detected")
        self.assertEqual(strep["species_typical_genes"], "")
        # genes that belong to no one organism are on every row of the sample, never
        # dropped; a held (review) gene is not shown at all
        for row in rows:
            self.assertEqual(row["unlinked_resistance"],
                             "Glycopeptide: vanA | Aminoglycoside: aac(6')-Ii (possible)")
            self.assertEqual(row["statement"],
                             "Bacteria detected: Enterococcus faecium, Streptococcus sp.")

    def test_a_sample_without_a_reported_organism_gets_one_statement_row(self):
        ids = [{"sample": "N", "status": "no_call"},
               {"sample": "H", "organism": "Massilia", "status": "review_override"}]
        amr = [gene("H", "tet(M)", "Tetracycline", "no_organism")]
        rows = {row["sample"]: row for row in self.build(ids, amr)}
        self.assertEqual(rows["N"]["statement"], "No bacteria present above reporting thresholds")
        self.assertEqual((rows["N"]["organism"], rows["N"]["unlinked_resistance"]), ("", ""))
        self.assertEqual(rows["H"]["outcome"], "held_for_review")
        self.assertEqual(rows["H"]["organism"], "")    # a held organism is not named
        self.assertEqual(rows["H"]["unlinked_resistance"], "Tetracycline: tet(M)")

    def test_controls_and_unrun_samples_are_left_out(self):
        rows = self.build([{"sample": "NEG_1", "status": "no_call"}], [])
        self.assertEqual(rows, [])

    def test_by_class_keeps_first_seen_order_and_drops_duplicates(self):
        calls = [gene("S", "tet(M)", "Tetracycline", "x"), gene("S", "erm(B)", "", "x"),
                 gene("S", "tet(M)", "Tetracycline", "x"), gene("S", "tet(L)", "Tetracycline", "x"),
                 gene("S", "tet(39)", "Tetracyclines", "x")]
        self.assertEqual(by_class(calls),
                         "Tetracycline: tet(M), tet(L), tet(39) | Unclassified: erm(B)")
        self.assertEqual(gene_label({"gene": "x", "phenotype_uncertain": "True"}), "x (possible)")
        # a species-typical gene's spectrum is not a prediction, so no "(possible)"
        self.assertEqual(gene_label({"gene": "x", "host_verdict": "intrinsic",
                                     "phenotype_uncertain": "True"}), "x")


class FinaliseClinicalTests(unittest.TestCase):
    def test_finalise_writes_the_clinical_summary_from_real_amr_output(self):
        from bridge.finalise import finalise_run
        from bridge.settings import load_settings
        from test_finalise import RMLST_HEADER, S16_HEADER, rmlst_row, s16_row
        tmp = Path(tempfile.mkdtemp())
        root = tmp / "results"
        amr_settings = settings_with_db(tmp)
        folder = root / "A"
        (folder / "rmlst").mkdir(parents=True)
        (folder / "16s" / "quant").mkdir(parents=True)
        (folder / "rmlst" / "species_calls.csv").write_text(
            RMLST_HEADER + rmlst_row("A", "Escherichia coli"))
        (folder / "16s" / "marker16s_calls.csv").write_text(
            S16_HEADER + s16_row("A", "Escherichia coli", "Escherichia", "coli", "species_strong"))
        (folder / "16s" / "quant" / "run_info.json").write_text(json.dumps({"n_processed": 1_000_000}))
        amr_sample(root, "A", [("blaCTX-M-15_1_AY044436", 100.0, 100.0, 80.0, 900)],
                   1_000_000, amr_settings)
        with patch("bridge.amr.settings.from_config", return_value=amr_settings):
            result = finalise_run(load_settings(), root, "RUN1")
        self.assertIsNone(result.outcome_error)
        outcome = read_csv(root / "sample_outcomes_all.csv")[0]
        self.assertEqual((outcome["outcome"], outcome["amr_outcome"]),
                         ("detected", "genes_reported"))
        row = read_csv(root / "clinical_summary_all.csv")[0]
        self.assertEqual(row["organism"], "Escherichia coli")
        self.assertEqual(row["resistance_genes"], "Beta-lactam: blaCTX-M-15")
        self.assertEqual(row["species_typical_genes"], "")
        self.assertEqual(row["unlinked_resistance"], "")


if __name__ == "__main__":
    unittest.main()
