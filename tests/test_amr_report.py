"""Host table, EUCAST table, the resistance gate, the report, and the CLI/scheduler wiring."""

from __future__ import annotations

import json
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

from bridge.amr.eucast import expected_for, load_expected
from bridge.amr.host import load_host_context, resolve
from bridge.amr.report import attribute, finalise_amr
from bridge.amr.settings import DATA_DIR

from amr_fixtures import amr_sample, id_summary, read_csv, settings_with_db

HOST = load_host_context(DATA_DIR / "gene_host_context.tsv")
TEM = ("blaTEM-1B_1_AY458016", 100.0, 100.0, 50.0, 600)


class HostTableTests(unittest.TestCase):
    def test_shipped_table_loads_and_answers(self):
        self.assertEqual(resolve("blaACT", "Enterobacter hormaechei", HOST)[0], "intrinsic")
        self.assertEqual(resolve("blaTEM", "Escherichia coli", HOST)[0], "acquired")
        self.assertEqual(resolve("blaTEM", "Staphylococcus aureus", HOST)[0],
                         "not_expected_in_taxon")
        self.assertEqual(resolve("blaTEM", "Nosuchgenus", HOST)[0], "undetermined")
        self.assertEqual(resolve(None, "Escherichia coli", HOST)[3], "no_family")

    def test_compound_genus_tries_each_part(self):
        self.assertEqual(resolve("blaTEM", "Escherichia-Shigella", HOST)[0], "acquired")

    def test_an_uncited_verdict_is_refused(self):
        path = Path(tempfile.mkdtemp()) / "t.tsv"
        path.write_text("gene_family\ttaxon\ttaxon_rank\tverdict\tcitation\tcitation_type\t"
                        "note\tsource\tcurated_on\nblaX\tE\tgenus\tintrinsic\t\t\t\t\t\n")
        with self.assertRaises(ValueError):
            load_host_context(path)


class EucastTests(unittest.TestCase):
    def setUp(self):
        self.rows = load_expected(DATA_DIR / "eucast_expected_resistant.tsv")

    def test_klebsiella_pneumoniae_ampicillin_and_the_enterobacterales_block(self):
        agents = {agent for agent, _ in expected_for("Klebsiella pneumoniae", self.rows)}
        self.assertIn("ampicillin", agents)
        self.assertIn("vancomycin", agents)
        self.assertNotIn("cefoxitin", agents)

    def test_enterobacter_complex_member_gets_ampc_agents(self):
        agents = {agent for agent, _ in expected_for("Enterobacter hormaechei", self.rows)}
        self.assertIn("cefoxitin", agents)

    def test_every_shipped_row_has_been_verified(self):
        # Every shipped row is checked against EUCAST Expected Resistant
        # Phenotypes v1.2 (Jan 2023). A row added without checking fails here
        # until it is verified.
        self.assertTrue(self.rows)
        self.assertTrue(all(row.verified for row in self.rows))
        self.assertEqual(expected_for("", self.rows), [])


class AttributionTests(unittest.TestCase):
    TEM_HIT = {"family": "blaTEM", "depth": "40"}
    ACT_HIT = {"family": "blaACT", "depth": "40"}

    def test_tiers(self):
        self.assertEqual(attribute(self.TEM_HIT, None, HOST)["attribution"], "no_id")
        self.assertEqual(attribute(self.TEM_HIT, [], HOST)["attribution"], "no_organism")
        sole = attribute(self.TEM_HIT, [{"organism": "Escherichia coli", "basis": "rmlst",
                                        "median_projected_depth": "20"}], HOST)
        self.assertEqual((sole["attribution"], sole["host_organism"], sole["host_verdict"]),
                         ("sole_organism", "Escherichia coli", "acquired"))
        self.assertEqual(sole["depth_ratio_to_candidates"], "Escherichia coli: 2")

        two = [{"organism": "Escherichia coli", "basis": "rmlst", "median_projected_depth": "20"},
               {"organism": "Enterobacter hormaechei", "basis": "rmlst",
                "median_projected_depth": "10"}]
        act = attribute(self.ACT_HIT, two, HOST)
        self.assertEqual((act["attribution"], act["host_organism"]),
                         ("intrinsic_host", "Enterobacter hormaechei"))
        tem = attribute(self.TEM_HIT, two, HOST)
        self.assertEqual(tem["attribution"], "unassigned")
        self.assertEqual(tem["host_organism"], "")
        self.assertIn("Escherichia coli (acquired)", tem["candidates"])

    def test_informational_id_rows_are_not_candidates(self):
        rows = [{"organism": "Escherichia coli", "basis": "rmlst"},
                {"organism": "Halomonas", "basis": "16s_only"}]
        self.assertEqual(attribute(self.TEM_HIT, rows, HOST)["attribution"], "sole_organism")

    def test_organisms_held_for_review_are_not_candidates(self):
        # Three 16S-fallback calls, all at reagent background
        rows = [{"organism": name, "basis": "16s_fallback", "status": "review"}
                for name in ("Escherichia ruysiae", "Salmonella enterica", "Shigella sonnei")]
        out = attribute(self.TEM_HIT, rows, HOST)
        self.assertEqual(out["attribution"], "no_organism")
        self.assertIn("Salmonella enterica", out["host_note"])
        rows.append({"organism": "Klebsiella pneumoniae", "basis": "rmlst", "status": "report"})
        out = attribute(self.TEM_HIT, rows, HOST)
        self.assertEqual((out["attribution"], out["host_organism"]),
                         ("sole_organism", "Klebsiella pneumoniae"))

    def test_several_16s_fallback_calls_are_not_candidates(self):
        # Three fallback calls, no control background for any
        rows = [{"organism": name, "basis": "16s_fallback", "status": "review_override"}
                for name in ("Moraxella tetraodonis", "Paracoccus", "Rivihabitans pingtungensis")]
        out = attribute(self.TEM_HIT, rows, HOST)
        self.assertEqual(out["attribution"], "no_organism")
        self.assertIn("Paracoccus", out["host_note"])

    def test_phenotype_table_lists_only_reported_organisms(self):
        from bridge.amr.report import phenotype_rows
        rows = [{"organism": "Salmonella enterica", "basis": "16s_fallback", "status": "review"},
                {"organism": "Klebsiella pneumoniae", "basis": "rmlst", "status": "report"}]
        pheno = phenotype_rows("S", "RUN1", [], rows, load_expected(DATA_DIR / "eucast_expected_resistant.tsv"))
        self.assertEqual({row["organism"] for row in pheno}, {"Klebsiella pneumoniae"})


class FinaliseAmrTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.root = self.tmp / "results"
        self.settings = settings_with_db(self.tmp)

    def test_nearest_allele_with_family_in_blank_is_held(self):
        # The blank carries blaTEM-1A. A sample's EXACT blaTEM-1B keeps its own
        # per-allele gate and is reported; a NEAREST blaTEM-1B (KMA's closest
        # template, not a match) may be the blank's reagent TEM under another
        # name, so it is held for review.
        amr_sample(self.root, "Neg_1", [("blaTEM-1A_1_HM749966", 100.0, 100.0, 5.0, 40)],
                   20_000, self.settings)
        amr_sample(self.root, "E", [("blaTEM-1B_1_AY458016", 100.0, 100.0, 50.0, 600)],
                   1_000_000, self.settings)
        amr_sample(self.root, "N", [("blaTEM-1B_1_AY458016", 97.5, 99.0, 50.0, 600)],
                   1_000_000, self.settings)
        for sample in ("E", "N"):
            id_summary(self.root, sample, [("Escherichia coli", "rmlst", 30.0)])

        finalise_amr(self.settings, self.root, "RUN1", controls=["Neg_1"])
        by = {row["sample"]: row for row in read_csv(self.root / "amr_summary_all.csv")
              if row["gene"] == "blaTEM-1B"}

        exact, nearest = by["E"], by["N"]
        self.assertEqual(exact["allele_level"], "exact")
        self.assertEqual(exact["status"], "report")
        self.assertIn("family_in_control", exact["flags"])
        self.assertNotIn("nearest_family_in_control", exact["flags"])

        self.assertEqual(nearest["allele_level"], "nearest")
        self.assertEqual(nearest["gate_band"], "no_control_background")
        self.assertEqual(nearest["status"], "review")
        self.assertIn("nearest_family_in_control", nearest["flags"])

    def test_background_attribution_and_phenotype_across_a_run(self):
        # The blank carries blaTEM at reagent level (no sample carries it credibly).
        amr_sample(self.root, "Neg_1", [("blaTEM-1A_1_HM749966", 100.0, 100.0, 5.0, 40)],
                   20_000, self.settings)
        # A: E. coli alone; the blank's own TEM-1A at a level near the blank's,
        # a different TEM allele (1B), and an exact CTX-M-15
        amr_sample(self.root, "A", [("blaTEM-1A_1_HM749966", 100.0, 100.0, 5.0, 60),
                                    ("blaTEM-1B_1_AY458016", 100.0, 100.0, 50.0, 600),
                                    ("blaCTX-M-15_1_AY044436", 100.0, 100.0, 80.0, 900)],
                   1_000_000, self.settings)
        id_summary(self.root, "A", [("Escherichia coli", "rmlst", 30.0)])
        # B: two organisms; blaACT is Enterobacter's own, tet(M) belongs to neither
        amr_sample(self.root, "B", [("blaACT-1_1_U58495", 100.0, 100.0, 30.0, 400),
                                    ("tet(M)_1_X92947", 100.0, 100.0, 30.0, 400)],
                   1_000_000, self.settings)
        id_summary(self.root, "B", [("Escherichia coli", "rmlst", 30.0),
                                    ("Enterobacter hormaechei", "rmlst", 15.0)])
        # C: AMR only, no identification; a nearest SHV-1
        amr_sample(self.root, "C", [("blaSHV-1_1_AF148850", 99.65, 100.0, 30.0, 400)],
                   1_000_000, self.settings)
        # D: nothing above the floor
        amr_sample(self.root, "D", [("tet(M)_1_X92947", 60.0, 20.0, 1.0, 2)],
                   1_000_000, self.settings)

        result = finalise_amr(self.settings, self.root, "RUN1", controls=["Neg_1"])
        self.assertEqual(result.failed, [])
        self.assertEqual(sorted(result.reported), ["A", "B", "C", "D"])
        self.assertTrue((self.root / "amr_baselines.csv").is_file())

        calls = read_csv(self.root / "amr_summary_all.csv")
        by = {(row["sample"], row["gene"]): row for row in calls}

        # the blank's own allele, at the blank's level, is background...
        tem1a = by[("A", "blaTEM-1A")]
        self.assertEqual(tem1a["gate_band"], "background")
        self.assertEqual(tem1a["status"], "review")
        self.assertIn("amr_below_background", tem1a["flags"])
        # ...but a different allele is gated on its own level (none in the blank),
        # with the family's presence in the blank flagged, not held against it
        tem1b = by[("A", "blaTEM-1B")]
        self.assertEqual(tem1b["gate_band"], "no_control_background")
        self.assertEqual(tem1b["status"], "report")
        self.assertIn("family_in_control", tem1b["flags"])
        self.assertIn("family_multiple_alleles", tem1b["flags"])
        ctx = by[("A", "blaCTX-M-15")]
        self.assertEqual((ctx["status"], ctx["attribution"], ctx["host_verdict"]),
                         ("report", "sole_organism", "acquired"))

        act = by[("B", "blaACT-1")]
        self.assertEqual((act["attribution"], act["host_organism"], act["host_verdict"]),
                         ("intrinsic_host", "Enterobacter hormaechei", "intrinsic"))
        self.assertEqual(by[("B", "tet(M)")]["attribution"], "unassigned")

        shv = by[("C", "blaSHV-1")]
        self.assertEqual((shv["attribution"], shv["allele_level"], shv["phenotype_uncertain"]),
                         ("no_id", "nearest", "True"))
        self.assertIn("cefotaxime", shv["possible_antibiotics"])

        self.assertEqual(by[("D", "")]["status"], "no_call")
        self.assertNotIn("Neg_1", {row["sample"] for row in calls})

        pheno = read_csv(self.root / "predicted_phenotype_all.csv")
        self.assertTrue(pheno)
        self.assertFalse(any(row["prediction"].lower().startswith("s") for row in pheno))
        index = {(r["sample"], r["organism"], r["agent"]): r for r in pheno}
        cefotaxime = index[("A", "Escherichia coli", "cefotaxime")]
        self.assertEqual(cefotaxime["prediction"],
                         "non-susceptible predicted (acquired determinant)")
        self.assertEqual(cefotaxime["acquired_evidence"], "blaCTX-M-15")
        # the background TEM-1A contributes nothing; the TEM-1B does
        evidence = " ".join(r["acquired_evidence"] for r in pheno)
        self.assertNotIn("blaTEM-1A", evidence)
        self.assertIn("blaTEM-1B", evidence)
        # E. coli's EUCAST block is present and marked verified
        vanco = index[("A", "Escherichia coli", "vancomycin")]
        self.assertEqual((vanco["eucast_expected"], vanco["eucast_verified"]), ("True", "True"))
        # blaACT on its own host is an intrinsic determinant
        cefoxitin = index[("B", "Enterobacter hormaechei", "cefoxitin")]
        self.assertEqual(cefoxitin["prediction"], "resistant (EUCAST expected phenotype)")
        self.assertEqual(cefoxitin["intrinsic_gene_evidence"], "blaACT-1")
        # the row says every basis it rests on, not just the first written
        self.assertEqual(cefoxitin["attribution"], "eucast_expected; intrinsic_host")
        self.assertEqual(vanco["attribution"], "eucast_expected")
        self.assertEqual(cefotaxime["attribution"], "sole_organism")
        # unassigned tet(M) sits at sample level with its candidates
        tet = index[("B", "", "tetracycline")]
        self.assertEqual(tet["attribution"], "unassigned")
        self.assertIn("candidates:", tet["notes"])
        # nearest SHV-1 is only ever "possible"
        shv_rows = [r for r in pheno if r["sample"] == "C"]
        self.assertTrue(shv_rows)
        self.assertTrue(all(r["prediction"].startswith("possible") for r in shv_rows))
        # every shipped EUCAST row is verified, so finalise counts none as unverified
        self.assertEqual(result.unverified_eucast_rows_used, 0)

    def test_an_intrinsic_gene_never_makes_a_prediction_of_its_own(self):
        # blaSHV-1 is intrinsic to K. pneumoniae. ResFinder lists piperacillin
        # for it; EUCAST expects only ampicillin/amoxicillin/ticarcillin. The
        # gene annotates EUCAST's rows and adds none: the same
        # rule keeps blaPDC from predicting ceftazidime in P. aeruginosa.
        amr_sample(self.root, "K", [("blaSHV-1_1_AF148850", 100.0, 100.0, 80.0, 900)],
                   1_000_000, self.settings)
        id_summary(self.root, "K", [("Klebsiella pneumoniae", "rmlst", 30.0)])
        finalise_amr(self.settings, self.root, "RUN1", controls=[])
        call = [r for r in read_csv(self.root / "amr_summary_all.csv") if r["sample"] == "K"][0]
        self.assertEqual((call["status"], call["host_verdict"]), ("report", "intrinsic"))
        pheno = {r["agent"]: r for r in read_csv(self.root / "predicted_phenotype_all.csv")
                 if r["sample"] == "K"}
        self.assertNotIn("piperacillin", pheno)
        self.assertEqual(pheno["ampicillin"]["prediction"], "resistant (EUCAST expected phenotype)")
        self.assertEqual(pheno["ampicillin"]["intrinsic_gene_evidence"], "blaSHV-1")
        self.assertFalse(any("intrinsic determinant" in r["prediction"] for r in pheno.values()))

    def test_families_come_from_todays_overrides_not_the_arm_run(self):
        # A results tree written with catB7 under the family "catB" must pick
        # up a later catB7 override through amr-report alone.
        amr_sample(self.root, "P", [("catB7_1_AF036933", 100.0, 100.0, 60.0, 800)],
                   1_000_000, self.settings)
        hits = self.root / "P" / "amr" / "amr_hits.csv"
        rows = read_csv(hits)
        for row in rows:
            row["family"] = "catB"            # as written before the override existed
        import csv as _csv
        with open(hits, "w", newline="", encoding="utf-8") as handle:
            writer = _csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        id_summary(self.root, "P", [("Pseudomonas aeruginosa", "rmlst", 30.0)])
        finalise_amr(self.settings, self.root, "RUN1", controls=[])
        call = [r for r in read_csv(self.root / "amr_summary_all.csv") if r["sample"] == "P"][0]
        self.assertEqual((call["family"], call["host_verdict"]), ("catB7", "intrinsic"))

    def test_a_split_gene_is_one_call_with_all_its_reads(self):
        # P. aeruginosa's aph(3')-IIb, filed by KMA under two submissions
        amr_sample(self.root, "P", [("aph(3')-IIb_1_X90856", 99.5, 100.0, 300.0, 90_852),
                                    ("aph(3')-IIb_2_AY000001", 98.0, 100.0, 80.0, 23_753)],
                   1_000_000, self.settings)
        finalise_amr(self.settings, self.root, "RUN1", controls=[])
        rows = [r for r in read_csv(self.root / "amr_summary_all.csv") if r["sample"] == "P"]
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["template"], "aph(3')-IIb_1_X90856")
        self.assertEqual(rows[0]["fragments"], "114605")
        self.assertIn("templates_merged:2", rows[0]["flags"])
        self.assertNotIn("family_multiple_alleles", rows[0]["flags"])

    def test_without_a_control_every_call_is_ungated(self):
        amr_sample(self.root, "A", [TEM], 1_000_000, self.settings)
        result = finalise_amr(self.settings, self.root, "RUN1", controls=[])
        self.assertIsNone(result.baseline_path)
        row = read_csv(self.root / "amr_summary_all.csv")[0]
        self.assertEqual(row["gate_band"], "ungated")
        self.assertIn("amr_gate_not_applied", row["flags"])
        self.assertEqual(row["status"], "report")

    def test_a_cross_contaminated_gene_is_gated_on_absolute_reads(self):
        # A carries CTX-M-15 in bulk; the blank's CTX-M-15 came from A.
        amr_sample(self.root, "Neg_1", [("blaCTX-M-15_1_AY044436", 100.0, 100.0, 3.0, 30)],
                   20_000, self.settings)
        amr_sample(self.root, "A", [("blaCTX-M-15_1_AY044436", 100.0, 100.0, 80.0, 90_000)],
                   2_000_000, self.settings)
        finalise_amr(self.settings, self.root, "RUN1", controls=["Neg_1"])
        base = read_csv(self.root / "amr_baselines.csv")
        ctx = [row for row in base if row["organism"] == "blaCTX-M-15"][0]
        self.assertEqual((ctx["mechanism"], ctx["source_sample"]), ("cross_contamination", "A"))
        row = [r for r in read_csv(self.root / "amr_summary_all.csv") if r["sample"] == "A"][0]
        self.assertEqual(row["gate_mechanism"], "cross_contamination")
        self.assertEqual(row["status"], "report")


class WiringTests(unittest.TestCase):
    def test_arms_parsing_and_rendering(self):
        from bridge.settings import ARMS, ID_ARMS, arms_arg, arms_from
        self.assertEqual(arms_from("all"), ARMS)
        self.assertEqual(arms_from(""), ARMS)
        self.assertEqual(arms_from("both"), ID_ARMS)
        self.assertEqual(arms_from("amr"), ("amr",))
        self.assertEqual(arms_from("amr,rmlst"), ("rmlst", "amr"))
        with self.assertRaises(ValueError):
            arms_from("rmlst,blast")
        for value in ("all", "both", "amr", "16s", "rmlst,amr"):
            self.assertEqual(arms_from(arms_arg(arms_from(value))), arms_from(value))
        self.assertEqual(arms_arg(ARMS), "all")
        self.assertEqual(arms_arg(("16s", "amr")), "16s,amr")

    def test_amr_only_submission_leaves_the_id_arms_alone(self):
        from bridge import pbs
        from bridge.reads import ReadPair
        from bridge.settings import load_settings
        tmp = Path(tempfile.mkdtemp())
        out = tmp / "results"
        for rel in ("rmlst/species_calls.csv", "16s/marker16s_calls.csv", "amr/stale.txt"):
            path = out / "S1" / rel
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("x")
        r1, r2 = tmp / "S1_R1_001.fastq.gz", tmp / "S1_R2_001.fastq.gz"
        r1.write_bytes(b"a")
        r2.write_bytes(b"b")
        settings = load_settings()
        with patch.object(pbs, "qsub", return_value="1.pbs"), \
             patch.object(pbs, "active_jobs", return_value=set()):
            subs = pbs.submit_samples(settings, [ReadPair("S1", r1, r2)], out, "RUN1",
                                      arms=("amr",))
        self.assertEqual(subs[0].action, "submitted")
        self.assertIn("--arms amr", (out / "oe" / "scripts" / "S1.pbs").read_text())
        self.assertTrue((out / "S1" / "rmlst" / "species_calls.csv").is_file())
        self.assertFalse((out / "S1" / "amr").exists())
        # and a sample whose AMR output exists is complete for --arms amr
        (out / "S1" / "amr").mkdir()
        (out / "S1" / "amr" / "amr_hits.csv").write_text("x")
        self.assertTrue(pbs.arms_complete(out / "S1", ("amr",)))

    def test_run_sample_runs_the_amr_arm_and_isolates_its_failure(self):
        from bridge.reads import ReadPair
        from bridge.sample import ArmFailed, run_sample
        from bridge.settings import load_settings
        tmp = Path(tempfile.mkdtemp())
        pair = ReadPair("S1", tmp / "R1.fq.gz", tmp / "R2.fq.gz")
        settings = load_settings()
        with patch("bridge.amr.kma.run_arm", return_value=tmp / "amr_hits.csv") as run:
            produced = run_sample(settings, pair, tmp / "out", arms=("amr",), threads=3)
        self.assertEqual(produced["amr"], tmp / "amr_hits.csv")
        self.assertEqual(run.call_args.kwargs["threads"], 3)
        self.assertEqual(run.call_args.args[4], tmp / "out" / "S1" / "amr")
        with patch("bridge.amr.kma.run_arm", side_effect=RuntimeError("kma died")):
            with self.assertRaises(ArmFailed) as caught:
                run_sample(settings, pair, tmp / "out", arms=("amr",))
        self.assertIn("amr arm", str(caught.exception))

    def test_finalise_runs_the_amr_report_after_reconciliation(self):
        from bridge.finalise import describe, finalise_run
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
        # an AMR-only sample and an AMR-only control must not break the ID steps
        amr_sample(root, "B", [TEM], 1_000_000, amr_settings)
        amr_sample(root, "Neg_1", [TEM], 20_000, amr_settings)

        with patch("bridge.amr.settings.from_config", return_value=amr_settings):
            result = finalise_run(load_settings(), root, "RUN1")
        self.assertEqual(result.reconciled, ["A"])
        self.assertIn(("B", "AMR arm only (no identification)"), result.skipped)
        self.assertIsNone(result.amr_error)
        self.assertEqual(sorted(result.amr.reported), ["A", "B"])
        ctx = [r for r in read_csv(root / "amr_summary_all.csv") if r["sample"] == "A"][0]
        # attribution read the organism the reconciliation step had just written
        self.assertEqual((ctx["attribution"], ctx["host_organism"]),
                         ("sole_organism", "Escherichia coli"))
        self.assertIn("amr_summary_all.csv", describe(result, root))

    def test_an_amr_failure_does_not_undo_identification(self):
        from bridge.finalise import finalise_run
        from bridge.settings import load_settings
        tmp = Path(tempfile.mkdtemp())
        root = tmp / "results"
        amr_sample(root, "B", [TEM], 1_000_000, settings_with_db(tmp))
        with patch("bridge.amr.report.finalise_amr", side_effect=RuntimeError("boom")):
            result = finalise_run(load_settings(), root, "RUN1")
        self.assertIn("boom", result.amr_error)
        self.assertTrue((root / "identification_summary_all.csv").is_file())

    def test_check_reports_the_kma_index_prefix(self):
        import io
        from contextlib import redirect_stdout
        from bridge.check import run_check
        from bridge.settings import load_settings
        tmp = Path(tempfile.mkdtemp())
        settings = load_settings()
        prefix = tmp / "kma" / "resfinder"
        prefix.parent.mkdir()
        phenotypes = tmp / "phenotypes.txt"
        phenotypes.write_text("x")
        settings = replace(settings, references=replace(
            settings.references, amr_db=prefix, amr_phenotypes=phenotypes))
        buffer = io.StringIO()
        with redirect_stdout(buffer):
            run_check(settings, arms=("amr",), scheduler="local")
        self.assertIn("no KMA index", buffer.getvalue())
        Path(f"{prefix}.name").write_text("x")
        buffer = io.StringIO()
        with redirect_stdout(buffer):
            run_check(settings, arms=("amr",), scheduler="local")
        text = buffer.getvalue()
        self.assertRegex(text, r"ok\s+amr_db")
        self.assertIn("kma", text)


if __name__ == "__main__":
    unittest.main()
