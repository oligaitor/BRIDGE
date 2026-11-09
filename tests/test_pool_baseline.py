"""bridge pool-baseline, and finalise --pooled-baseline (bridge.pool / finalise)."""

from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path

from bridge.finalise import finalise_run, load_pooled_rows
from bridge.pool import main as pool_main, merge_rows, pool_baselines
from bridge.settings import load_settings

from test_finalise import RMLST_HEADER, S16_HEADER, rmlst_row, s16_row


def make_sample(root: Path, name: str, rmlst_rows, s16_rows, total=1_000_000) -> None:
    folder = root / name
    (folder / "16s" / "quant").mkdir(parents=True)
    (folder / "rmlst").mkdir(parents=True)
    (folder / "rmlst" / "species_calls.csv").write_text(RMLST_HEADER + "".join(rmlst_rows))
    (folder / "16s" / "marker16s_calls.csv").write_text(S16_HEADER + "".join(s16_rows))
    (folder / "16s" / "quant" / "run_info.json").write_text(json.dumps({"n_processed": total}))


def read_csv(path: Path) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


class PoolTests(unittest.TestCase):
    def setUp(self):
        self.settings = load_settings()
        base = Path(tempfile.mkdtemp())
        self.run_a, self.run_b = base / "A", base / "B"
        self.out = base / "pooled.csv"

        # Run A: a deep blank with reagent Moraxella, and E. coli carried over
        # from a strong positive on the same run.
        make_sample(self.run_a, "POS", [rmlst_row("POS", "Escherichia coli", est=500_000.0)],
                    [s16_row("POS", "Escherichia coli", "Escherichia", "coli",
                             "species_strong", raw=700_000.0)])
        make_sample(self.run_a, "NEGMATRIX_1",
                    [rmlst_row("NEGMATRIX_1", "Escherichia coli", reported=False, est=20.0)],
                    [s16_row("NEGMATRIX_1", "Moraxella", "Moraxella", "", "genus_only", raw=500.0),
                     s16_row("NEGMATRIX_1", "Escherichia", "Escherichia", "", "genus_only",
                             raw=70.0)],
                    total=20_000)
        # Run B: a blank far too shallow to use, and a pure fungal dilution
        # that carries reagent Bradyrhizobium and Moraxella at a lower level.
        make_sample(self.run_b, "Neg_B",
                    [rmlst_row("Neg_B", "Escherichia coli", reported=False, est=2.0)],
                    [s16_row("Neg_B", "Moraxella", "Moraxella", "", "genus_only", raw=40.0)],
                    total=150)
        make_sample(self.run_b, "Cglab_D1",
                    [rmlst_row("Cglab_D1", "Escherichia coli", reported=False, est=2.0)],
                    [s16_row("Cglab_D1", "Bradyrhizobium", "Bradyrhizobium", "", "genus_only",
                             raw=4_000.0),
                     s16_row("Cglab_D1", "Moraxella", "Moraxella", "", "genus_only", raw=100.0)],
                    total=2_000_000)

    def pool(self):
        return pool_baselines(self.settings, [("RUNA", self.run_a), ("RUNB", self.run_b)],
                              extra_patterns=["cglab*"])

    def test_reagent_rows_are_pooled_at_their_maximum(self):
        result = self.pool()
        rows = {(row["arm"], row["organism"]): row for row in result.rows}
        self.assertTrue(all(row["run_id"] == "*" for row in result.rows))
        moraxella = rows[("marker16s", "Moraxella")]
        # 500 / 20,000 = 25,000 RPM on run A beats 100 / 2 M = 50 RPM on run B
        self.assertAlmostEqual(float(moraxella["rpm_ntc"]), 25_000.0)
        self.assertEqual(moraxella["control_ids"], "RUNA/NEGMATRIX_1; RUNB/Cglab_D1")
        self.assertEqual(moraxella["mechanism"], "reagent")
        self.assertIn("highest on RUNA", moraxella["note"])
        # the fungal dilution contributes its own reagent organism
        self.assertIn(("marker16s", "Bradyrhizobium"), rows)

    def test_cross_contamination_stays_with_its_run(self):
        result = self.pool()
        self.assertGreater(result.dropped_cross, 0)
        # Run A's blank E. coli came from POS on run A: not pooled. (Run B's
        # fungal dilution holds a little E. coli with no source on its run, so
        # that one is reagent and may be.)
        self.assertFalse([row for row in result.rows
                          if row["organism"].lower().startswith("escherichia")
                          and "RUNA/" in row["control_ids"]])

    def test_shallow_controls_are_listed_but_not_used(self):
        result = self.pool()
        used = {(row["run_id"], row["control"]): row["used"] for row in result.controls}
        self.assertEqual(used[("RUNA", "NEGMATRIX_1")], "yes")
        self.assertEqual(used[("RUNB", "Neg_B")], "no (too shallow)")
        self.assertEqual(used[("RUNB", "Cglab_D1")], "yes")
        self.assertEqual(result.runs_without_usable_control, [])

    def test_command_writes_the_baseline_and_the_control_list(self):
        code = pool_main(self.settings, ["--run", f"RUNA={self.run_a}",
                                         "--run", f"RUNB={self.run_b}",
                                         "--extra-control", "Cglab*", "-o", str(self.out)])
        self.assertEqual(code, 0)
        self.assertTrue(read_csv(self.out))
        controls = read_csv(self.out.with_name("pooled_controls.csv"))
        self.assertEqual({row["control"] for row in controls},
                         {"NEGMATRIX_1", "Neg_B", "Cglab_D1"})

    def test_merge_keeps_the_highest_row(self):
        low = {"arm": "amr", "organism": "blaTEM-1A", "rpm_ntc": "10", "reads_ntc": "1",
               "total_reads_ntc": "100000", "control_ids": "N1"}
        high = {**low, "rpm_ntc": "900", "reads_ntc": "45", "total_reads_ntc": "50000",
                "control_ids": "N2"}
        [row] = merge_rows([("R1", low), ("R2", high)])
        self.assertEqual((row["rpm_ntc"], row["reads_ntc"]), ("900", "45"))
        self.assertEqual(row["n_controls"], 2)


class FinaliseWithPoolTests(unittest.TestCase):
    def setUp(self):
        self.settings = load_settings()
        self.base = Path(tempfile.mkdtemp())
        self.pooled = self.base / "pooled.csv"
        with open(self.pooled, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["run_id", "arm", "organism", "rpm_ntc", "reads_ntc",
                             "total_reads_ntc", "n_controls", "control_ids", "source",
                             "mechanism", "source_sample", "source_reads",
                             "ntc_over_source_frac", "note"])
            writer.writerow(["*", "marker16s", "Moraxella", "25000", "500", "20000", "1",
                             "RUNA/NEG_1", "pooled_max_across_runs", "reagent", "", "", "", ""])
            writer.writerow(["*", "marker16s", "Escherichia", "50000", "1000", "20000", "1",
                             "RUNA/NEG_1", "pooled_max_across_runs", "reagent", "", "", "", ""])

    def summary(self, result) -> dict[str, dict[str, str]]:
        return {row["organism"]: row for row in read_csv(result.combined_path)}

    def test_a_run_without_a_blank_is_gated_on_the_pool(self):
        root = self.base / "noblank"
        empty = [rmlst_row("S", "Escherichia coli", reported=False, est=20.0)]
        make_sample(root, "S", empty,
                    [s16_row("S", "Moraxella tetraodonis", "Moraxella", "tetraodonis",
                             "species_strong", raw=100.0)])
        make_sample(root, "K", [rmlst_row("K", "Klebsiella pneumoniae")],
                    [s16_row("K", "Klebsiella pneumoniae", "Klebsiella", "pneumoniae",
                             "species_strong")])
        result = finalise_run(self.settings, root, "RUNX", pooled_baseline=self.pooled)
        self.assertEqual(result.pooled_rows, 2)
        rows = self.summary(result)
        moraxella = rows["Moraxella tetraodonis"]
        self.assertEqual((moraxella["rpm_band"], moraxella["status"]), ("background", "review"))
        self.assertIn("rpm_pooled_baseline", moraxella["flags"])
        klebsiella = rows["Klebsiella pneumoniae"]
        self.assertEqual(klebsiella["rpm_band"], "no_control_background")
        self.assertNotIn("rpm_pooled_baseline", klebsiella["flags"])
        written = read_csv(root / "rpm_baselines.csv")
        self.assertEqual({row["run_id"] for row in written}, {"*"})

    def test_the_runs_own_blank_wins_over_the_pool(self):
        root = self.base / "withblank"
        make_sample(root, "A", [rmlst_row("A", "Escherichia coli")],
                    [s16_row("A", "Escherichia coli", "Escherichia", "coli", "species_strong")])
        # The run's own blank holds a little Escherichia (3,500 RPM); the
        # pool's 50,000 RPM would have put A's 7,000 RPM call under background.
        make_sample(root, "NEG_1",
                    [rmlst_row("NEG_1", "Escherichia coli", reported=False, est=2.0)],
                    [s16_row("NEG_1", "Escherichia", "Escherichia", "", "genus_only", raw=70.0)],
                    total=20_000)
        result = finalise_run(self.settings, root, "RUNY", pooled_baseline=self.pooled)
        row = self.summary(result)["Escherichia coli"]
        self.assertNotIn("rpm_pooled_baseline", row["flags"])
        self.assertEqual({r["run_id"] for r in read_csv(root / "rpm_baselines.csv")},
                         {"RUNY", "*"})

    def test_a_per_run_file_is_refused_as_a_pool(self):
        bad = self.base / "per_run.csv"
        bad.write_text("run_id,arm,organism,rpm_ntc\nRUN1,marker16s,Moraxella,10\n")
        with self.assertRaises(ValueError):
            load_pooled_rows(bad)


if __name__ == "__main__":
    unittest.main()
