"""Running one sample's two arms (bridge.sample).

The arms are independent, and reconciliation reports a sample from 16S alone,
so a failure in one must not cost the other: a crash in the rMLST arm must not
stop the 16S arm, or the sample is absent from the run table rather than
reported from 16S.
"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from bridge import sample as sample_mod
from bridge.reads import ReadPair
from bridge.sample import ArmFailed, run_sample
from bridge.settings import ID_ARMS, load_settings


class ArmIndependenceTests(unittest.TestCase):

    def setUp(self):
        self.settings = load_settings()
        self.tmp = Path(tempfile.mkdtemp())
        r1, r2 = self.tmp / "S_R1.fq", self.tmp / "S_R2.fq"
        r1.write_text("@r\nACGT\n+\nIIII\n")
        r2.write_text("@r\nACGT\n+\nIIII\n")
        self.pair = ReadPair("S", r1, r2)
        self.out = self.tmp / "out"

    def test_the_16s_arm_still_runs_when_rmlst_fails(self):
        marker16s = []
        with patch.object(sample_mod, "run_pipeline",
                          side_effect=RuntimeError("kallisto index died with SIGSEGV")), \
             patch.object(sample_mod, "marker16s_main",
                          side_effect=lambda argv: marker16s.append(argv) or 0):
            with self.assertRaises(ArmFailed) as raised:
                run_sample(self.settings, self.pair, self.out)
        # the 16S arm was attempted despite the rMLST failure...
        self.assertEqual(len(marker16s), 1)
        self.assertIn("--sample-id", marker16s[0])
        # ...and the job still reports the failure, so PBS marks it failed and
        # the submitter will pick the sample up again.
        self.assertIn("rmlst arm", str(raised.exception))
        self.assertIn("SIGSEGV", str(raised.exception))

    def test_a_16s_failure_does_not_hide_a_successful_rmlst_run(self):
        with patch.object(sample_mod, "run_pipeline"), \
             patch.object(sample_mod, "marker16s_main", return_value=3):
            with self.assertRaises(ArmFailed) as raised:
                run_sample(self.settings, self.pair, self.out)
        self.assertIn("16s arm", str(raised.exception))

    def test_both_arms_failing_names_both(self):
        with patch.object(sample_mod, "run_pipeline", side_effect=RuntimeError("a")), \
             patch.object(sample_mod, "marker16s_main", side_effect=RuntimeError("b")):
            with self.assertRaises(ArmFailed) as raised:
                run_sample(self.settings, self.pair, self.out)
        message = str(raised.exception)
        self.assertIn("rmlst arm", message)
        self.assertIn("16s arm", message)

    def test_a_clean_run_returns_both_paths_and_raises_nothing(self):
        with patch.object(sample_mod, "run_pipeline"), \
             patch.object(sample_mod, "marker16s_main", return_value=0):
            produced = run_sample(self.settings, self.pair, self.out, arms=ID_ARMS)
        self.assertEqual(sorted(produced), ["16s", "rmlst"])
        self.assertEqual(produced["16s"], self.out / "S" / "16s" / "marker16s_calls.csv")

    def test_a_clean_default_run_includes_the_amr_arm(self):
        amr_hits = self.out / "S" / "amr" / "amr_hits.csv"
        with patch.object(sample_mod, "run_pipeline"), \
             patch.object(sample_mod, "marker16s_main", return_value=0), \
             patch("bridge.amr.kma.run_arm", return_value=amr_hits):
            produced = run_sample(self.settings, self.pair, self.out)
        self.assertEqual(sorted(produced), ["16s", "amr", "rmlst"])
        self.assertEqual(produced["amr"], amr_hits)

    def test_an_amr_failure_costs_neither_identification_arm(self):
        with patch.object(sample_mod, "run_pipeline") as rmlst, \
             patch.object(sample_mod, "marker16s_main", return_value=0) as s16, \
             patch("bridge.amr.kma.run_arm", side_effect=RuntimeError("kma died")):
            with self.assertRaises(ArmFailed) as raised:
                run_sample(self.settings, self.pair, self.out)
        rmlst.assert_called_once()
        s16.assert_called_once()
        self.assertIn("amr arm", str(raised.exception))
        self.assertNotIn("rmlst arm", str(raised.exception))

    def test_one_arm_only_does_not_attempt_the_other(self):
        with patch.object(sample_mod, "run_pipeline") as rmlst, \
             patch.object(sample_mod, "marker16s_main", return_value=0) as s16:
            run_sample(self.settings, self.pair, self.out, arms=("16s",))
        rmlst.assert_not_called()
        s16.assert_called_once()


if __name__ == "__main__":
    unittest.main()
