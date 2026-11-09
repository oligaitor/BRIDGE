"""config.yaml -> Settings (bridge.settings)."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from bridge.settings import (
    DEFAULT_CONFIG, DEFAULT_CONTROL_PREFIXES, DEFAULT_TIERS, arms_from,
    is_control, load_settings, tier_for,
)

# The smallest config bridge.reconcile.config.load_config accepts: nothing from the
# references / tools / hpc blocks.
LEGACY_CONFIG = """\
taxonomy:
  genus_species_separator: "~"
marker16s:
  genus_min_frac: 0.01
  species_dominance: 0.70
  species_margin: 0.20
rpm:
  enabled: true
  rpm_r_min: 10.0
  rpm_min: 0.05
"""


class RepoConfigTests(unittest.TestCase):
    def setUp(self):
        self.settings = load_settings()

    def test_all_three_new_blocks_are_read(self):
        s = self.settings
        self.assertEqual(s.path, DEFAULT_CONFIG)
        self.assertTrue(s.hpc.project)
        self.assertEqual(s.hpc.python_module, s.cfg.get("hpc.python_module") or None)
        self.assertEqual(s.hpc.controls_prefixes, DEFAULT_CONTROL_PREFIXES)
        self.assertEqual([t.name for t in s.hpc.tiers], ["small", "medium", "large", "xlarge"])
        self.assertIsNone(s.hpc.tiers[-1].max_bytes)
        self.assertEqual(s.hpc.cliquesnv_time_limit_seconds, 600)
        self.assertEqual(s.tools.tools_env, DEFAULT_CONFIG.parent / "tools.env")
        self.assertTrue(s.tools.tools_env.is_file())
        self.assertIsNotNone(s.references.rmlst_alleles)
        self.assertIsNotNone(s.references.s16_copy_numbers)
        self.assertEqual(s.references.required_for(("16s",)), ["s16_reference", "s16_index"])

    def test_tier_boundaries_are_strict_upper_bounds(self):
        hpc = self.settings.hpc
        self.assertEqual(tier_for(hpc, 0).name, "small")
        self.assertEqual(tier_for(hpc, 9_999_999).name, "small")
        self.assertEqual(tier_for(hpc, 10_000_000).name, "medium")
        self.assertEqual(tier_for(hpc, 999_999_999).name, "large")
        self.assertEqual(tier_for(hpc, 1_000_000_000).name, "xlarge")
        self.assertEqual(tier_for(hpc, 50_000_000_000).name, "xlarge")


class LegacyConfigTests(unittest.TestCase):
    def test_config_without_the_new_blocks_still_loads(self):
        path = Path(tempfile.mkdtemp()) / "config.yaml"
        path.write_text(LEGACY_CONFIG, encoding="utf-8")
        s = load_settings(path)
        self.assertEqual(s.hpc.tiers, DEFAULT_TIERS)
        self.assertEqual(s.hpc.controls_prefixes, DEFAULT_CONTROL_PREFIXES)
        self.assertIsNone(s.hpc.python_module)
        self.assertEqual(s.hpc.project, "PROJECT")
        self.assertEqual(s.hpc.storage, "")
        self.assertIsNone(s.tools.tools_env)
        self.assertTrue(all(value is None for _, value in s.references.items()))
        self.assertEqual(s.cfg.get("rpm.rpm_r_min"), 10.0)


class HelperTests(unittest.TestCase):
    def test_arms_from(self):
        self.assertEqual(arms_from("both"), ("rmlst", "16s"))
        self.assertEqual(arms_from("all"), ("rmlst", "16s", "amr"))
        self.assertEqual(arms_from("16s"), ("16s",))
        self.assertEqual(arms_from("rmlst"), ("rmlst",))
        self.assertEqual(arms_from("amr"), ("amr",))
        with self.assertRaises(ValueError):
            arms_from("blast")

    def test_is_control(self):
        prefixes = DEFAULT_CONTROL_PREFIXES
        self.assertTrue(is_control("Neg_20260416", prefixes))
        self.assertTrue(is_control("NTC_1", prefixes))
        self.assertFalse(is_control("SAMPLE-001", prefixes))
        # Case-insensitive, and not tied to an underscore.
        self.assertTrue(is_control("NEGMATRIX_20260129", prefixes))
        self.assertTrue(is_control("neg_lower", prefixes))
        self.assertTrue(is_control("blank_2", prefixes))
        self.assertFalse(is_control("NAI01G-2026", prefixes))
        self.assertFalse(is_control("LIB-26-P067-370", prefixes))


if __name__ == "__main__":
    unittest.main()
