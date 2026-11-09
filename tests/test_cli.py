"""Argument parsing for the bridge command.

The PBS scripts write ``bridge sample --config ... --sample-id ...``, with
the option AFTER the subcommand. argparse accepts a top-level option only
BEFORE it, which would fail every job with "unrecognized arguments:
--config". These tests pin both positions.
"""

from __future__ import annotations

import unittest
from pathlib import Path

from bridge.cli import _config_path, build_parser
from bridge.settings import DEFAULT_CONFIG

CONFIG = Path("/some/where/config.yaml")


class ConfigPositionTests(unittest.TestCase):

    def parse(self, argv):
        return build_parser().parse_args(argv)

    def test_config_after_each_subcommand(self):
        cases = {
            "sample": ["sample", "--config", str(CONFIG), "--sample-id", "S",
                       "--reads", "r1.fq", "r2.fq", "--out", "o"],
            "run": ["run", "--config", str(CONFIG), "--reads-list", "l",
                    "--out", "o", "--run-id", "R"],
            "finalise": ["finalise", "--config", str(CONFIG), "--out", "o", "--run-id", "R"],
            "check": ["check", "--config", str(CONFIG)],
        }
        for name, argv in cases.items():
            with self.subTest(subcommand=name):
                args = self.parse(argv)
                self.assertEqual(_config_path(args), CONFIG)

    def test_config_before_the_subcommand_still_works(self):
        args = self.parse(["--config", str(CONFIG), "finalise", "--out", "o", "--run-id", "R"])
        self.assertEqual(_config_path(args), CONFIG)

    def test_the_subcommand_default_does_not_erase_the_earlier_one(self):
        """A subparser copies its defaults over the parent's namespace, so a
        shared dest would replace a path given in front of the subcommand."""
        args = self.parse(["--config", str(CONFIG), "sample", "--sample-id", "S",
                           "--reads", "a", "b", "--out", "o"])
        self.assertIsNone(args.config_after)
        self.assertEqual(_config_path(args), CONFIG)

    def test_no_config_falls_back_to_the_repository_copy(self):
        args = self.parse(["check"])
        self.assertIsNone(_config_path(args))
        from bridge.settings import load_settings
        self.assertEqual(load_settings(_config_path(args)).path, DEFAULT_CONFIG)

    def test_the_rendered_job_script_uses_a_form_that_parses(self):
        """Whatever the template writes must survive this parser."""
        import re
        from importlib import resources
        template = resources.files("bridge").joinpath("templates/sample.pbs").read_text(encoding="utf-8")
        body = template.split("{bridge_command}", 1)[1]
        flags = re.findall(r"--[a-z-]+", body)
        self.assertEqual(flags[0], "--config")          # the case that broke
        args = self.parse(["sample", "--config", str(CONFIG), "--sample-id", "S",
                           "--reads", "a", "b", "--out", "o", "--arms", "both",
                           "--threads", "16", "--memory-gb", "64", "--tempdir", "t"])
        self.assertEqual(_config_path(args), CONFIG)
        self.assertEqual(args.arms, "both")


class PreflightTests(unittest.TestCase):
    """`run` checks the environment before submitting, not once per job."""

    def test_a_failing_check_submits_nothing(self):
        import tempfile
        from unittest.mock import patch
        from bridge import cli
        from bridge.reads import ReadPair

        tmp = Path(tempfile.mkdtemp())
        r1, r2 = tmp / "S_R1.fq", tmp / "S_R2.fq"
        r1.write_text("x")
        r2.write_text("x")
        pair = ReadPair("S", r1, r2)

        args = build_parser().parse_args(
            ["run", "--reads-list", str(tmp / "l"), "--out", str(tmp / "o"), "--run-id", "R"])
        self.assertFalse(args.skip_check)
        with patch("bridge.reads.parse_read_list", return_value=([pair], [])), \
             patch("bridge.check.run_check", return_value=1) as check, \
             patch("bridge.pbs.submit_samples") as submit, \
             patch("bridge.pbs.prepare_blast_cache") as blast:
            code = cli.cmd_run(args)
        self.assertEqual(code, 1)
        check.assert_called_once()
        submit.assert_not_called()
        blast.assert_not_called()

    def test_skip_check_is_available_as_an_override(self):
        args = build_parser().parse_args(
            ["run", "--reads-list", "l", "--out", "o", "--run-id", "R", "--skip-check"])
        self.assertTrue(args.skip_check)


if __name__ == "__main__":
    unittest.main()
