from contextlib import redirect_stderr, redirect_stdout
from io import StringIO
from pathlib import Path
import os
import sys
import unittest
from unittest.mock import patch

from bridge.rmlst.cli import main, parser
from bridge.rmlst.dependencies import DependencyError, check_dependencies


PYTHON_VERSIONS = {
    "pandas": "2.2.0",
    "numpy": "1.26.0",
    "matplotlib": "3.8.0",
    "seaborn": "0.13.2",
}


class DependencyTests(unittest.TestCase):
    def test_check_deps_does_not_require_run_inputs(self):
        args = parser().parse_args(["--check-deps"])
        self.assertTrue(args.check_deps)
        self.assertIsNone(args.allele_target)
        self.assertIsNone(args.output_folder)

    def test_cliquesnv_default_can_be_set_to_a_jar_by_environment(self):
        jar = "/data/apps/clique-snv/2.0.3/clique-snv.jar"
        with patch.dict(os.environ, {"BRIDGE_CLIQUESNV": jar}):
            args = parser().parse_args(["--check-deps"])
        self.assertEqual(args.cliquesnv, jar)

    def test_default_route_checks_minimap2_cliquesnv_and_java(self):
        def resolve(command: str) -> str:
            return f"/tools/{Path(command).name}"

        with (
            patch(
                "bridge.rmlst.dependencies.metadata.version",
                side_effect=lambda package: PYTHON_VERSIONS[package],
            ),
            patch("bridge.rmlst.dependencies.executable", side_effect=resolve) as resolver,
            patch(
                "bridge.rmlst.dependencies.companion_executable",
                return_value="/tools/makeblastdb",
            ),
        ):
            statuses = check_dependencies()

        checked = [call.args[0] for call in resolver.call_args_list]
        self.assertIn("minimap2", checked)
        self.assertNotIn("bwa-mem2", checked)
        self.assertIn("cliquesnv", checked)
        self.assertIn("java", checked)
        self.assertIn("makeblastdb", {status.name for status in statuses})

    def test_viralconsensus_only_does_not_require_cliquesnv_or_java(self):
        with (
            patch(
                "bridge.rmlst.dependencies.metadata.version",
                side_effect=lambda package: PYTHON_VERSIONS[package],
            ),
            patch(
                "bridge.rmlst.dependencies.executable",
                side_effect=lambda command: f"/tools/{command}",
            ) as resolver,
            patch(
                "bridge.rmlst.dependencies.companion_executable",
                return_value="/tools/makeblastdb",
            ),
        ):
            check_dependencies(viralconsensus_only=True)

        checked = [call.args[0] for call in resolver.call_args_list]
        self.assertNotIn("cliquesnv", checked)
        self.assertNotIn("java", checked)

    def test_missing_dependencies_are_reported_together(self):
        def resolve(command: str) -> str:
            if command in {"mash", "java"}:
                raise FileNotFoundError(f"Executable not found: {command}")
            return f"/tools/{command}"

        with (
            patch(
                "bridge.rmlst.dependencies.metadata.version",
                side_effect=lambda package: PYTHON_VERSIONS[package],
            ),
            patch("bridge.rmlst.dependencies.executable", side_effect=resolve),
            patch(
                "bridge.rmlst.dependencies.companion_executable",
                side_effect=FileNotFoundError("Executable not found: makeblastdb"),
            ),
        ):
            with self.assertRaises(DependencyError) as captured:
                check_dependencies(cliquesnv="/tools/clique-snv.jar")

        message = str(captured.exception)
        self.assertIn("mash", message)
        self.assertIn("java", message)
        self.assertIn("makeblastdb", message)

    def test_check_deps_cli_prints_success(self):
        output = StringIO()
        with (
            patch.object(sys, "argv", ["bridge.rmlst", "--check-deps"]),
            patch("bridge.rmlst.cli.check_dependencies", return_value=()),
            redirect_stdout(output),
        ):
            main()
        self.assertIn("Dependency check passed", output.getvalue())

    def test_check_deps_cli_exits_nonzero_on_failure(self):
        error = DependencyError(("mash (mash): missing",))
        stderr = StringIO()
        with (
            patch.object(sys, "argv", ["bridge.rmlst", "--check-deps"]),
            patch("bridge.rmlst.cli.check_dependencies", side_effect=error),
            redirect_stderr(stderr),
            self.assertRaises(SystemExit) as captured,
        ):
            main()
        self.assertEqual(captured.exception.code, 1)
        self.assertIn("mash", stderr.getvalue())


if __name__ == "__main__":
    unittest.main()
