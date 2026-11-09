"""PBS submission logic (bridge.pbs), with qsub/qselect patched out."""

from __future__ import annotations

import re
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

from bridge import pbs
from bridge.reads import ReadPair
from bridge.settings import ARMS, load_settings, tier_for


class PbsTestCase(unittest.TestCase):
    def setUp(self):
        settings = load_settings()
        # Fixed scheduler values, so the assertions do not depend on config.yaml.
        self.settings = replace(settings, hpc=replace(
            settings.hpc, project="PROJECT", python_module="python3/3.12.1",
            storage="gdata/PROJECT+scratch/PROJECT"))
        self.tmp = Path(tempfile.mkdtemp())
        self.out = self.tmp / "results"
        self.out.mkdir()

    def pair(self, sample: str, size: int = 1024) -> ReadPair:
        folder = self.tmp / "reads" / sample
        folder.mkdir(parents=True, exist_ok=True)
        r1 = folder / f"{sample}_S1_L001_R1_001.fastq.gz"
        r2 = folder / f"{sample}_S1_L001_R2_001.fastq.gz"
        r1.write_bytes(b"a" * size)
        r2.write_bytes(b"b" * size)
        return ReadPair(sample, r1, r2)


class RenderTests(PbsTestCase):
    def test_sample_script_carries_resources_environment_and_command(self):
        pair = self.pair("S1")
        tier = tier_for(self.settings.hpc, 2048)
        text = pbs.render_sample_script(self.settings, pair, self.out, "RUN1", tier,
                                        arms=ARMS, size_bytes=2048)
        for needle in ("#PBS -N S1_RUN1", "#PBS -P PROJECT", "#PBS -l storage=gdata/PROJECT+scratch/PROJECT", "#PBS -q normal",
                       f"#PBS -l ncpus={tier.ncpus}", f"#PBS -l mem={tier.mem_gb}GB",
                       f"#PBS -l jobfs={tier.jobfs_gb}GB",
                       f"#PBS -l walltime={tier.walltime}",
                       "module purge", "module load python3/3.12.1",
                       "tools.env", '--sample-id "S1"', "--arms all",
                       f'--threads "${{PBS_NCPUS:-{tier.ncpus}}}"', "sample \\"):
            self.assertIn(needle, text)
        self.assertTrue("bridge" in text)
        self.assertIn("/bridge_S1", text)
        # every {placeholder} filled; ${PBS_NCPUS:-4} is shell, not a placeholder
        self.assertIsNone(re.search(r"\{[a-z_]+\}", text), text)

    def test_single_arm_and_no_tools_env(self):
        from dataclasses import replace
        settings = replace(self.settings, tools=replace(self.settings.tools, tools_env=None))
        text = pbs.render_sample_script(settings, self.pair("S1"), self.out, "RUN1",
                                        self.settings.hpc.tiers[0], arms=("16s",), size_bytes=1)
        self.assertIn("--arms 16s", text)
        self.assertIn("binaries expected on PATH", text)

    def test_finalise_script_waits_for_every_sample_job(self):
        text = pbs.render_finalise_script(self.settings, self.out, "RUN1",
                                          depends_on=["1.pbs", "2.pbs"],
                                          controls=["Neg_A", "Neg_B"])
        self.assertIn("#PBS -W depend=afterany:1.pbs:2.pbs", text)
        self.assertIn('--control "Neg_A"', text)
        self.assertIn('--control "Neg_B"', text)
        self.assertIn('--run-id "RUN1"', text)
        self.assertIn("#PBS -l ncpus=1", text)
        without = pbs.render_finalise_script(self.settings, self.out, "RUN1",
                                             depends_on=[], controls=[])
        self.assertNotIn("depend=", without)


class SubmitTests(PbsTestCase):
    def test_dry_run_renders_scripts_and_submits_nothing(self):
        with patch.object(pbs, "qsub") as qsub:
            subs = pbs.submit_samples(self.settings, [self.pair("S1")], self.out, "RUN1",
                                      dry_run=True)
        qsub.assert_not_called()
        self.assertEqual(subs[0].action, "dry_run")
        self.assertTrue((self.out / "oe" / "scripts" / "S1.pbs").is_file())
        script, job = pbs.submit_finalise(self.settings, self.out, "RUN1",
                                          depends_on=[], controls=[], dry_run=True)
        self.assertIsNone(job)
        self.assertTrue(script.is_file())

    def test_submission_records_the_job_id_and_clears_a_partial_arm(self):
        pair = self.pair("S1")
        stale = self.out / "S1" / "rmlst" / "partial.txt"
        stale.parent.mkdir(parents=True)
        stale.write_text("half done")
        with patch.object(pbs, "qsub", return_value="123.pbs") as qsub, \
             patch.object(pbs, "active_jobs", return_value=set()):
            subs = pbs.submit_samples(self.settings, [pair], self.out, "RUN1")
        self.assertEqual(qsub.call_count, 1)
        self.assertEqual(subs[0].action, "submitted")
        self.assertEqual(subs[0].job_id, "123.pbs")
        self.assertEqual((self.out / "oe" / ".jobids" / "S1").read_text().strip(), "123.pbs")
        self.assertFalse(stale.exists())

    def test_complete_inflight_deferred_and_castanet_cases(self):
        done, running, extra, castanet = (self.pair(n) for n in ("DONE", "RUN", "EXTRA", "CAST"))
        for rel in ("rmlst/species_calls.csv", "16s/marker16s_calls.csv", "amr/amr_hits.csv"):
            path = self.out / "DONE" / rel
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("x")
        jid = self.out / "oe" / ".jobids"
        jid.mkdir(parents=True)
        (jid / "RUN").write_text("77.pbs\n")
        (self.out / "CAST" / "graphcalls").mkdir(parents=True)
        with patch.object(pbs, "qsub", return_value="9.pbs"), \
             patch.object(pbs, "active_jobs", return_value={"77.pbs"}):
            subs = pbs.submit_samples(self.settings, [done, running, castanet, extra],
                                      self.out, "RUN1", max_submit=1)
        actions = {s.sample: s.action for s in subs}
        self.assertEqual(actions["DONE"], "skipped_complete")
        self.assertEqual(actions["RUN"], "skipped_inflight")
        self.assertEqual(actions["CAST"], "refused")
        self.assertEqual(actions["EXTRA"], "submitted")
        # our own bacterial_ID/ folder is not a Castanet marker
        self.assertIsNone(pbs.castanet_marker(self.out / "DONE"))
        (self.out / "DONE" / "bacterial_ID").mkdir()
        self.assertIsNone(pbs.castanet_marker(self.out / "DONE"))

    def test_max_submit_defers_the_rest(self):
        pairs = [self.pair("A"), self.pair("B"), self.pair("C")]
        with patch.object(pbs, "qsub", side_effect=["1.g", "2.g"]), \
             patch.object(pbs, "active_jobs", return_value=set()):
            subs = pbs.submit_samples(self.settings, pairs, self.out, "RUN1", max_submit=2)
        self.assertEqual([s.action for s in subs], ["submitted", "submitted", "deferred"])

    def test_force_resubmits_a_complete_sample_but_not_a_running_one(self):
        done, running = self.pair("DONE"), self.pair("RUN")
        for name in ("DONE", "RUN"):
            for rel in ("rmlst/species_calls.csv", "16s/marker16s_calls.csv"):
                path = self.out / name / rel
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text("x")
        jid = self.out / "oe" / ".jobids"
        jid.mkdir(parents=True)
        (jid / "RUN").write_text("5.pbs\n")
        with patch.object(pbs, "qsub", return_value="9.pbs"), \
             patch.object(pbs, "active_jobs", return_value={"5.pbs"}):
            subs = pbs.submit_samples(self.settings, [done, running], self.out, "RUN1",
                                      arms=("16s",), force=True)
        actions = {s.sample: s.action for s in subs}
        self.assertEqual(actions, {"DONE": "submitted", "RUN": "skipped_inflight"})
        # only the requested arm's folder was cleared
        self.assertFalse((self.out / "DONE" / "16s").exists())
        self.assertTrue((self.out / "DONE" / "rmlst" / "species_calls.csv").is_file())

    def test_scripts_pin_pythonpath_to_this_repository(self):
        text = pbs.render_sample_script(self.settings, self.pair("S1"), self.out, "RUN1",
                                        self.settings.hpc.tiers[0], arms=ARMS, size_bytes=1)
        self.assertIn("export PYTHONPATH=", text)
        self.assertIn("src", text.split("export PYTHONPATH=")[1].splitlines()[0])
        fin = pbs.render_finalise_script(self.settings, self.out, "RUN1", depends_on=[], controls=[])
        self.assertIn("export PYTHONPATH=", fin)


class BlastCacheTests(PbsTestCase):

    def test_cache_is_built_once_up_front_when_makeblastdb_exists(self):
        from dataclasses import replace
        reference = self.tmp / "alleles.fasta"
        reference.write_text(">a\nACGT\n")
        settings = replace(self.settings, references=replace(self.settings.references,
                                                             rmlst_alleles=reference))
        with patch.object(pbs, "companion_executable", return_value="/tools/makeblastdb"), \
             patch("bridge.rmlst.consensus._ensure_blast_database",
                   return_value=self.out / ".bridge_blastdb" / "reference_abc") as build:
            status = pbs.prepare_blast_cache(settings, self.out)
        self.assertTrue(status.startswith("ready:"))
        build.assert_called_once_with(reference, self.out.resolve(), "/tools/makeblastdb")

    def test_missing_makeblastdb_is_a_warning_not_an_error(self):
        from dataclasses import replace
        reference = self.tmp / "alleles.fasta"
        reference.write_text(">a\nACGT\n")
        settings = replace(self.settings, references=replace(self.settings.references,
                                                             rmlst_alleles=reference))
        with patch.object(pbs, "companion_executable", side_effect=FileNotFoundError("no")):
            status = pbs.prepare_blast_cache(settings, self.out)
        self.assertTrue(status.startswith("WARNING"))

    def test_missing_reference_is_reported(self):
        from dataclasses import replace
        settings = replace(self.settings, references=replace(self.settings.references,
                                                             rmlst_alleles=self.tmp / "nope.fasta"))
        self.assertTrue(pbs.prepare_blast_cache(settings, self.out).startswith("skipped"))


if __name__ == "__main__":
    unittest.main()
