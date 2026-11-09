"""Building the candidate-allele index (bridge.rmlst.kallisto).

`kallisto index` can die with SIGSEGV on malformed input, and the input
usually lives in per-job scratch that is deleted with the job. These tests pin
three behaviours: unusable sequences are dropped before kallisto sees them, a
crash is retried once, and a repeat crash keeps its input.
"""

from __future__ import annotations

import subprocess
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from bridge.rmlst import kallisto
from bridge.rmlst.errors import InsufficientMaterial
from bridge.rmlst.models import Config, FastaRecord


def allele(name: str, length: int) -> FastaRecord:
    return FastaRecord(name, name, "ACGT" * (length // 4) + "A" * (length % 4))


class CandidateIndexTests(unittest.TestCase):

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.config = Config(
            allele_target=self.tmp / "ref.fasta", seqdef=self.tmp / "seqdef.txt",
            reads=(self.tmp / "r1.fq", self.tmp / "r2.fq"), sample_id="S",
            output=self.tmp / "out", tempdir=self.tmp / "work", threads=4,
            kallisto="kallisto")
        self.config.output.mkdir(parents=True, exist_ok=True)

    def frame(self, names):
        return pd.DataFrame({"Template": names})

    def test_sequences_shorter_than_k_are_dropped_before_kallisto_sees_them(self):
        records = {"BACT000001_1": allele("BACT000001_1", 400),
                   "BACT000001_2": allele("BACT000001_2", 0),    # bare header
                   "BACT000001_3": allele("BACT000001_3", 12)}   # shorter than k
        written = {}

        def fake_run(command, tool, *, stdout_path=None):
            if command[1] == "index":
                written["fasta"] = Path(command[-1]).read_text()
            elif stdout_path is not None:
                Path(stdout_path).parent.mkdir(parents=True, exist_ok=True)
                Path(stdout_path).write_bytes(b"BAM")

        with patch.object(kallisto, "executable", return_value="kallisto"), \
             patch.object(kallisto, "run", side_effect=fake_run), \
             patch.object(kallisto, "_bam_has_eof", return_value=True), \
             self.assertLogs("kallisto", level="WARNING") as logs:
            kallisto.run_candidate_pseudobam(
                self.config, self.frame(list(records)), records, index_kmer_size=31)

        self.assertIn("BACT000001_1", written["fasta"])
        self.assertNotIn("BACT000001_2", written["fasta"])
        self.assertNotIn("BACT000001_3", written["fasta"])
        self.assertIn("candidate_alleles_status=dropped", "".join(logs.output))

    def test_nothing_usable_is_insufficient_material_not_a_crash(self):
        records = {"BACT000001_1": allele("BACT000001_1", 0)}
        with patch.object(kallisto, "executable", return_value="kallisto"), \
             patch.object(kallisto, "run") as run, \
             self.assertRaises(InsufficientMaterial):
            kallisto.run_candidate_pseudobam(
                self.config, self.frame(list(records)), records, index_kmer_size=31)
        run.assert_not_called()          # kallisto is never handed an empty file

    def test_a_crash_is_retried_once_and_then_succeeds(self):
        records = {"BACT000001_1": allele("BACT000001_1", 400)}
        attempts = []

        def flaky(command, tool, *, stdout_path=None):
            if command[1] == "index":
                attempts.append(1)
                if len(attempts) == 1:
                    raise subprocess.CalledProcessError(-11, command)
            elif stdout_path is not None:
                Path(stdout_path).parent.mkdir(parents=True, exist_ok=True)
                Path(stdout_path).write_bytes(b"BAM")

        with patch.object(kallisto, "executable", return_value="kallisto"), \
             patch.object(kallisto, "run", side_effect=flaky), \
             patch.object(kallisto, "_bam_has_eof", return_value=True), \
             self.assertLogs("kallisto", level="WARNING"):
            kallisto.run_candidate_pseudobam(
                self.config, self.frame(list(records)), records, index_kmer_size=31)
        self.assertEqual(len(attempts), 2)

    def test_a_repeated_crash_keeps_the_input_and_says_where(self):
        records = {"BACT000001_1": allele("BACT000001_1", 400)}

        def always_segfault(command, tool, *, stdout_path=None):
            raise subprocess.CalledProcessError(-11, command)

        with patch.object(kallisto, "executable", return_value="kallisto"), \
             patch.object(kallisto, "run", side_effect=always_segfault), \
             self.assertLogs("kallisto", level="WARNING"):
            with self.assertRaises(RuntimeError) as raised:
                kallisto.run_candidate_pseudobam(
                    self.config, self.frame(list(records)), records, index_kmer_size=31)

        kept = self.config.output / "failed_candidate_alleles.fasta"
        self.assertTrue(kept.is_file())
        self.assertIn("BACT000001_1", kept.read_text())
        self.assertIn(str(kept), str(raised.exception))
        self.assertIn("1 candidate alleles", str(raised.exception))
        # ...and it is NOT reported as a negative result, which would turn a
        # sample with real signal into a silent false negative.
        self.assertNotIsInstance(raised.exception, InsufficientMaterial)


if __name__ == "__main__":
    unittest.main()
