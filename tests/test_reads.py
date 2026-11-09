"""Read lists and read-pair discovery (bridge.reads)."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from bridge.reads import (
    ReadPair, ReadProblem, discover_reads, mate_basename, parse_read_list,
)


class MateNameTests(unittest.TestCase):
    def test_both_naming_conventions(self):
        self.assertEqual(mate_basename("SAMPLE-086X_S59_L001_R1_001.fastq.gz"),
                         "SAMPLE-086X_S59_L001_R2_001.fastq.gz")
        self.assertEqual(mate_basename("MIX-170x292_80-20_R1.fastq.gz"),
                         "MIX-170x292_80-20_R2.fastq.gz")
        self.assertIsNone(mate_basename("sample_1.fastq.gz"))

    def test_only_the_first_token_is_replaced(self):
        self.assertEqual(mate_basename("A_R1_x_R1_.fq.gz"), "A_R2_x_R1_.fq.gz")


class ReadListTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def make(self, sample: str, r1_name: str, with_mate: bool = True,
             nested: bool = True, size: int = 10) -> Path:
        folder = self.tmp / sample if nested else self.tmp
        folder.mkdir(parents=True, exist_ok=True)
        r1 = folder / r1_name
        r1.write_bytes(b"x" * size)
        if with_mate:
            (folder / mate_basename(r1_name)).write_bytes(b"y" * size)
        return r1

    def test_list_with_comments_crlf_duplicates_and_missing_mates(self):
        a = self.make("A", "A_S1_L001_R1_001.fastq.gz")
        b = self.make("B", "B_R1.fastq.gz")
        no_mate = self.make("C", "C_R1.fastq.gz", with_mate=False)
        no_token = self.make("D", "D_1.fastq.gz", with_mate=False)
        listing = self.tmp / "sample.list"
        listing.write_bytes(
            (f"# a comment\r\n{a}\r\n\r\n{b}\r\n{a}\r\n{no_mate}\r\n{no_token}\r\n"
             f"{self.tmp / 'E' / 'E_R1.fastq.gz'}\r\n").encode())
        pairs, problems = parse_read_list(listing)
        self.assertEqual([p.sample for p in pairs], ["A", "B"])
        self.assertEqual(pairs[0], ReadPair("A", a, a.with_name("A_S1_L001_R2_001.fastq.gz")))
        self.assertEqual(pairs[0].size_bytes(), 20)
        reasons = {p.sample: p.reason for p in problems}
        self.assertIn("more than once", reasons["A"])
        self.assertIn("R2 mate not found", reasons["C"])
        self.assertIn("no '_R1_' or '_R1.'", reasons["D"])
        self.assertIn("R1 not found", reasons["E"])
        self.assertTrue(all(isinstance(p, ReadProblem) for p in problems))

    def test_discover_nested_layout(self):
        self.make("S2", "S2_S2_L001_R1_001.fastq.gz")
        self.make("S1", "S1_R1.fq.gz")
        pairs, problems = discover_reads(self.tmp)
        self.assertEqual([p.sample for p in pairs], ["S1", "S2"])
        self.assertEqual(problems, [])

    def test_discover_flat_layout_names_by_basename(self):
        self.make("", "X_S9_L001_R1_001.fastq.gz", nested=False)
        self.make("", "Y_R1.fastq.gz", nested=False)
        self.make("", "MIX-170x292_80-20_R1.fastq.gz", nested=False)
        pairs, _ = discover_reads(self.tmp)
        self.assertEqual([p.sample for p in pairs], ["MIX-170x292_80-20", "X", "Y"])

    def test_flat_sample_name_strips_illumina_suffixes(self):
        from bridge.reads import flat_sample_name
        self.assertEqual(flat_sample_name("SAMPLE-086X_S59_L001_R1_001.fastq.gz"), "SAMPLE-086X")
        self.assertEqual(flat_sample_name("A_S1_R1.fq.gz"), "A")
        self.assertEqual(flat_sample_name("B_L001_R1_001.fq.gz"), "B")


if __name__ == "__main__":
    unittest.main()
