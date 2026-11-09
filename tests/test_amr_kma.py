"""KMA parsing, the command, and one whole arm run with KMA faked (bridge.amr.kma)."""

from __future__ import annotations

import gzip
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from bridge.amr.kma import (
    HIT_COLUMNS, count_read_pairs, kma_command, merge_hits, parse_mapstat, parse_res,
    read_hits, run_arm,
)

from amr_fixtures import settings_with_db, write_kma_output


class ParseTests(unittest.TestCase):
    def setUp(self):
        self.prefix = Path(tempfile.mkdtemp()) / "kma"
        write_kma_output(self.prefix, [
            ("blaTEM-1B_1_AY458016", 100.0, 100.0, 51.2, 600),
            ("blaSHV-12_1_AF462395", 99.77, 100.0, 20.0, 240),
        ])

    def test_res_rows_keyed_by_trimmed_template(self):
        rows = parse_res(Path(f"{self.prefix}.res"))
        self.assertEqual(set(rows), {"blaTEM-1B_1_AY458016", "blaSHV-12_1_AF462395"})
        self.assertEqual(rows["blaSHV-12_1_AF462395"]["template_identity"], "99.77")
        self.assertEqual(rows["blaTEM-1B_1_AY458016"]["depth"], "51.20")

    def test_mapstat_metadata_and_counts(self):
        meta, rows = parse_mapstat(Path(f"{self.prefix}.mapstat"))
        self.assertEqual(meta["version"], "1.6.6")
        self.assertEqual(rows["blaTEM-1B_1_AY458016"]["fragments"], "600")
        self.assertEqual(rows["blaTEM-1B_1_AY458016"]["reads"], "1200")

    def test_merge_parses_names_and_families(self):
        hits = merge_hits("S1", parse_res(Path(f"{self.prefix}.res")),
                          parse_mapstat(Path(f"{self.prefix}.mapstat"))[1])
        by_gene = {hit["gene"]: hit for hit in hits}
        self.assertEqual(by_gene["blaSHV-12"]["family"], "blaSHV")
        self.assertEqual(by_gene["blaTEM-1B"]["fragments"], "600")
        self.assertEqual(set(hits[0]), set(HIT_COLUMNS))

    def test_missing_files_are_empty_not_errors(self):
        self.assertEqual(parse_res(Path("nope.res")), {})
        self.assertEqual(parse_mapstat(Path("nope.mapstat")), ({}, {}))


class CommandTests(unittest.TestCase):
    def test_command_carries_paired_reads_db_threads_and_args(self):
        command = kma_command("kma", Path("R1.fq.gz"), Path("R2.fq.gz"), Path("/db/resfinder"),
                              Path("/out/kma"), 8, "-1t1 -ef -pm f")
        self.assertEqual(command[:4], ["kma", "-ipe", "R1.fq.gz", "R2.fq.gz"])
        self.assertIn("-t_db", command)
        self.assertEqual(command[command.index("-t") + 1], "8")
        self.assertEqual(command[-4:], ["-1t1", "-ef", "-pm", "f"])

    def test_count_read_pairs_gzip_and_plain(self):
        folder = Path(tempfile.mkdtemp())
        record = b"@r\nACGT\n+\nIIII\n"
        with gzip.open(folder / "R1.fastq.gz", "wb") as handle:
            handle.write(record * 7)
        (folder / "R1.fastq").write_bytes(record * 3)
        self.assertEqual(count_read_pairs(folder / "R1.fastq.gz"), 7)
        self.assertEqual(count_read_pairs(folder / "R1.fastq"), 3)


class RunArmTests(unittest.TestCase):
    def test_a_whole_arm_run_with_kma_faked(self):
        tmp = Path(tempfile.mkdtemp())
        settings = settings_with_db(tmp)
        Path(f"{settings.db}.name").parent.mkdir(parents=True, exist_ok=True)
        Path(f"{settings.db}.name").write_text("blaTEM-1B_1_AY458016\n")
        reads = tmp / "reads"
        reads.mkdir()
        with gzip.open(reads / "S1_R1.fastq.gz", "wb") as handle:
            handle.write(b"@r\nACGT\n+\nIIII\n" * 100)
        (reads / "S1_R2.fastq.gz").write_bytes(b"")
        out = tmp / "results" / "S1" / "amr"

        def fake_run(command, tool, **_):
            prefix = Path(command[command.index("-o") + 1])
            write_kma_output(prefix, [("blaTEM-1B_1_AY458016", 100.0, 100.0, 51.2, 60),
                                      ("tet(M)_1_X92947", 70.0, 30.0, 0.4, 2)])

        with patch("bridge.tools.executable", return_value="/fake/kma"), \
             patch("bridge.tools.run", side_effect=fake_run) as run, \
             patch("bridge.amr.kma.kma_version", return_value="KMA-1.6.6"):
            path = run_arm(settings, "S1", reads / "S1_R1.fastq.gz", reads / "S1_R2.fastq.gz",
                           out, threads=2)

        self.assertEqual(run.call_args[0][1], "kma")
        hits = {hit["gene"]: hit for hit in read_hits(path)}
        self.assertEqual(hits["blaTEM-1B"]["tier"], "confident")
        self.assertEqual(hits["blaTEM-1B"]["allele_level"], "exact")
        self.assertIn("Ampicillin", hits["blaTEM-1B"]["antibiotics"])
        self.assertEqual(hits["tet(M)"]["tier"], "below")
        manifest = json.loads((out / "amr_manifest.json").read_text())
        self.assertEqual(manifest["total_read_pairs"], 100)
        self.assertEqual(manifest["kma_version"], "KMA-1.6.6")
        self.assertEqual(manifest["templates_above_floor"], 1)
        self.assertTrue((out / "amr.log").is_file())

    def test_missing_index_is_a_clear_error(self):
        tmp = Path(tempfile.mkdtemp())
        settings = settings_with_db(tmp)
        with self.assertRaises(FileNotFoundError) as caught:
            run_arm(settings, "S1", tmp / "R1", tmp / "R2", tmp / "out")
        self.assertIn("build-amr-db", str(caught.exception))


if __name__ == "__main__":
    unittest.main()
