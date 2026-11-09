"""ResFinder names, families, phenotypes and the index build (bridge.amr.database)."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from bridge.amr.database import (
    build_db, concatenate_fastas, family_of, load_family_overrides, load_phenotypes,
    normalise_agent, parse_template,
)
from bridge.amr.settings import DATA_DIR

from amr_fixtures import write_phenotypes

OVERRIDES = load_family_overrides(DATA_DIR / "family_overrides.tsv")


class TemplateNameTests(unittest.TestCase):
    def test_gene_variant_accession(self):
        name = parse_template("blaCTX-M-15_1_AY044436")
        self.assertEqual((name.gene, name.variant, name.accession),
                         ("blaCTX-M-15", "1", "AY044436"))

    def test_accession_with_underscore_and_primes(self):
        name = parse_template("aac(6')-Ib-cr_1_NG_047471  ")
        self.assertEqual(name.gene, "aac(6')-Ib-cr")
        self.assertEqual(name.accession, "NG_047471")

    def test_unparseable_name_is_its_own_gene(self):
        self.assertEqual(parse_template("mystery").gene, "mystery")


class FamilyTests(unittest.TestCase):
    CASES = {
        "blaCTX-M-15": "blaCTX-M", "blaTEM-1B": "blaTEM", "blaOXA-48": "blaOXA",
        "blaSHV-12": "blaSHV", "blaOXY-1-1": "blaOXY", "blaOKP-B-1": "blaOKP-B",
        "blaZ": "blaZ", "blaACT-1": "blaACT",
        "aac(6')-Ib-cr": "aac(6')", "aph(3')-Ia": "aph(3')", "aph(6)-Id": "aph(6)",
        "qnrB19": "qnrB", "dfrA17": "dfrA", "catA1": "catA", "mcr-1.1": "mcr",
        "tet(X4)": "tet(X)", "tet(M)": "tet(M)", "sul1": "sul1", "erm(B)": "erm(B)",
        "oqxA": "oqxA", "tmexC2": "tmexC",
    }
    OVERRIDDEN = {
        "aac(6')-aph(2'')": "aac(6')-aph(2'')", "ant(3'')-Ia": "ant(3'')/aadA",
        "aadA1": "ant(3'')/aadA", "blaCKO-1": "blaCKO/MAL", "blaMAL-1": "blaCKO/MAL",
        "blaSRT-2": "blaSRT/SST", "blaPDC-1": "blaPAO/PDC", "fosA": "fosA",
        "fosA7": "fosA7", "fosA3": "fosA3", "tet(W/32/O)": "tet(W)", "vga(A)LC": "vga(A)",
        # Chromosomal genes kept apart from their acquired class, and
        # ResFinder's capitalised OqxA/OqxB mapped onto the host table's keys
        "aph(3')-IIb": "aph(3')-IIb", "aac(6')-Ii": "aac(6')-Ii",
        "OqxA": "oqxA", "OqxB": "oqxB",
        # Chromosomal alleles that share a family with acquired ones
        "catB7": "catB7", "aac(6')-Ic": "aac(6')-Ic", "aac(6')-Iaa": "aac(6')-Iaa",
    }

    def test_new_families_reach_their_host_rows(self):
        """Each chromosomal-gene override must land on a curated host row."""
        from bridge.amr.host import load_host_context, resolve
        table = load_host_context(DATA_DIR / "gene_host_context.tsv")
        for gene, organism, verdict in (("aph(3')-IIb", "Pseudomonas aeruginosa", "intrinsic"),
                                        ("aph(3')-Ia", "Pseudomonas aeruginosa", "acquired"),
                                        ("aac(6')-Ii", "Enterococcus faecium", "intrinsic"),
                                        ("OqxB", "Klebsiella pneumoniae", "intrinsic"),
                                        # intrinsic only where the trait lives
                                        ("catB7", "Pseudomonas aeruginosa", "intrinsic"),
                                        ("catB8", "Pseudomonas aeruginosa", "acquired"),
                                        ("blaPDC-1", "Pseudomonas aeruginosa", "intrinsic"),
                                        ("blaPDC-1", "Pseudomonas putida", "undetermined"),
                                        ("lsa(A)", "Enterococcus faecalis", "intrinsic"),
                                        ("lsa(A)", "Enterococcus faecium", "undetermined"),
                                        ("aac(6')-Ic", "Serratia marcescens", "intrinsic"),
                                        ("aac(6')-Ib", "Serratia marcescens", "undetermined"),
                                        ("aac(6')-Iaa", "Salmonella enterica", "intrinsic"),
                                        ("blaCKO-1", "Citrobacter koseri", "intrinsic"),
                                        ("blaCKO-1", "Citrobacter freundii", "undetermined"),
                                        ("tet(Z)", "Corynebacterium striatum", "undetermined")):
            self.assertEqual(resolve(family_of(gene, OVERRIDES), organism, table)[0], verdict,
                             f"{gene} on {organism}")

    def test_rules(self):
        for gene, family in self.CASES.items():
            self.assertEqual(family_of(gene, OVERRIDES), family, gene)

    def test_overrides_beat_rules_and_keep_fosA7_apart(self):
        for gene, family in self.OVERRIDDEN.items():
            self.assertEqual(family_of(gene, OVERRIDES), family, gene)

    def test_every_host_table_family_is_reachable(self):
        """The host table is keyed on gene families; the ones ResFinder names
        by a rule or override must land exactly on those keys."""
        from bridge.amr.host import load_host_context
        keys = {family for family, _ in load_host_context(DATA_DIR / "gene_host_context.tsv")}
        produced = set(self.CASES.values()) | set(self.OVERRIDDEN.values())
        for family in ("blaACT", "blaCTX-M", "blaTEM", "blaOXA", "aac(6')", "aph(3')",
                       "ant(3'')/aadA", "blaSRT/SST", "blaCKO/MAL", "fosA", "fosA7",
                       "qnrB", "dfrA", "mcr", "tet(M)", "tet(W)", "tet(X)"):
            self.assertIn(family, keys)
            self.assertIn(family, produced)


class PhenotypeTests(unittest.TestCase):
    def setUp(self):
        self.db = load_phenotypes(write_phenotypes(Path(tempfile.mkdtemp())), OVERRIDES)

    def test_lookup_by_template_then_gene(self):
        self.assertEqual(self.db.lookup("blaSHV-12_1_AF462395").drug_class, "Beta-lactam")
        # another submission of the same allele falls back to the gene
        self.assertEqual(self.db.lookup("blaSHV-12_2_XX000001").gene, "blaSHV-12")
        self.assertIsNone(self.db.lookup("blaNOPE-1_1_X"))

    def test_family_spectrum_varies_for_shv_not_tem(self):
        union, inter = self.db.family_spectrum("blaSHV")
        self.assertIn("cefotaxime", union)
        self.assertNotIn("cefotaxime", inter)
        union, inter = self.db.family_spectrum("blaTEM")
        self.assertEqual(union, inter)

    def test_agent_names_normalise_across_conventions(self):
        self.assertEqual(normalise_agent("Amoxicillin+Clavulanic acid"),
                         normalise_agent("amoxicillin-clavulanic acid"))
        self.assertEqual(normalise_agent("Trimethoprim/Sulfamethoxazole"),
                         "trimethoprim-sulfamethoxazole")


class BuildTests(unittest.TestCase):
    def test_concatenate_keeps_first_copy_and_skips_all_fsa(self):
        folder = Path(tempfile.mkdtemp())
        (folder / "beta-lactam.fsa").write_text(">blaTEM-1B_1_A\nACGT\nACGT\n>blaSHV-1_1_B\nTTTT\n")
        (folder / "other.fsa").write_text(">blaTEM-1B_1_A\nACGTACGT\n>tet(M)_1_C\nGGGG\n")
        (folder / "all.fsa").write_text(">junk_1_X\nAAAA\n")
        sources, count = concatenate_fastas(folder, folder / "kma" / "resfinder.fsa")
        self.assertEqual(count, 3)
        self.assertEqual([p.name for p in sources], ["beta-lactam.fsa", "other.fsa"])
        text = (folder / "kma" / "resfinder.fsa").read_text()
        self.assertEqual(text.count(">blaTEM-1B_1_A"), 1)
        self.assertNotIn("junk", text)

    def test_dry_run_writes_the_manifest(self):
        folder = Path(tempfile.mkdtemp())
        write_phenotypes(folder)
        (folder / "beta-lactam.fsa").write_text(">blaTEM-1B_1_A\nACGT\n")
        manifest = build_db(folder, kma="kma", kma_index="kma_index", dry_run=True)
        self.assertEqual(manifest["templates"], 1)
        self.assertEqual(manifest["command"][0], "kma_index")
        written = json.loads((folder / "kma" / "db_manifest.json").read_text())
        self.assertEqual(written["index_prefix"], str((folder / "kma" / "resfinder").resolve()))

    def test_missing_phenotypes_is_refused(self):
        folder = Path(tempfile.mkdtemp())
        (folder / "beta-lactam.fsa").write_text(">a_1_b\nACGT\n")
        with self.assertRaises(FileNotFoundError):
            build_db(folder, dry_run=True)


if __name__ == "__main__":
    unittest.main()
