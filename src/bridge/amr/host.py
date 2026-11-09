"""Is a resistance gene INTRINSIC or ACQUIRED for the organism it was found with?

The curated table is ``data/gene_host_context.tsv``, every asserted verdict
cited. Its keys are gene families as derived from the ResFinder gene name
(``database.family_of``).

ResFinder is a catalogue of *acquired* genes, but a large share of confident
calls are chromosomal genes their host carries as a matter of course -
blaACT in Enterobacter, fosA in Klebsiella, blaSRT in Serratia - listed
because they are mobilised in *other* genera. Calling those acquired
overstates the finding; this table is what stops that.

SAFETY PROPERTY.``intrinsic``, ``acquired`` and
``not_expected_in_taxon`` are only ever asserted by an explicit curated row
carrying a citation. An uncurated combination is ``undetermined``: a to-do,
never a silent wrong verdict.
"""

from __future__ import annotations

import re
from pathlib import Path

VERDICTS = ("intrinsic", "acquired", "not_expected_in_taxon", "undetermined")
FIELDS = ("gene_family", "taxon", "taxon_rank", "verdict", "citation", "citation_type",
          "note", "source", "curated_on")

# Epithets that are not species, and genus renames that are the same taxon.
NON_SPECIES = {"species", "sp", "sp.", "spp", "spp.", "complex", "group", "unidentified"}
GENUS_SYNONYMS = {
    "Metamycoplasma": "Mycoplasma",
    "Mycoplasmoides": "Mycoplasma",
}


def _clean(name) -> str:
    """Strip provisional-name brackets: '[Pantoea]' -> 'Pantoea'."""
    return re.sub(r"[\[\]]", "", str(name or "")).strip()


def genus_of(name) -> str:
    tokens = _clean(name).split()
    if not tokens:
        return ""
    genus = tokens[0]
    return GENUS_SYNONYMS.get(genus, genus)


def species_of(name) -> str | None:
    tokens = _clean(name).split()
    if len(tokens) < 2:
        return None
    epithet = tokens[1].lower()
    return None if epithet in NON_SPECIES else epithet


def _taxon_key(name) -> str:
    name = str(name).strip()
    if name in ("*", ""):
        return "*"
    genus, species = genus_of(name), species_of(name)
    return f"{genus} {species}" if species else genus


def load_host_context(path: Path) -> dict[tuple[str, str], dict]:
    """Parse the curated TSV -> {(family, taxon_key): record}. First entry wins."""
    path = Path(path)
    table: dict[tuple[str, str], dict] = {}
    header = None
    with open(path, encoding="utf-8") as handle:
        for lineno, raw in enumerate(handle, 1):
            line = raw.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            cells = line.split("\t")
            if header is None:
                header = cells
                if tuple(header) != FIELDS:
                    raise ValueError(f"{path.name}:{lineno}: header is {header}, "
                                     f"expected {list(FIELDS)}")
                continue
            if len(cells) != len(FIELDS):
                raise ValueError(f"{path.name}:{lineno}: {len(cells)} cells, expected "
                                 f"{len(FIELDS)}. Line: {line[:90]!r}")
            record = dict(zip(FIELDS, (cell.strip() for cell in cells)))
            if record["verdict"] not in VERDICTS:
                raise ValueError(f"{path.name}:{lineno}: verdict {record['verdict']!r} "
                                 f"not one of {VERDICTS}")
            if record["verdict"] != "undetermined" and not record["citation"]:
                raise ValueError(f"{path.name}:{lineno}: verdict {record['verdict']!r} for "
                                 f"{record['gene_family']}/{record['taxon']} has no citation")
            table.setdefault((record["gene_family"], _taxon_key(record["taxon"])), record)
    if header is None:
        raise ValueError(f"{path.name}: no header row found")
    return table


def resolve(family: str | None, organism: str | None,
            table: dict[tuple[str, str], dict]) -> tuple[str, str, str, str]:
    """(verdict, citation, note, match_rank). Species row > genus row > family default."""
    if not family:
        return "undetermined", "", "gene not mapped to a family", "no_family"
    if not organism:
        return "undetermined", "", "no organism to judge against", "no_taxon"
    genus, species = genus_of(organism), species_of(organism)
    candidates = []
    if species:
        candidates.append((f"{genus} {species}", "species"))
    candidates.append((genus, "genus"))
    # 'Escherichia-Shigella' style compound genera: try each part.
    if "-" in genus:
        candidates += [(part, "genus") for part in genus.split("-") if part]
    candidates.append(("*", "family_default"))
    for key, rank in candidates:
        record = table.get((family, key))
        if record:
            return record["verdict"], record["citation"], record["note"], rank
    return "undetermined", "", f"no curated row for {family} on {organism}", "uncurated"
