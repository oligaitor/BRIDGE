"""Read taxonomy out of MIMt reference headers.

The database is ``MIMt-16S_M2c_25_04_tax_bacteria_tax_in_id.fna``. Every header
carries a complete seven-rank lineage, semicolon-delimited, with no whitespace:

    >NR_199355.1;K__Bacteria;P__Cyanobacteriota;C__Cyanophyceae;O__Oculatellales;
     F__Oculatellaceae;G__Eurychoronema;S__Eurychoronema_bolivianum_strain_LEGE_231228

That last property is why this variant of the file is the one to use. kallisto
takes a target's name as the header up to the first whitespace, so the entire
lineage survives into ``abundance.tsv`` as the ``target_id``. The arm therefore
needs no mapping table, no taxonomy dump and no database lookups - the answer is
already attached to every row kallisto emits.

Two hazards in the real data that this parser exists to absorb:

``G__`` and ``S__`` can disagree
    ``G__Agrilutibacter;S__Lysobacter_terrestris_strain_II4`` - a reclassified
    organism whose species string still carries the old genus. The ``G__`` field
    is the curated one and wins; the species epithet is taken positionally from
    the ``S__`` field, so a stale genus there does no harm.

Strain and isolate suffixes
    ``S__`` is ``<Genus>_<epithet>_<anything>`` - ``_strain_MI-28``,
    ``_isolate_S23``, ``_DSM_43160``, or a bare code like ``_Hxd3``. Only the
    epithet is wanted. Taking the second underscore-token rather than trying to
    strip suffixes handles all of these shapes without a list of exceptions.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

# Rank prefixes in the order MIMt writes them, after the accession.
RANK_PREFIXES = ("K__", "P__", "C__", "O__", "F__", "G__", "S__")
LINEAGE_FIELDS = 1 + len(RANK_PREFIXES)  # accession + seven ranks

# A species epithet is lowercase alphabetic. Anything else in that position
# (a strain code, a digit, an authority) means the epithet is not resolved.
EPITHET_RE = re.compile(r"^[a-z]+$")

# Genus names occasionally carry a nomenclatural authority:
# ``Methylacidiphilum_(ex_Ratnadevi_et_al._2023)``. Keep the name, drop the rest.
GENUS_AUTHORITY_RE = re.compile(r"_?\(.*$")

# Placeholders that mean "not identified to this rank" rather than a real name.
UNRESOLVED = {"", "sp", "spp", "uncultured", "unclassified", "bacterium",
              "unidentified", "metagenome"}


@dataclass(frozen=True)
class ReferenceTaxon:
    """The taxonomy of one reference sequence, parsed from its header."""

    target_id: str
    accession: str
    genus: str | None
    species: str | None
    family: str | None = None
    order: str | None = None

    @property
    def species_full(self) -> str | None:
        """``"Genus species"``, or None where the epithet is unresolved."""
        if self.genus and self.species:
            return f"{self.genus} {self.species}"
        return None

    @property
    def label(self) -> str:
        """Best available name: species if resolved, else genus, else accession."""
        return self.species_full or self.genus or self.accession


def _clean(value: str, prefix: str) -> str:
    """Strip a rank prefix and surrounding whitespace."""
    value = value.strip()
    return value[len(prefix):].strip() if value.startswith(prefix) else value


def parse_lineage_header(header: str) -> ReferenceTaxon:
    """Parse one MIMt header into a :class:`ReferenceTaxon`.

    Never raises on a malformed header. A header this cannot read yields a
    taxon with ``genus=None``, which the aggregation step reports as
    unassignable rather than silently dropping - a reference the arm cannot
    interpret is a fact worth surfacing, not an error worth crashing on.
    """
    header = header.lstrip(">").split()[0] if header.strip() else ""
    parts = header.split(";")
    accession = parts[0] if parts else ""

    if len(parts) < LINEAGE_FIELDS:
        return ReferenceTaxon(target_id=header, accession=accession,
                              genus=None, species=None)

    genus_raw = _clean(parts[6], "G__")
    species_raw = _clean(parts[7], "S__")
    family = _clean(parts[5], "F__") or None
    order = _clean(parts[4], "O__") or None

    genus = GENUS_AUTHORITY_RE.sub("", genus_raw).strip("_") or None
    if genus and genus.lower() in UNRESOLVED:
        genus = None

    # Positional: <Genus>_<epithet>_<strain junk>. The genus token here may
    # disagree with G__ (reclassified organisms); the epithet is still correct.
    species = None
    tokens = species_raw.split("_")
    if len(tokens) >= 2:
        candidate = tokens[1].lower()
        if EPITHET_RE.match(candidate) and candidate not in UNRESOLVED:
            species = candidate

    return ReferenceTaxon(target_id=header, accession=accession, genus=genus,
                          species=species, family=family, order=order)


def load_reference(fasta: Path) -> dict[str, ReferenceTaxon]:
    """Parse every header in a MIMt FASTA into ``{target_id: ReferenceTaxon}``.

    Reads headers only - the sequences are never held in memory, so this costs
    a few seconds and a few MB on the 29,090-sequence database.
    """
    taxa: dict[str, ReferenceTaxon] = {}
    with open(fasta, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            taxon = parse_lineage_header(line)
            if taxon.target_id:
                taxa[taxon.target_id] = taxon
    if not taxa:
        raise ValueError(f"No FASTA headers found in {fasta}")
    return taxa


def find_index(reference: Path, index: Path | None = None) -> Path:
    """Locate the kallisto index for a reference FASTA.

    An explicit path wins. Otherwise look for ``<reference stem>.idx`` beside
    the FASTA - which is how the MIMt distribution ships, with a prebuilt index
    next to the sequences.
    """
    if index is not None:
        if not Path(index).is_file():
            raise FileNotFoundError(f"kallisto index not found: {index}")
        return Path(index)

    candidate = reference.with_suffix(".idx")
    if candidate.is_file():
        return candidate
    raise FileNotFoundError(
        f"No kallisto index beside {reference}. Expected {candidate}; "
        f"build one with `kallisto index -i {candidate} {reference}` or pass "
        f"--kallisto-index explicitly."
    )
