"""Taxonomy normalisation — pure, offline, deterministic.

Parses a taxon string into genus + species epithet (handling the
``Genus~species`` convention and plain ``Genus species``) and applies the
optional manual synonym maps from config. Network-free, so strength ranking and
reconciliation are fully testable offline.
"""

from __future__ import annotations

from dataclasses import dataclass

from .config import Config


@dataclass(frozen=True)
class TaxonName:
    """A parsed/normalised taxon.

    ``genus`` is ``None`` only when the input had no usable name. ``species`` is
    the lowercased epithet (e.g. ``"mucilaginosa"``) or ``None`` for a
    genus-only name.
    """

    raw: str
    genus: str | None
    species: str | None

    @property
    def species_full(self) -> str | None:
        if self.genus and self.species:
            return f"{self.genus} {self.species}"
        return None

    @property
    def rank(self) -> str:
        if self.genus and self.species:
            return "species"
        if self.genus:
            return "genus"
        return "unknown"


def parse_taxon(raw: str | None, cfg: Config) -> TaxonName:
    """Parse a raw taxon string into genus + species epithet (no synonyms applied)."""
    if raw is None:
        return TaxonName(raw="", genus=None, species=None)

    sep = cfg.require("taxonomy.genus_species_separator")
    text = raw.replace(sep, " ").strip()
    tokens = text.split()
    if not tokens:
        return TaxonName(raw=raw, genus=None, species=None)

    genus = tokens[0]
    # Skip nomenclatural prefixes that are not the genus.
    if genus.lower() in {"candidatus", "uncultured"} and len(tokens) > 1:
        tokens = tokens[1:]
        genus = tokens[0]

    species = None
    if len(tokens) >= 2 and tokens[1].isalpha():
        species = tokens[1].lower()
    return TaxonName(raw=raw, genus=genus, species=species)


def normalise(raw: str | None, cfg: Config) -> TaxonName:
    """Parse and apply the manual synonym maps from config.

    Order: species-level overrides first (they may rewrite both genus and
    species), then genus-level overrides. Pure and offline.
    """
    name = parse_taxon(raw, cfg)
    if name.genus is None:
        return name

    sep = cfg.require("taxonomy.genus_species_separator")
    species_map: dict = cfg.get("taxonomy.species_synonym_map") or {}
    genus_map: dict = cfg.get("taxonomy.synonym_map") or {}

    if name.species is not None:
        key = f"{name.genus}{sep}{name.species}"
        if key in species_map:
            name = parse_taxon(species_map[key], cfg)
            if name.genus is None:
                return name

    if name.genus in genus_map:
        canonical_genus = genus_map[name.genus]
        return TaxonName(raw=name.raw, genus=canonical_genus, species=name.species)

    return name


# ---------------------------------------------------------------------------
# Compound genus names
# ---------------------------------------------------------------------------


def genus_variants(genus: str | None, cfg: Config) -> set[str]:
    """The comparable forms of a genus name (lowercased).

    SILVA does not split clades it cannot resolve, naming them with a hyphen —
    ``Escherichia-Shigella``. rMLST names the same organisms
    ``Escherichia``. The compound form only reaches the pipeline when SILVA's
    organism field is an uninformative placeholder and ``silva_organism`` falls
    back to the genus rank, but when it does, treating the two spellings as
    different genera splits one organism into two calls that then read as a
    metric conflict. Including each component as a variant lets them join.

    Disabled by ``taxonomy.merge_compound_genus: false``, which reduces this to
    a plain case-insensitive comparison.
    """
    key = (genus or "").strip().lower()
    if not key:
        return set()
    variants = {key}
    if cfg.get("taxonomy.merge_compound_genus", True):
        sep = cfg.get("taxonomy.compound_genus_separator", "-")
        if sep and sep in key:
            variants.update(part for part in key.split(sep) if part)
    return variants


def genus_matches(a: str | None, b: str | None, cfg: Config) -> bool:
    """True if two genus strings name the same genus (see :func:`genus_variants`)."""
    if not a or not b:
        return False
    return bool(genus_variants(a, cfg) & genus_variants(b, cfg))


def is_compound_genus(genus: str | None, cfg: Config) -> bool:
    """True if ``genus`` is a hyphenated SILVA-style compound name."""
    if not genus or not cfg.get("taxonomy.merge_compound_genus", True):
        return False
    sep = cfg.get("taxonomy.compound_genus_separator", "-")
    return bool(sep) and sep in genus


def same_genus(a: TaxonName, b: TaxonName) -> bool:
    return a.genus is not None and a.genus == b.genus


def same_species(a: TaxonName, b: TaxonName) -> bool:
    return (
        a.genus is not None
        and a.genus == b.genus
        and a.species is not None
        and a.species == b.species
    )
