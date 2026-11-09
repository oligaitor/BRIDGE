"""Allele -> antibiotics, from ResFinder's ``phenotypes.txt``.

Written into ``amr_hits.csv`` by the arm, so the report never needs the
database. Three columns carry it:

* ``antibiotics`` - the matched allele's list, verbatim ResFinder names;
* ``family_spectrum_varies`` - whether the family's alleles disagree about
  which antibiotics they affect (SHV, TEM, OXA do; tet(M) does not);
* ``family_possible_antibiotics`` - antibiotics some OTHER allele of the
  family affects but this one does not.

The report combines the last two with the allele level: an ``exact`` allele
is taken at its word, while a ``nearest`` allele in a family whose spectrum
varies is ``phenotype_uncertain``, and the extra antibiotics are reported as
*possible*, never as predicted.
"""

from __future__ import annotations

from .database import PhenotypeDB, normalise_agent


def annotate(hit: dict, db: PhenotypeDB | None) -> dict:
    """Fill the phenotype columns of one hit, in place."""
    entry = db.lookup(hit.get("template", "")) if db is not None else None
    if entry is None:
        hit.update(drug_class="", antibiotics="", mechanism="", phenotype_notes="",
                   family_spectrum_varies="", family_possible_antibiotics="")
        if db is not None:
            hit["flags"] = "; ".join(filter(None, [hit.get("flags", ""), "no_phenotype_entry"]))
        return hit
    own = {normalise_agent(agent) for agent in entry.antibiotics}
    union, inter = db.family_spectrum(hit.get("family", ""))
    extra = sorted(union - own)
    hit.update(
        drug_class=entry.drug_class,
        antibiotics="; ".join(entry.antibiotics),
        mechanism=entry.mechanism,
        phenotype_notes="; ".join(filter(None, [entry.notes, entry.required_gene and
                                               f"requires {entry.required_gene}"])),
        family_spectrum_varies=str(union != inter),
        family_possible_antibiotics="; ".join(extra),
    )
    return hit


def is_uncertain(hit: dict) -> bool:
    """Named antibiotics rest on an allele that was not established."""
    return (str(hit.get("allele_level")) != "exact"
            and str(hit.get("family_spectrum_varies")).strip().lower() == "true")
