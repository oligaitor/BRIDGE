"""The two arms' verdict vocabulary.

Each arm reduces its own evidence to one of three strengths, and the reconciler
reads nothing else. That is what lets the two arms measure completely different
things - allele profiles across ~53 ribosomal protein loci on one side, 16S
abundance against MIMt on the other - and still be compared.

The ranking functions themselves live with the data they read:

    species_calls.csv    -> adapters.rmlst.rmlst_results
    marker16s_calls.csv  -> rank_marker16s.rank_marker16s
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum

from .taxonomy import TaxonName


# ---------------------------------------------------------------------------
# rMLST
# ---------------------------------------------------------------------------


class RmlstStrength(str, Enum):
    SPECIES_SUPPORTED = "species_supported"
    GENUS_ONLY = "genus_only"
    SILENT = "silent"


@dataclass(frozen=True)
class RmlstStrengthResult:
    strength: RmlstStrength
    taxon: TaxonName | None
    support: float | None
    n_loci: int | None
    flags: list[str] = field(default_factory=list)
    reason: str = ""
    coverage: float | None = None
    # This result's row number in species_calls.csv. It is how the report gets
    # back to the arm's own per-organism numbers without matching on a name.
    group_index: int | None = None


# ---------------------------------------------------------------------------
# 16S
# ---------------------------------------------------------------------------


class BlastStrength(str, Enum):
    SPECIES_STRONG = "species_strong"
    GENUS_ONLY = "genus_only"
    SILENT_WEAK = "silent_weak"


@dataclass(frozen=True)
class BlastStrengthResult:
    strength: BlastStrength
    top_taxon: TaxonName | None
    top_pident: float | None
    top_qcovs: float | None
    margin: float | None
    candidate_taxa: list[TaxonName] = field(default_factory=list)
    flags: list[str] = field(default_factory=list)
    reason: str = ""
    query: str = ""
    coverage: float | None = None
    coverage_matched_by: str = ""
    # This result's row number in marker16s_calls.csv; see group_index above.
    row_index: int | None = None
