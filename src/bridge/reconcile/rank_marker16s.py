"""Turn a 16S marker call into the reconciler's ``BlastStrengthResult``.

The result type is named for BLAST, the classic way of calling a 16S read set.
This arm is not BLAST-based: it decides on kallisto abundance against MIMt,
not on percent identity to a top hit. What carries across is the *shape of the
judgement*, which is what lets one reconciler read both:

    BLAST                              marker16s
    -----                              ---------
    percent identity to the top hit -> that genus's share of the profile
    identity gap to the runner-up   -> the top species' lead within its genus
    SPECIES_STRONG / GENUS_ONLY /      the same three, decided in calls.py
      SILENT_WEAK

Both answer "is the top candidate far enough ahead of its nearest rival to be
named?" - one in units of sequence identity, the other in units of abundance.
The reconciler only ever reads ``strength``, the genus and the species name, so
the substitution is invisible to it.

Pseudoalignment produces no alignment to score, so after the call is made the
arm aligns the reads to the called genus's references and BLASTs a consensus
against them. ``top_pident`` and ``top_qcovs`` carry that consensus's percent
identity and query coverage - the same two numbers a Sanger 16S report gives -
and ``coverage`` carries the breadth of the top reference. They are evidence
columns; the strength decision is still the abundance one. A calls file
without these columns leaves them ``None``.
"""

from __future__ import annotations

from dataclasses import dataclass

from .config import Config
from .strength import BlastStrength, BlastStrengthResult
from .taxonomy import TaxonName, parse_taxon

# marker16s_calls.csv strength strings -> the reconciler's enum.
_STRENGTH = {
    "species_strong": BlastStrength.SPECIES_STRONG,
    "genus_only": BlastStrength.GENUS_ONLY,
    "silent_weak": BlastStrength.SILENT_WEAK,
}


@dataclass(frozen=True)
class Marker16SRow:
    """One row of ``marker16s_calls.csv``, as the adapter reads it."""

    sample_id: str
    taxon: str
    genus: str
    species: str
    rank: str
    strength: str
    abundance_frac: float
    genus_frac: float
    reads_assigned: float
    raw_reads: float
    copy_number: float
    within_genus_share: float
    margin: float
    n_references: int
    competing_species: str
    reason: str
    # Alignment metrics. None where the file lacks them or
    # the stage did not run for this genus - never 0, which would assert a
    # measurement that was not made.
    n_aligned_reads: int | None = None
    mean_read_identity: float | None = None
    ref_breadth: float | None = None
    mean_depth: float | None = None
    consensus_length: int | None = None
    consensus_pident: float | None = None
    consensus_qcov: float | None = None
    consensus_top_hit: str = ""
    alignment_status: str = ""


def rank_marker16s(
    row: Marker16SRow, cfg: Config, *, row_index: int | None = None
) -> BlastStrengthResult:
    """Map one marker16s call onto a ``BlastStrengthResult``.

    The arm has already made the strength decision in
    :mod:`bridge.marker16s.calls`, where the abundance data lives. This does
    not re-decide it - re-deriving the verdict from summary columns would let
    the two copies of the rule drift apart. It translates.

    ``row_index`` is this row's position in ``marker16s_calls.csv``. It is how
    the reconciled report gets back to the arm's own per-organism numbers
    without matching on a name - see ``reconcile_cli._s16_evidence``. Matching
    by name silently misses whenever the reconciler reports an organism under a
    different label than the arm used.
    """
    strength = _STRENGTH.get(row.strength.strip().lower())
    if strength is None:
        raise ValueError(
            f"Unknown marker16s strength {row.strength!r} for {row.taxon!r}. "
            f"Expected one of {sorted(_STRENGTH)}."
        )

    # Name the organism at the rank the arm was willing to commit to. A
    # GENUS_ONLY call must not carry a species name into the reconciler, or the
    # matrix would treat a declined call as a confident one.
    label = row.taxon if strength is BlastStrength.SPECIES_STRONG else (row.genus or row.taxon)
    taxon: TaxonName | None = parse_taxon(label, cfg) if label else None
    if taxon is not None and strength is not BlastStrength.SPECIES_STRONG:
        taxon = TaxonName(raw=taxon.raw, genus=taxon.genus, species=None)

    candidates = [
        parse_taxon(name.strip(), cfg)
        for name in (row.competing_species or "").split(";")
        if name.strip()
    ]
    # candidate_taxa is what the reconciler reads when matching this result to
    # an rMLST call's genus (``reconcile._partner_16s``), so an empty list makes
    # the arm's opinion unmatchable. That happens for a genus whose species are
    # all unresolved, which is precisely a case worth reporting at genus level.
    # Guarantee at least the taxon this call is about.
    if not candidates and taxon is not None:
        candidates = [taxon]

    flags: list[str] = []
    if strength is BlastStrength.GENUS_ONLY and len(candidates) > 1:
        # Records *why* the species was withheld, so a reviewer sees that the
        # marker was ambiguous rather than absent.
        flags.append("marker16s_within_genus_ambiguous")
    if row.copy_number and row.copy_number > 1.0:
        flags.append(f"copy_number_corrected_{row.copy_number:.1f}x")

    return BlastStrengthResult(
        strength=strength,
        top_taxon=taxon,
        # From the alignment stage; None where absent.
        # The reconciler never branches on them - they are evidence.
        top_pident=row.consensus_pident,
        top_qcovs=row.consensus_qcov,
        margin=row.margin,
        candidate_taxa=candidates,
        flags=flags,
        reason=row.reason,
        query=f"marker16s:{row.genus or row.taxon}",
        coverage=row.ref_breadth,
        coverage_matched_by="16s_alignment" if row.ref_breadth is not None else "",
        row_index=row_index,
    )
