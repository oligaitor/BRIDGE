"""Reconciliation: rMLST governs, 16S is the fallback.

* **rMLST called something** - the call set is rMLST's, one row per species it
  resolved. 16S contributes no calls.
* **rMLST called nothing at all** - 16S produces the call set instead.

A call carried by tens of ribosomal-protein loci does not need a single 16S
gene to vouch for it, and requiring joint agreement sends clean rMLST species
calls to review whenever 16S names an unrelated genus. 16S is still run on
every sample, and its numbers still travel to the report; it just does not
vote.

Nothing disappears. The rMLST reference is PANDAv7 *target genera*, a restricted
genus set, so 16S can legitimately see an organism rMLST structurally cannot.
Those appear as ``16s_only`` rows: present, labelled, and not part of the
reported call set. The same goes for a demoted near neighbour, which appears as
``rmlst_alternative``.

There is no confidence tier and no detection gate. What stands in their place
is measurement: the RPM-r band against the run's own negative control,
per-organism read counts, and per-organism coverage of the loci behind each
call. A reader weighs those directly instead of reading a letter grade derived
from them.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum

from .config import Config
from .strength import (
    BlastStrength,
    BlastStrengthResult,
    RmlstStrength,
    RmlstStrengthResult,
)
from .taxonomy import genus_matches


class ReportedRank(str, Enum):
    SPECIES = "species"
    GENUS = "genus"


# The two values that mean "this row is a reported organism". Every "is this a
# call?" question in the package routes through this set, so adding an
# informational row type cannot accidentally make it count as a finding.
CALL_BASES = frozenset({"rmlst", "16s_fallback"})

#: ``basis`` vocabulary, in full.
#:
#: ``rmlst``              rMLST reported this species. A reported call.
#: ``16s_fallback``       rMLST resolved nothing anywhere in this sample, so the
#:                        16S arm produced the call set. A reported call.
#: ``rmlst_alternative``  A demoted near neighbour: it held up as an allele
#:                        profile but fell short on exclusive loci, and names
#:                        the organism it is a second name *for*. Informational.
#: ``16s_only``           rMLST governed this sample and did not call this
#:                        genus, but 16S saw it. Informational - it may be a
#:                        genuine off-target organism or reagent background, and
#:                        the RPM columns are how you tell.
BASIS_VALUES = ("rmlst", "16s_fallback", "rmlst_alternative", "16s_only")

#: Flag set on the report-eligible ``16s_fallback`` calls of a sample whose
#: calls are all 16S fallback and hold at least
#: ``reconcile.multi_fallback_review_min`` such eligible calls. See
#: :func:`flag_multi_fallback`.
MULTI_FALLBACK_FLAG = "multiple_16s_fallback"
DEFAULT_MULTI_FALLBACK_REVIEW_MIN = 2


def below_background(flags) -> bool:
    """Did the RPM gate put this call at or below the run's reagent background?

    The gate flags the arm result (``rpm_below_background``) and the reconciler
    prefixes it with the arm (``16s:rpm_below_background``), hence ``endswith``.
    Shared with summary_csv._status so the two cannot disagree about which
    calls were eligible for ``report``.
    """
    return any(flag.endswith("rpm_below_background") for flag in flags)


@dataclass
class OrganismCall:
    """One organism (one reported row)."""

    organism: str
    genus: str | None
    reported_rank: ReportedRank
    basis: str
    # Where this organism's evidence sits in each arm's own CSV. This is how the
    # report fetches per-organism numbers back out of the arm files. Deliberately
    # an index and not a name: the reported name need not match the arm's row
    # (a genus-level name, a compound genus), and a failed name lookup would
    # silently lose every number behind the row.
    rmlst_row_index: int | None = None
    s16_row_index: int | None = None
    rmlst_taxon: str | None = None
    s16_top_hit: str | None = None
    s16_margin: float | None = None
    flags: list[str] = field(default_factory=list)
    notes: str = ""

    @property
    def is_call(self) -> bool:
        return self.basis in CALL_BASES


def _named(taxon) -> tuple[ReportedRank, str]:
    """The finest name this taxon commits to, and the rank that name sits at."""
    if taxon is not None and taxon.species_full:
        return ReportedRank.SPECIES, taxon.species_full
    return ReportedRank.GENUS, (taxon.genus if taxon else "") or ""


def _arm_flags(prefix: str, result) -> list[str]:
    return [f"{prefix}:{flag}" for flag in (result.flags if result else [])]


def _partner_16s(
    rmlst_result: RmlstStrengthResult,
    s16_called: list[BlastStrengthResult],
    cfg: Config,
) -> BlastStrengthResult | None:
    """The 16S result covering the same genus as this rMLST call, if any.

    The genus is the join key between the arms. 16S reference names include
    compound spellings, so ``Escherichia-Shigella`` has to find
    ``Escherichia`` or it stands up as a spurious ``16s_only`` row for an
    organism rMLST already reported. :func:`taxonomy.genus_matches` is
    compound-aware, which is why it rather than string equality is used.

    One 16S row may partner several rMLST rows, and that is correct. A genus-
    level 16S call in a sample where rMLST resolved two species of that genus
    corroborates both equally - 16S saw the genus and could not split it, which
    supports every species in it, not one picked arbitrarily.
    """
    genus = rmlst_result.taxon.genus if rmlst_result.taxon else None
    if not genus:
        return None
    for result in s16_called:
        candidates = [
            candidate.genus for candidate in
            ([result.top_taxon] if result.top_taxon else []) + list(result.candidate_taxa)
            if candidate is not None and candidate.genus
        ]
        if any(genus_matches(candidate, genus, cfg) for candidate in candidates):
            return result
    return None


def _rmlst_governed(
    rmlst_called: list[RmlstStrengthResult],
    s16_called: list[BlastStrengthResult],
    cfg: Config,
) -> list[OrganismCall]:
    """The call set when rMLST resolved anything at all."""
    calls: list[OrganismCall] = []
    partnered: set[int] = set()

    for result in rmlst_called:
        partner = _partner_16s(result, s16_called, cfg)
        if partner is not None:
            partnered.add(id(partner))

        rank, organism = _named(result.taxon)
        # SPECIES_SUPPORTED is a reported call. GENUS_ONLY from this arm only
        # ever means an `alternative_to` row - a second acceptable name for an
        # organism already called, not a second organism. Counting it as a call
        # would make every near-neighbour pair read as polymicrobial.
        basis = ("rmlst" if result.strength is RmlstStrength.SPECIES_SUPPORTED
                 else "rmlst_alternative")

        flags = _arm_flags("rmlst", result) + _arm_flags("16s", partner)
        notes = ""
        if partner is None:
            # Informational only: the second marker had nothing to say here.
            flags.append("no_16s_corroboration")

        calls.append(OrganismCall(
            organism=organism,
            genus=result.taxon.genus if result.taxon else None,
            reported_rank=rank,
            basis=basis,
            rmlst_row_index=result.group_index,
            s16_row_index=partner.row_index if partner is not None else None,
            rmlst_taxon=organism or None,
            s16_top_hit=_named(partner.top_taxon)[1] or None if partner else None,
            s16_margin=partner.margin if partner else None,
            flags=flags,
            notes=notes,
        ))

    # Anything 16S saw in a genus rMLST did not call. Not a reported call - but
    # it must not vanish either, because the rMLST reference covers a restricted
    # set of target genera and cannot see outside it.
    for result in s16_called:
        if id(result) in partnered:
            continue
        rank, organism = _named(result.top_taxon)
        if not organism:
            continue
        calls.append(OrganismCall(
            organism=organism,
            genus=result.top_taxon.genus if result.top_taxon else None,
            reported_rank=rank,
            basis="16s_only",
            s16_row_index=result.row_index,
            s16_top_hit=organism,
            s16_margin=result.margin,
            flags=_arm_flags("16s", result),
            notes="seen by the 16S arm only; rMLST resolved other organisms in "
                  "this sample but did not call this genus, which its "
                  "target-genus reference may not cover",
        ))
    return calls


def _s16_fallback(s16_called: list[BlastStrengthResult]) -> list[OrganismCall]:
    """The call set when rMLST resolved nothing anywhere in the sample."""
    calls: list[OrganismCall] = []
    for result in s16_called:
        rank, organism = _named(result.top_taxon)
        if not organism:
            continue
        calls.append(OrganismCall(
            organism=organism,
            genus=result.top_taxon.genus if result.top_taxon else None,
            reported_rank=rank,
            basis="16s_fallback",
            s16_row_index=result.row_index,
            s16_top_hit=organism,
            s16_margin=result.margin,
            flags=_arm_flags("16s", result),
            notes="16S fallback: the rMLST arm resolved nothing in this sample",
        ))
    return calls


def reconcile_sample(
    rmlst_results: list[RmlstStrengthResult],
    s16_results: list[BlastStrengthResult],
    cfg: Config,
) -> list[OrganismCall]:
    """Full reconciliation for one sample -> one OrganismCall per organism."""
    rmlst_called = [
        result for result in rmlst_results
        if result.strength is not RmlstStrength.SILENT
        and result.taxon is not None and result.taxon.genus
    ]
    s16_called = [
        result for result in s16_results
        if result.strength is not BlastStrength.SILENT_WEAK
        and result.top_taxon is not None and result.top_taxon.genus
    ]

    calls = (_rmlst_governed(rmlst_called, s16_called, cfg) if rmlst_called
             else _s16_fallback(s16_called))

    # Deterministic: reported calls first, then by genus and organism.
    return sorted(
        calls,
        key=lambda call: (0 if call.is_call else 1, call.genus or "~", call.organism),
    )


def flag_multi_fallback(calls: list[OrganismCall], threshold: int) -> int:
    """Hold a sample's report-eligible 16S fallback calls when there are several.

    Fires only when BOTH hold:

    * every reported call in the sample is ``16s_fallback`` - rMLST resolved
      nothing, so the sample carried little target-genus material;
    * at least ``threshold`` of those calls are still eligible for ``report``,
      i.e. the RPM gate did not already put them at or below background.

    Those eligible calls - and only those - are flagged, and ``status`` reads
    the flag (summary_csv._status) to give ``review_override`` in place of
    ``report``. A fallback call the gate already failed stays ``review`` and is
    not counted: it is already held, on a measured verdict.

    Why: when 16S names several organisms at once with no rMLST signal, the
    likeliest explanation is reagent or environmental background the run's
    control did not model (Paracoccus, Sphingomonas, Micromonospora and the
    like), not a genuine polymicrobial infection.

    Flagged, not removed: the rows stay in the table. A threshold of 0 or less
    disables the rule. Returns the number of eligible fallback calls seen.

    Mutates ``flags``, so call it before serialising.
    """
    reported = [call for call in calls if call.is_call]
    if any(call.basis != "16s_fallback" for call in reported):
        return 0
    eligible = [call for call in reported if not below_background(call.flags)]
    if threshold > 0 and len(eligible) >= threshold:
        for call in eligible:
            call.flags.append(MULTI_FALLBACK_FLAG)
            call.notes = "; ".join(filter(None, [
                call.notes,
                f"one of {len(eligible)} report-eligible 16S fallback calls in this "
                f"sample - held for review as possible unmodelled "
                f"reagent/environmental background",
            ]))
    return len(eligible)


def sample_outcome(calls: list[OrganismCall]) -> str:
    """The sample's headline result, and the polymicrobial flag that goes with it.

    Only reported calls count. An ``rmlst_alternative`` row is a second name for
    an organism already in the list, and a ``16s_only`` row is an observation the
    governing arm did not make - neither is a second organism, and letting either
    reach this count would report a single mitis-group isolate as a mixed
    infection.

    Mutates ``flags``, so call it before serialising.
    """
    reported = [call for call in calls if call.is_call]
    if not reported:
        return "no_call"
    if len(reported) >= 2:
        for call in reported:
            call.flags.append("polymicrobial")
        return "polymicrobial"
    return "reported"
