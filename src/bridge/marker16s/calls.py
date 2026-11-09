"""Turn per-taxon abundances into calls, with an honest resolution limit.

The rule this module implements is the arm's one piece of real judgement, so it
is worth stating plainly.

16S cannot separate some organisms. *Escherichia coli* and *Shigella* are one
genomospecies by this marker; the *Streptococcus mitis* group (mitis, oralis,
pneumoniae, pseudopneumoniae) differs by a handful of bases; the *Bacillus
cereus* group likewise. The usual failure mode of a 16S pipeline is to pick the
top hit anyway and report a species it cannot actually distinguish, such as
an *S. pneumoniae* over-call.

Rather than curate a list of "difficult" taxa, this arm measures the ambiguity
directly and lets the data declare it. Within a genus that has enough signal to
be real, look at how that genus's mass is distributed across its species:

* one species holds a dominant share  -> SPECIES_STRONG, name the species
* the mass is spread across several   -> GENUS_ONLY, name only the genus
* the genus has too little signal     -> SILENT_WEAK, report nothing

The **margin** - the gap between the top species' share of its genus and the
runner-up's - is the measurement. It is the direct analogue of BLAST's percent
identity margin, and it makes the arm self-diagnose: *E. coli*/*Shigella*, the
mitis group and the *B. cereus* group all land in GENUS_ONLY by construction,
because their mass genuinely does split. No exception list to maintain, and an
organism nobody anticipated gets the same honest treatment.

MIMt makes this sharper than it would be on a larger database. It is curated
rather than exhaustive, and thin coverage of the mitis group (3 *S. mitis*
references, 5 *S. oralis*, 3 *S. pneumoniae*) means kallisto's EM splits reads
near-arbitrarily between them - which the margin rule correctly reads as "cannot
resolve" rather than as a species call.
"""

from __future__ import annotations

import csv
import logging
from dataclasses import dataclass
from pathlib import Path

from .aggregate import TaxonAbundance, genus_totals

# --- Defaults -------------------------------------------------------------
#
# A genus must hold this share of the assigned 16S signal before it is called at
# all. Below it, the arm stays quiet: trace-level genera in a 16S profile are
# overwhelmingly reagent background, and the reconciler treats SILENT_WEAK as
# "no opinion" rather than "absent", so nothing is lost by declining.
DEFAULT_GENUS_MIN_FRAC = 0.01

# The top species must hold this share OF ITS GENUS to be named. At 0.70 a
# genus split 60/40 between two species reports as GENUS_ONLY, which is the
# right answer for the mitis group.
DEFAULT_SPECIES_DOMINANCE = 0.70

# ...and must beat the runner-up by this much of the genus's mass. The
# dominance test alone would name a species at 0.71 against a runner-up at 0.29;
# requiring a margin as well means a genuinely close call stays at genus.
DEFAULT_SPECIES_MARGIN = 0.20

# An absolute floor in reads, independent of proportion. A genus can hold 100%
# of a near-empty profile; that is not a finding.
DEFAULT_MIN_GENUS_READS = 10.0

from .alignment_metrics import METRIC_COLUMNS, GenusAlignmentMetrics, metric_cells

CALLS_COLUMNS = [
    "sample_id", "taxon", "genus", "species", "rank", "strength",
    "abundance_frac", "genus_frac", "reads_assigned", "raw_reads",
    "copy_number", "within_genus_share", "margin", "n_references",
    "competing_species", "reason",
    # Identity and coverage from a real alignment, per genus. Blank for silent
    # genera and where the alignment stage did not run.
    *METRIC_COLUMNS,
]


@dataclass(frozen=True)
class Marker16SCall:
    """One genus-level verdict from the 16S arm."""

    sample_id: str
    taxon: str
    genus: str | None
    species: str | None
    rank: str                  # "species" | "genus" | "none"
    strength: str              # species_strong | genus_only | silent_weak
    abundance_frac: float      # share of the whole corrected profile
    genus_frac: float          # this genus's share of the profile
    reads_assigned: float      # corrected reads for the named taxon
    raw_reads: float           # before copy-number correction
    copy_number: float
    within_genus_share: float  # top species' share of its genus
    margin: float              # gap to the runner-up species, as genus share
    n_references: int
    competing_species: list[str]
    reason: str


def marker16s_calls(
    sample_id: str,
    taxa: list[TaxonAbundance],
    *,
    genus_min_frac: float = DEFAULT_GENUS_MIN_FRAC,
    species_dominance: float = DEFAULT_SPECIES_DOMINANCE,
    species_margin: float = DEFAULT_SPECIES_MARGIN,
    min_genus_reads: float = DEFAULT_MIN_GENUS_READS,
) -> list[Marker16SCall]:
    """Rank each genus present into species_strong / genus_only / silent_weak.

    One call per genus. The reconciler buckets evidence by genus and collapses
    to the best result per genus, so emitting a genus-level verdict is what it
    is built to consume.
    """
    logger = logging.getLogger("marker16s")
    totals = genus_totals(taxa)
    grand_total = sum(item.reads for item in taxa) or 1.0
    calls: list[Marker16SCall] = []

    for genus, genus_reads in sorted(totals.items(), key=lambda kv: -kv[1]):
        members = sorted(
            (item for item in taxa if item.genus == genus),
            key=lambda item: (-item.reads, item.taxon),
        )
        genus_frac = genus_reads / grand_total
        named = [item for item in members if item.species]
        top = named[0] if named else None
        runner = named[1] if len(named) > 1 else None

        share = (top.reads / genus_reads) if (top and genus_reads) else 0.0
        runner_share = (runner.reads / genus_reads) if (runner and genus_reads) else 0.0
        margin = share - runner_share

        # Below either floor the arm declines to comment. Proportion alone is
        # not enough: a genus can hold all of a near-empty profile.
        #
        # The read floor applies to COPY-NUMBER-CORRECTED reads, so in raw reads
        # it scales with copy number (~42 at the 4.2 fallback, 10 for a
        # single-copy organism). Corrected reads are always <= raw, so moving
        # the floor to raw reads (``genus_totals(taxa, raw=True)``) would need
        # the value re-picked.
        if genus_frac < genus_min_frac or genus_reads < min_genus_reads:
            strength, rank, species = "silent_weak", "none", None
            taxon_label = genus
            # Name only the clause(s) that actually failed.
            failed = []
            if genus_frac < genus_min_frac:
                failed.append(f"frac {genus_frac:.4f} < {genus_min_frac}")
            if genus_reads < min_genus_reads:
                failed.append(f"corrected reads {genus_reads:.1f} < {min_genus_reads}")
            reason = f"genus below reporting floor ({', '.join(failed)})"
        elif top is not None and share >= species_dominance and margin >= species_margin:
            strength, rank, species = "species_strong", "species", top.species
            # TaxonAbundance.taxon is already "Genus species" where the epithet
            # resolved, so it is the label; fall back to the genus otherwise.
            taxon_label = top.taxon or genus
            reason = (
                f"one species holds {share:.0%} of the genus, "
                f"{margin:.0%} clear of the next"
            )
        else:
            strength, rank, species = "genus_only", "genus", None
            taxon_label = genus
            if top is None:
                reason = "no species-level assignment within the genus"
            else:
                reason = (
                    f"genus mass split across species "
                    f"(top {share:.0%}, margin {margin:.0%}); "
                    f"16S cannot resolve these"
                )

        calls.append(Marker16SCall(
            sample_id=sample_id,
            taxon=taxon_label,
            genus=genus,
            species=species,
            rank=rank,
            strength=strength,
            abundance_frac=(top.reads / grand_total) if (top and rank == "species") else genus_frac,
            genus_frac=genus_frac,
            reads_assigned=top.reads if (top and rank == "species") else genus_reads,
            raw_reads=top.raw_reads if (top and rank == "species") else sum(m.raw_reads for m in members),
            copy_number=top.copy_number if top else 0.0,
            within_genus_share=share,
            margin=margin,
            n_references=top.n_references if top else sum(m.n_references for m in members),
            competing_species=[m.taxon for m in named[:5]],
            reason=reason,
        ))

    reported = sum(1 for c in calls if c.strength != "silent_weak")
    logger.info(
        "marker16s_calls_status=complete genera=%d reported=%d species_strong=%d",
        len(calls), reported,
        sum(1 for c in calls if c.strength == "species_strong"),
    )
    return calls


def write_marker16s_calls(
    calls: list[Marker16SCall], path: Path,
    metrics: dict[str, GenusAlignmentMetrics] | None = None,
) -> Path:
    """Write ``marker16s_calls.csv`` - the arm's output, and the adapter's input.

    ``metrics`` is the alignment stage's result keyed by genus; a genus it does
    not cover gets blank metric cells.
    """
    metrics = metrics or {}
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CALLS_COLUMNS)
        writer.writeheader()
        for call in calls:
            writer.writerow({
                **metric_cells(metrics.get(call.genus) if call.genus else None),
                "sample_id": call.sample_id,
                "taxon": call.taxon,
                "genus": call.genus or "",
                "species": call.species or "",
                "rank": call.rank,
                "strength": call.strength,
                "abundance_frac": f"{call.abundance_frac:.6g}",
                "genus_frac": f"{call.genus_frac:.6g}",
                "reads_assigned": f"{call.reads_assigned:.3f}",
                "raw_reads": f"{call.raw_reads:.3f}",
                "copy_number": f"{call.copy_number:.3f}",
                "within_genus_share": f"{call.within_genus_share:.6g}",
                "margin": f"{call.margin:.6g}",
                "n_references": call.n_references,
                "competing_species": "; ".join(call.competing_species),
                "reason": call.reason,
            })
    return path
