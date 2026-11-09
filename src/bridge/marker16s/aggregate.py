"""Roll per-reference abundances up to per-taxon abundances.

kallisto reports one row per reference sequence. MIMt holds several sequences
for most organisms, so the first job is to sum them per taxon. The second is to
correct a bias that is specific to 16S and would otherwise corrupt every
comparison the arm makes.

The 16S copy-number problem
---------------------------
The 16S rRNA gene sits in an *rrn operon*, and bacteria carry between 1 and
about 15 copies of that operon depending on species. A genome with 7 copies
contributes roughly 7 times the 16S reads of a genome with 1 copy at the same
cell density, so raw 16S abundance overstates multi-copy organisms in direct
proportion to their copy number.

This is the standard correction in the amplicon literature (rrnDB catalogues
copy number; CopyRighter and similar tools apply it): weight each taxon by the
inverse of its estimated copy number. Copy number is phylogenetically
conserved, so a genus-level estimate is a reasonable fallback where the species
is unknown, and a database-wide mean is the last resort.

It matters here for a specific reason. The rMLST arm targets ~53 ribosomal
*protein* genes, which are single-copy. Comparing an uncorrected 16S abundance
against an rMLST abundance is therefore biased by exactly the copy number - so
any cross-arm magnitude check (see ``rank_marker16s`` and the RPM-r gate) is
meaningless without this step.

The correction is applied transparently: both the raw and corrected figures are
carried through, so a reader can always see what it did.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import pandas as pd

from .reference import ReferenceTaxon

# Fallback when neither the species nor its genus is in the copy-number table.
# The median across sequenced bacteria is a little over 4; 4.2 is the commonly
# cited rrnDB mean and is what tools default to.
DEFAULT_COPY_NUMBER = 4.2


@dataclass
class TaxonAbundance:
    """One organism's share of the 16S signal."""

    taxon: str
    genus: str | None
    species: str | None
    reads: float                 # copy-number-corrected, summed over references
    raw_reads: float             # what kallisto actually assigned
    copy_number: float           # divisor used
    n_references: int            # MIMt sequences contributing
    abundance_frac: float = 0.0  # share of the corrected total
    references: list[str] = field(default_factory=list)
    # How the divisor was arrived at: "species", "genus" or "default". Without
    # this a table that matched nothing is indistinguishable from no table at
    # all - both silently divide everything by DEFAULT_COPY_NUMBER, which
    # cancels out of every ratio and leaves the bias fully intact while the
    # output still says "copy_number_corrected".
    copy_number_source: str = "default"


def _copy_number(
    taxon: ReferenceTaxon, table: dict[str, float] | None
) -> tuple[float, str]:
    """Estimated rrn copy number: species, else genus, else the global mean.

    Returns the divisor and which level supplied it.
    """
    if not table:
        return DEFAULT_COPY_NUMBER, "default"
    if taxon.species_full:
        value = table.get(taxon.species_full.lower())
        if value:
            return float(value), "species"
    if taxon.genus:
        value = table.get(taxon.genus.lower())
        if value:
            return float(value), "genus"
    return DEFAULT_COPY_NUMBER, "default"


def aggregate_taxa(
    abundance: pd.DataFrame,
    reference: dict[str, ReferenceTaxon],
    *,
    copy_numbers: dict[str, float] | None = None,
    min_reads: float = 1.0,
) -> list[TaxonAbundance]:
    """Sum kallisto's per-reference counts into per-taxon abundances.

    Taxa are keyed on the finest resolved name - ``"Genus species"`` where the
    epithet parsed, otherwise the genus alone. References whose header could not
    be parsed are collected under ``"unassigned"`` rather than dropped, so the
    total always reconciles against what kallisto reported.

    ``min_reads`` discards taxa below a floor before fractions are computed.
    kallisto's EM spreads a trace of mass across many references, and without a
    floor those fragments dominate the taxon count while contributing nothing.
    """
    logger = logging.getLogger("marker16s")
    buckets: dict[str, dict] = {}

    for row in abundance.itertuples(index=False):
        counts = float(row.est_counts)
        if counts <= 0:
            continue
        taxon = reference.get(str(row.target_id))
        if taxon is None:
            key, genus, species = "unassigned", None, None
            divisor, source = DEFAULT_COPY_NUMBER, "default"
        else:
            key = taxon.species_full or taxon.genus or "unassigned"
            genus, species = taxon.genus, taxon.species
            divisor, source = _copy_number(taxon, copy_numbers)

        bucket = buckets.setdefault(key, dict(
            taxon=key, genus=genus, species=species, raw=0.0, corrected=0.0,
            divisor_weighted=0.0, n=0, refs=[], source=source,
        ))
        bucket["raw"] += counts
        bucket["corrected"] += counts / divisor
        bucket["divisor_weighted"] += divisor * counts
        bucket["n"] += 1
        if len(bucket["refs"]) < 5:
            bucket["refs"].append(str(row.target_id))

    results: list[TaxonAbundance] = []
    for bucket in buckets.values():
        if bucket["raw"] < min_reads:
            continue
        results.append(TaxonAbundance(
            taxon=bucket["taxon"],
            genus=bucket["genus"],
            species=bucket["species"],
            reads=bucket["corrected"],
            raw_reads=bucket["raw"],
            # Read-weighted mean divisor, so the reported copy number reflects
            # the references that actually carried the signal.
            copy_number=bucket["divisor_weighted"] / bucket["raw"],
            n_references=bucket["n"],
            references=bucket["refs"],
            copy_number_source=bucket["source"],
        ))

    total = sum(item.reads for item in results)
    for item in results:
        item.abundance_frac = (item.reads / total) if total else 0.0

    results.sort(key=lambda item: (-item.reads, item.taxon))
    logger.info(
        "aggregate_status=complete taxa=%d assigned_reads=%.0f corrected_total=%.1f",
        len(results), sum(i.raw_reads for i in results), total,
    )
    breakdown = copy_number_coverage(results)
    # Logged every run, not only when a table is supplied: "everything fell back
    # to the default" is exactly the state that must not pass unnoticed, because
    # a uniform divisor cancels out of every downstream ratio.
    logger.info(
        "copy_number_status=%s table=%s species=%d genus=%d default=%d "
        "reads_corrected_frac=%.3f",
        "applied" if copy_numbers else "absent",
        len(copy_numbers) if copy_numbers else 0,
        breakdown["species"], breakdown["genus"], breakdown["default"],
        breakdown["reads_resolved_frac"],
    )
    if copy_numbers and breakdown["reads_resolved_frac"] < 0.5:
        logger.warning(
            "copy_number_status=mostly_unmatched - a table was supplied but "
            "%.0f%% of assigned reads still fell back to %.1f. Check the table's "
            "names against the reference's.",
            100 * (1 - breakdown["reads_resolved_frac"]), DEFAULT_COPY_NUMBER,
        )
    return results


def copy_number_coverage(taxa: list[TaxonAbundance]) -> dict[str, float]:
    """How much of the signal got a real copy number, and at what rank.

    ``reads_resolved_frac`` weights by reads rather than by taxon: a table that
    resolves one dominant organism and misses fifty trace ones has done most of
    the job, and a taxon count would say the opposite.
    """
    counts = {"species": 0, "genus": 0, "default": 0}
    resolved_reads = 0.0
    for item in taxa:
        counts[item.copy_number_source] = counts.get(item.copy_number_source, 0) + 1
        if item.copy_number_source != "default":
            resolved_reads += item.raw_reads
    total_reads = sum(item.raw_reads for item in taxa)
    return {
        **counts,
        "reads_resolved_frac": (resolved_reads / total_reads) if total_reads else 0.0,
    }


def genus_totals(taxa: list[TaxonAbundance], *, raw: bool = False) -> dict[str, float]:
    """Reads summed per genus.

    Corrected by default - that is the denominator for the margin rule, which
    is a composition question and so wants genome equivalents.

    ``raw=True`` gives what kallisto actually assigned, for the questions that
    are about *sequencing depth* rather than composition. A detection floor is
    one of those: "enough reads for this estimate to be stable" is a property of
    the run, not of the organism's genome, and a floor applied to corrected
    reads silently varies 15-fold across taxa with copy number.
    """
    totals: dict[str, float] = {}
    for item in taxa:
        if item.genus:
            value = item.raw_reads if raw else item.reads
            totals[item.genus] = totals.get(item.genus, 0.0) + value
    return totals
