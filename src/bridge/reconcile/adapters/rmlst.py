"""``species_calls.csv`` (the rMLST arm) -> ``list[RmlstStrengthResult]``.

This is where the rMLST verdict is decided. The rMLST arm reports one row per
species with locus counts, a call method and an abundance - no rank vocabulary
and no support percentage - and it can report several species for one sample.

How the three strengths are decided
-----------------------------------
The mapping keys on the two fields the reproducibility work showed to be stable
across runs, and avoids the ones that are not:

``reported`` is True
    The species cleared the unambiguous-locus bar. SPECIES_SUPPORTED.

``alternative_to`` is set
    A demoted near neighbour: it held up as an allele profile but fell short on
    exclusive loci, and names the organism it is a second name *for*. That is
    precisely "the marker cannot separate these", so it becomes GENUS_ONLY -
    corroborating the genus without the reconciler treating it as a second
    organism. Reporting it as SILENT would throw away a real observation, and as
    SPECIES_SUPPORTED would over-call.
    The species name is kept on the result and reaches the report as an
    ``rmlst_alternative`` row; it is the name a reader needs to see.

everything else
    EM-fallback background. SILENT.

This is deliberate. Duplicate runs of the same samples showed species calls and
``alternative_to`` labelling identical in every case, while the underlying locus
counts moved in 8 of 48 calls and abundances in 42 of 48. Keying the adapter on
the stable layer means a borderline sample cannot flip verdict between runs
because a count wobbled by one.
"""

from __future__ import annotations

import csv
import logging
from dataclasses import dataclass
from pathlib import Path

from ..config import Config
from ..strength import RmlstStrength, RmlstStrengthResult
from ..taxonomy import parse_taxon

REQUIRED_COLUMNS = {"sample_id", "species", "reported", "call_method"}


@dataclass(frozen=True)
class RmlstRow:
    """One row of the rMLST arm's ``species_calls.csv``."""

    sample_id: str
    species: str
    reported: bool
    call_method: str
    pi: float | None
    exact_profile_loci: int
    exact_profile_unambiguous_loci: int
    unambiguous_loci: int
    ambiguous_loci: int
    em_evidence_unambiguous_loci: int | None
    em_evidence_ambiguous_loci: int | None
    unambiguous_est_counts: float | None
    ambiguous_est_counts: float | None
    # Coverage of the loci that name this species ALONE - the same loci
    # unambiguous_loci counts. Blank where the species has no unambiguous
    # evidence, which is every em_fallback row. See
    # bridge.rmlst.species.species_coverage. The 16S arm's own breadth/depth
    # (one gene, one reference) live on Marker16SRow and are not comparable.
    mean_callable_breadth: float | None
    min_callable_breadth: float | None
    median_projected_depth: float | None
    min_projected_depth: float | None
    total_callable_bases: float | None
    total_allele_length: float | None
    supporting_fragments: float | None
    alternative_to: str

    @property
    def total_est_counts(self) -> float | None:
        """Reads behind this call, for the RPM-r gate. None where unknown.

        Both components count: an allele shared between near neighbours is
        still material in the tube, and excluding it would understate an
        organism purely because it has close relatives in the database.

        **None is not zero.** The rMLST arm leaves these cells blank for a
        species called from its allele profile whose consensus groups never
        reached EM, so no abundance was estimated. Reading that as zero reads
        would fail the RPM-r floor and silence a real call. Returning None
        routes it to the gate's "no read count known" path, which passes the
        evidence through ungated rather than destroying it.
        """
        if self.unambiguous_est_counts is None and self.ambiguous_est_counts is None:
            return None
        return (self.unambiguous_est_counts or 0.0) + (self.ambiguous_est_counts or 0.0)


def _number(row: dict, key: str, cast, default):
    try:
        value = row.get(key, "")
        return cast(value) if str(value).strip() != "" else default
    except (TypeError, ValueError):
        return default


def load_rmlst_calls(path: Path) -> list[RmlstRow]:
    """Parse ``species_calls.csv``.

    A file without ``exact_profile_unambiguous_loci`` and ``alternative_to``
    has no demoted rows to recognise, so it degrades to "reported or silent"
    rather than failing.
    """
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = REQUIRED_COLUMNS - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"{path} is not an rMLST species_calls.csv: "
                f"missing {sorted(missing)}"
            )
        rows = []
        for row in reader:
            rows.append(RmlstRow(
                sample_id=(row.get("sample_id") or "").strip(),
                species=(row.get("species") or "").strip(),
                reported=str(row.get("reported", "")).strip().lower() == "true",
                call_method=(row.get("call_method") or "").strip(),
                # Default None, not 0.0: a blank abundance cell means "not
                # estimated", and the two must stay distinguishable all the way
                # to the gates. See RmlstRow.total_est_counts.
                pi=_number(row, "pi", float, None),
                exact_profile_loci=_number(row, "exact_profile_loci", int, 0),
                exact_profile_unambiguous_loci=_number(
                    row, "exact_profile_unambiguous_loci", int, 0),
                unambiguous_loci=_number(
                    row, "allele_evidence_unambiguous_loci", int, 0),
                ambiguous_loci=_number(
                    row, "allele_evidence_ambiguous_loci", int, 0),
                em_evidence_unambiguous_loci=_number(
                    row, "em_evidence_unambiguous_loci", int, None),
                em_evidence_ambiguous_loci=_number(
                    row, "em_evidence_ambiguous_loci", int, None),
                unambiguous_est_counts=_number(
                    row, "unambiguous_est_counts", float, None),
                ambiguous_est_counts=_number(
                    row, "ambiguous_est_counts", float, None),
                # Default None throughout, so a species_calls.csv without the
                # coverage columns degrades to blank cells, not fabricated zeros.
                mean_callable_breadth=_number(
                    row, "mean_callable_breadth", float, None),
                min_callable_breadth=_number(
                    row, "min_callable_breadth", float, None),
                median_projected_depth=_number(
                    row, "median_projected_depth", float, None),
                min_projected_depth=_number(
                    row, "min_projected_depth", float, None),
                total_callable_bases=_number(
                    row, "total_callable_bases", float, None),
                total_allele_length=_number(
                    row, "total_allele_length", float, None),
                supporting_fragments=_number(
                    row, "supporting_fragments", float, None),
                alternative_to=(row.get("alternative_to") or "").strip(),
            ))
    return rows


def rmlst_results(
    rows: list[RmlstRow],
    cfg: Config,
    *,
    include_demoted: bool = True,
) -> list[RmlstStrengthResult]:
    """Build the reconciler's rMLST side.

    One result per species, not per sample. The reconciler groups by genus and
    collapses to the best result per genus, so several species from one sample
    is the shape it expects - and a polymicrobial sample must be able to carry
    one result per organism.
    """
    logger = logging.getLogger("bridge.reconcile.adapters.rmlst")
    results: list[RmlstStrengthResult] = []

    for index, row in enumerate(rows):
        if not row.species:
            continue

        if row.reported:
            strength = RmlstStrength.SPECIES_SUPPORTED
            taxon = parse_taxon(row.species, cfg)
            flags: list[str] = []
            reason = (f"{row.call_method}: {row.unambiguous_loci} unambiguous loci "
                      f"({row.exact_profile_loci} compatible)")
        elif row.alternative_to and include_demoted:
            # A second acceptable name for an organism already called. The
            # species is kept: `basis=rmlst_alternative` already says the row is
            # not a reported call, and a reader needs to see "Streptococcus
            # oralis, an alternative name for Streptococcus mitis", not
            # "Streptococcus". Keeping it also gives the RPM gate the species key it files read
            # counts under, so a demoted name is measured against background on
            # its own counts rather than on its genus's total.
            strength = RmlstStrength.GENUS_ONLY
            taxon = parse_taxon(row.species, cfg)
            flags = ["rmlst_alternative_name"]
            reason = (f"demoted near neighbour of {row.alternative_to}; "
                      f"only {row.unambiguous_loci} unambiguous loci - the "
                      f"marker cannot separate these")
        else:
            continue  # EM-fallback background: no opinion, not a negative.

        results.append(RmlstStrengthResult(
            strength=strength,
            taxon=taxon,
            # The rMLST arm has no support percentage. pi is a mixture proportion
            # and is carried as a percentage for the report only - the
            # reconciler never reads it. None where no abundance was estimated.
            support=None if row.pi is None else row.pi * 100.0,
            n_loci=row.unambiguous_loci,
            flags=flags,
            reason=reason,
            coverage=None,
            group_index=index,
        ))

    logger.info(
        "adapter_status=complete rows=%d results=%d species_supported=%d genus_only=%d",
        len(rows), len(results),
        sum(1 for r in results if r.strength is RmlstStrength.SPECIES_SUPPORTED),
        sum(1 for r in results if r.strength is RmlstStrength.GENUS_ONLY),
    )
    return results
