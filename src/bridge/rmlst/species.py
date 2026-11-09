from __future__ import annotations

from collections import defaultdict
import re

import numpy as np
import pandas as pd

from .models import ConsensusCall
from .seqdef import base_species_label


ACCEPTED_ASSIGNMENTS = {
    "complete_exact_match",
    "unique_partial_match",
    "ambiguous_partial_match",
}
GENUS_PLACEHOLDER_RE = re.compile(r"^([^\s]+)\s+sp\.$", re.IGNORECASE)
# Unambiguous loci a species needs before it may be reported. Mirrors
# models.Config.min_unambiguous_loci; kept here so the calling convention has a
# default without importing Config into the scoring code.
DEFAULT_MIN_UNAMBIGUOUS_LOCI = 3

# Whether a locus must carry an ACCEPTED (exact) consensus match before it may
# count as evidence for a species profile.
#
# False - the default - counts any locus whose consensus resolves to exactly one
# species, whether or not the assembled sequence exactly matched a catalogued
# allele. The strict setting silently penalises low-abundance organisms, which
# typically show broad locus support but imperfect assembly (many loci naming
# the organism alone, few exact allele matches).
#
# Scored against a 450-sample synthetic mixture panel at a bar of 3: 0 of 87
# spurious near-neighbour calls readmitted, 0 of 900 true calls lost. (At a bar
# of 2 one spurious call returns, which is why the bar stays at 3.)
#
# Kept switchable because this is a scientific threshold.
DEFAULT_REQUIRE_ACCEPTED_CONSENSUS = False
EM_COLUMNS = [
    "species", "pi",
    "allele_evidence_unambiguous_loci", "allele_evidence_ambiguous_loci",
    "em_evidence_unambiguous_loci", "em_evidence_ambiguous_loci",
    "unambiguous_est_counts", "ambiguous_est_counts", "detected",
]
# How well assembled the evidence behind a call actually was, summarised across
# the loci that name the species ALONE - the same loci
# allele_evidence_unambiguous_loci counts. See species_coverage.
#
# Every one of these is blank, never zero, for a species with no unambiguous
# evidence (an em_fallback row, or a demoted near neighbour that fell short of
# the bar). Blank means "not measured"; zero would mean "measured, and there was
# nothing there". Same rule as pi and the est_counts columns.
COVERAGE_COLUMNS = [
    "mean_callable_breadth", "min_callable_breadth",
    "median_projected_depth", "min_projected_depth",
    "total_callable_bases", "total_allele_length", "supporting_fragments",
]
FINAL_COLUMNS = [
    "sample_id", "species", "reported", "call_method", "exact_profile_loci",
    # How many of exact_profile_loci map to this species ALONE. The two differ
    # sharply for a near neighbour: 27 compatible loci against 1 exclusive one.
    "exact_profile_unambiguous_loci",
    "pi", "allele_evidence_unambiguous_loci",
    "allele_evidence_ambiguous_loci",
    "em_evidence_unambiguous_loci", "em_evidence_ambiguous_loci",
    "unambiguous_est_counts", "ambiguous_est_counts",
    *COVERAGE_COLUMNS,
    # Set on an unreported row that is a second acceptable name for a reported
    # organism, rather than a separate organism that failed the bar.
    "alternative_to",
]


def _placeholder_genus(species: str) -> str | None:
    match = GENUS_PLACEHOLDER_RE.fullmatch(base_species_label(species))
    return match.group(1).casefold() if match else None


def _definitive_genus(species: str) -> str:
    return base_species_label(species).split(None, 1)[0].casefold()


def merge_consensus_calls(calls: list[ConsensusCall]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for locus, locus_calls in _group_by(calls, lambda call: call.locus_id).items():
        groups: list[list[ConsensusCall]] = []
        for call in sorted(locus_calls, key=lambda value: value.runtime_cluster_id):
            destination = None
            for group in groups:
                representative = group[0]
                same_sequence = bool(call.consensus) and call.consensus == representative.consensus
                if same_sequence:
                    destination = group
                    break
            if destination is None:
                groups.append([call])
            else:
                destination.append(call)
        for ordinal, group in enumerate(groups, 1):
            # Identical consensuses represent one sequence-level observation.
            # Keep the assignment from the strongest runtime cluster and
            # reassign the duplicate clusters' quantitative evidence to it.
            winner = sorted(
                group,
                key=lambda call: (
                    -call.est_counts,
                    -call.exclusive_fragments,
                    -call.assigned_fragments,
                    call.runtime_cluster_id,
                ),
            )[0]
            compatible = set(winner.compatible_alleles)
            species = set(winner.compatible_species)
            rows.append({
                "sample_id": winner.sample_id,
                "locus_id": locus,
                "consensus_group_id": f"{locus}_CG{ordinal:04d}",
                "source_subconsensuses": ";".join(
                    call.subconsensus_id or call.runtime_cluster_id for call in group
                ),
                "source_runtime_clusters": ";".join(call.runtime_cluster_id for call in group),
                "source_representatives": ";".join(sorted({call.representative for call in group})),
                "merged_cluster_count": len(group),
                "combined_est_counts": sum(call.est_counts for call in group),
                "combined_assigned_fragments": sum(call.assigned_fragments for call in group),
                "exclusive_fragments": sum(call.exclusive_fragments for call in group),
                "ambiguous_fragments": sum(call.ambiguous_fragments for call in group),
                "consensus_status": winner.status,
                "consensus_method": winner.consensus_method,
                "consensus": winner.consensus,
                "callable_bases": winner.callable_bases,
                "allele_length": winner.allele_length,
                "callable_breadth": winner.callable_breadth,
                "winner_haplotype_proportion": winner.haplotype_proportion,
                "supporting_fragments": sum(call.supporting_fragments for call in group),
                "projected_median_depth": winner.projected_median_depth,
                "projected_breadth": winner.projected_breadth,
                "variant_evidence_passed": winner.variant_evidence_passed,
                "composite_from_members": any(
                    call.composite_from_members for call in group
                ),
                "compatible_alleles": ";".join(sorted(compatible)),
                "compatible_allele_count": len(compatible),
                "compatible_species": ";".join(sorted(species)),
                "compatible_species_count": len(species),
                "allele_search_scope": winner.allele_search_scope,
            })
    return pd.DataFrame(rows)


def _group_by(values, key):
    result = defaultdict(list)
    for value in values:
        result[key(value)].append(value)
    return result


def _species_set(value: object) -> set[str]:
    return {item for item in str(value).split(";") if item and item != "nan"}


def _resolved_evidence(groups: pd.DataFrame) -> pd.DataFrame:
    """Resolve genus placeholders without making generic evidence exclusive."""
    evidence = groups.copy()
    evidence["species_set"] = evidence["compatible_species"].map(_species_set)
    evidence["species_set"] = evidence["species_set"].map(
        lambda values: {
            value for value in values if _placeholder_genus(value) is None
        } or set(values)
    )

    definitive_by_genus: dict[str, set[str]] = defaultdict(set)
    accepted = evidence[evidence["consensus_status"].isin(ACCEPTED_ASSIGNMENTS)]
    for values in accepted["species_set"]:
        if len(values) == 1:
            species = next(iter(values))
            if _placeholder_genus(species) is None:
                definitive_by_genus[_definitive_genus(species)].add(species)

    resolved: list[set[str]] = []
    forced_ambiguous: list[bool] = []
    for values in evidence["species_set"]:
        genera = [_placeholder_genus(value) for value in values]
        if genera and all(genus is not None for genus in genera):
            candidates = {
                species
                for genus in genera
                for species in definitive_by_genus.get(str(genus), set())
            }
            if candidates:
                resolved.append(candidates)
                forced_ambiguous.append(True)
                continue
        resolved.append(set(values))
        forced_ambiguous.append(False)
    evidence["species_set"] = resolved
    composite = (
        evidence["composite_from_members"].fillna(False).astype(bool)
        if "composite_from_members" in evidence
        else pd.Series(False, index=evidence.index)
    )
    evidence["forced_ambiguous"] = [
        placeholder or is_composite
        for placeholder, is_composite in zip(forced_ambiguous, composite)
    ]
    return evidence


def direct_profile_evidence(
    groups: pd.DataFrame,
    minimum_loci: int,
    min_unambiguous_loci: int = DEFAULT_MIN_UNAMBIGUOUS_LOCI,
    require_accepted_consensus: bool = DEFAULT_REQUIRE_ACCEPTED_CONSENSUS,
) -> dict[str, dict[str, int]]:
    """Species callable directly from their allele profile, with the evidence.

    Two counts matter and they are not the same number:

    ``unambiguous_loci``
        loci whose consensus maps to this species *alone*. This is the real
        evidence for the call.
    ``profile_loci``
        loci whose consensus is merely *compatible* with this species - it may
        equally support several near neighbours.

    Only requiring one of the former let a near neighbour ride in on the
    latter. The profile-wide intersection collapses to ``{species}`` as soon as
    a single locus is exclusive to it, because
    ``{X} & {X, Y, Z} & ... == {X}`` - so one exclusive locus plus
    ``minimum_loci - 1`` loci of mere compatibility was enough to report a
    species holding a thousandth of the sample. Requiring ``min_unambiguous_loci``
    exclusive loci is what separates an organism that is present from one whose
    name simply appears in other organisms' shared evidence.

    Note ``unambiguous_loci <= len(per_locus)`` always, so a
    ``min_unambiguous_loci`` at or above ``minimum_loci`` subsumes the profile
    breadth test rather than competing with it.
    """
    return {
        species: counts
        for species, counts in _profile_candidates(
            groups, minimum_loci, require_accepted_consensus).items()
        if counts["unambiguous_loci"] >= min_unambiguous_loci
    }


def _profile_candidates(
    groups: pd.DataFrame,
    minimum_loci: int,
    require_accepted_consensus: bool = DEFAULT_REQUIRE_ACCEPTED_CONSENSUS,
) -> dict[str, dict[str, int]]:
    """Every species whose allele profile holds up, before the bar is applied.

    Split out of :func:`direct_profile_evidence` so the species that fail *only*
    the unambiguous-locus bar can be recovered rather than discarded - those are
    the second acceptable names for an organism already called, and
    :func:`demoted_profile_species` is what turns them into ``alternative_to``.

    A species with no exclusive locus anywhere is skipped outright, which keeps
    the expensive per-locus loop short and is sound: such a species has no
    evidence naming it alone, so it is neither callable nor a second name for
    anything.
    """
    evidence = _resolved_evidence(groups)
    usable = ~evidence["forced_ambiguous"]
    # A locus that carried no reads is not evidence of anything (the same filter
    # as ``em_species``). Consensus can be assembled for a locus the EM then
    # assigns no abundance; without this a negative control can report a
    # species on zero-count loci when accepted consensus is not required.
    usable &= evidence["combined_est_counts"].gt(0)
    if require_accepted_consensus:
        usable &= evidence["consensus_status"].isin(ACCEPTED_ASSIGNMENTS)
    accepted = evidence[usable].copy()
    candidates = sorted({species for values in accepted["species_set"] for species in values})
    passed: dict[str, dict[str, int]] = {}
    for species in candidates:
        unambiguous_loci = accepted.loc[
            accepted["species_set"].map(lambda values: values == {species}),
            "locus_id",
        ].nunique()
        if unambiguous_loci < 1:
            continue
        per_locus: list[set[str]] = []
        for _, locus in accepted.groupby("locus_id", sort=False):
            supporting = [values for values in locus["species_set"] if species in values]
            if supporting:
                per_locus.append(set.union(*supporting))
        if len(per_locus) < minimum_loci:
            continue
        compatible_across_profile = set.intersection(*per_locus)
        if compatible_across_profile == {species}:
            passed[species] = {
                "profile_loci": len(per_locus),
                "unambiguous_loci": int(unambiguous_loci),
            }
    return passed


def demoted_profile_species(
    groups: pd.DataFrame,
    minimum_loci: int,
    min_unambiguous_loci: int = DEFAULT_MIN_UNAMBIGUOUS_LOCI,
    require_accepted_consensus: bool = DEFAULT_REQUIRE_ACCEPTED_CONSENSUS,
) -> dict[str, dict[str, int]]:
    """Species that hold up as a profile but fall short of the reporting bar.

    The complement of :func:`direct_profile_evidence` over the same candidate
    set: the profile geometry passed, and the *only* thing missing is exclusive
    evidence. That is what a near neighbour looks like - on the synthetic
    mixture panel every one of these carried exactly one exclusive locus, while
    the 16,790 EM-fallback rows carried none at all.

    These are not reported, but they are not noise either: they are the same
    organism under a name the evidence cannot exclude, which is why they earn an
    ``alternative_to`` label and a plain EM-fallback row does not.
    """
    return {
        species: counts
        for species, counts in _profile_candidates(
            groups, minimum_loci, require_accepted_consensus).items()
        if counts["unambiguous_loci"] < min_unambiguous_loci
    }


def direct_profile_species(
    groups: pd.DataFrame,
    minimum_loci: int,
    min_unambiguous_loci: int = DEFAULT_MIN_UNAMBIGUOUS_LOCI,
) -> dict[str, int]:
    """``{species: profile_loci}`` - the compatibility-breadth count only.

    Use :func:`direct_profile_evidence` where the unambiguous count matters,
    because this one overstates how well supported a call is.
    """
    return {
        species: counts["profile_loci"]
        for species, counts in direct_profile_evidence(
            groups, minimum_loci, min_unambiguous_loci
        ).items()
    }


def _scored_evidence(groups: pd.DataFrame) -> pd.DataFrame:
    """Consensus groups that count as evidence, with the ambiguity flags set.

    One definition, shared by :func:`em_species` and :func:`species_coverage`,
    so the coverage figures and ``allele_evidence_unambiguous_loci`` provably
    describe the SAME loci.

    ``combined_est_counts > 0`` is part of the definition, not an optimisation:
    a locus with no estimated abundance is not evidence. See
    :func:`_profile_candidates`.
    """
    evidence = _resolved_evidence(groups)
    evidence = evidence[
        evidence["species_set"].map(bool)
        & evidence["combined_est_counts"].gt(0)
    ].copy()
    if evidence.empty:
        return evidence
    evidence["allele_evidence"] = evidence["consensus_status"].isin(
        ACCEPTED_ASSIGNMENTS
    )
    evidence["source_unambiguous"] = (
        ~evidence["forced_ambiguous"]
        & evidence["species_set"].map(len).eq(1)
    )
    evidence["source_ambiguous"] = ~evidence["source_unambiguous"]
    # "Unambiguous" evidence for a species: an accepted consensus that names
    # that species and nothing else. This is what unambiguous_est_counts sums,
    # what allele_evidence_unambiguous_loci counts, and what species_coverage
    # measures breadth and depth over.
    evidence["unambiguous"] = (
        evidence["allele_evidence"] & evidence["source_unambiguous"]
    )
    evidence["ambiguous"] = ~evidence["unambiguous"]
    return evidence


def species_coverage(groups: pd.DataFrame) -> dict[str, dict[str, float]]:
    """How well assembled each species' unambiguous evidence was.

    The per-locus breadth and depth in ``consensus_groups.csv``, summarised per
    species and carried into ``species_calls.csv`` as a quality signal.

    Aggregated over the loci that name the species ALONE, which is the evidence
    the call actually rests on - loci shared with a near neighbour would let a
    well-covered relative inflate the figure.

    **One value per locus, not per consensus group.** A locus that resolved into
    several haplotypes contributes several groups; the best-assembled one speaks
    for the locus, so the loci aggregated here are exactly the loci
    ``allele_evidence_unambiguous_loci`` counts. Aggregating over groups instead
    would quietly double-count a locus in a within-species mixture.

    Returns ``{species: {metric: value}}``; a species with no unambiguous
    evidence is simply absent, which renders as a blank cell rather than a zero.
    """
    evidence = _scored_evidence(groups)
    if evidence.empty:
        return {}
    unambiguous = evidence[evidence["unambiguous"]].copy()
    if unambiguous.empty:
        return {}

    # Safe because "unambiguous" means len(species_set) == 1.
    unambiguous["_species"] = unambiguous["species_set"].map(
        lambda values: next(iter(values))
    )
    numeric = ("callable_breadth", "projected_median_depth", "callable_bases",
               "allele_length", "supporting_fragments")
    for column in numeric:
        unambiguous[column] = (
            pd.to_numeric(unambiguous.get(column), errors="coerce")
            if column in unambiguous else float("nan")
        )
    # The locus's best-assembled unambiguous consensus speaks for that locus.
    # NaN breadth sorts last, so a locus is only represented by an unmeasured
    # group when every group at that locus was unmeasured.
    per_locus = (
        unambiguous.sort_values("callable_breadth", ascending=False,
                                na_position="last")
        .drop_duplicates(subset=["_species", "locus_id"], keep="first")
    )

    def _stat(values: pd.Series, how: str) -> float:
        # projected_median_depth is NaN for some consensus methods, so drop
        # missing values rather than poisoning the whole aggregate; all-missing
        # stays NaN, which is the honest answer.
        present = values.dropna()
        if present.empty:
            return float("nan")
        return float(getattr(present, how)())

    coverage: dict[str, dict[str, float]] = {}
    for species, rows in per_locus.groupby("_species"):
        coverage[str(species)] = {
            "mean_callable_breadth": _stat(rows["callable_breadth"], "mean"),
            "min_callable_breadth": _stat(rows["callable_breadth"], "min"),
            "median_projected_depth": _stat(rows["projected_median_depth"], "median"),
            "min_projected_depth": _stat(rows["projected_median_depth"], "min"),
            "total_callable_bases": _stat(rows["callable_bases"], "sum"),
            "total_allele_length": _stat(rows["allele_length"], "sum"),
            "supporting_fragments": _stat(rows["supporting_fragments"], "sum"),
        }
    return coverage


def _attach_coverage(
    final: pd.DataFrame, coverage: dict[str, dict[str, float]]
) -> pd.DataFrame:
    """Fill the coverage columns from :func:`species_coverage`, blank where absent."""
    final = final.copy()
    for column in COVERAGE_COLUMNS:
        final[column] = [
            coverage.get(species, {}).get(column, float("nan"))
            for species in final["species"]
        ]
    return final


def em_species(
    groups: pd.DataFrame,
    minimum_loci: int,
    min_unambiguous_loci: int = DEFAULT_MIN_UNAMBIGUOUS_LOCI,
) -> pd.DataFrame:
    evidence = _scored_evidence(groups)
    if evidence.empty:
        return pd.DataFrame(columns=EM_COLUMNS)
    totals = evidence.groupby("locus_id")["combined_est_counts"].transform("sum")
    evidence["weight"] = evidence["combined_est_counts"] / totals
    species = sorted({value for values in evidence["species_set"] for value in values})
    index = {value: number for number, value in enumerate(species)}
    pi = np.ones(len(species), dtype=float) / len(species)
    for _ in range(500):
        previous = pi.copy()
        attributed = np.zeros(len(species), dtype=float)
        for row in evidence.itertuples(index=False):
            indices = [index[value] for value in row.species_set]
            mass = pi[indices]
            fractions = mass / mass.sum() if mass.sum() else np.ones(len(indices)) / len(indices)
            attributed[indices] += fractions * float(row.weight)
        pi = attributed / attributed.sum()
        if np.abs(pi - previous).sum() < 1e-6:
            break
    allele_unambiguous: dict[str, int] = {}
    allele_ambiguous: dict[str, int] = {}
    em_unambiguous: dict[str, int] = {}
    em_ambiguous: dict[str, int] = {}
    unambiguous_counts: dict[str, float] = {}
    ambiguous_counts: dict[str, float] = {}
    for value in species:
        includes_species = evidence["species_set"].map(
            lambda values: value in values
        )
        allele_unambiguous[value] = evidence.loc[
            evidence["allele_evidence"]
            & evidence["source_unambiguous"]
            & includes_species,
            "locus_id",
        ].nunique()
        allele_ambiguous[value] = evidence.loc[
            evidence["allele_evidence"]
            & evidence["source_ambiguous"]
            & includes_species,
            "locus_id",
        ].nunique()
        em_unambiguous[value] = evidence.loc[
            evidence["source_unambiguous"] & includes_species,
            "locus_id",
        ].nunique()
        em_ambiguous[value] = evidence.loc[
            evidence["source_ambiguous"] & includes_species,
            "locus_id",
        ].nunique()
        unambiguous_counts[value] = evidence.loc[
            evidence["unambiguous"] & includes_species, "combined_est_counts"
        ].sum()
        ambiguous_counts[value] = evidence.loc[
            evidence["ambiguous"] & includes_species, "combined_est_counts"
        ].sum()
    return pd.DataFrame({
        "species": species,
        "pi": [float(pi[index[value]]) for value in species],
        "allele_evidence_unambiguous_loci": [
            allele_unambiguous[value] for value in species
        ],
        "allele_evidence_ambiguous_loci": [
            allele_ambiguous[value] for value in species
        ],
        "em_evidence_unambiguous_loci": [
            em_unambiguous[value] for value in species
        ],
        "em_evidence_ambiguous_loci": [
            em_ambiguous[value] for value in species
        ],
        "unambiguous_est_counts": [unambiguous_counts[value] for value in species],
        "ambiguous_est_counts": [ambiguous_counts[value] for value in species],
        "detected": [
            allele_unambiguous[value] >= min_unambiguous_loci for value in species
        ],
    }).sort_values("pi", ascending=False).reset_index(drop=True)


def final_species_calls(
    sample_id: str,
    groups: pd.DataFrame,
    minimum_loci: int,
    min_unambiguous_loci: int = DEFAULT_MIN_UNAMBIGUOUS_LOCI,
    require_accepted_consensus: bool = DEFAULT_REQUIRE_ACCEPTED_CONSENSUS,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    direct = direct_profile_evidence(groups, minimum_loci, min_unambiguous_loci,
                                     require_accepted_consensus)
    demoted = demoted_profile_species(groups, minimum_loci, min_unambiguous_loci,
                                      require_accepted_consensus)
    em = em_species(groups, minimum_loci, min_unambiguous_loci)
    if em.empty:
        em = pd.DataFrame.from_records([
            {
                # These rows exist only because every consensus group of this
                # species had zero estimated abundance and so never reached EM.
                # The locus counts carry the direct call's own evidence, which
                # is what the reporting bar reads. Abundance is NaN (an empty
                # cell), not 0.0: it was not estimated, and downstream gates
                # would read 0.0 as "measured and absent" and silence the call.
                "species": species, "pi": float("nan"),
                "allele_evidence_unambiguous_loci": direct[species]["unambiguous_loci"],
                "allele_evidence_ambiguous_loci": 0,
                "em_evidence_unambiguous_loci": 0,
                "em_evidence_ambiguous_loci": 0,
                "unambiguous_est_counts": float("nan"),
                "ambiguous_est_counts": float("nan"), "detected": False,
            }
            for species in sorted(direct)
        ], columns=EM_COLUMNS)
    rows = []
    for row in em.itertuples(index=False):
        is_direct = row.species in direct
        # Both branches already carry the unambiguous-locus bar: a direct call
        # cleared it in direct_profile_evidence, and `detected` is that same
        # comparison.
        reported = is_direct or bool(row.detected)
        rows.append({
            "sample_id": sample_id,
            "species": row.species,
            "reported": reported,
            "call_method": "exact_allele_profile" if is_direct else "em_fallback",
            "exact_profile_loci": (
                direct[row.species]["profile_loci"] if is_direct else 0
            ),
            "exact_profile_unambiguous_loci": (
                direct[row.species]["unambiguous_loci"] if is_direct else 0
            ),
            "pi": float(row.pi),
            "allele_evidence_unambiguous_loci": int(
                row.allele_evidence_unambiguous_loci
            ),
            "allele_evidence_ambiguous_loci": int(
                row.allele_evidence_ambiguous_loci
            ),
            "em_evidence_unambiguous_loci": int(
                row.em_evidence_unambiguous_loci
            ),
            "em_evidence_ambiguous_loci": int(
                row.em_evidence_ambiguous_loci
            ),
            "unambiguous_est_counts": float(row.unambiguous_est_counts),
            "ambiguous_est_counts": float(row.ambiguous_est_counts),
        })
    missing = set(direct) - set(em["species"])
    rows.extend({
        "sample_id": sample_id,
        "species": species,
        "reported": True,
        "call_method": "exact_allele_profile",
        "exact_profile_loci": direct[species]["profile_loci"],
        "exact_profile_unambiguous_loci": direct[species]["unambiguous_loci"],
        # Same reasoning as the em.empty branch above: this species was called
        # from its allele profile but never reached EM, so its abundance was not
        # estimated. NaN, not 0.0.
        "pi": float("nan"),
        "allele_evidence_unambiguous_loci": direct[species]["unambiguous_loci"],
        "allele_evidence_ambiguous_loci": 0,
        "em_evidence_unambiguous_loci": 0,
        "em_evidence_ambiguous_loci": 0,
        "unambiguous_est_counts": float("nan"),
        "ambiguous_est_counts": float("nan"),
    } for species in sorted(missing))
    final = pd.DataFrame.from_records(rows, columns=FINAL_COLUMNS)
    if not final.empty:
        final = _attach_coverage(final, species_coverage(groups))
        final = _label_alternative_names(final, demoted)
        final = final.sort_values(["reported", "pi"], ascending=[False, False])
    return final, em


def _label_alternative_names(
    final: pd.DataFrame,
    demoted: dict[str, dict[str, int]],
) -> pd.DataFrame:
    """Say which unreported species are a second name for a reported one.

    A culture whose 16S/rMLST evidence is shared across near neighbours draws
    more than one species name - a mitis-group isolate will produce both
    ``Streptococcus mitis`` and ``Streptococcus oralis``. Only the
    best-supported name is reported, but the others are not noise: they are the
    same organism under a name the evidence cannot exclude. Without this column
    a reader cannot tell that from a genuinely separate organism that failed the
    bar, and two names in one sample read as polymicrobial.

    Two conditions, and the evidence one is what makes the column mean anything:

    * the species is in ``demoted`` - it held up as an allele profile in its own
      right and fell short only on exclusive loci. A same-genus match alone is
      not enough: it would also tag the large EM-fallback background.
    * its genus has a reported species to be an alternative *to*. A demoted
      species standing alone in its genus is a weak call, not a second name.
    """
    reported = final[final["reported"].astype(bool)]
    best_by_genus: dict[str, str] = {}
    for row in reported.sort_values("pi", ascending=False).itertuples(index=False):
        genus = _definitive_genus(row.species)
        best_by_genus.setdefault(genus, row.species)

    labels = []
    for row in final.itertuples(index=False):
        if bool(row.reported) or row.species not in demoted:
            labels.append("")
            continue
        sibling = best_by_genus.get(_definitive_genus(row.species), "")
        labels.append(sibling if sibling != row.species else "")
    final = final.copy()
    final["alternative_to"] = labels
    return final
