"""From KMA's numbers to a call: detected or not, how firmly, and which allele.

Three tiers (the same shape as abritAMR's exact / ``*`` / ``^``):

* ``confident`` - identity >= ``min_identity``, coverage >= ``confident_coverage``,
  and enough material under it (``min_depth`` and ``min_fragments``);
* ``partial``   - identity and coverage over the detection floor
  (``min_identity``, ``min_coverage``) but not the rest;
* ``below``     - everything else. Kept in ``amr_hits.csv``, never reported.

Identity means identity WHERE THE GENE IS COVERED - KMA's Query_Identity -
not Template_Identity, which counts every uncovered base as a mismatch and so
charges a coverage gap twice: a gene matching at 99.77% over 89% of the
template reads as 88.75% Template_Identity and would fall under the 90% floor.
Real divergence still fails: chromosomal oqxAB homologues in Enterobacter and
Serratia sit at 85-88% where covered.

Allele level, for anything above ``below``:

* ``exact``   - the whole template, every base matching (identity AND
  coverage >= ``exact_identity``). The allele is named.
* ``nearest`` - KMA's best template, but not a perfect match: the allele is
  "closest to", with an estimate of how many bases differ.

Why the allele level is separate from the tier: capture reads can cover a
gene completely at 99.8% identity - a confident *family* call whose allele
is not established. For blaSHV that is the difference between SHV-1 (a
narrow-spectrum penicillinase) and SHV-12 (an ESBL), so the phenotype step
reads the level, not the tier.

Every threshold is provisional until the blood-culture AST validation.
"""

from __future__ import annotations

from .settings import Thresholds

TIERS = ("confident", "partial", "below")
ALLELE_LEVELS = ("exact", "nearest")


def _num(value, default: float = 0.0) -> float:
    try:
        text = str(value).strip()
        return float(text) if text else default
    except (TypeError, ValueError):
        return default


def covered_identity(hit: dict) -> float:
    """Identity over the covered part (Query_Identity; Template_Identity if absent)."""
    return _num(hit.get("query_identity"), _num(hit.get("template_identity")))


def tier(hit: dict, th: Thresholds) -> str:
    identity = covered_identity(hit)
    coverage = _num(hit.get("template_coverage"))
    if identity < th.min_identity or coverage < th.min_coverage:
        return "below"
    if (coverage >= th.confident_coverage
            and _num(hit.get("depth")) >= th.min_depth
            and _num(hit.get("fragments")) >= th.min_fragments):
        return "confident"
    return "partial"


def allele_level(hit: dict, th: Thresholds) -> str:
    identity = _num(hit.get("template_identity"))
    coverage = _num(hit.get("template_coverage"))
    query_identity = _num(hit.get("query_identity"), identity)
    # KMA reports percentages to two decimals; 99.995 rounds to 100.00.
    floor = th.exact_identity - 0.005
    if identity >= floor and coverage >= floor and query_identity >= floor:
        return "exact"
    return "nearest"


def est_differences(hit: dict) -> int:
    """Bases estimated to differ from the template over the part that was covered.

    Template_Identity counts matching bases over the whole template, so an
    uncovered stretch would read as mismatch; Query_Identity is over the
    aligned part only, which is the honest denominator for "how different".
    """
    length = _num(hit.get("template_length"))
    coverage = min(_num(hit.get("template_coverage")), 100.0) / 100.0
    query_identity = _num(hit.get("query_identity"), 100.0) / 100.0
    return int(round(length * coverage * (1.0 - query_identity)))


def _set_flag(hit: dict, flag: str, on: bool) -> None:
    flags = [f.strip() for f in str(hit.get("flags") or "").split(";")
             if f.strip() and not (f.strip() == flag or f.strip().startswith(flag + ":"))]
    if on:
        flags.append(flag)
    hit["flags"] = "; ".join(flags)


def flag_multiple_alleles(hits: list[dict]) -> list[dict]:
    """``family_multiple_alleles``: two or more DIFFERENT genes of one family
    above ``below`` in one sample.

    Counted on distinct gene names, so two ResFinder submissions of the same
    allele (blaPAO_1 and blaPAO_2) are not "multiple alleles". With KMA's -1t1
    each read goes to one template, so the flag means either two genuinely
    different alleles (two organisms, or two copies) or reads split between
    close alleles. Either way the allele call for that family needs a second
    look.
    """
    genes: dict[str, set[str]] = {}
    for hit in hits:
        if hit.get("tier") != "below":
            genes.setdefault(hit.get("family", ""), set()).add(hit.get("gene", ""))
    for hit in hits:
        _set_flag(hit, "family_multiple_alleles",
                  hit.get("tier") != "below" and len(genes.get(hit.get("family", ""), ())) > 1)
    return hits


def classify(hits: list[dict], th: Thresholds) -> list[dict]:
    """Add tier, allele_level, est_differences and the multi-allele flag, in place."""
    for hit in hits:
        hit["tier"] = tier(hit, th)
        hit["allele_level"] = allele_level(hit, th) if hit["tier"] != "below" else ""
        hit["est_differences"] = est_differences(hit)
    return flag_multiple_alleles(hits)


def _rank(hit: dict) -> tuple:
    return (hit.get("allele_level") == "exact",
            min(_num(hit.get("template_coverage")), 100.0),
            _num(hit.get("template_identity")),
            _num(hit.get("depth")))


def collapse_genes(hits: list[dict], th: Thresholds) -> list[dict]:
    """One call per gene name: the best template, carrying every template's reads.

    ResFinder often holds several submissions of one allele (``tet(M)_1_X92947``,
    ``tet(M)_2_...``), and KMA divides a gene's reads between them. Reported
    separately they double-count the gene and split its evidence. So, among
    templates above ``below`` that share a gene name, the best-matching one
    (exact first, then coverage, identity, depth) is the call, and its
    ``fragments`` and ``depth`` become the sum over all of them - the gene's
    reads, wherever KMA filed them. ``templates_merged:<n>`` in ``flags``
    records that it happened. Templates below the floor are returned
    untouched. Works on copies; the input is not modified.
    """
    classified = classify([dict(hit) for hit in hits], th)
    groups: dict[str, list[dict]] = {}
    out: list[dict] = []
    for hit in classified:
        if hit.get("tier") == "below":
            out.append(hit)
        else:
            groups.setdefault(hit.get("gene", ""), []).append(hit)
    for members in groups.values():
        best = max(members, key=_rank)
        if len(members) > 1:
            best["fragments"] = f"{sum(_num(m.get('fragments')) for m in members):g}"
            best["depth"] = f"{sum(_num(m.get('depth')) for m in members):.2f}"
            best["flags"] = "; ".join(filter(None, [best.get("flags", ""),
                                                    f"templates_merged:{len(members)}"]))
            best["tier"] = tier(best, th)
        out.append(best)
    return flag_multiple_alleles(out)
