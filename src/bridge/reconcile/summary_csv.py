"""The reported result as a CSV: one row per organism, both arms side by side.

``reconciled_calls.json`` is the complete record and stays the source of truth,
but nobody reads JSON to answer "what grew in these samples". This writes the
same information as a flat table that opens in Excel and concatenates across
samples with ``cat``.

One row per **reported organism**, not per sample: a polymicrobial sample gets a
row each (each carrying the ``polymicrobial`` flag), and a sample with nothing
reportable still gets one row with ``status=no_call``, so no sample silently
disappears from the table. The sample-level outcome itself (``reported`` /
``polymicrobial`` / ``no_call``) lives in ``reconciled_calls.json`` as
``outcome``; as a column it would say nothing ``status`` and ``flags`` do not
already say.

Negative controls are **not** in this table. They are not patient results, no
organism is called from them, and mixing them in makes a results table that has
to be filtered before it can be read. What the blank contained is recorded in
``rpm_baselines.csv``, per run and per arm, which is also the file the gate
actually divides by.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any, Iterable

from .reconcile import CALL_BASES, MULTI_FALLBACK_FLAG, below_background

# Ordered for reading left to right: what was found, on what evidence, how much
# of it there was, and whether it cleared reagent background.
#
# There is no `confidence` column: a letter grade derived from these numbers
# would obscure them. `basis` says which arm produced the row and whether it is
# a reported call; the evidence columns say how good it is.
SUMMARY_COLUMNS = [
    # --- identity -----------------------------------------------------------
    # `basis` is one of: rmlst | 16s_fallback (reported calls), or
    # rmlst_alternative | 16s_only (informational rows). See reconcile.py.
    # `status` is the one-word reading of it, for filtering. See _status.
    "sample", "run_id", "organism", "genus", "reported_rank", "basis",
    "status",
    # --- rMLST arm, PER ORGANISM, verbatim from species_calls.csv -----------
    # A blank here means "not estimated", never zero. The rMLST arm leaves pi
    # and the est_counts empty for a species called from its allele profile
    # whose consensus groups never reached EM.
    "rmlst_taxon", "call_method", "pi",
    "exact_profile_unambiguous_loci", "exact_profile_loci",
    "allele_evidence_unambiguous_loci", "allele_evidence_ambiguous_loci",
    "em_evidence_unambiguous_loci", "em_evidence_ambiguous_loci",
    "unambiguous_est_counts", "ambiguous_est_counts", "alternative_to",
    # --- rMLST coverage, over those same unambiguous loci -------------------
    # rMLST only; the 16S arm's equivalents are the 16s_* columns further down,
    # measured on one gene against one reference and not comparable to these.
    "mean_callable_breadth", "min_callable_breadth",
    "median_projected_depth", "min_projected_depth",
    "total_callable_bases", "total_allele_length", "supporting_fragments",
    # --- 16S arm, PER ORGANISM, verbatim from marker16s_calls.csv ----------
    # Populated wherever the 16S arm named a matching genus, including on rows
    # rMLST governed. The arm runs on every sample; it just does not decide
    # what gets reported.
    #
    # `16s_abundance_frac` is a share of the COPY-NUMBER-CORRECTED profile, and
    # its denominator includes taxa below the reporting floor - so the fractions
    # in one sample do not sum to 1. `16s_raw_reads` is what kallisto actually
    # assigned, before correction; `16s_reads_assigned` is after it.
    "16s_taxon", "16s_margin",
    "16s_abundance_frac", "16s_reads_assigned", "16s_raw_reads",
    # Identity and coverage from an alignment of the
    # reads to the called genus's MIMt references. `16s_pident` / `16s_qcov`
    # are the consensus's BLAST percent identity and query coverage against
    # the top reference - the two numbers a Sanger 16S report gives;
    # `16s_ref_breadth` / `16s_mean_depth` describe the top reference's
    # coverage; `16s_read_identity` is the mean per-read identity and
    # `16s_aligned_reads` how many reads it rests on. Blank on a calls file
    # without these columns.
    "16s_pident", "16s_qcov", "16s_ref_breadth", "16s_mean_depth",
    "16s_read_identity", "16s_aligned_reads",
    # --- contamination gate -------------------------------------------------
    # The verdict, the arm it came from, then BOTH arms' figures.
    #
    # There is deliberately NO single ``rpm`` column. An organism has one RPM
    # per arm and they are not comparable, let alone addable: the rMLST figure
    # counts fragments over ~53 protein-coding loci, the 16S figure counts
    # copy-number-corrected reads over one gene, and the two can differ
    # several-fold for one organism. A column called ``rpm`` carrying one of
    # those reads as "the organism's RPM" and hides which arm
    # measured it. ``rpm_band`` is the decision; ``rpm_band_arm`` says who made
    # it; the columns after it are the evidence.
    #
    # ``abs_r_*`` is sample reads over the
    # blank's reads, per arm; ``bg_mechanism`` is the governing arm's
    # attribution of the organism's presence in the blank (reagent, gated on
    # RPM-r; cross_contamination, gated on abs-r), which tells the reader
    # which of the two ratios the band came from. Band ``ungated`` means no
    # baseline at all: RPM computed, nothing tested.
    "rpm_band", "rpm_band_arm",
    "rpm_rmlst", "rpm_r_rmlst", "rpm_16s", "rpm_r_16s",
    "abs_r_rmlst", "abs_r_16s", "bg_mechanism",
    # --- provenance ---------------------------------------------------------
    "flags", "notes",
]

# Worst to best, among verdicts that actually measured something.
_BAND_ORDER = {
    "background": 0,
    "low": 1,
    "no_control_background": 2,
    "ungated": 2,
    "moderate": 3,
    "high": 4,
}


def _rank(detail: dict) -> tuple[int, int]:
    """Sort key deciding which arm's verdict the headline columns carry.

    Three tiers, in order:

    0. **A failed gate governs unconditionally.** One arm calling an organism
       reagent background is the finding, whatever the other arm says.
    1. **Then any arm that had a denominator.** ``no_control_background`` is not
       a clean result, it is an absent measurement: the organism had no control
       level to divide by, so it was held to the 0.05 RPM absolute floor, which
       at clinical depth is under a fifth of one read. An arm that did measure a
       ratio is the only one saying anything about contamination, so it governs
       - even when its verdict is the more comfortable of the two.
    2. **Unmeasured passes last**, shown only when nothing else measured.

    Within a tier, the worse band governs.
    """
    band = detail.get("band", "")
    if band == "background":
        return (0, 0)
    rule = detail.get("rule")
    if rule:
        measured = rule in ("rpm_r", "abs_r")
    else:  # a record written before `rule` existed
        measured = detail.get("rpm_r") is not None
    return (1 if measured else 2, _BAND_ORDER.get(band, 99))


def _fmt(value: Any, digits: int = 4) -> str:
    """Render a value for a CSV cell.

    Significant-figure formatting, so a trace abundance keeps its precision
    rather than rounding to 0.0000. Very small values come out in scientific
    notation (9.58e-05), which spreadsheets read back as numbers correctly.
    """
    if value is None or value == "":
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int,)):
        return str(value)
    if isinstance(value, float):
        return f"{value:.{digits}g}"
    return str(value)


def _lookup(detail: dict, names: Iterable[str]) -> dict:
    """First of ``names`` present in one arm's RPM detail, case-insensitively."""
    wanted = [name for name in names if name]
    for name in wanted:
        if name in detail:
            return detail[name]
    lowered = {key.lower(): value for key, value in detail.items()}
    for name in wanted:
        if name.lower() in lowered:
            return lowered[name.lower()]
    return {}


def _rpm_by_arm(record: dict, call: dict) -> dict[str, dict]:
    """Each arm's RPM verdict for one call, keyed by arm.

    The lookup has to try the **arm's own taxon name** first, because that is
    the key the gate wrote under: it stores results by
    ``taxon.species_full or taxon.genus``, which is the arm's finest resolution,
    not the reconciler's final organism name. A genus-level call named
    ``Clostridium`` therefore misses an entry filed under
    ``Clostridium nigeriense``, and a name-only lookup would show a blank RPM
    for a figure that had in fact been computed.
    """
    organism = call.get("organism") or ""
    genus = call.get("genus") or ""
    per_arm: dict[str, dict] = {}
    for arm, arm_taxon in (("rmlst", call.get("rmlst_taxon")),
                           ("marker16s", call.get("16s_taxon"))):
        detail = (record.get("rpm") or {}).get(arm) or {}
        found = _lookup(detail, (arm_taxon, organism, genus))
        if not found and genus:
            # Last resort: any entry in this genus. A 16S arm that resolved a
            # species the reconciler declined to name still measured the same
            # organism's background.
            prefix = f"{genus.lower()} "
            for key, value in detail.items():
                if key.lower() == genus.lower() or key.lower().startswith(prefix):
                    found = value
                    break
        if found:
            per_arm[arm] = found
    return per_arm


def _governing_band(per_arm: dict[str, dict]) -> tuple[str, dict]:
    """The arm whose verdict the headline columns should carry. See :func:`_rank`.

    Not simply the rMLST arm: where rMLST has no baseline for an organism it
    passes on the absolute floor as ``no_control_background``, while the 16S
    arm may score the same organism barely clear of reagent background.
    Showing the rMLST arm would turn a marginal gate result into a silent
    pass.

    Ties go to rMLST, which is the discriminating arm.
    """
    if not per_arm:
        return "", {}
    ordered = sorted(
        per_arm.items(),
        key=lambda item: (*_rank(item[1]), 0 if item[0] == "rmlst" else 1),
    )
    return ordered[0]


# Per-organism evidence carried verbatim from the two arm CSVs, and the number
# of significant figures each is worth. Six for the quantities a trace organism
# lives in - a pi of 9.58e-05 must not round to 0.0001 - and integers pass
# through _fmt untouched.
_EVIDENCE_DIGITS = {
    "pi": 6, "unambiguous_est_counts": 6, "ambiguous_est_counts": 6,
    "mean_callable_breadth": 4, "min_callable_breadth": 4,
    "median_projected_depth": 4, "min_projected_depth": 4,
    "total_callable_bases": 6, "total_allele_length": 6,
    "supporting_fragments": 6,
    "16s_abundance_frac": 6, "16s_reads_assigned": 6, "16s_raw_reads": 6,
    "16s_pident": 5, "16s_qcov": 4, "16s_ref_breadth": 4, "16s_mean_depth": 4,
    "16s_read_identity": 4, "16s_aligned_reads": 6,
}
_EVIDENCE_COLUMNS = (
    "call_method", "exact_profile_loci", "exact_profile_unambiguous_loci",
    "allele_evidence_unambiguous_loci", "allele_evidence_ambiguous_loci",
    "em_evidence_unambiguous_loci", "em_evidence_ambiguous_loci",
    "alternative_to", *_EVIDENCE_DIGITS,
)


def _status(call: dict, flags: list[str]) -> str:
    """``report``, ``review`` or ``review_override``, for filtering a run.

    Derived, not decided here - it is a one-word reading of ``basis``, the
    contamination verdict and the reconciler's flags, so it cannot disagree
    with the columns beside it. A row is ``review`` when any of:

    * it is not a reported call at all (``16s_only``, ``rmlst_alternative``) -
      real observations, but not organisms this pipeline is asserting;
    * the RPM-r gate put it at or below the run's own reagent background.

    Otherwise it is ``review_override`` when it would have been ``report`` but
    is one of several such 16S-fallback calls in a sample with no other
    reported call (the ``multiple_16s_fallback`` flag, set in
    reconcile.flag_multi_fallback): an ambiguous 16S-only call that may be
    background the control did not model. A separate value rather than
    ``review`` so those cases can be pulled out and looked at on their own.
    The reconciler never flags a below-background call, so the two do not
    overlap; the order here is a guard, not a tie-break.

    Deliberately NOT a revived confidence tier. It reads no locus count, no
    coverage and no abundance. Where there is no negative control - as on the
    synthetic panel - the RPM half cannot fire. Anything other than ``report``
    keeps an organism out of AMR attribution (bridge.amr.report).
    """
    if call.get("basis") not in CALL_BASES:
        return "review"
    if below_background(flags):
        return "review"
    if MULTI_FALLBACK_FLAG in flags:
        return "review_override"
    return "report"


def summary_rows(record: dict) -> list[dict[str, str]]:
    """Flatten one ``reconciled_calls.json`` record into CSV rows."""
    # Only two things are genuinely sample-level. Sample read totals and the
    # sample outcome stay in the JSON record: repeated on every organism row
    # they would mislead, and the outcome is `status` + the `polymicrobial`
    # flag by another name.
    shared = {
        "sample": record.get("sample", ""),
        "run_id": record.get("run_id", ""),
    }

    calls = record.get("calls") or []
    if not calls:
        # A sample with no reportable organism still gets a row. Dropping it
        # would make "not in the table" ambiguous between "nothing found" and
        # "never ran", which is exactly the distinction a negative result needs.
        return [{
            **{column: "" for column in SUMMARY_COLUMNS}, **shared,
            "organism": "", "reported_rank": "none", "basis": "",
            "status": "no_call",
            "flags": "; ".join(record.get("flags") or []),
        }]

    rows = []
    for call in calls:
        per_arm = _rpm_by_arm(record, call)
        arm, governing = _governing_band(per_arm)
        rmlst_rpm = per_arm.get("rmlst", {})
        s16_rpm = per_arm.get("marker16s", {})

        # Sample-level flags (today only `rpm_gate_not_applied`) reach every
        # organism row, not just the no-calls row. Order kept, duplicates
        # dropped.
        flags = list(dict.fromkeys([*(record.get("flags") or []),
                                    *(call.get("flags") or [])]))
        # Name the other arm's control evidence in the row itself, not only in
        # the JSON: "no_control_background" is otherwise read as "clean". The
        # gate already flags the arm result it applies to, and the reconciler
        # prefixes those with the arm name ("rmlst:..."), so only add one here
        # where that has not already happened - a flag listed twice under two
        # spellings reads as two findings.
        for detail in per_arm.values():
            for other in detail.get("other_arm_background") or []:
                flag = f"rpm_control_background_in_{other}"
                if not any(existing.endswith(flag) for existing in flags):
                    flags.append(flag)

        rows.append({
            **{column: "" for column in SUMMARY_COLUMNS}, **shared,
            "organism": call.get("organism", ""),
            "genus": call.get("genus") or "",
            "reported_rank": call.get("reported_rank", ""),
            "basis": call.get("basis", ""),
            "status": _status(call, flags),
            "rmlst_taxon": call.get("rmlst_taxon") or "",
            # Each arm's own numbers for THIS organism, straight through. A key
            # the record does not carry stays blank rather than becoming zero.
            **{
                column: _fmt(call[column], _EVIDENCE_DIGITS.get(column, 4))
                for column in _EVIDENCE_COLUMNS if column in call
            },
            "16s_taxon": call.get("16s_taxon") or "",
            "16s_margin": _fmt(call.get("16s_margin"), 3),
            "rpm_band": governing.get("band", ""),
            "rpm_band_arm": arm,
            "rpm_rmlst": _fmt(rmlst_rpm.get("rpm"), 4),
            "rpm_r_rmlst": _fmt(rmlst_rpm.get("rpm_r"), 4),
            "rpm_16s": _fmt(s16_rpm.get("rpm"), 4),
            "rpm_r_16s": _fmt(s16_rpm.get("rpm_r"), 4),
            "abs_r_rmlst": _fmt(rmlst_rpm.get("abs_r"), 4),
            "abs_r_16s": _fmt(s16_rpm.get("abs_r"), 4),
            "bg_mechanism": governing.get("mechanism") or "",
            "flags": "; ".join(flags),
            "notes": call.get("notes", "") or "",
        })
    return rows


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))



def write_summary_csv(record: dict, path: Path) -> Path:
    """Write one sample's reported result to ``identification_summary.csv``."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_COLUMNS)
        writer.writeheader()
        writer.writerows(summary_rows(record))
    return path


def combine_summaries(results_dir: Path, out_path: Path) -> tuple[Path, int]:
    """Concatenate every sample's summary into one run-level table.

    Negative controls are absent by construction: the reconciler skips them, so
    they have no per-sample summary to pick up. See the module docstring.

    Returns the path written and the number of samples included.
    """
    rows: list[dict[str, str]] = []
    samples = 0
    for csv_path in sorted(results_dir.glob("*/bacterial_ID/identification_summary.csv")):
        found = _read_csv(csv_path)
        if not found:
            continue
        rows.extend(found)
        samples += 1

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_COLUMNS,
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return out_path, samples
