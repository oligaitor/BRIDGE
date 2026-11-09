"""``bridge reconcile`` - join the two arms into one report.

This is the whole of the identification decision layer, and it is
deliberately small. Each arm has already done its own work and written a CSV;
this reads both, measures each organism against reagent background, hands the
two typed lists to the reconciler, and writes the result.

    rmlst/species_calls.csv  --adapters.rmlst-----> RmlstStrengthResult[] --+
                                                                            |--> reconcile_sample
    16s/marker16s_calls.csv  --adapters.marker16s-> BlastStrengthResult[] --+

Running it as a separate step from the two arms is intentional. The arms are
expensive and run on a compute node; reconciliation is seconds of arithmetic and
no external tools. Keeping them apart means a threshold can be changed and every
sample re-reconciled in a minute, without re-running any sequencing analysis.
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

from .adapters.marker16s import load_marker16s_calls, marker16s_results
from .adapters.rmlst import load_rmlst_calls, rmlst_results
from .config import Config, load_config
from .reconcile import (
    DEFAULT_MULTI_FALLBACK_REVIEW_MIN,
    flag_multi_fallback,
    reconcile_sample,
    sample_outcome,
)
from .rpm import (
    DEFAULT_ABS_R_MIN, DEFAULT_MIN_SAMPLE_READS, DEFAULT_RPM_MIN, DEFAULT_RPM_R_MIN,
    MIN_NTC_RPM, arms_with_background, evaluate, load_baselines, lookup_baseline,
    rpm_detail, ungated,
)
from .summary_csv import write_summary_csv

LOGGER = logging.getLogger("bridge.reconcile.reconcile_cli")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="bridge reconcile",
        description="Reconcile the rMLST and 16S marker arms into one call set.",
    )
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--rmlst-calls", type=Path, required=True,
                        help="the rMLST arm's species_calls.csv")
    parser.add_argument("--marker16s-calls", type=Path, default=None,
                        help="marker16s_calls.csv; omit to run the rMLST arm alone")
    parser.add_argument("-c", "--config", type=Path, required=True)
    parser.add_argument("-o", "--output", type=Path, required=True)
    parser.add_argument("--run-id", default="*",
                        help="sequencing run, used to pick the negative-control "
                             "baseline; reagent background is run-specific")
    parser.add_argument("--rpm-baselines", type=Path, default=None,
                        help="negative-control baseline CSV (see rpm.py). "
                             "Without it the RPM-r gate is skipped and every "
                             "call is flagged as ungated.")
    parser.add_argument("--total-reads", type=int, default=None,
                        help="sample's total raw read count, the RPM denominator")
    parser.add_argument("-v", "--verbose", action="store_true")
    return parser


def _rmlst_reads_by_organism(rmlst_rows) -> dict[str, float]:
    """Read counts keyed by species, with genus totals as a fallback.

    The genus entries matter. ``_gate_results`` falls back to the genus when it
    cannot find an organism by name - e.g. a call named only at genus level,
    which a species-keyed dict can never contain. Without them such a row would
    reach the report flagged ``rpm_gate_no_read_count``, never checked against
    background at all.

    A genus total is at least as large as any species within it, so as a
    denominator-free numerator it makes the gate stricter, never looser. The
    species key always wins where it exists (``setdefault``).
    """
    reads = {row.species.lower(): row.total_est_counts
             for row in rmlst_rows if row.species}

    genus_totals: dict[str, float] = {}
    for row in rmlst_rows:
        if not row.species or row.total_est_counts is None:
            continue
        genus = row.species.split()[0].lower()
        genus_totals[genus] = genus_totals.get(genus, 0.0) + row.total_est_counts
    for genus, total in genus_totals.items():
        reads.setdefault(genus, total)
    return reads


def _gate_results(results, arm, reads_by_organism, baselines, run_id,
                  total_reads, cfg):
    """Measure every arm result against its reagent background.

    **This annotates; it does not remove.** The RPM band *is* the verdict a
    reader acts on, so an organism that failed is more useful visible and
    marked than absent - "not in the table" would otherwise mean both
    "nothing there" and "found, and judged to be kitome".

    The knock-on is deliberate: because a failed organism stays in the call set,
    the 16S fallback opens only when the rMLST arm genuinely resolved
    nothing, not when the gate happened to clear the board.

    With no baseline (``baselines is None``) nothing is tested, but RPM is
    still computed for every organism and the sample is flagged
    ``rpm_gate_not_applied``; without a read total nothing can be computed at
    all and only the flag is returned.
    """
    if total_reads is None:
        return results, {}, ["rpm_gate_not_applied"]

    gated, detail = [], {}
    for result in results:
        # The two result types name their taxon field differently:
        # RmlstStrengthResult.taxon, BlastStrengthResult.top_taxon.
        taxon = getattr(result, "taxon", None) or getattr(result, "top_taxon", None)
        name = (taxon.species_full or taxon.genus) if taxon else None
        if not name:
            gated.append(result)
            continue
        reads = reads_by_organism.get(name.lower())
        if reads is None and taxon.genus:
            reads = reads_by_organism.get(taxon.genus.lower())
        if reads is None:
            # No read count for this organism, so there is nothing to compare
            # against its background. Pass the evidence through rather than
            # destroying it, but say so - an ungated call must not be mistaken
            # for one that cleared the gate.
            #
            # This is reached when the rMLST arm called a species from its
            # allele profile but never estimated its abundance (blank
            # est_counts); see RmlstRow.total_est_counts.
            gated.append(type(result)(**{
                **result.__dict__,
                "flags": list(result.flags) + ["rpm_gate_no_read_count"],
            }))
            LOGGER.info("rpm_gate arm=%s organism=%s skipped: no read count", arm, name)
            continue

        genus = taxon.genus if taxon else None
        others: list[str] = []
        if baselines is None:
            outcome = ungated(name, arm, reads, total_reads)
        else:
            outcome = evaluate(
                name, arm, reads, total_reads,
                lookup_baseline(baselines, run_id, arm, name, genus=genus),
                rpm_r_min=float(cfg.get("rpm.rpm_r_min", DEFAULT_RPM_R_MIN)),
                rpm_min=float(cfg.get("rpm.rpm_min", DEFAULT_RPM_MIN)),
                min_ntc_rpm=float(cfg.get("rpm.min_ntc_rpm", MIN_NTC_RPM)),
                abs_r_min=float(cfg.get("rpm.abs_r_min", DEFAULT_ABS_R_MIN)),
                min_sample_reads=float(cfg.get("rpm.min_sample_reads",
                                               DEFAULT_MIN_SAMPLE_READS)),
            )
            # "No background" is an arm-local statement. Where the OTHER arm's
            # control saw this organism, say so rather than letting the 0.05
            # RPM floor - which nothing can fail - read as a clean result.
            if outcome.rule == "floor":
                others = [a for a in arms_with_background(baselines, run_id, name, genus)
                          if a != arm]
        detail[name] = rpm_detail(outcome, others)

        flags = list(result.flags)
        if outcome.pooled:
            # Gated against another run's blank (bridge pool-baseline), not
            # this run's: weaker evidence, so the row says which it was.
            flags.append("rpm_pooled_baseline")
        if others:
            flags.append(f"rpm_control_background_in_{'_'.join(others)}")
            LOGGER.info("rpm_gate arm=%s organism=%s ungated here but present "
                        "in the control via %s", arm, name, ", ".join(others))
        if not outcome.passed:
            flags.append("rpm_below_background")
            LOGGER.info("rpm_gate arm=%s organism=%s band=%s %s",
                        arm, name, outcome.band, outcome.reason)

        gated.append(type(result)(**{
            **result.__dict__,
            "flags": flags,
            "reason": (result.reason if outcome.passed
                       else f"{result.reason}; {outcome.reason}"),
        }))
    return gated, detail, ([] if baselines is not None else ["rpm_gate_not_applied"])


def _rmlst_evidence(row) -> dict:
    """This organism's own numbers, verbatim from ``species_calls.csv``.

    Carried per organism rather than per sample: in a polymicrobial sample the
    read count is exactly the figure that has to differ between organisms for
    the result to mean anything.
    """
    return {
        "call_method": row.call_method,
        "pi": row.pi,
        "exact_profile_loci": row.exact_profile_loci,
        "exact_profile_unambiguous_loci": row.exact_profile_unambiguous_loci,
        "allele_evidence_unambiguous_loci": row.unambiguous_loci,
        "allele_evidence_ambiguous_loci": row.ambiguous_loci,
        "em_evidence_unambiguous_loci": row.em_evidence_unambiguous_loci,
        "em_evidence_ambiguous_loci": row.em_evidence_ambiguous_loci,
        "unambiguous_est_counts": row.unambiguous_est_counts,
        "ambiguous_est_counts": row.ambiguous_est_counts,
        "alternative_to": row.alternative_to,
        # Coverage of the loci that name this species alone - the same loci
        # allele_evidence_unambiguous_loci counts. The 16S arm's own coverage
        # figures (one gene, one reference) are the 16s_* columns below.
        "mean_callable_breadth": row.mean_callable_breadth,
        "min_callable_breadth": row.min_callable_breadth,
        "median_projected_depth": row.median_projected_depth,
        "min_projected_depth": row.min_projected_depth,
        "total_callable_bases": row.total_callable_bases,
        "total_allele_length": row.total_allele_length,
        "supporting_fragments": row.supporting_fragments,
    }


def _s16_evidence(row) -> dict:
    """This organism's own numbers, verbatim from ``marker16s_calls.csv``.

    Keyed ``16s_*`` to match the column names a reader sees. The dataclass
    fields behind them stay ``s16_*`` only because a Python identifier cannot
    begin with a digit; nothing outside the source needs to know that.
    """
    return {
        # Share of the whole COPY-NUMBER-CORRECTED profile, not of raw reads,
        # and the denominator includes taxa below the reporting floor - so a
        # sample's reported fractions do not sum to 1.
        "16s_abundance_frac": row.abundance_frac,
        "16s_reads_assigned": row.reads_assigned,
        "16s_raw_reads": row.raw_reads,
        # From an alignment of the reads to the called genus's references.
        # None on a calls file without these columns.
        "16s_pident": row.consensus_pident,
        "16s_qcov": row.consensus_qcov,
        "16s_ref_breadth": row.ref_breadth,
        "16s_mean_depth": row.mean_depth,
        "16s_read_identity": row.mean_read_identity,
        "16s_aligned_reads": row.n_aligned_reads,
    }


def _evidence_for(call, rmlst_rows, s16_rows) -> dict:
    """Both arms' per-organism numbers for one call, by row index.

    By index, not by name: the call's name need not match its row (a genus-level
    name, a compound genus), and a failed name lookup would quietly lose every
    number behind the row.
    """
    evidence: dict = {}
    index = call.rmlst_row_index
    if index is not None and 0 <= index < len(rmlst_rows):
        evidence |= _rmlst_evidence(rmlst_rows[index])
    index = call.s16_row_index
    if index is not None and 0 <= index < len(s16_rows):
        evidence |= _s16_evidence(s16_rows[index])
    return evidence


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING,
                        format="%(levelname)s %(name)s %(message)s")
    cfg: Config = load_config(args.config)
    args.output.mkdir(parents=True, exist_ok=True)

    rmlst_rows = load_rmlst_calls(args.rmlst_calls)
    rmlst = rmlst_results(rmlst_rows, cfg)
    rmlst_reads = _rmlst_reads_by_organism(rmlst_rows)

    s16, s16_rows, s16_reads = [], [], {}
    if args.marker16s_calls:
        s16_rows = load_marker16s_calls(args.marker16s_calls)
        s16 = marker16s_results(s16_rows, cfg)
        for row in s16_rows:
            # RAW reads, not the copy-number-corrected figure, because these
            # feed the RPM-r gate - see rpm.py's "Copy number must not reach the
            # ratio". The corrected value stays in marker16s_calls.csv, where it
            # answers a different question.
            if row.taxon:
                s16_reads[row.taxon.lower()] = row.raw_reads
            if row.genus:
                s16_reads.setdefault(row.genus.lower(), row.raw_reads)

    # The baseline file: named on the command line, else `rpm.baselines` in the
    # config. A path given explicitly and missing is an error; one that came
    # from the config is a warning, because the config describes a usual
    # layout and this run may simply have no blank.
    baseline_path = args.rpm_baselines
    if baseline_path is None and cfg.get("rpm.baselines"):
        baseline_path = Path(str(cfg.get("rpm.baselines"))).expanduser()
    baselines = None
    if baseline_path is not None:
        if baseline_path.is_file():
            baselines = load_baselines(baseline_path)
        elif args.rpm_baselines is not None:
            build_parser().error(f"--rpm-baselines: no such file: {baseline_path}")
        else:
            LOGGER.warning("rpm.baselines %s not found; the gate is not applied",
                           baseline_path)
    rpm_flags: list[str] = []
    if bool(cfg.get("rpm.enabled", True)):
        rmlst, rmlst_rpm, flags = _gate_results(
            rmlst, "rmlst", rmlst_reads, baselines, args.run_id,
            args.total_reads, cfg)
        rpm_flags += flags
        s16, s16_rpm, flags = _gate_results(
            s16, "marker16s", s16_reads, baselines, args.run_id,
            args.total_reads, cfg)
        rpm_flags += flags
    else:
        rmlst_rpm = s16_rpm = {}

    calls = reconcile_sample(rmlst, s16, cfg)
    flag_multi_fallback(calls, int(cfg.get("reconcile.multi_fallback_review_min",
                                           DEFAULT_MULTI_FALLBACK_REVIEW_MIN)))
    outcome = sample_outcome(calls)

    # Whole-sample totals. They are context in the JSON record, not columns in
    # the CSV, where a sample total repeated against every organism in a
    # polymicrobial result would mislead.
    #
    # Summed over ROWS, not over the reads dicts. Those dicts deliberately file
    # one 16S row under both its taxon and its genus, so summing their values
    # counts the first row of every genus twice. None means "not estimated", not
    # zero, so it contributes nothing rather than dragging the total down.
    rmlst_total = sum(row.total_est_counts for row in rmlst_rows
                      if row.total_est_counts is not None)
    s16_total = sum(row.raw_reads for row in s16_rows
                    if row.raw_reads is not None)

    record = {
        "sample": args.sample_id,
        "run_id": args.run_id,
        "outcome": outcome,
        "rmlst_reads": rmlst_total,
        "reads_16s": s16_total,
        "total_raw_reads": args.total_reads,
        "flags": sorted(set(rpm_flags)),
        "calls": [
            {
                "organism": call.organism,
                "genus": call.genus,
                "reported_rank": getattr(call.reported_rank, "value", str(call.reported_rank)),
                "basis": call.basis,
                # Per-arm evidence. Carried so the record shows WHY a call was
                # made, not just what it was - a reviewer needs the locus count,
                # the coverage and the read count behind it.
                "rmlst_taxon": call.rmlst_taxon,
                "16s_taxon": call.s16_top_hit,
                "16s_margin": call.s16_margin,
                **_evidence_for(call, rmlst_rows, s16_rows),
                "flags": list(call.flags),
                "notes": call.notes,
            }
            for call in calls
        ],
        # Both arms' verdicts, in full. The summary CSV picks a governing one to
        # show in its headline columns, but the record keeps both, because the
        # two arms can disagree.
        "rpm": {
            arm: dict(detail)
            for arm, detail in (("rmlst", rmlst_rpm), ("marker16s", s16_rpm))
        },
    }
    path = args.output / "reconciled_calls.json"
    # newline="\n" so the record is byte-identical whichever platform wrote it.
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(record, indent=2) + "\n")
    # The same result as a flat table, for reading and for concatenating across
    # samples. The JSON remains the complete record.
    csv_path = write_summary_csv(record, args.output / "identification_summary.csv")
    print(f"{args.sample_id}: {outcome}, {len(calls)} row(s) -> {path}")
    print(f"{' ' * len(args.sample_id)}  summary CSV -> {csv_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
