"""``bridge evaluate`` - re-score a finished run under two contamination-gate rules.

    bridge evaluate --summary identification_summary_all.csv \\
        --baselines rpm_baselines.csv --ntc-total-reads 22539 \\
        --truth truth_table.csv

Reads a run's summary table, re-applies the contamination gate to every
reported call two ways - reagent-only RPM-r for everything ("old") and the
mechanism-attributed rule ("new") - and, given a truth table, shows what each
rule gets right. It runs on the CSVs alone, so a threshold can be checked
before any re-run.

What it needs and where it gets it:

* the per-organism read counts per arm - in the summary already
  (``unambiguous_est_counts + ambiguous_est_counts``, ``16s_raw_reads``);
* the sample's depth - not in the summary, so it is derived from a read count
  and its RPM (``16s_raw_reads / rpm_16s * 1e6``);
* the blank's read counts - the baseline's ``reads_ntc`` column; for a file
  without it pass ``--ntc-total-reads`` and they are derived from ``rpm_ntc``;
* the mechanism per organism - the baseline's ``mechanism`` column; for a file
  without it, it is derived from the summary itself: an organism (or its genus) with a
  credible ``basis=rmlst`` row somewhere on the run is cross-contamination,
  the same test ``baseline.scan_cohort`` applies to the arm outputs.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, replace
from pathlib import Path

from bridge.reconcile.reconcile import CALL_BASES
from bridge.reconcile.rpm import (
    DEFAULT_ABS_R_MIN, DEFAULT_MIN_SAMPLE_READS, DEFAULT_RPM_MIN, DEFAULT_RPM_R_MIN,
    MECHANISM_CROSS, MECHANISM_REAGENT, MIN_NTC_RPM, Baseline, evaluate,
    load_baselines, lookup_baseline, rpm_detail,
)
from bridge.reconcile.summary_csv import _governing_band

TRUTH_NONE = {"", "n/a", "na", "none", "nan"}


def _num(value) -> float | None:
    try:
        text = str(value).strip()
        return float(text) if text else None
    except (TypeError, ValueError):
        return None


@dataclass
class Scored:
    sample: str
    organism: str
    genus: str
    basis: str
    depth: float
    truth: str                 # species | genus | not_in_truth | no_truth
    known: str
    old_band: str
    new_band: str
    new_rule: str
    new_arm: str
    rpm_r: float | None
    abs_r: float | None

    @property
    def old_verdict(self) -> str:
        return "background" if self.old_band == "background" else "report"

    @property
    def new_verdict(self) -> str:
        return "background" if self.new_band == "background" else "report"


def load_truth(path: Path, sample_col: str, organism_col: str) -> dict[str, str]:
    truth: dict[str, set[str]] = defaultdict(set)
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            sample = (row.get(sample_col) or "").strip()
            known = (row.get(organism_col) or "").strip()
            if sample and known.lower() not in TRUTH_NONE:
                truth[sample].add(known)
    return {sample: "; ".join(sorted(names)) for sample, names in truth.items()}


def truth_for(sample: str, truth: dict[str, str]) -> str | None:
    if sample in truth:
        return truth[sample]
    stripped = re.sub(r"[RX]$", "", sample)
    if stripped in truth:
        return truth[stripped]
    for key, value in truth.items():
        if re.sub(r"_nofilt$", "", key) == sample or re.sub(r"_nofilt$", "", key) == stripped:
            return value
    return None


def truth_match(organism: str, genus: str, known: str | None) -> str:
    if known is None:
        return "no_truth"
    lowered = known.lower()
    if organism and organism.lower() in lowered:
        return "species"
    if genus and genus.lower() in lowered:
        return "genus"
    return "not_in_truth"


def derive_depth(row: dict) -> float | None:
    raw, rpm = _num(row.get("16s_raw_reads")), _num(row.get("rpm_16s"))
    if raw and rpm:
        return raw / rpm * 1e6
    est = (_num(row.get("unambiguous_est_counts")) or 0) + (_num(row.get("ambiguous_est_counts")) or 0)
    rpm = _num(row.get("rpm_rmlst"))
    if est and rpm:
        return est / rpm * 1e6
    return None


def arm_reads(row: dict) -> dict[str, float | None]:
    unamb, amb = _num(row.get("unambiguous_est_counts")), _num(row.get("ambiguous_est_counts"))
    rmlst = None if unamb is None and amb is None else (unamb or 0.0) + (amb or 0.0)
    return {"rmlst": rmlst, "marker16s": _num(row.get("16s_raw_reads"))}


def credible_sources(rows: list[dict], rpm_r_min: float) -> tuple[set[str], set[str]]:
    """Species and genera with a credible rMLST-reported row on the run.

    Mirrors ``baseline.scan_cohort``: the row must be ``basis=rmlst`` and clear
    the reagent rule on at least one arm - a blank ``rpm_r`` beside a filled
    ``rpm`` means no control entry (floor pass); otherwise RPM-r >= threshold.
    """
    species: set[str] = set()
    genera: set[str] = set()
    for row in rows:
        if row.get("basis") != "rmlst":
            continue
        credible = False
        for arm in ("rmlst", "16s"):
            if _num(row.get(f"rpm_{arm}")) is None:
                continue
            ratio = _num(row.get(f"rpm_r_{arm}"))
            if ratio is None or ratio >= rpm_r_min:
                credible = True
        if not credible:
            continue
        organism = (row.get("organism") or "").strip().lower()
        genus = (row.get("genus") or "").strip().lower() or organism.split()[0]
        if organism:
            species.add(organism)
        if genus:
            genera.add(genus)
    return species, genera


def attribute_legacy(baselines: dict, species: set[str], genera: set[str],
                     ntc_total_reads: float | None) -> dict:
    """Fill mechanism (and reads_ntc, when depth is known) into a legacy baseline."""
    out = {}
    for key, baseline in baselines.items():
        organism = key[2]
        genus = organism.split()[0] if organism else ""
        cross = organism in species or organism in genera or genus in genera
        reads_ntc = baseline.reads_ntc
        if reads_ntc is None and ntc_total_reads:
            reads_ntc = baseline.rpm_ntc * ntc_total_reads / 1e6
        out[key] = replace(baseline, reads_ntc=reads_ntc,
                           mechanism=MECHANISM_CROSS if cross else MECHANISM_REAGENT)
    return out


def score_rows(rows: list[dict], baselines: dict, truth: dict[str, str], *,
               rpm_r_min: float, rpm_min: float, min_ntc_rpm: float,
               abs_r_min: float, min_sample_reads: float) -> tuple[list[Scored], int]:
    scored: list[Scored] = []
    skipped = 0
    for row in rows:
        if row.get("basis") not in CALL_BASES:
            continue
        depth = derive_depth(row)
        if not depth:
            skipped += 1
            continue
        organism = (row.get("organism") or "").strip()
        genus = (row.get("genus") or "").strip()
        run_id = (row.get("run_id") or "*").strip()
        old_detail, new_detail = {}, {}
        rpm_r = abs_r = None
        for arm, reads in arm_reads(row).items():
            if reads is None:
                continue
            baseline = lookup_baseline(baselines, run_id, arm, organism, genus=genus)
            reagent_only = (Baseline(rpm_ntc=baseline.rpm_ntc)
                            if isinstance(baseline, Baseline) else baseline)
            old = evaluate(organism, arm, reads, depth, reagent_only,
                           rpm_r_min=rpm_r_min, rpm_min=rpm_min, min_ntc_rpm=min_ntc_rpm)
            new = evaluate(organism, arm, reads, depth, baseline,
                           rpm_r_min=rpm_r_min, rpm_min=rpm_min, min_ntc_rpm=min_ntc_rpm,
                           abs_r_min=abs_r_min, min_sample_reads=min_sample_reads)
            old_detail[arm] = rpm_detail(old)
            new_detail[arm] = rpm_detail(new)
        old_arm, old_gov = _governing_band(old_detail)
        new_arm, new_gov = _governing_band(new_detail)
        if new_gov:
            rpm_r, abs_r = new_gov.get("rpm_r"), new_gov.get("abs_r")
        known = truth_for(row.get("sample", ""), truth) if truth else None
        scored.append(Scored(
            sample=row.get("sample", ""), organism=organism, genus=genus,
            basis=row.get("basis", ""), depth=depth,
            truth=truth_match(organism, genus, known) if truth else "no_truth",
            known=known or "", old_band=old_gov.get("band", ""),
            new_band=new_gov.get("band", ""), new_rule=str(new_gov.get("rule", "")),
            new_arm=new_arm, rpm_r=rpm_r, abs_r=abs_r))
    return scored, skipped


def confusion(scored: list[Scored], which: str) -> str:
    counts: Counter = Counter()
    for item in scored:
        verdict = item.old_verdict if which == "old" else item.new_verdict
        counts[(item.truth, verdict)] += 1
    lines = [f"{'truth':<14}{'report':>8}{'background':>12}"]
    for truth in ("species", "genus", "not_in_truth", "no_truth"):
        lines.append(f"{truth:<14}{counts[(truth, 'report')]:>8}{counts[(truth, 'background')]:>12}")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="bridge evaluate", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--baselines", type=Path, required=True)
    parser.add_argument("--truth", type=Path, default=None,
                        help="CSV with a sample column and a known-organism column")
    parser.add_argument("--truth-sample-col", default="sample")
    parser.add_argument("--truth-organism-col", default="known_organism")
    parser.add_argument("--ntc-total-reads", type=float, default=None,
                        help="the blank's n_processed, for a baseline file without reads_ntc")
    parser.add_argument("--rpm-r-min", type=float, default=DEFAULT_RPM_R_MIN)
    parser.add_argument("--rpm-min", type=float, default=DEFAULT_RPM_MIN)
    parser.add_argument("--min-ntc-rpm", type=float, default=MIN_NTC_RPM)
    parser.add_argument("--abs-r-min", type=float, default=DEFAULT_ABS_R_MIN)
    parser.add_argument("--min-sample-reads", type=float, default=DEFAULT_MIN_SAMPLE_READS)
    parser.add_argument("--sweep-abs-r", default="",
                        help="comma-separated abs_r_min values to compare, e.g. 5,10,20")
    parser.add_argument("--out", type=Path, default=None,
                        help="write every scored row (old and new verdicts) to this CSV")
    args = parser.parse_args(argv)

    with open(args.summary, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    baselines = load_baselines(args.baselines)
    legacy = all(b.reads_ntc is None for b in baselines.values())
    if legacy:
        if args.ntc_total_reads is None:
            print("ERROR: this baseline file has no reads_ntc column; pass --ntc-total-reads "
                  "(the blank's n_processed) so they can be derived from rpm_ntc.",
                  file=sys.stderr)
            return 1
        species, genera = credible_sources(rows, args.rpm_r_min)
        baselines = attribute_legacy(baselines, species, genera, args.ntc_total_reads)
        print(f"Legacy baseline: mechanism derived from the summary "
              f"({len(species)} credible species, {len(genera)} genera); "
              f"reads_ntc from rpm_ntc x {args.ntc_total_reads:,.0f}")
    cross = sum(1 for b in baselines.values() if b.mechanism == MECHANISM_CROSS)
    print(f"Baseline rows: {len(baselines)} ({cross} cross_contamination, "
          f"{len(baselines) - cross} reagent)")

    truth = (load_truth(args.truth, args.truth_sample_col, args.truth_organism_col)
             if args.truth else {})
    common = dict(rpm_r_min=args.rpm_r_min, rpm_min=args.rpm_min,
                  min_ntc_rpm=args.min_ntc_rpm, min_sample_reads=args.min_sample_reads)
    scored, skipped = score_rows(rows, baselines, truth, abs_r_min=args.abs_r_min, **common)
    print(f"Scored {len(scored)} reported calls ({skipped} skipped: depth not derivable); "
          f"truth for {sum(1 for s in scored if s.truth != 'no_truth')} of them")

    print("\n=== old rule: reagent RPM-r for everything ===")
    print(confusion(scored, "old"))
    print(f"\n=== new rule: mechanism-attributed (abs_r_min {args.abs_r_min:g}) ===")
    print(confusion(scored, "new"))

    changed = [s for s in scored if s.old_verdict != s.new_verdict]
    print(f"\n=== verdict changed on {len(changed)} row(s) ===")
    if changed:
        print(f"{'sample':<16}{'organism':<30}{'basis':<13}{'truth':<13}"
              f"{'old':<12}{'new':<12}{'rule':<7}{'rpm_r':>8}{'abs_r':>9}")
        for s in sorted(changed, key=lambda s: (s.truth, s.sample)):
            print(f"{s.sample:<16}{s.organism[:29]:<30}{s.basis:<13}{s.truth:<13}"
                  f"{s.old_band:<12}{s.new_band:<12}{s.new_rule:<7}"
                  f"{(f'{s.rpm_r:.3g}' if s.rpm_r is not None else ''):>8}"
                  f"{(f'{s.abs_r:.3g}' if s.abs_r is not None else ''):>9}")

    if args.sweep_abs_r:
        print("\n=== abs_r_min sweep: rows flipped to report, by truth ===")
        print(f"{'abs_r_min':>10}{'species':>9}{'genus':>7}{'not_in_truth':>14}{'no_truth':>10}")
        for value in [float(v) for v in args.sweep_abs_r.split(",") if v.strip()]:
            swept, _ = score_rows(rows, baselines, truth, abs_r_min=value, **common)
            flips = Counter(s.truth for s in swept
                            if s.old_verdict == "background" and s.new_verdict == "report")
            print(f"{value:>10g}{flips['species']:>9}{flips['genus']:>7}"
                  f"{flips['not_in_truth']:>14}{flips['no_truth']:>10}")

    if args.out:
        with open(args.out, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["sample", "organism", "genus", "basis", "depth", "truth", "known",
                             "old_band", "new_band", "new_rule", "new_arm", "rpm_r", "abs_r",
                             "changed"])
            for s in scored:
                writer.writerow([s.sample, s.organism, s.genus, s.basis, f"{s.depth:.0f}",
                                 s.truth, s.known, s.old_band, s.new_band, s.new_rule, s.new_arm,
                                 "" if s.rpm_r is None else f"{s.rpm_r:.4g}",
                                 "" if s.abs_r is None else f"{s.abs_r:.4g}",
                                 s.old_verdict != s.new_verdict])
        print(f"\nRows written to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
