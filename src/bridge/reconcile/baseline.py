"""Turn a negative control's output into per-organism baselines for the gate.

    bridge build-baseline -i <results_dir> -n NEG_1 -r RUN1 -o <results_dir>/rpm_baselines.csv

(``bridge finalise`` calls :func:`build_run_baselines` for you.)

For every organism the negative control saw, this records how high it got there
- in reads per million AND in reads - and which of two mechanisms put it there.
Reconciliation then reports an organism from a real sample only if it is well
clear of that level: ten times its RPM for reagent background, ten times its
read count for cross-contamination (see ``rpm.py``, "Two mechanisms").

Why per run
-----------
Reagent lots change between sequencing runs and the contamination they carry
changes with them, so a control is only properly comparable to the samples
extracted alongside it. ``--run-id`` tags the baseline; give reconciliation the
same value. ``*`` makes it a pooled fallback usable by any run.

The mechanism is a run-level fact too: it comes from what the *other* samples
on the run carry, so the baseline must be rebuilt for each run, after every
arm has finished.

The denominator
---------------
RPM needs a total read count. This reads ``n_processed`` out of kallisto's
``run_info.json`` - the number of read pairs the quantification actually saw.
Taking it from the run itself means the sample and the control are measured the
same way, which matters more than which convention is used.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

from .rpm import (
    DEFAULT_RPM_MIN, DEFAULT_RPM_R_MIN, MAX_NTC_SOURCE_FRAC, MECHANISM_CROSS,
    MIN_CONTROL_READS, MIN_NTC_RPM, Baseline, CohortEvidence, baselines_from_rows,
    build_baselines, evaluate, lookup_baseline, write_baselines,
)

# Negative controls are recognised by a CASE-INSENSITIVE name prefix. "NEG"
# rather than "NEG_", so a blank named e.g. NEGMATRIX_... is not reconciled as
# a patient sample, which would leave the run with no baseline at all.
DEFAULT_CONTROL_PREFIXES = ("NEG", "NTC", "Blank")
DEFAULT_EXCLUDE_PREFIXES = DEFAULT_CONTROL_PREFIXES
ARMS = ("rmlst", "marker16s")


def is_control_name(name: str, prefixes: tuple[str, ...] = DEFAULT_CONTROL_PREFIXES) -> bool:
    """Whether a sample name starts with any control prefix, ignoring case."""
    folded = name.casefold()
    return any(folded.startswith(prefix.casefold()) for prefix in prefixes if prefix)


def total_reads(sample_dir: Path) -> int | None:
    """Read pairs processed, from either arm's kallisto run_info.json.

    The last two entries are an older flat layout, still accepted so an
    existing results tree yields a baseline without being restructured.
    """
    sample_dir = Path(sample_dir)
    for candidate in (sample_dir / "16s" / "quant" / "run_info.json",
                      sample_dir / "rmlst" / "quant" / "run_info.json",
                      sample_dir / "marker16s" / "quant" / "run_info.json",
                      sample_dir / "quant" / "run_info.json"):
        if candidate.is_file():
            try:
                return int(json.loads(candidate.read_text())["n_processed"])
            except (ValueError, KeyError, TypeError):
                continue
    return None


def _arm_file(sample_dir: Path, *candidates: tuple[str, ...]) -> Path:
    """The first of ``candidates`` that exists under ``sample_dir``.

    An older flat layout is still accepted so an existing results tree yields
    a baseline without being restructured. Falls back to the first candidate, so a
    genuinely missing file is reported against the current layout.
    """
    sample_dir = Path(sample_dir)
    for parts in candidates:
        path = sample_dir.joinpath(*parts)
        if path.is_file():
            return path
    return sample_dir.joinpath(*candidates[0])


def rmlst_file(sample_dir: Path) -> Path:
    return _arm_file(sample_dir, ("rmlst", "species_calls.csv"), ("species_calls.csv",))


def marker16s_file(sample_dir: Path) -> Path:
    return _arm_file(sample_dir, ("16s", "marker16s_calls.csv"),
                     ("marker16s", "marker16s_calls.csv"))


def _est_counts(row: dict) -> float | None:
    try:
        return (float(row.get("unambiguous_est_counts") or 0)
                + float(row.get("ambiguous_est_counts") or 0))
    except ValueError:
        return None


def rmlst_profile(path: Path) -> dict[str, float]:
    """Per-species read counts from species_calls.csv.

    Every row counts, reported or not. A baseline is a picture of what the
    reagents carry, and an organism sitting below the reporting bar in the
    control is exactly the background we want to measure - filtering to reported
    calls would record a blank as clean when it was not.
    """
    profile: dict[str, float] = {}
    if not path.is_file():
        return profile
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            species = (row.get("species") or "").strip()
            reads = _est_counts(row)
            if species and reads is not None and reads > 0:
                profile[species] = profile.get(species, 0.0) + reads
    return profile


def reported_species(path: Path) -> dict[str, float]:
    """Species the rMLST arm REPORTED in a sample, with their est_counts.

    This is the credible-source side of the mechanism question: only an
    organism rMLST named from its allele profile counts as something a sample
    carries in bulk. EM background rows and demoted near neighbours do not.
    """
    reported: dict[str, float] = {}
    if not path.is_file():
        return reported
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if str(row.get("reported", "")).strip().lower() != "true":
                continue
            species = (row.get("species") or "").strip()
            if not species:
                continue
            reads = _est_counts(row)
            reported[species] = reported.get(species, 0.0) + (reads or 0.0)
    return reported


def marker16s_profile(path: Path) -> dict[str, float]:
    """Per-taxon RAW read counts from marker16s_calls.csv, at both ranks.

    Recorded under the genus as well as the full name, because an organism the
    control saw only at genus level still has to gate a sample call made at
    species level, and vice versa.

    ``raw_reads``, not ``reads_assigned``: the copy-number-corrected column
    must not reach a ratio comparing one organism to itself. This has to match
    what ``reconcile_cli`` feeds the numerator - if the two ever disagree, every
    RPM-r in the run is wrong by the ratio of two copy numbers. See the "Copy
    number must not reach the ratio" section of ``rpm.py``.
    """
    profile: dict[str, float] = {}
    if not path.is_file():
        return profile
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            try:
                reads = float(row.get("raw_reads") or 0)
            except ValueError:
                continue
            if reads <= 0:
                continue
            for key in ((row.get("taxon") or "").strip(),
                        (row.get("genus") or "").strip()):
                if key:
                    profile[key] = max(profile.get(key, 0.0), reads)
    return profile


ControlProfiles = dict[str, dict[str, tuple[float, int]]]


def control_profiles(in_dir: Path, controls: list[str]) -> tuple[
        ControlProfiles, ControlProfiles, dict[str, int]]:
    """Both arms' profiles for each named control, plus each control's depth."""
    rmlst_controls: ControlProfiles = {}
    s16_controls: ControlProfiles = {}
    depths: dict[str, int] = {}
    for name in controls:
        sample_dir = Path(in_dir) / name
        if not sample_dir.is_dir():
            raise FileNotFoundError(f"no such sample folder: {sample_dir}")
        depth = total_reads(sample_dir)
        if depth is None:
            raise FileNotFoundError(
                f"no run_info.json under {sample_dir}; did both arms finish?")
        depths[name] = depth
        rmlst = rmlst_profile(rmlst_file(sample_dir))
        s16 = marker16s_profile(marker16s_file(sample_dir))
        rmlst_controls[name] = {k: (v, depth) for k, v in rmlst.items()}
        s16_controls[name] = {k: (v, depth) for k, v in s16.items()}
    return rmlst_controls, s16_controls, depths


def sample_dirs(in_dir: Path, controls: list[str],
                exclude_prefixes: tuple[str, ...] = DEFAULT_EXCLUDE_PREFIXES) -> list[Path]:
    """The run's non-control sample folders: those with an rMLST result file."""
    excluded = set(controls)
    found = []
    for path in sorted(Path(in_dir).iterdir()):
        name = path.name
        if not path.is_dir() or name in excluded or name.startswith("."):
            continue
        if is_control_name(name, exclude_prefixes):
            continue
        if rmlst_file(path).is_file():
            found.append(path)
    return found


def scan_cohort(
    in_dir: Path,
    controls: list[str],
    reagent_baselines: dict[tuple[str, str, str], Baseline],
    run_id: str,
    *,
    exclude_prefixes: tuple[str, ...] = DEFAULT_EXCLUDE_PREFIXES,
    rpm_r_min: float = DEFAULT_RPM_R_MIN,
    rpm_min: float = DEFAULT_RPM_MIN,
    min_ntc_rpm: float = MIN_NTC_RPM,
) -> dict[str, CohortEvidence]:
    """What the run's samples credibly carry, per arm.

    A sample is a **credible source** for an organism when the rMLST arm
    reported that species there AND, measured against the reagent-only baseline
    on that sample's own depth, the organism clears the gate on at least one
    arm (no control entry on that arm, or RPM-r >= ``rpm_r_min``). The second
    half is what stops a marginal call in a low-biomass sample - a
    *Corynebacterium* on the reporting bar in a blank-like specimen - from
    flipping the whole run's mechanism for that genus.

    Returns the strongest credible source per organism and per arm: rMLST
    est_counts for the rMLST arm, 16S raw reads for the 16S arm, keyed by
    species and by genus.
    """
    species: dict[str, dict[str, tuple[str, float]]] = {arm: {} for arm in ARMS}
    genera: dict[str, dict[str, tuple[str, float]]] = {arm: {} for arm in ARMS}

    for sample_dir in sample_dirs(in_dir, controls, exclude_prefixes):
        reported = reported_species(rmlst_file(sample_dir))
        if not reported:
            continue
        depth = total_reads(sample_dir)
        if not depth:
            continue
        s16 = {k.lower(): v for k, v in marker16s_profile(marker16s_file(sample_dir)).items()}
        for name, est in reported.items():
            genus = name.split()[0]
            per_arm = {
                "rmlst": est,
                "marker16s": s16.get(name.lower()) or s16.get(genus.lower()),
            }
            credible = False
            for arm, reads in per_arm.items():
                if reads is None:
                    continue
                verdict = evaluate(
                    name, arm, reads, depth,
                    lookup_baseline(reagent_baselines, run_id, arm, name, genus=genus),
                    rpm_r_min=rpm_r_min, rpm_min=rpm_min, min_ntc_rpm=min_ntc_rpm)
                if verdict.passed:
                    credible = True
                    break
            if not credible:
                continue
            for arm, reads in per_arm.items():
                if reads is None or reads <= 0:
                    continue
                for key, table in ((name.lower(), species[arm]), (genus.lower(), genera[arm])):
                    current = table.get(key)
                    if current is None or reads > current[1]:
                        table[key] = (sample_dir.name, reads)

    return {arm: CohortEvidence(species[arm], genera[arm]) for arm in ARMS}


def build_run_baselines(
    in_dir: Path,
    controls: list[str],
    run_id: str,
    *,
    min_control_reads: int = MIN_CONTROL_READS,
    cohort: bool = True,
    exclude_prefixes: tuple[str, ...] = DEFAULT_EXCLUDE_PREFIXES,
    rpm_r_min: float = DEFAULT_RPM_R_MIN,
    rpm_min: float = DEFAULT_RPM_MIN,
    min_ntc_rpm: float = MIN_NTC_RPM,
    max_ntc_source_frac: float = MAX_NTC_SOURCE_FRAC,
) -> list[dict[str, object]]:
    """Baseline rows for a run, from its named negative controls.

    Two passes when ``cohort`` is on: reagent-only rows first, which are what
    the credible-source test in :func:`scan_cohort` measures against; then the
    same rows with each organism's mechanism filled in.
    """
    rmlst_controls, s16_controls, _ = control_profiles(in_dir, controls)

    def rows_with(evidence: dict[str, CohortEvidence] | None) -> list[dict[str, object]]:
        return (build_baselines(rmlst_controls, run_id, "rmlst",
                                min_control_reads=min_control_reads,
                                cohort=evidence["rmlst"] if evidence else None,
                                max_ntc_source_frac=max_ntc_source_frac)
                + build_baselines(s16_controls, run_id, "marker16s",
                                  min_control_reads=min_control_reads,
                                  cohort=evidence["marker16s"] if evidence else None,
                                  max_ntc_source_frac=max_ntc_source_frac))

    reagent_rows = rows_with(None)
    if not cohort or not reagent_rows:
        return reagent_rows
    evidence = scan_cohort(
        in_dir, controls, baselines_from_rows(reagent_rows), run_id,
        exclude_prefixes=exclude_prefixes, rpm_r_min=rpm_r_min,
        rpm_min=rpm_min, min_ntc_rpm=min_ntc_rpm)
    return rows_with(evidence)


def describe_rows(rows: list[dict[str, object]]) -> str:
    """A short human summary of a baseline: counts, the highest levels, the notes."""
    cross = [row for row in rows if row.get("mechanism") == MECHANISM_CROSS]
    lines = [f"Wrote {len(rows)} baseline rows: {len(cross)} cross-contamination, "
             f"{len(rows) - len(cross)} reagent"]
    if cross:
        lines.append("Cross-contamination (gated on absolute reads; source sample in brackets):")
        for row in sorted(cross, key=lambda r: -float(r["reads_ntc"]))[:12]:
            lines.append(f"  {row['arm']:10s} {str(row['organism'])[:34]:36s} "
                         f"{float(row['reads_ntc']):8.1f} reads  [{row['source_sample']}]")
    lines.append("Highest reagent levels (these are the organisms hardest to report):")
    reagent = [row for row in rows if row.get("mechanism") != MECHANISM_CROSS]
    for row in sorted(reagent, key=lambda r: -float(r["rpm_ntc"]))[:8]:
        lines.append(f"  {row['arm']:10s} {str(row['organism'])[:34]:36s} "
                     f"{float(row['rpm_ntc']):12.1f} RPM")
    noted = [row for row in rows if row.get("note")]
    if noted:
        lines.append("Notes:")
        for row in noted:
            lines.append(f"  {row['arm']:10s} {row['organism']}: {row['note']}")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="bridge build-baseline", description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--in-dir", type=Path, required=True,
                        help="results root holding each sample's folder")
    parser.add_argument("-n", "--negative", action="append", required=True,
                        help="negative control sample name; repeat for several")
    parser.add_argument("-r", "--run-id", default="*",
                        help="sequencing run id to tag the baseline with "
                             "('*' = pooled fallback for any run)")
    parser.add_argument("-o", "--out", type=Path, required=True)
    parser.add_argument("--min-control-reads", type=int, default=MIN_CONTROL_READS,
                        help="controls shallower than this are excluded: one "
                             "stray read dominates their RPM and would block "
                             "genuine organisms")
    parser.add_argument("--no-cohort", action="store_true",
                        help="skip the cohort scan: every organism is reagent "
                             "background and gated on RPM-r")
    parser.add_argument("--exclude-prefix", action="append", default=None,
                        help="sample-name prefix to leave out of the cohort scan "
                             f"(repeatable; default {' '.join(DEFAULT_EXCLUDE_PREFIXES)})")
    parser.add_argument("--max-ntc-source-frac", type=float, default=MAX_NTC_SOURCE_FRAC,
                        help="note a row whose control reads exceed this fraction "
                             "of its source sample's reads")
    parser.add_argument("--rpm-r-min", type=float, default=DEFAULT_RPM_R_MIN)
    parser.add_argument("--min-ntc-rpm", type=float, default=MIN_NTC_RPM)
    args = parser.parse_args(argv)

    try:
        _, _, depths = control_profiles(args.in_dir, args.negative)
    except FileNotFoundError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    for name in args.negative:
        depth = depths[name]
        print(f"{name}: {depth:,} read pairs processed")
        if depth < args.min_control_reads:
            print(f"  WARNING: below --min-control-reads ({args.min_control_reads:,}); "
                  f"this control will be excluded", file=sys.stderr)

    rows = build_run_baselines(
        args.in_dir, args.negative, args.run_id,
        min_control_reads=args.min_control_reads, cohort=not args.no_cohort,
        exclude_prefixes=tuple(args.exclude_prefix or DEFAULT_EXCLUDE_PREFIXES),
        rpm_r_min=args.rpm_r_min, min_ntc_rpm=args.min_ntc_rpm,
        max_ntc_source_frac=args.max_ntc_source_frac)
    if not rows:
        print("\nERROR: no baseline rows produced. Either every control was too "
              "shallow, or neither arm found anything in them.", file=sys.stderr)
        return 1

    write_baselines(rows, args.out)
    print()
    print(describe_rows(rows))
    print(f"\n-> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
