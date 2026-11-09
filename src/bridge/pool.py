"""``bridge pool-baseline`` - one reagent baseline from many runs' blanks.

    bridge pool-baseline -o /scratch/project/runs/pooled_baselines.csv \\
        --run RUN1=/scratch/project/runs/RUN1 \\
        --run RUN2=/scratch/project/runs/RUN2 \\
        --extra-control 'Fungal*'

Why
---
A blank only gates a run if it was sequenced deeply enough to measure
anything (``rpm.MIN_CONTROL_READS``). A run whose every blank falls below that
has no baseline at all: every call reads ``no_control_background`` and reagent
genes that deeper blanks carry go ungated in the resistance arm.

This builds a run_id ``*`` baseline that ``bridge finalise
--pooled-baseline`` adds beside each run's own. The lookup reads a pooled row
only where the run's own blank never saw the organism, so it is a fallback,
never an override.

What goes in
------------
Each run's controls (recognised by ``hpc.controls_prefixes``), plus any sample
matching ``--extra-control`` - a glob, case-insensitive - for specimens that
should hold no bacteria at all, such as pure fungal-culture dilutions. Each
run's baseline is built exactly as ``finalise`` builds it (same depth cut,
same cohort scan), and only its **reagent** rows are pooled. Cross-
contamination needs a source on the same run; another run's carry-over says
nothing about this one's, and pooling it would hold genuine organisms
everywhere else.

How rows combine
----------------
Per arm and organism (or gene), the highest RPM across every contributing
control - the maximum, as within one run, because the question a baseline
answers is "how high does this get in a blank". ``control_ids`` lists every
control that saw it as ``RUN/CONTROL``; ``note`` names the run that set the
level.

Alongside the baseline, ``<out stem>_controls.csv`` lists every control that
was considered, its depth on each arm and whether it was deep enough to use.
"""

from __future__ import annotations

import argparse
import csv
import fnmatch
import sys
from dataclasses import dataclass, field, replace
from pathlib import Path

from bridge.reconcile import baseline
from bridge.reconcile.rpm import BASELINE_COLUMNS, MECHANISM_REAGENT, write_baselines

from .finalise import find_arm_outputs, sample_dirs
from .settings import Settings, is_control

POOLED_SOURCE = "pooled_max_across_runs"
CONTROL_COLUMNS = ["run_id", "control", "kind", "id_read_pairs", "amr_read_pairs", "used"]


@dataclass
class PoolResult:
    rows: list[dict] = field(default_factory=list)
    controls: list[dict] = field(default_factory=list)   # one per (run, control)
    dropped_cross: int = 0                               # cross-contamination rows left out
    runs_without_usable_control: list[str] = field(default_factory=list)


def parse_run(value: str) -> tuple[str, Path]:
    """``RUN_ID=DIR`` -> (run id, results root)."""
    run_id, sep, path = value.partition("=")
    if not sep or not run_id.strip() or not path.strip():
        raise ValueError(f"--run must be RUN_ID=DIR, got {value!r}")
    return run_id.strip(), Path(path.strip()).expanduser().resolve()


def matches_any(name: str, patterns: list[str]) -> bool:
    folded = name.casefold()
    return any(fnmatch.fnmatchcase(folded, pattern.casefold()) for pattern in patterns)


def pool_controls(out_root: Path, prefixes: tuple[str, ...],
                  extra_patterns: list[str]) -> list[tuple[str, str]]:
    """(name, kind) of every folder that enters the pool from one run."""
    found = []
    for path in sample_dirs(out_root):
        if is_control(path.name, prefixes):
            found.append((path.name, "negative_control"))
        elif matches_any(path.name, extra_patterns):
            found.append((path.name, "extra_control"))
    return found


def _run_rows(settings: Settings, run_id: str, out_root: Path,
              controls: list[str], *, min_control_reads: int | None,
              cohort: bool) -> tuple[list[dict], dict[str, tuple[int | None, int | None]]]:
    """One run's baseline rows on all three arms, and each control's depths."""
    from bridge.amr.gate import build_amr_baselines, hits_file, total_read_pairs
    from bridge.amr.report import amr_sample_dirs
    from bridge.amr.settings import from_config

    cfg = settings.cfg
    prefixes = settings.hpc.controls_prefixes
    depths: dict[str, tuple[int | None, int | None]] = {}
    for name in controls:
        folder = out_root / name
        depths[name] = (baseline.total_reads(folder), total_read_pairs(folder))

    rows: list[dict] = []
    id_controls = [name for name in controls
                   if depths[name][0] is not None
                   and any(path is not None for path in find_arm_outputs(out_root / name))]
    if id_controls:
        rows += baseline.build_run_baselines(
            out_root, id_controls, run_id,
            min_control_reads=(min_control_reads if min_control_reads is not None
                               else baseline.MIN_CONTROL_READS),
            cohort=cohort, exclude_prefixes=prefixes,
            rpm_r_min=float(cfg.get("rpm.rpm_r_min", baseline.DEFAULT_RPM_R_MIN)),
            rpm_min=float(cfg.get("rpm.rpm_min", baseline.DEFAULT_RPM_MIN)),
            min_ntc_rpm=float(cfg.get("rpm.min_ntc_rpm", baseline.MIN_NTC_RPM)),
            max_ntc_source_frac=float(cfg.get("rpm.max_ntc_source_frac",
                                              baseline.MAX_NTC_SOURCE_FRAC)))

    amr_controls = [name for name in controls if hits_file(out_root / name).is_file()]
    if amr_controls:
        amr = from_config(cfg, settings.path)
        gate = (amr.gate if min_control_reads is None
                else replace(amr.gate, min_control_reads=min_control_reads))
        excluded = set(controls)
        samples = [folder for folder in amr_sample_dirs(out_root) if folder.name not in excluded]
        rows += build_amr_baselines(out_root, amr_controls, samples, run_id, gate,
                                    amr.thresholds, cohort=cohort)
    return rows, depths


def merge_rows(per_run: list[tuple[str, dict]]) -> list[dict]:
    """Reagent rows from several runs -> one ``*`` row per (arm, organism): the maximum."""
    groups: dict[tuple[str, str], list[tuple[str, dict]]] = {}
    for run_id, row in per_run:
        key = (str(row["arm"]), str(row["organism"]).strip().lower())
        groups.setdefault(key, []).append((run_id, row))
    merged = []
    for (arm, _), items in sorted(groups.items()):
        top_run, top = max(items, key=lambda item: float(item[1]["rpm_ntc"]))
        contributors = sorted({f"{run_id}/{name.strip()}"
                               for run_id, row in items
                               for name in str(row.get("control_ids") or "").split(";")
                               if name.strip()})
        runs = sorted({run_id for run_id, _ in items})
        merged.append({
            **{column: "" for column in BASELINE_COLUMNS},
            "run_id": "*",
            "arm": arm,
            "organism": top["organism"],
            "rpm_ntc": top["rpm_ntc"],
            "reads_ntc": top["reads_ntc"],
            "total_reads_ntc": top["total_reads_ntc"],
            "n_controls": len(contributors),
            "control_ids": "; ".join(contributors),
            "source": POOLED_SOURCE,
            "mechanism": MECHANISM_REAGENT,
            "note": (f"highest on {top_run}; seen in {len(contributors)} control(s) "
                     f"on {len(runs)} run(s)"),
        })
    return merged


def pool_baselines(settings: Settings, runs: list[tuple[str, Path]], *,
                   extra_patterns: list[str] | None = None,
                   min_control_reads: int | None = None,
                   cohort: bool = True) -> PoolResult:
    result = PoolResult()
    extra_patterns = list(extra_patterns or [])
    id_floor = (min_control_reads if min_control_reads is not None
                else baseline.MIN_CONTROL_READS)
    per_run: list[tuple[str, dict]] = []
    for run_id, out_root in runs:
        if not out_root.is_dir():
            raise FileNotFoundError(f"no such results directory for {run_id}: {out_root}")
        found = pool_controls(out_root, settings.hpc.controls_prefixes, extra_patterns)
        names = [name for name, _ in found]
        rows, depths = _run_rows(settings, run_id, out_root, names,
                                 min_control_reads=min_control_reads, cohort=cohort)
        usable_here = False
        for name, kind in found:
            id_depth, amr_depth = depths[name]
            used = any(depth is not None and depth >= id_floor
                       for depth in (id_depth, amr_depth))
            usable_here = usable_here or used
            result.controls.append({
                "run_id": run_id, "control": name, "kind": kind,
                "id_read_pairs": "" if id_depth is None else id_depth,
                "amr_read_pairs": "" if amr_depth is None else amr_depth,
                "used": "yes" if used else "no (too shallow)",
            })
        if not usable_here:
            result.runs_without_usable_control.append(run_id)
        for row in rows:
            if row.get("mechanism") == MECHANISM_REAGENT:
                per_run.append((run_id, row))
            else:
                result.dropped_cross += 1
    result.rows = merge_rows(per_run)
    return result


def write_controls(controls: list[dict], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CONTROL_COLUMNS)
        writer.writeheader()
        writer.writerows(controls)
    return path


def describe(result: PoolResult) -> str:
    lines = ["Controls considered:"]
    for row in result.controls:
        lines.append(f"  {row['run_id']:12s} {row['control']:28s} {row['kind']:17s} "
                     f"id={row['id_read_pairs'] or '-':>9} amr={row['amr_read_pairs'] or '-':>9}  "
                     f"{row['used']}")
    if result.runs_without_usable_control:
        lines.append("Runs with no usable control of their own (the pool is their only "
                     f"background): {', '.join(result.runs_without_usable_control)}")
    by_arm: dict[str, int] = {}
    for row in result.rows:
        by_arm[row["arm"]] = by_arm.get(row["arm"], 0) + 1
    lines.append(f"Pooled rows: {len(result.rows)} "
                 f"({', '.join(f'{arm} {n}' for arm, n in sorted(by_arm.items())) or 'none'}); "
                 f"{result.dropped_cross} cross-contamination row(s) left out (run-specific)")
    for arm in sorted(by_arm):
        top = sorted((row for row in result.rows if row["arm"] == arm),
                     key=lambda row: -float(row["rpm_ntc"]))[:8]
        lines.append(f"Highest {arm}:")
        for row in top:
            lines.append(f"  {str(row['organism'])[:36]:38s} {float(row['rpm_ntc']):12.1f} RPM  "
                         f"({row['note']})")
    return "\n".join(lines)


def main(settings: Settings, argv: list[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="bridge pool-baseline", description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run", action="append", required=True, metavar="RUN_ID=DIR",
                        help="a finished results root and the run id it was run with "
                             "(repeatable; one per run)")
    parser.add_argument("-o", "--out", type=Path, required=True,
                        help="pooled baseline CSV to write")
    parser.add_argument("--extra-control", action="append", default=[], metavar="GLOB",
                        help="also pool samples matching this name pattern (case-insensitive; "
                             "repeatable) - specimens that should hold no bacteria, e.g. 'Fungal*'")
    parser.add_argument("--min-control-reads", type=int, default=None,
                        help="depth a control needs to be used (default: as finalise)")
    parser.add_argument("--no-cohort", action="store_true",
                        help="pool every row as reagent, skipping each run's "
                             "cross-contamination attribution")
    args = parser.parse_args(argv)
    try:
        runs = [parse_run(value) for value in args.run]
        result = pool_baselines(settings, runs, extra_patterns=args.extra_control,
                                min_control_reads=args.min_control_reads,
                                cohort=not args.no_cohort)
    except (ValueError, FileNotFoundError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    print(describe(result))
    controls_path = write_controls(result.controls,
                                   args.out.with_name(f"{args.out.stem}_controls.csv"))
    if not result.rows:
        print("\nERROR: no pooled rows - every control was too shallow, or held nothing.",
              file=sys.stderr)
        return 1
    write_baselines(result.rows, args.out)
    print(f"\n-> {args.out}\n-> {controls_path}")
    return 0
