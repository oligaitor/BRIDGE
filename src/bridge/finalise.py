"""``bridge finalise`` - from arm outputs to the run-level table.

Five steps, seconds each, no external tools, so this runs on a login node and
is what you re-run after changing a threshold in ``config.yaml``:

1. build ``rpm_baselines.csv`` from the run's negative control(s);
2. reconcile every non-control sample into ``<SAMPLE>/bacterial_ID/``;
3. concatenate the per-sample summaries into ``identification_summary_all.csv``;
4. where the AMR arm ran: gate, attribute and interpret the resistance genes
   (``bridge.amr.report``) into ``<SAMPLE>/amr_report/`` and the run tables
   ``amr_summary_all.csv`` / ``predicted_phenotype_all.csv``. After step 2 on
   purpose, because attribution reads the reconciled organism calls. A
   failure here is reported but cannot undo steps 1-3;
5. one outcome per sample on ``<out>/sample_bacteria.list`` (``--sample-list``), negatives
   and failures included, into ``sample_outcomes_all.csv``, and the short
   organism + resistance view into ``clinical_summary_all.csv`` (``outcome.py``,
   ``clinical.py``). Like step 4, a failure here leaves steps 1-4 standing.

Controls are named with ``--control`` or recognised by the prefixes in
``hpc.controls_prefixes``.

A pooled baseline (``--pooled-baseline``, or ``rpm.pooled_baseline`` in the
config; built by ``bridge pool-baseline``) adds run_id ``*`` rows beside the
run's own in both baseline files. The lookup reads them only where this run's
blank never saw the organism or gene, so a run with a good blank is gated as
before for everything that blank measured.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path

from bridge.reconcile import baseline
from bridge.reconcile.reconcile_cli import main as reconcile_main
from bridge.reconcile.summary_csv import combine_summaries

from .settings import Settings, is_control

SKIP_DIRS = {"oe"}
ID_BASELINE_ARMS = ("rmlst", "marker16s")
AMR_BASELINE_ARM = "amr"


@dataclass
class FinaliseResult:
    controls: list[str] = field(default_factory=list)
    baseline_path: Path | None = None
    baseline_rows: int = 0
    pooled_path: Path | None = None
    pooled_rows: int = 0                 # identification rows taken from the pool
    reconciled: list[str] = field(default_factory=list)
    skipped: list[tuple[str, str]] = field(default_factory=list)
    failed: list[tuple[str, str]] = field(default_factory=list)
    combined_path: Path | None = None
    combined_samples: int = 0
    amr: object | None = None            # bridge.amr.report.AmrResult, when the arm ran
    amr_error: str | None = None
    outcome: object | None = None        # outcome.OutcomeResult
    outcome_error: str | None = None


def find_arm_outputs(sample_dir: Path) -> tuple[Path | None, Path | None]:
    """The two arm result files, accepting the older flat layout too."""
    rmlst = baseline._arm_file(sample_dir, ("rmlst", "species_calls.csv"),
                               ("species_calls.csv",))
    s16 = baseline._arm_file(sample_dir, ("16s", "marker16s_calls.csv"),
                             ("marker16s", "marker16s_calls.csv"))
    return (rmlst if rmlst.is_file() else None, s16 if s16.is_file() else None)


def sample_dirs(out_root: Path) -> list[Path]:
    return sorted(
        path for path in Path(out_root).iterdir()
        if path.is_dir() and path.name not in SKIP_DIRS and not path.name.startswith(".")
    )


def detect_controls(out_root: Path, prefixes: tuple[str, ...]) -> list[str]:
    return [path.name for path in sample_dirs(out_root) if is_control(path.name, prefixes)]


def pooled_baseline_path(settings: Settings, explicit: Path | None = None) -> Path | None:
    """``--pooled-baseline`` if given, else ``rpm.pooled_baseline`` (relative to the config)."""
    if explicit is not None:
        return Path(explicit).expanduser().resolve()
    value = settings.cfg.get("rpm.pooled_baseline")
    if value in (None, "", False):
        return None
    path = Path(str(value)).expanduser()
    return path if path.is_absolute() else (settings.path.parent / path).resolve()


def load_pooled_rows(path: Path) -> list[dict[str, str]]:
    """The rows of a pooled baseline file, refusing one that is not pooled.

    A per-run file passed by mistake would carry another run's id, match
    nothing here and silently gate nothing - so it is an error, not a no-op.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"no such pooled baseline: {path}")
    with open(path, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    other = sorted({str(row.get("run_id", "")).strip() for row in rows} - {"*"})
    if other:
        raise ValueError(f"{path} is not a pooled baseline: it has run_id "
                         f"{', '.join(other)} (expected only '*'; build one with "
                         "`bridge pool-baseline`)")
    return rows


def finalise_run(settings: Settings, out_root: Path, run_id: str, *,
                 controls: list[str] | None = None,
                 min_control_reads: int | None = None,
                 cohort: bool = True,
                 pooled_baseline: Path | None = None,
                 sample_list: Path | None = None,
                 recurrence: Path | None = None,
                 verbose: bool = False) -> FinaliseResult:
    out_root = Path(out_root).resolve()
    if not out_root.is_dir():
        raise FileNotFoundError(f"no such results directory: {out_root}")
    result = FinaliseResult()
    prefixes = settings.hpc.controls_prefixes
    cfg = settings.cfg
    result.controls = list(controls) if controls else detect_controls(out_root, prefixes)

    pooled: list[dict[str, str]] = []
    result.pooled_path = pooled_baseline_path(settings, pooled_baseline)
    if result.pooled_path is not None:
        pooled = load_pooled_rows(result.pooled_path)
    pooled_id = [row for row in pooled if row.get("arm") in ID_BASELINE_ARMS]
    pooled_amr = [row for row in pooled if row.get("arm") == AMR_BASELINE_ARM]

    # 1. baseline - with the cohort scan that attributes each organism in the
    #    blank to reagent background or cross-contamination from a run positive.
    baseline_path = out_root / "rpm_baselines.csv"
    id_controls = [name for name in result.controls
                   if any(path is not None for path in find_arm_outputs(out_root / name))]
    if result.controls and not id_controls:
        print("NOTE: no negative control has identification-arm output (an AMR-only "
              "run?); the identification baseline is skipped.")
    rows: list[dict] = []
    if id_controls:
        rows = baseline.build_run_baselines(
            out_root, id_controls, run_id,
            min_control_reads=(min_control_reads if min_control_reads is not None
                               else baseline.MIN_CONTROL_READS),
            cohort=cohort, exclude_prefixes=prefixes,
            rpm_r_min=float(cfg.get("rpm.rpm_r_min", baseline.DEFAULT_RPM_R_MIN)),
            rpm_min=float(cfg.get("rpm.rpm_min", baseline.DEFAULT_RPM_MIN)),
            min_ntc_rpm=float(cfg.get("rpm.min_ntc_rpm", baseline.MIN_NTC_RPM)),
            max_ntc_source_frac=float(cfg.get("rpm.max_ntc_source_frac",
                                              baseline.MAX_NTC_SOURCE_FRAC)))
        print(f"Baseline   : {baseline_path}  (from {', '.join(id_controls)})")
        print(baseline.describe_rows(rows))
    if id_controls or pooled_id:
        baseline.write_baselines(rows + pooled_id, baseline_path)
        result.baseline_path = baseline_path
        result.baseline_rows = len(rows) + len(pooled_id)
        result.pooled_rows = len(pooled_id)
        if pooled_id:
            print(f"Pooled     : {len(pooled_id)} identification rows from {result.pooled_path}"
                  + ("" if rows else "  (this run has no usable blank of its own)"))
    elif not result.controls:
        print("WARNING: no negative control found (name one with --control, or use a "
              f"name starting with {' '.join(prefixes)}). The contamination gate is "
              "skipped and every call is flagged rpm_gate_not_applied.")
        if baseline_path.is_file():
            print(f"           (an older {baseline_path.name} is present and will be used)")
            result.baseline_path = baseline_path

    # 2. reconcile
    control_set = set(result.controls)
    for sample_dir in sample_dirs(out_root):
        name = sample_dir.name
        if name in control_set:
            result.skipped.append((name, "negative control"))
            continue
        rmlst, s16 = find_arm_outputs(sample_dir)
        if rmlst is None:
            result.skipped.append((name, "no rmlst/species_calls.csv (arm incomplete)"
                                   if s16 is not None or not has_amr(sample_dir)
                                   else "AMR arm only (no identification)"))
            continue
        argv = ["--sample-id", name, "--rmlst-calls", str(rmlst),
                "-c", str(settings.path), "-o", str(sample_dir / "bacterial_ID"),
                "--run-id", run_id]
        if s16 is not None:
            argv += ["--marker16s-calls", str(s16)]
        total = baseline.total_reads(sample_dir)
        if total is not None:
            argv += ["--total-reads", str(total)]
        if result.baseline_path is not None and result.baseline_path.is_file():
            argv += ["--rpm-baselines", str(result.baseline_path)]
        if verbose:
            argv.append("-v")
        try:
            code = reconcile_main(argv)
        except Exception as error:  # noqa: BLE001 - one sample must not stop the run
            result.failed.append((name, f"{type(error).__name__}: {error}"))
            continue
        if code:
            result.failed.append((name, f"exit {code}"))
        else:
            result.reconciled.append(name)

    # 3. combine
    combined, n = combine_summaries(out_root, out_root / "identification_summary_all.csv")
    result.combined_path, result.combined_samples = combined, n

    # 4. resistance, where the AMR arm ran
    if any(has_amr(path) for path in sample_dirs(out_root)):
        result.amr, result.amr_error = run_amr_report(settings, out_root, run_id,
                                                      result.controls, cohort=cohort,
                                                      pooled_rows=pooled_amr)

    # 5. one outcome per sample, from what this finalise produced. Only the
    #    samples reconciled now count: combine_summaries also picks up a
    #    summary left by an earlier finalise for a sample that failed today.
    try:
        from .outcome import outcome_settings, outcomes_from_tree, read_csv
        reconciled = set(result.reconciled)
        failures = {name: why for name, why in result.skipped + result.failed
                    if why != "negative control"}
        result.outcome = outcomes_from_tree(
            out_root, run_id, outcome_settings(cfg, recurrence), controls=result.controls,
            control_prefixes=prefixes, sample_list=sample_list, failures=failures,
            id_rows=[row for row in read_csv(result.combined_path)
                     if row.get("sample") in reconciled],
            use_amr=result.amr is not None)
    except Exception as error:  # noqa: BLE001 - must not undo steps 1-4
        result.outcome_error = f"{type(error).__name__}: {error}"
    return result


def has_amr(sample_dir: Path) -> bool:
    return (Path(sample_dir) / "amr" / "amr_hits.csv").is_file()


def run_amr_report(settings: Settings, out_root: Path, run_id: str, controls: list[str],
                   *, cohort: bool = True, pooled_rows: list[dict] | None = None):
    """Step 4 on its own (``bridge amr-report``). Returns (AmrResult | None, error)."""
    try:
        from bridge.amr.report import finalise_amr
        from bridge.amr.settings import from_config
        amr = finalise_amr(from_config(settings.cfg, settings.path), out_root, run_id,
                           controls=controls, cohort=cohort, pooled_rows=pooled_rows)
    except Exception as error:  # noqa: BLE001 - must not undo the identification steps
        return None, f"{type(error).__name__}: {error}"
    return amr, None


def describe(result: FinaliseResult, out_root: Path) -> str:
    lines = [
        "",
        f"Reconciled : {len(result.reconciled)}",
        f"Skipped    : {len(result.skipped)}",
    ]
    for name, why in result.skipped:
        if why != "negative control":
            lines.append(f"             {name}: {why}")
    if result.failed:
        lines.append(f"FAILED     : {len(result.failed)}")
        for name, why in result.failed:
            lines.append(f"             {name}: {why}")
    lines += [
        f"Combined   : {result.combined_path}  ({result.combined_samples} samples)",
        "",
        f"Results    : {out_root}/<SAMPLE>/bacterial_ID/reconciled_calls.json",
        f"             {out_root}/<SAMPLE>/bacterial_ID/identification_summary.csv",
        f"             {out_root}/identification_summary_all.csv   <- the run in one table",
        f"             {out_root}/rpm_baselines.csv                <- what the control held",
    ]
    if result.pooled_path is not None:
        lines.append(f"Pooled     : {result.pooled_path}  ({result.pooled_rows} identification "
                     "rows; calls gated on them are flagged rpm_pooled_baseline)")
    if result.amr_error:
        lines += ["", f"AMR FAILED : {result.amr_error}",
                  "             (identification results above are unaffected; fix and run "
                  "bridge amr-report)"]
    elif result.amr is not None:
        from bridge.amr.report import describe as describe_amr
        lines.append(describe_amr(result.amr, out_root))
    if result.outcome_error:
        lines += ["", f"OUTCOMES FAILED: {result.outcome_error}",
                  "             (tables above are unaffected; fix and run bridge outcomes)"]
    elif result.outcome is not None:
        from .outcome import describe as describe_outcomes
        lines += [""] + describe_outcomes(result.outcome)
    return "\n".join(lines)
