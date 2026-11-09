"""``bridge`` - one command for the whole pipeline.

    bridge check
    bridge run --reads-list sample.list --out <dir> --run-id RUN1
    bridge finalise --out <dir> --run-id RUN1
    bridge outcomes --out <dir> --run-id RUN1     # outcome per sample, alone

``run`` submits one PBS job per sample and a finalise job that waits for them;
``--scheduler local`` does the same work in this session. Everything else is a
building block those two use, exposed for re-running one step by hand.

The resistance arm runs with the others by default (``--arms all``), or alone:

    bridge run --arms amr --reads-list sample.list --out <existing dir> --run-id RUN1
    bridge amr-report --out <dir> --run-id RUN1       # re-interpret only

Intermediates are removed as each arm finishes (``--keep-temp`` keeps them);
``bridge clean --out <dir>`` measures, and with ``--apply`` removes, any that
remain in a results tree.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

from . import __version__
from .settings import DEFAULT_CONFIG, arms_from, is_control, load_settings


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="bridge",
        description="BRIDGE: dual-marker (rMLST + 16S) bacterial identification plus ResFinder "
                    "resistance (AMR): run, gate and reconcile a sequencing run.",
    )
    arms_help = ("all (default) | both (= rmlst,16s) | rmlst | 16s | amr | a comma list "
                 "such as rmlst,amr")
    parser.add_argument("--version", action="version", version=f"bridge {__version__}")
    parser.add_argument("--config", type=Path, default=None,
                        help=f"config.yaml (default: {DEFAULT_CONFIG})")
    sub = parser.add_subparsers(dest="command", required=True)

    def add_config(subparser: argparse.ArgumentParser) -> argparse.ArgumentParser:
        """Accept ``--config`` after the subcommand as well as before it.

        argparse only accepts a top-level option before the subcommand, which
        is not where anyone writes it - the PBS scripts included. A separate
        dest keeps the two from overwriting each other: a subparser's default
        is copied over the parent's value, so reusing ``config`` would erase a
        path given in front of the subcommand.
        """
        subparser.add_argument("--config", dest="config_after", type=Path, default=None,
                               help=argparse.SUPPRESS)
        return subparser

    check = add_config(sub.add_parser("check", help="verify Python, binaries, references and the scheduler"))
    check.add_argument("--arms", default="all", help=arms_help)
    check.add_argument("--scheduler", choices=("pbs", "local"), default="pbs")
    check.add_argument("--out", type=Path, default=None, help="results directory to test for writability")

    run = add_config(sub.add_parser("run", help="run every sample of a read list, then finalise"))
    source = run.add_mutually_exclusive_group(required=True)
    source.add_argument("--reads-list", type=Path,
                        help="one absolute R1 path per line; sample = parent directory")
    source.add_argument("--reads-dir", type=Path,
                        help="directory of <SAMPLE>/*_R1*.fastq.gz (or flat *_R1*.fastq.gz)")
    run.add_argument("--out", type=Path, required=True,
                     help="results root; each sample lands in <out>/<SAMPLE>/")
    run.add_argument("--run-id", required=True,
                     help="sequencing run id; tags the baseline and every result")
    run.add_argument("--arms", default="all",
                     help=arms_help + ". A sample is skipped when every requested arm's "
                          "output already exists, so --arms amr on an existing results "
                          "tree adds the AMR arm without touching rmlst/ or 16s/")
    run.add_argument("--scheduler", choices=("pbs", "local"), default="pbs")
    run.add_argument("--control", action="append", default=[],
                     help="negative control sample name (repeatable); otherwise recognised by prefix")
    run.add_argument("--max-submit", type=int, default=0,
                     help="submit at most N jobs this time (0 = no limit)")
    run.add_argument("--dry-run", action="store_true", help="render the job scripts, submit nothing")
    run.add_argument("--force", action="store_true",
                     help="re-run the requested arm(s) even where their output already exists "
                          "(e.g. --arms 16s --force to refresh the 16S arm on an old tree)")
    run.add_argument("--no-finalise", action="store_true", help="submit/run the arms only")
    run.add_argument("--skip-check", action="store_true",
                     help="submit without verifying the tools and reference files first "
                          "(the check is there so a missing program fails once, here, "
                          "rather than once per sample on a compute node)")
    run.add_argument("--threads", type=int, default=4, help="local scheduler only")
    run.add_argument("--memory-gb", type=float, default=None, help="local scheduler only")
    run.add_argument("--tempdir", type=Path, default=None, help="local scheduler only")
    run.add_argument("--keep-temp", action="store_true",
                     help="keep each arm's intermediates (16S alignments, KMA .aln/.mat.gz, "
                          "the rMLST work folder); by default they are removed once the arm "
                          "succeeds. See `bridge clean`")
    run.add_argument("-v", "--verbose", action="store_true")

    sample = add_config(sub.add_parser("sample", help="the requested arms for one sample, in this process (what a PBS job runs)"))
    sample.add_argument("--sample-id", required=True)
    sample.add_argument("--reads", nargs=2, type=Path, required=True, metavar=("R1", "R2"))
    sample.add_argument("--out", type=Path, required=True)
    sample.add_argument("--arms", default="all", help=arms_help)
    sample.add_argument("--threads", type=int, default=4)
    sample.add_argument("--memory-gb", type=float, default=None)
    sample.add_argument("--tempdir", type=Path, default=None)
    sample.add_argument("--keep-temp", action="store_true",
                        help="keep each arm's intermediates (removed by default)")
    sample.add_argument("-v", "--verbose", action="store_true")

    clean = add_config(sub.add_parser(
        "clean", help="measure, and with --apply remove, the intermediates of a finished "
                      "results tree (what finalise and amr-report never read)"))
    clean.add_argument("--out", type=Path, required=True)
    clean.add_argument("--apply", action="store_true",
                       help="delete them (without it, only report what would go and its size)")
    clean.add_argument("--blast-cache", action="store_true",
                       help="also remove the run's shared rMLST BLAST database "
                            "(.bridge_blastdb; `bridge run` rebuilds it if needed)")

    fin = add_config(sub.add_parser("finalise", aliases=["finalize"],
                         help="baseline + reconcile + run table (login node, seconds)"))
    fin.add_argument("--out", type=Path, required=True)
    fin.add_argument("--run-id", required=True)
    fin.add_argument("--control", action="append", default=[])
    fin.add_argument("--min-control-reads", type=int, default=None)
    fin.add_argument("--no-cohort", action="store_true",
                     help="build a reagent-only baseline (every organism gated on RPM-r), "
                          "skipping the run-level cross-contamination attribution")
    fin.add_argument("--pooled-baseline", type=Path, default=None,
                     help="cross-run baseline from `bridge pool-baseline`, used where this "
                          "run's blank never saw an organism or gene (default: "
                          "rpm.pooled_baseline in the config)")
    fin.add_argument("--sample-list", type=Path, default=None,
                     help="the run's read list, so a sample that never ran still gets an "
                          "outcome row (default: <out>/sample_bacteria.list, written by "
                          "bridge run; else <out>/sample.list, read only)")
    fin.add_argument("--recurrence", type=Path, default=None,
                     help="genus recurrence table from `bridge recurrence`, used to split "
                          "held samples (default: outcome.held_split.recurrence_table)")
    fin.add_argument("-v", "--verbose", action="store_true")

    outc = add_config(sub.add_parser(
        "outcomes", help="one outcome per sample (negatives and failures included) and the "
                         "clinician summary, without re-reconciling (seconds)"))
    where = outc.add_mutually_exclusive_group(required=True)
    where.add_argument("--out", type=Path, help="a results tree (as finalise sees it)")
    where.add_argument("--tables", type=Path,
                       help="a folder of run-level *_all.csv tables + a sample list "
                            "(e.g. a local copy of a run's results)")
    outc.add_argument("--run-id", default=None,
                      help="required with --out; with --tables, read from the tables")
    outc.add_argument("--sample-list", type=Path, default=None,
                      help="default: sample_bacteria.list, else sample.list, beside the results")
    outc.add_argument("--control", action="append", default=[])
    outc.add_argument("--depth-csv", type=Path, default=None,
                      help="--tables only: sample,total_read_pairs (otherwise depth is unknown)")
    outc.add_argument("--dest", type=Path, default=None,
                      help="--tables only: write the two tables here instead of into the folder")
    outc.add_argument("--recurrence", type=Path, default=None,
                      help="genus recurrence table (default: outcome.held_split.recurrence_table)")

    rec = add_config(sub.add_parser(
        "recurrence", help="genus prevalence per run, from runs' identification_summary_all.csv: "
                           "the cross-run evidence that splits held samples (seconds)"))
    rec.add_argument("runs", nargs="+", type=Path,
                     help="run folders (results trees or copies of their tables)")
    rec.add_argument("--out", type=Path, required=True, help="the table to write (CSV)")

    amr = add_config(sub.add_parser(
        "amr-report", help="re-interpret the AMR arm only: resistance gate, organism "
                           "attribution, predicted phenotype (login node, seconds)"))
    amr.add_argument("--out", type=Path, required=True)
    amr.add_argument("--run-id", required=True)
    amr.add_argument("--control", action="append", default=[])
    amr.add_argument("--no-cohort", action="store_true",
                     help="reagent-only resistance baseline (no cross-contamination attribution)")
    amr.add_argument("--pooled-baseline", type=Path, default=None,
                     help="cross-run baseline from `bridge pool-baseline` (default: "
                          "rpm.pooled_baseline in the config)")

    sub.add_parser("build-amr-db", add_help=False,
                   help="index a ResFinder checkout for the AMR arm (its own --help)")
    sub.add_parser("build-baseline", add_help=False,
                   help="rpm_baselines.csv from a control (its own --help)")
    sub.add_parser("pool-baseline", add_help=False,
                   help="one cross-run reagent baseline (run_id *) from several runs' "
                        "controls, for finalise --pooled-baseline (its own --help)")
    sub.add_parser("reconcile", add_help=False,
                   help="reconcile one sample (its own --help; -c defaults to --config)")
    sub.add_parser("evaluate", add_help=False,
                   help="re-score a run's identification_summary_all.csv under the reagent-only and "
                        "mechanism-attributed gate rules against a truth table (its own --help)")
    sub.add_parser("build-copy-numbers", add_help=False,
                   help="16S rrn copy-number table from an rrnDB download")
    return parser


def _config_path(args) -> Path | None:
    """``--config`` from either side of the subcommand."""
    return getattr(args, "config_after", None) or args.config


def _load(args):
    return load_settings(_config_path(args))


def cmd_check(args) -> int:
    from .check import run_check
    return run_check(_load(args), arms=arms_from(args.arms),
                     scheduler=args.scheduler, out=args.out)


def cmd_sample(args) -> int:
    from .reads import ReadPair
    from .sample import run_sample
    settings = _load(args)
    pair = ReadPair(args.sample_id, args.reads[0].resolve(), args.reads[1].resolve())
    run_sample(settings, pair, args.out, arms=arms_from(args.arms),
               threads=args.threads, memory_gb=args.memory_gb, tempdir=args.tempdir,
               verbose=args.verbose, keep_temp=args.keep_temp)
    return 0


def cmd_finalise(args) -> int:
    from .finalise import describe, finalise_run
    settings = _load(args)
    try:
        result = finalise_run(settings, args.out, args.run_id, controls=args.control or None,
                              min_control_reads=args.min_control_reads,
                              cohort=not args.no_cohort, pooled_baseline=args.pooled_baseline,
                              sample_list=args.sample_list, recurrence=args.recurrence,
                              verbose=args.verbose)
    except (FileNotFoundError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    print(describe(result, args.out.resolve()))
    return 1 if _finalise_failed(result) else 0


def _finalise_failed(result) -> bool:
    return bool(result.failed or result.amr_error or result.outcome_error
                or (result.amr is not None and result.amr.failed))


def cmd_outcomes(args) -> int:
    from .outcome import describe, outcome_settings, outcomes_from_tables, outcomes_from_tree
    settings = _load(args)
    if args.recurrence is not None and not args.recurrence.is_file():
        print(f"ERROR: no such recurrence table: {args.recurrence}", file=sys.stderr)
        return 1
    chosen = outcome_settings(settings.cfg, args.recurrence)
    prefixes = settings.hpc.controls_prefixes
    try:
        if args.out is not None:
            if not args.run_id:
                print("ERROR: --run-id is required with --out", file=sys.stderr)
                return 1
            out_root = args.out.resolve()
            if not out_root.is_dir():
                raise FileNotFoundError(f"no such results directory: {out_root}")
            from .finalise import detect_controls
            controls = list(args.control) or detect_controls(out_root, prefixes)
            result = outcomes_from_tree(out_root, args.run_id, chosen, controls=controls,
                                        control_prefixes=prefixes,
                                        sample_list=args.sample_list)
        else:
            if not args.tables.is_dir():
                raise FileNotFoundError(f"no such folder: {args.tables}")
            result = outcomes_from_tables(args.tables, chosen, run_id=args.run_id,
                                          sample_list=args.sample_list,
                                          depth_csv=args.depth_csv, dest=args.dest,
                                          controls=args.control, control_prefixes=prefixes)
    except (FileNotFoundError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    print("\n".join(describe(result)))
    split = chosen.split
    print(f"Held split : {'on' if split.enabled else 'off'}"
          + (f", recurrence from {len({r.get('run_id') for r in split.recurrence})} run(s)"
             if split.recurrence else ", no recurrence table (run-level evidence only)"))
    return 0


def cmd_recurrence(args) -> int:
    from .held_split import build_recurrence
    settings = _load(args)
    missing = [str(path) for path in args.runs
               if not (path / "identification_summary_all.csv").is_file()]
    for path in missing:
        print(f"SKIP {path}: no identification_summary_all.csv", file=sys.stderr)
    out, runs = build_recurrence(args.runs, args.out, settings.hpc.controls_prefixes)
    print(f"Recurrence : {out}  ({runs} run(s))")
    return 0 if runs else 1


def cmd_amr_report(args) -> int:
    from bridge.amr.report import amr_sample_dirs, describe as describe_amr
    from .finalise import AMR_BASELINE_ARM, load_pooled_rows, pooled_baseline_path, run_amr_report
    settings = _load(args)
    out_root = args.out.resolve()
    controls = list(args.control) or [
        folder.name for folder in amr_sample_dirs(out_root)
        if is_control(folder.name, settings.hpc.controls_prefixes)]
    pooled_path = pooled_baseline_path(settings, args.pooled_baseline)
    try:
        pooled = load_pooled_rows(pooled_path) if pooled_path is not None else []
    except (FileNotFoundError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    amr, error = run_amr_report(settings, out_root, args.run_id, controls,
                                cohort=not args.no_cohort,
                                pooled_rows=[row for row in pooled
                                             if row.get("arm") == AMR_BASELINE_ARM])
    if error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    print(describe_amr(amr, out_root))
    return 1 if amr.failed else 0


def cmd_clean(args) -> int:
    from .clean import clean_run, describe as describe_clean
    out_root = args.out.resolve()
    try:
        report = clean_run(out_root, apply=args.apply, blast_cache=args.blast_cache)
    except FileNotFoundError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    print(describe_clean(report, out_root))
    return 0


def cmd_run(args) -> int:
    from .finalise import describe, finalise_run
    from .pbs import (
        arms_complete, prepare_blast_cache,
        submit_finalise, submit_samples,
    )
    from .reads import discover_reads, parse_read_list
    from .sample import run_sample
    from .settings import is_control

    settings = _load(args)
    arms = arms_from(args.arms)
    pairs, problems = (parse_read_list(args.reads_list) if args.reads_list
                       else discover_reads(args.reads_dir))
    for problem in problems:
        print(f"SKIP {problem.sample}: {problem.reason}", file=sys.stderr)
    if not pairs:
        print("ERROR: no usable read pairs found", file=sys.stderr)
        return 1
    out_root = args.out.resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    from .outcome import record_sample_list
    list_path, added = record_sample_list(out_root, pairs)
    if added:
        print(f"Sample list: {list_path}  (+{added})")
    controls = list(args.control) or [
        pair.sample for pair in pairs if is_control(pair.sample, settings.hpc.controls_prefixes)]
    print(f"Samples   : {len(pairs)}  (controls: {', '.join(controls) or 'none found'})")
    print(f"Results   : {out_root}")

    if args.scheduler == "local":
        done = 0
        for pair in pairs:
            if not args.force and arms_complete(out_root / pair.sample, arms):
                print(f"SKIP {pair.sample}: already complete")
                continue
            run_sample(settings, pair, out_root, arms=arms, threads=args.threads,
                       memory_gb=args.memory_gb, tempdir=args.tempdir, verbose=args.verbose,
                       keep_temp=args.keep_temp)
            done += 1
        print(f"\nRan {done} sample(s) locally.")
        if args.no_finalise:
            return 0
        result = finalise_run(settings, out_root, args.run_id, controls=controls or None)
        print(describe(result, out_root))
        return 1 if _finalise_failed(result) else 0

    # Preflight. Submitting hundreds of jobs that each discover the same
    # missing program wastes a queue slot per sample and a day of waiting.
    if args.scheduler == "pbs" and not args.dry_run and not args.skip_check:
        from .check import run_check
        if run_check(settings, arms=arms, scheduler="pbs", out=out_root):
            print("\nNothing was submitted. Fix the MISSING lines above and run the same "
                  "command again, or pass --skip-check to submit regardless.",
                  file=sys.stderr)
            return 1
        print()

    if "rmlst" in arms and not args.dry_run:
        # Once, up front, on the login node: concurrent jobs must not each
        # build (and possibly corrupt) the shared BLAST database.
        print(f"BLAST cache: {prepare_blast_cache(settings, out_root)}")
    submissions = submit_samples(settings, pairs, out_root, args.run_id, arms=arms,
                                 max_submit=args.max_submit, dry_run=args.dry_run,
                                 force=args.force, keep_temp=args.keep_temp)
    counts: dict[str, int] = {}
    for item in submissions:
        counts[item.action] = counts.get(item.action, 0) + 1
        tag = {"submitted": "->", "dry_run": "DRY", "failed": "ERROR"}.get(item.action, "SKIP")
        job = f" {item.job_id}" if item.job_id else ""
        print(f"{tag:<5} {item.sample}: {item.action}{job}  {item.detail}")
    job_ids = [item.job_id for item in submissions if item.action == "submitted" and item.job_id]
    print()
    for action, n in sorted(counts.items()):
        print(f"{action:<17}: {n}")
    if args.no_finalise:
        return 0
    manual = (f"bridge finalise --out {out_root} --run-id {args.run_id}"
              + "".join(f" --control {name}" for name in controls))
    try:
        script, job_id = submit_finalise(settings, out_root, args.run_id, depends_on=job_ids,
                                         controls=controls, dry_run=args.dry_run)
    except (subprocess.CalledProcessError, OSError) as error:
        detail = (getattr(error, "stderr", "") or str(error)).strip()
        print(f"\nWARNING: the finalise job could not be submitted ({detail}).", file=sys.stderr)
        print("The sample jobs are unaffected. When they have all finished, run this on "
              "the login node instead:", file=sys.stderr)
        print(f"  {manual}", file=sys.stderr)
        return 1 if counts.get("failed") else 0
    if args.dry_run:
        print(f"\nDry run. Job scripts under {out_root / 'oe' / 'scripts'} (finalise: {script.name}).")
    else:
        print(f"\nFinalise job {job_id} waits for {len(job_ids)} sample job(s); watch with: qstat -u $USER")
        print(f"When it has run: {out_root / 'identification_summary_all.csv'}")
        if "amr" in arms:
            print(f"                 {out_root / 'amr_summary_all.csv'}")
        print(f"(or run it yourself at any time: {manual})")
    return 1 if counts.get("failed") else 0


def cmd_build_baseline(args, rest: list[str]) -> int:
    from bridge.reconcile.baseline import main as baseline_main
    return baseline_main(rest)


def cmd_pool_baseline(args, rest: list[str]) -> int:
    from .pool import main as pool_main
    return pool_main(_load(args), rest)


def cmd_reconcile(args, rest: list[str]) -> int:
    from bridge.reconcile.reconcile_cli import main as reconcile_main
    if "-c" not in rest and "--config" not in rest:
        rest = ["-c", str(_config_path(args) or DEFAULT_CONFIG)] + rest
    return reconcile_main(rest)


def cmd_build_copy_numbers(args, rest: list[str]) -> int:
    from .copy_numbers import main as copy_numbers_main
    return copy_numbers_main(rest)


def cmd_evaluate(args, rest: list[str]) -> int:
    from .evaluate import main as evaluate_main
    return evaluate_main(rest)


def cmd_build_amr_db(args, rest: list[str]) -> int:
    from bridge.amr.cli import main as amr_main
    if "-c" not in rest and "--config" not in rest:
        rest = ["-c", str(_config_path(args) or DEFAULT_CONFIG)] + rest
    return amr_main(["build-db"] + rest)


PASSTHROUGH = {
    "build-baseline": cmd_build_baseline,
    "pool-baseline": cmd_pool_baseline,
    "reconcile": cmd_reconcile,
    "evaluate": cmd_evaluate,
    "build-copy-numbers": cmd_build_copy_numbers,
    "build-amr-db": cmd_build_amr_db,
}


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = build_parser()
    # Pass-through subcommands hand their own arguments to the underlying tool.
    for name in PASSTHROUGH:
        if name in argv:
            index = argv.index(name)
            args = parser.parse_args(argv[:index + 1])
            return PASSTHROUGH[name](args, argv[index + 1:])
    args = parser.parse_args(argv)
    handler = {
        "check": cmd_check, "run": cmd_run, "sample": cmd_sample,
        "finalise": cmd_finalise, "finalize": cmd_finalise,
        "amr-report": cmd_amr_report, "clean": cmd_clean, "outcomes": cmd_outcomes,
        "recurrence": cmd_recurrence,
    }[args.command]
    return handler(args)


if __name__ == "__main__":
    raise SystemExit(main())
