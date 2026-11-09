"""``python -m bridge.amr.cli`` - the resistance arm on its own.

    python -m bridge.amr.cli build-db --db-dir /data/db/resfinder_db
    python -m bridge.amr.cli sample   --sample-id S1 --reads R1.fastq.gz R2.fastq.gz --out <results>
    python -m bridge.amr.cli report   --out <results> --run-id RUN1

Inside a full run you never need this: ``bridge run`` runs the arm per
sample (``--arms amr`` for this arm alone, on a new or an existing results
tree) and ``bridge finalise`` runs the report. ``bridge build-amr-db``
and ``bridge amr-report`` hand their arguments straight to the matching
subcommand here.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from bridge import __version__
from .settings import DEFAULT_CONFIG, load_amr_settings


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="python -m bridge.amr.cli", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", action="version", version=f"bridge {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    def common(p: argparse.ArgumentParser) -> argparse.ArgumentParser:
        p.add_argument("-c", "--config", type=Path, default=None,
                       help=f"config.yaml (default: {DEFAULT_CONFIG})")
        return p

    db = common(sub.add_parser("build-db", help="index a ResFinder checkout for KMA"))
    db.add_argument("--db-dir", type=Path, required=True,
                    help="resfinder_db checkout (holds the *.fsa files and phenotypes.txt)")
    db.add_argument("--prefix", type=Path, default=None,
                    help="index prefix (default: <db-dir>/kma/resfinder); "
                         "put this in references.amr_db")
    db.add_argument("--kma", default=None, help="kma binary (default: tools.kma)")
    db.add_argument("--kma-index", default=None,
                    help="kma_index binary (default: beside kma, else 'kma index')")
    db.add_argument("--dry-run", action="store_true",
                    help="concatenate and write the manifest, skip indexing")

    sample = common(sub.add_parser("sample", help="run the arm for one sample"))
    sample.add_argument("--sample-id", required=True)
    sample.add_argument("--reads", nargs=2, type=Path, required=True, metavar=("R1", "R2"))
    sample.add_argument("--out", type=Path, required=True,
                        help="results root; output lands in <out>/<SAMPLE>/amr/")
    sample.add_argument("-t", "--threads", type=int, default=4)

    report = common(sub.add_parser("report", help="gate, attribute and interpret a results tree"))
    report.add_argument("--out", type=Path, required=True)
    report.add_argument("--run-id", required=True)
    report.add_argument("--control", action="append", default=[],
                        help="negative control sample name (repeatable); otherwise by prefix")
    report.add_argument("--no-cohort", action="store_true",
                        help="reagent-only baseline: skip the cross-contamination attribution")
    return parser


def cmd_build_db(args) -> int:
    from .database import build_db
    settings = load_amr_settings(args.config)
    manifest = build_db(args.db_dir, args.prefix, kma=args.kma or settings.kma,
                        kma_index=args.kma_index, dry_run=args.dry_run)
    print(f"Templates : {manifest['templates']}")
    print(f"Commit    : {manifest['git_commit'] or 'not a git checkout'}"
          + (f"  ({manifest['git_date']})" if manifest["git_date"] else ""))
    print(f"Index     : {manifest['index_prefix']}")
    print("\nIn config.yaml:")
    print(f"  references.amr_db:         {manifest['index_prefix']}")
    print(f"  references.amr_phenotypes: {manifest['phenotypes']}")
    return 0


def cmd_sample(args) -> int:
    from .kma import run_arm
    settings = load_amr_settings(args.config)
    out = Path(args.out) / args.sample_id / "amr"
    path = run_arm(settings, args.sample_id, args.reads[0].resolve(), args.reads[1].resolve(),
                   out, threads=args.threads)
    print(f"-> {path}")
    return 0


def cmd_report(args) -> int:
    from bridge.reconcile.baseline import is_control_name
    from .report import amr_sample_dirs, describe, finalise_amr
    settings = load_amr_settings(args.config)
    out = args.out.resolve()
    controls = list(args.control) or [
        folder.name for folder in amr_sample_dirs(out)
        if is_control_name(folder.name, settings.control_prefixes)]
    result = finalise_amr(settings, out, args.run_id, controls=controls,
                          cohort=not args.no_cohort)
    print(describe(result, out))
    return 1 if result.failed else 0


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    handler = {"build-db": cmd_build_db, "sample": cmd_sample, "report": cmd_report}[args.command]
    try:
        return handler(args)
    except (FileNotFoundError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
