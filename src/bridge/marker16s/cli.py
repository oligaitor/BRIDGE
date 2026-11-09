"""Run BRIDGE's 16S marker arm on one sample (``python -m bridge.marker16s.cli``)."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

from bridge import __version__

from bridge.tools import (
    close_logging, configure_logging, reset_runtime_accounting,
    runtime_accounting, stage,
)
from .aggregate import DEFAULT_COPY_NUMBER, aggregate_taxa, copy_number_coverage
from .alignment_metrics import alignment_metrics_for_calls, resolve_tools
from .calls import (
    DEFAULT_GENUS_MIN_FRAC, DEFAULT_MIN_GENUS_READS, DEFAULT_SPECIES_DOMINANCE,
    DEFAULT_SPECIES_MARGIN, marker16s_calls, write_marker16s_calls,
)
from .quant import quant_16s, read_run_info
from .reference import find_index, load_reference

# Default MIMt location for standalone use; bridge always passes --reference.
DEFAULT_DB_DIR = Path("/data/db/16S")
DEFAULT_REFERENCE = DEFAULT_DB_DIR / "MIMt-16S_M2c_25_04_tax_bacteria_tax_in_id.fna"

# rrn copy numbers, ``name<TAB>copy_number``, built from rrnDB by
# build_copy_number_table.py. Picked up automatically when it exists: without a
# table every taxon is divided by the same 4.2, which corrects nothing.
DEFAULT_COPY_NUMBER_TABLE = DEFAULT_DB_DIR / "rrnDB_copy_numbers.tsv"


def positive_float(value: str) -> float:
    number = float(value)
    if number <= 0:
        raise argparse.ArgumentTypeError(f"must be greater than zero: {value}")
    return number


def fraction(value: str) -> float:
    number = float(value)
    if not 0.0 <= number <= 1.0:
        raise argparse.ArgumentTypeError(f"must be between 0 and 1: {value}")
    return number


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m bridge.marker16s.cli",
        description=(
            "BRIDGE 16S marker arm: abundance and taxonomic assignment against MIMt. "
            "Emits marker16s_calls.csv for the dual-marker reconciler."
        ),
    )
    parser.add_argument("--reads", nargs=2, required=True, type=Path,
                        metavar=("R1", "R2"), help="paired FASTQ files")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("-o", "--output", required=True, type=Path)
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE,
                        help="MIMt FASTA with taxonomy in the header "
                             "(use the _tax_in_id variant)")
    parser.add_argument("--kallisto-index", type=Path, default=None,
                        help="prebuilt index; defaults to <reference>.idx")
    parser.add_argument("--kallisto", default="kallisto")
    parser.add_argument("-t", "--threads", type=int, default=8)
    parser.add_argument("--copy-number-table", type=Path, default=None,
                        help="TSV of name<TAB>copy_number (rrnDB-style) used to "
                             f"correct 16S multi-copy bias. Defaults to "
                             f"{DEFAULT_COPY_NUMBER_TABLE} when that file "
                             "exists; without a table a database-wide mean is "
                             "used for every taxon, which corrects nothing")
    parser.add_argument("--no-copy-number-table", action="store_true",
                        help="run uncorrected even if the default table exists "
                             "- for measuring what the correction changes")
    parser.add_argument("--genus-min-frac", type=fraction,
                        default=DEFAULT_GENUS_MIN_FRAC,
                        help="genus must hold this share of the profile")
    parser.add_argument("--min-genus-reads", type=positive_float,
                        default=DEFAULT_MIN_GENUS_READS,
                        help="absolute floor on a genus's COPY-NUMBER-CORRECTED "
                             "reads, the same figure the proportion tests use "
                             "(about 42 raw reads at the 4.2 default copy number, "
                             "10 for a single-copy organism)")
    parser.add_argument("--species-dominance", type=fraction,
                        default=DEFAULT_SPECIES_DOMINANCE,
                        help="top species' share of its genus needed to name it")
    parser.add_argument("--species-margin", type=fraction,
                        default=DEFAULT_SPECIES_MARGIN,
                        help="gap to the runner-up species needed to name it")
    parser.add_argument("--no-reuse", action="store_true",
                        help="re-run the quant even if abundance.tsv exists")
    parser.add_argument("--no-alignment-metrics", action="store_true",
                        help="skip the identity/coverage stage (minimap2 + samtools + "
                             "viral_consensus + blastn on the called genera)")
    parser.add_argument("--minimap2", default="minimap2")
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--viral-consensus", default="viral_consensus")
    parser.add_argument("--blastn", default="blastn")
    parser.add_argument("--makeblastdb", default="makeblastdb")
    parser.add_argument("-v", "--verbose", action="store_true")
    return parser


def load_copy_numbers(path: Path | None) -> dict[str, float] | None:
    """Read an rrnDB-style ``name<TAB>copy_number`` table, keyed lowercase.

    Accepts both species (``escherichia coli``) and genus (``escherichia``)
    rows; :func:`..aggregate._copy_number` prefers the more specific match.
    """
    if path is None:
        return None
    table: dict[str, float] = {}
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            try:
                table[parts[0].strip().lower()] = float(parts[1])
            except ValueError:
                continue
    if not table:
        raise ValueError(f"No usable rows in copy-number table {path}")
    return table


def resolve_copy_number_table(args) -> Path | None:
    """Which copy-number table to use: explicit, then the default, then none.

    An explicitly named table that does not exist is an error, not a silent
    fallback - if someone passed ``--copy-number-table`` they meant it, and a
    run that quietly corrects nothing is worse than one that stops.
    """
    if args.no_copy_number_table:
        return None
    if args.copy_number_table is not None:
        if not args.copy_number_table.is_file():
            raise FileNotFoundError(
                f"no copy-number table at {args.copy_number_table}")
        return args.copy_number_table
    return (DEFAULT_COPY_NUMBER_TABLE
            if DEFAULT_COPY_NUMBER_TABLE.is_file() else None)


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    args.output.mkdir(parents=True, exist_ok=True)
    handlers = configure_logging(args.output / "marker16s.log", args.verbose)
    reset_runtime_accounting()
    logger = logging.getLogger("marker16s")

    try:
        logger.info("marker16s_status=started sample=%s version=%s",
                    args.sample_id, __version__)

        index = find_index(args.reference, args.kallisto_index)
        table_path = resolve_copy_number_table(args)
        copy_numbers = load_copy_numbers(table_path)
        if copy_numbers is None:
            # Not fatal - the arm still works, and every within-sample ratio is
            # unaffected because the divisor is then a constant. What breaks is
            # any comparison BETWEEN organisms or between the two arms.
            logger.warning(
                "copy_number_status=no_table - every taxon will be divided by "
                "%.1f, so the correction is a constant and cross-organism and "
                "cross-arm magnitudes stay biased by copy number. Build one "
                "with build_copy_number_table.py.", DEFAULT_COPY_NUMBER)

        with stage("marker16s_reference"):
            reference = load_reference(args.reference)
            logger.info("reference_status=loaded targets=%d", len(reference))

        abundance = quant_16s(
            list(args.reads), index, args.output,
            threads=args.threads, kallisto=args.kallisto,
            reuse=not args.no_reuse,
        )

        with stage("marker16s_aggregate"):
            taxa = aggregate_taxa(abundance, reference, copy_numbers=copy_numbers)

        with stage("marker16s_calls"):
            calls = marker16s_calls(
                args.sample_id, taxa,
                genus_min_frac=args.genus_min_frac,
                species_dominance=args.species_dominance,
                species_margin=args.species_margin,
                min_genus_reads=args.min_genus_reads,
            )

        # Identity and coverage from a real alignment, for the called genera.
        # Optional at runtime: a missing tool skips it with one warning and the
        # calls file carries blanks in those columns.
        metrics = {}
        alignment_info: dict[str, object] = {"enabled": not args.no_alignment_metrics}
        if not args.no_alignment_metrics:
            with stage("marker16s_alignment_metrics"):
                tools = resolve_tools(
                    minimap2=args.minimap2, samtools=args.samtools,
                    viral_consensus=args.viral_consensus, blastn=args.blastn,
                    makeblastdb=args.makeblastdb)
                if tools is None:
                    alignment_info["skipped_reason"] = "missing tool (see log)"
                else:
                    metrics = alignment_metrics_for_calls(
                        calls, reference, args.reference, list(args.reads),
                        args.output / "alignment", tools=tools, threads=args.threads)
                    alignment_info["genera_attempted"] = len(metrics)
                    alignment_info["genera_scored"] = sum(
                        1 for item in metrics.values() if item.status == "ok")
        calls_path = write_marker16s_calls(
            calls, args.output / "marker16s_calls.csv", metrics)

        manifest = {
            "schema": "bridge.rmlst-16s",
            "schema_version": 2,
            "package": "bridge",
            "package_version": __version__,
            "sample_id": args.sample_id,
            "reference": str(args.reference),
            "kallisto_index": str(index),
            "reference_targets": len(reference),
            "reads": [str(path) for path in args.reads],
            "copy_number_table": str(table_path) if table_path else None,
            "copy_number_entries": len(copy_numbers) if copy_numbers else 0,
            # Recorded per run so a result can be read knowing how much of the
            # correction was real rather than the 4.2 fallback.
            "copy_number_coverage": copy_number_coverage(taxa),
            "parameters": {
                "threads": args.threads,
                "genus_min_frac": args.genus_min_frac,
                "min_genus_reads": args.min_genus_reads,
                "species_dominance": args.species_dominance,
                "species_margin": args.species_margin,
            },
            "taxa_detected": len(taxa),
            "genera_called": sum(1 for c in calls if c.strength != "silent_weak"),
            # n_processed is the RPM denominator.
            "read_counts": read_run_info(args.output),
            "alignment_metrics": alignment_info,
            "runtime_breakdown": runtime_accounting(),
        }
        manifest_path = args.output / "marker16s_manifest.json"
        with open(manifest_path, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(json.dumps(manifest, indent=2) + "\n")

        logger.info("marker16s_status=complete calls=%s manifest=%s",
                    calls_path, manifest_path)
        return 0
    except Exception as error:  # noqa: BLE001 - surfaced to the caller's exit code
        logger.error("marker16s_status=failed error=%s", error)
        raise
    finally:
        close_logging(handlers)


if __name__ == "__main__":
    raise SystemExit(main())
