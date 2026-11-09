from __future__ import annotations

import argparse
from math import isfinite
import os
from pathlib import Path
import sys

from .dependencies import DependencyError, check_dependencies
from .models import Config


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be at least one")
    return parsed


def positive_float(value: str) -> float:
    parsed = float(value)
    if not isfinite(parsed) or parsed <= 0:
        raise argparse.ArgumentTypeError("must be greater than zero")
    return parsed


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        prog="python -m bridge.rmlst",
        description=(
            "BRIDGE rMLST arm: post-kallisto runtime clustering, consensus, "
            "exact allele and mixture-aware rMLST calling"
        ),
    )
    result.add_argument(
        "--check-deps", action="store_true",
        help="check dependencies for the selected pipeline route and exit",
    )
    result.add_argument("--allele-target", type=Path)
    result.add_argument("--rmlstseqdef", type=Path)
    result.add_argument("--reads", type=Path, nargs=2, metavar=("READ1", "READ2"))
    result.add_argument("--sample-id", default="sample")
    result.add_argument("-o", "--output-folder", type=Path)
    result.add_argument("--tempdir", type=Path)
    result.add_argument(
        "--keeptmp", action="store_true",
        help=(
            "retain the work directory and intermediate diagnostic outputs; "
            "by default only final species outputs, manifest, and log are kept"
        ),
    )
    result.add_argument(
        "--verbose", action="store_true",
        help="show detailed pipeline messages and internally called tool output",
    )
    result.add_argument("-t", "--threads", type=positive_int, default=4)
    result.add_argument(
        "--memory-limit-gb", type=positive_float,
        help=(
            "overall memory budget in GB used to size concurrent CliqueSNV "
            "Java heaps; overrides detected system/cgroup memory"
        ),
    )
    result.add_argument("--kallisto", default="kallisto")
    result.add_argument(
        "--kallisto-index", type=Path,
        help=(
            "existing full-reference Kallisto index; bypass automatic index "
            "discovery and construction, and reuse its k-mer size for the "
            "candidate pseudobAM index"
        ),
    )
    result.add_argument("--mash", default="mash")
    result.add_argument(
        "--short-read-aligner", "--aligner",
        choices=("bwa-mem2", "minimap2"), default="minimap2",
        help=(
            "cluster-read alignment backend (default: minimap2, using "
            "validated short-read parameters and MAPQ 10)"
        ),
    )
    result.add_argument("--bwa-mem2", default="bwa-mem2")
    result.add_argument("--minimap2", default="minimap2")
    result.add_argument("--samtools", default="samtools")
    result.add_argument(
        "--cliquesnv",
        default=os.environ.get("BRIDGE_CLIQUESNV", "cliquesnv"),
        help=(
            "CliqueSNV executable or JAR (default: BRIDGE_CLIQUESNV, "
            "otherwise cliquesnv)"
        ),
    )
    result.add_argument(
        "--cliquesnv-time-limit-seconds", "--cliquesnv-runtime-limit-seconds",
        dest="cliquesnv_time_limit_seconds", type=positive_int,
        help=(
            "pass CliqueSNV's native -tl time limit for each invocation; on expiry "
            "CliqueSNV returns a consensus haplotype (default: CliqueSNV's 10800 seconds)"
        ),
    )
    result.add_argument(
        "--viral-consensus", default="viral_consensus",
        help="ViralConsensus executable (default: viral_consensus)",
    )
    result.add_argument(
        "--viralconsensus-only", action="store_true",
        help=(
            "always produce one ViralConsensus sequence per runtime cluster; "
            "disable CliqueSNV subconsensus reconstruction"
        ),
    )
    result.add_argument("--blastn", default="blastn")
    result.add_argument("--makeblastdb", default="makeblastdb")
    result.add_argument(
        "--blastdb-cache-dir", type=Path, default=None,
        help=(
            "directory holding the shared .bridge_blastdb cache "
            "(default: the parent of --output-folder). Name the run root here "
            "whenever --output-folder is nested deeper than one level inside "
            "it, so concurrent samples keep sharing one database"
        ),
    )
    result.add_argument("--kallisto-bootstraps", type=positive_int, default=100)
    result.add_argument(
        "--bootstrap-support-threshold", type=int, default=100,
        help=(
            "requested per-locus support threshold; an empty locus relaxes to "
            "its highest observed support (default: 100)"
        ),
    )
    result.add_argument("--kallisto-seed", type=int, default=1000)
    result.add_argument("--initial-abundance", action="store_true", help="use initial rather than bootstrap-mean est_counts")
    result.add_argument("--mash-k", type=positive_int, default=21)
    result.add_argument("--mash-threshold", type=float, default=0.01)
    result.add_argument(
        "--consensus-min-fragments", type=positive_int, default=1,
        help="minimum assigned read pairs needed to attempt consensus (default: 1)",
    )
    result.add_argument(
        "--subconsensus-min-proportion", type=float, default=0.01,
        help="retain haplotypes strictly above this proportion (default: 0.01)",
    )
    result.add_argument(
        "--subconsensus-min-fragments", type=positive_int, default=20,
        help="minimum estimated supporting read pairs (default: 20)",
    )
    result.add_argument(
        "--subconsensus-min-depth", type=float, default=10.0,
        help="minimum projected haplotype depth (default: 10)",
    )
    result.add_argument(
        "--subconsensus-min-breadth", type=float, default=0.80,
        help="fraction of the locus meeting projected haplotype depth (default: 0.80)",
    )
    result.add_argument(
        "--subconsensus-min-end-distance", type=int, default=5,
        help="bases from a read end required for interior variant support (default: 5)",
    )
    result.add_argument(
        "--minor-allele-min-depth", type=positive_int, default=5,
        help=(
            "second-highest quality-filtered base depth that triggers CliqueSNV "
            "(default: 5; no strand or frequency filter)"
        ),
    )
    result.add_argument("--minimum-profile-loci", type=positive_int, default=3)
    result.add_argument(
        "--require-accepted-consensus", action="store_true",
        help=(
            "count a locus as profile evidence only when its consensus exactly "
            "matches a catalogued allele. Off by "
            "default: the strict rule silently penalises low-abundance "
            "organisms, which typically have broad locus support but imperfect "
            "assembly."
        ),
    )
    result.add_argument(
        "--min-unambiguous-loci", type=positive_int, default=3,
        help=(
            "loci mapping to a species ALONE before it may be reported "
            "(default: 3). Loci that are merely compatible with the species do "
            "not count, so a near neighbour cannot be reported on evidence that "
            "belongs to the organism actually present"
        ),
    )
    result.add_argument(
        "--distinguish-subspecies", action="store_true",
        help=(
            "treat each species/subspecies combination in rmlstseqdef as a "
            "distinct taxon, including rows with a blank subspecies"
        ),
    )
    result.add_argument(
        "--cluster-read-assignment", choices=("compatible", "exclusive"),
        default="compatible",
        help="include a fragment in every compatible runtime cluster, or only uniquely compatible fragments",
    )
    result.add_argument(
        "--filter-composite-consensus",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "treat a consensus absent from its cluster but constructible from "
            "cluster-member bases as ambiguous evidence (default: enabled; "
            "disable with --no-filter-composite-consensus)"
        ),
    )
    return result


def main(argv: list[str] | None = None) -> None:
    argument_parser = parser()
    args = argument_parser.parse_args(argv)
    if args.check_deps:
        try:
            statuses = check_dependencies(
                kallisto=args.kallisto,
                mash=args.mash,
                short_read_aligner=args.short_read_aligner,
                bwa_mem2=args.bwa_mem2,
                minimap2=args.minimap2,
                samtools=args.samtools,
                cliquesnv=args.cliquesnv,
                viral_consensus=args.viral_consensus,
                viralconsensus_only=args.viralconsensus_only,
                blastn=args.blastn,
                makeblastdb=args.makeblastdb,
            )
        except DependencyError as error:
            print(error, file=sys.stderr)
            raise SystemExit(1) from error
        for status in statuses:
            print(
                f"available\t{status.name}\t{status.requirement}\t{status.resolved}"
            )
        print(f"Dependency check passed: {len(statuses)} dependencies available")
        return

    missing = [
        option for option, value in (
            ("--allele-target", args.allele_target),
            ("--rmlstseqdef", args.rmlstseqdef),
            ("--reads", args.reads),
            ("--output-folder", args.output_folder),
        ) if value is None
    ]
    if missing:
        argument_parser.error(
            "the following arguments are required unless --check-deps is used: "
            + ", ".join(missing)
        )

    assert args.output_folder is not None
    assert args.allele_target is not None
    assert args.rmlstseqdef is not None
    assert args.reads is not None
    output = args.output_folder.expanduser().resolve()
    from .pipeline import run_pipeline

    run_pipeline(Config(
        allele_target=args.allele_target.expanduser().resolve(),
        seqdef=args.rmlstseqdef.expanduser().resolve(),
        reads=tuple(path.expanduser().resolve() for path in args.reads),
        sample_id=args.sample_id,
        output=output,
        tempdir=(args.tempdir or output / "work").expanduser().resolve(),
        threads=args.threads,
        memory_limit_gb=args.memory_limit_gb,
        keep_temp=args.keeptmp,
        verbose=args.verbose,
        kallisto=args.kallisto,
        kallisto_index=(
            args.kallisto_index.expanduser().resolve()
            if args.kallisto_index is not None else None
        ),
        mash=args.mash,
        short_read_aligner=args.short_read_aligner,
        bwa_mem2=args.bwa_mem2,
        minimap2=args.minimap2,
        samtools=args.samtools,
        cliquesnv=args.cliquesnv,
        cliquesnv_time_limit_seconds=args.cliquesnv_time_limit_seconds,
        viral_consensus=args.viral_consensus,
        viralconsensus_only=args.viralconsensus_only,
        blastn=args.blastn,
        makeblastdb=args.makeblastdb,
        blastdb_cache_dir=(
            args.blastdb_cache_dir.expanduser().resolve()
            if args.blastdb_cache_dir is not None else None
        ),
        bootstraps=args.kallisto_bootstraps,
        support_threshold=args.bootstrap_support_threshold,
        seed=args.kallisto_seed,
        use_bootstrap_mean=not args.initial_abundance,
        mash_k=args.mash_k,
        mash_threshold=args.mash_threshold,
        min_fragments=args.consensus_min_fragments,
        subconsensus_min_proportion=args.subconsensus_min_proportion,
        subconsensus_min_fragments=args.subconsensus_min_fragments,
        subconsensus_min_depth=args.subconsensus_min_depth,
        subconsensus_min_breadth=args.subconsensus_min_breadth,
        subconsensus_min_end_distance=args.subconsensus_min_end_distance,
        minor_allele_min_depth=args.minor_allele_min_depth,
        min_profile_loci=args.minimum_profile_loci,
        min_unambiguous_loci=args.min_unambiguous_loci,
        require_accepted_consensus=args.require_accepted_consensus,
        read_assignment=args.cluster_read_assignment,
        filter_composite_consensus=args.filter_composite_consensus,
        distinguish_subspecies=args.distinguish_subspecies,
    ))
