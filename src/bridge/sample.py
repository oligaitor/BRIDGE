"""Run the requested arms for one sample, in this process.

This is what a PBS job executes (``bridge sample``), and what
``--scheduler local`` runs directly. Every path and threshold comes from
``config.yaml``.

Outputs land in ``<out>/<SAMPLE>/rmlst/``, ``<out>/<SAMPLE>/16s/`` and
``<out>/<SAMPLE>/amr/``, the layout ``bridge finalise`` reads. Each arm's
intermediates (alignments, KMA matrices) are removed once it succeeds, unless
``keep_temp``.
"""

from __future__ import annotations

import logging
from pathlib import Path

from bridge.marker16s.cli import main as marker16s_main
from bridge.rmlst.models import Config as RmlstConfig
from bridge.rmlst.pipeline import run_pipeline

from . import __version__
from .clean import human, remove_intermediates
from .reads import ReadPair
from .settings import ARMS, Settings

LOGGER = logging.getLogger("bridge.sample")


def code_versions() -> str:
    """Which code is actually running - printed at the top of every job log, so
    a result can be traced to the folder that produced it."""
    import bridge
    return f"bridge {__version__} ({Path(bridge.__file__).resolve().parent})"


def rmlst_config(settings: Settings, pair: ReadPair, out_root: Path, *,
                 threads: int, memory_gb: float | None, tempdir: Path | None,
                 verbose: bool, keep_temp: bool) -> RmlstConfig:
    refs, tools, cfg = settings.references, settings.tools, settings.cfg
    if refs.rmlst_alleles is None or refs.rmlst_seqdef is None:
        raise ValueError("references.rmlst_alleles and references.rmlst_seqdef must be set")
    output = (Path(out_root) / pair.sample / "rmlst").resolve()
    return RmlstConfig(
        allele_target=refs.rmlst_alleles.resolve(),
        seqdef=refs.rmlst_seqdef.resolve(),
        reads=(pair.r1.resolve(), pair.r2.resolve()),
        sample_id=pair.sample,
        output=output,
        tempdir=(Path(tempdir) if tempdir else output / "work").resolve(),
        threads=threads,
        memory_limit_gb=memory_gb,
        keep_temp=keep_temp,
        verbose=verbose,
        kallisto=tools.kallisto,
        kallisto_index=refs.rmlst_index.resolve() if refs.rmlst_index else None,
        mash=tools.mash,
        minimap2=tools.minimap2,
        bwa_mem2=tools.bwa_mem2,
        samtools=tools.samtools,
        cliquesnv=tools.cliquesnv,
        cliquesnv_time_limit_seconds=settings.hpc.cliquesnv_time_limit_seconds,
        viral_consensus=tools.viral_consensus,
        blastn=tools.blastn,
        makeblastdb=tools.makeblastdb,
        # The shared full-reference BLAST cache lives at the run root, so the
        # samples of a run share one database instead of each building its own.
        blastdb_cache_dir=Path(out_root).resolve(),
        min_unambiguous_loci=int(cfg.get("rmlst.min_unambiguous_loci", 3)),
        require_accepted_consensus=bool(cfg.get("rmlst.require_accepted_consensus", False)),
    )


def marker16s_argv(settings: Settings, pair: ReadPair, out_root: Path, *,
                   threads: int, verbose: bool) -> list[str]:
    refs, cfg = settings.references, settings.cfg
    if refs.s16_reference is None:
        raise ValueError("references.s16_reference must be set")
    tools = settings.tools
    argv = [
        "--sample-id", pair.sample,
        "--reads", str(pair.r1), str(pair.r2),
        "--reference", str(refs.s16_reference),
        "--kallisto", tools.kallisto,
        "--minimap2", tools.minimap2, "--samtools", tools.samtools,
        "--viral-consensus", tools.viral_consensus,
        "--blastn", tools.blastn, "--makeblastdb", tools.makeblastdb,
        "-t", str(threads),
        "-o", str(Path(out_root) / pair.sample / "16s"),
        "--genus-min-frac", str(cfg.get("marker16s.genus_min_frac", 0.01)),
        "--min-genus-reads", str(cfg.get("marker16s.min_genus_reads", 10.0)),
        "--species-dominance", str(cfg.get("marker16s.species_dominance", 0.70)),
        "--species-margin", str(cfg.get("marker16s.species_margin", 0.20)),
    ]
    if refs.s16_index is not None:
        argv += ["--kallisto-index", str(refs.s16_index)]
    # Explicit either way: an unset table must mean "uncorrected on purpose",
    # not whatever the arm's own default happens to point at.
    if refs.s16_copy_numbers is not None:
        argv += ["--copy-number-table", str(refs.s16_copy_numbers)]
    else:
        argv += ["--no-copy-number-table"]
    if not bool(cfg.get("marker16s.alignment_metrics", True)):
        argv.append("--no-alignment-metrics")
    if verbose:
        argv.append("-v")
    return argv


class ArmFailed(RuntimeError):
    """One or more arms failed, after every requested arm was attempted."""


def tidy_arm(sample_dir: Path, arm: str, keep_temp: bool) -> None:
    """Drop a finished arm's intermediates (see ``clean``), as the rMLST arm
    does for itself. Only after the arm succeeded: a failed arm keeps its
    files for diagnosis. Never fatal - the results are already written."""
    if keep_temp:
        return
    try:
        freed = remove_intermediates(sample_dir, (arm,))
    except OSError as error:
        LOGGER.warning("could not remove %s intermediates in %s: %s", arm, sample_dir, error)
        return
    if freed:
        print(f"    removed {arm} intermediates ({human(freed)}); keep them with --keep-temp",
              flush=True)


def run_sample(settings: Settings, pair: ReadPair, out_root: Path, *,
               arms: tuple[str, ...] = ARMS, threads: int = 4,
               memory_gb: float | None = None, tempdir: Path | None = None,
               verbose: bool = False, keep_temp: bool = False) -> dict[str, Path]:
    """Run the requested arms; return the per-arm result files.

    The arms are independent - different reference, different algorithm,
    no shared state - and reconciliation is built to report a sample from the
    16S arm alone (``16s_fallback``). So a failure in one arm must not stop the
    others: losing rMLST costs resolution, losing AMR costs the resistance
    report, while losing both ID arms costs the identification. Every
    requested arm is attempted, then the failure is raised, so the job still
    exits non-zero and the sample is resubmitted rather than quietly half-done.
    """
    out_root = Path(out_root)
    sample_dir = out_root / pair.sample
    sample_dir.mkdir(parents=True, exist_ok=True)
    produced: dict[str, Path] = {}
    failures: list[tuple[str, Exception]] = []
    print(f"=== {pair.sample}: {code_versions()} ===", flush=True)

    if "rmlst" in arms:
        print(f"=== rMLST arm: {pair.sample} ===", flush=True)
        try:
            config = rmlst_config(settings, pair, out_root, threads=threads,
                                  memory_gb=memory_gb, tempdir=tempdir,
                                  verbose=verbose, keep_temp=keep_temp)
            run_pipeline(config)
            produced["rmlst"] = config.output / "species_calls.csv"
        except Exception as error:  # noqa: BLE001 - the 16S arm still gets its turn
            failures.append(("rmlst", error))
            LOGGER.exception("rmlst arm failed for %s", pair.sample)
            print(f"!!! rMLST arm FAILED for {pair.sample}: "
                  f"{type(error).__name__}: {error}", flush=True)
            print("    continuing to the 16S arm; this sample can still be "
                  "reported from 16S alone.", flush=True)

    if "16s" in arms:
        print(f"=== 16S arm: {pair.sample} ===", flush=True)
        try:
            argv = marker16s_argv(settings, pair, out_root, threads=threads,
                                  verbose=verbose)
            code = marker16s_main(argv)
            if code:
                raise RuntimeError(f"16S arm exited {code}")
            produced["16s"] = sample_dir / "16s" / "marker16s_calls.csv"
            tidy_arm(sample_dir, "16s", keep_temp)
        except Exception as error:  # noqa: BLE001 - reported with the other arm
            failures.append(("16s", error))
            LOGGER.exception("16S arm failed for %s", pair.sample)
            print(f"!!! 16S arm FAILED for {pair.sample}: "
                  f"{type(error).__name__}: {error}", flush=True)

    if "amr" in arms:
        print(f"=== AMR arm: {pair.sample} ===", flush=True)
        try:
            from bridge.amr.kma import run_arm
            from bridge.amr.settings import from_config
            amr_settings = from_config(settings.cfg, settings.path)
            produced["amr"] = run_arm(amr_settings, pair.sample, pair.r1, pair.r2,
                                      sample_dir / "amr", threads=threads)
            tidy_arm(sample_dir, "amr", keep_temp)
        except Exception as error:  # noqa: BLE001 - reported with the other arms
            failures.append(("amr", error))
            LOGGER.exception("AMR arm failed for %s", pair.sample)
            print(f"!!! AMR arm FAILED for {pair.sample}: "
                  f"{type(error).__name__}: {error}", flush=True)

    print(f"=== done: {pair.sample} ===")
    for arm, path in produced.items():
        print(f"  {arm:<6}: {path}")
    if failures:
        print(f"  FAILED: {', '.join(arm for arm, _ in failures)}")
        raise ArmFailed(
            f"{pair.sample}: " + "; ".join(f"{arm} arm: {error}" for arm, error in failures))
    return produced
