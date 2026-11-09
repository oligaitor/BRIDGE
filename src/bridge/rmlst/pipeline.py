from __future__ import annotations

import json
import logging
from math import isfinite
import shutil
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from time import perf_counter

import pandas as pd

from bridge import __version__

from .clustering import cluster_abundance, cluster_candidates
from .consensus import run_cluster_consensus
from .dependencies import check_config_dependencies
from .errors import InsufficientMaterial
from .fasta import read_fasta
from .kallisto import derive_candidates, run_candidate_pseudobam, run_quant, sha256
from .models import Config, ConsensusCall
from .reads import assign_and_extract_reads
from .reporting import write_direct_species_alleles
from .seqdef import add_rst_matches, allele_species
from .species import FINAL_COLUMNS, final_species_calls, merge_consensus_calls
from bridge.tools import (
    close_logging,
    companion_executable,
    configure_logging,
    executable,
    reset_runtime_accounting,
    runtime_accounting,
    stage,
)


# Deleted on a normal run; retained by --keeptmp. consensus_groups.csv is kept
# deliberately: it is the per-locus detail behind the per-species coverage
# aggregates in species_calls.csv, and the only way to audit them without
# re-running the sample.
INTERMEDIATE_OUTPUT_NAMES = (
    "kallisto_candidates.csv",
    "runtime_clusters.csv",
    "runtime_cluster_abundance.csv",
    "runtime_cluster_read_counts.csv",
    "runtime_cluster_consensus.csv",
    "species_allele_fastas.csv",
    "em_fallback.csv",
    "cluster_consensus",
)


def run_pipeline(config: Config) -> dict[str, Path]:
    started = perf_counter()
    reset_runtime_accounting()
    config.output.mkdir(parents=True, exist_ok=True)
    log_path = config.output / "rmlst.log"
    log_handlers = configure_logging(log_path, config.verbose)
    logger = logging.getLogger("bridge.rmlst")
    tempdir_active = False
    try:
        _validate(config)
        with stage("dependency_check"):
            dependencies = check_config_dependencies(config)
            for dependency in dependencies:
                logger.info(
                    "dependency_status=available name=%s requirement=%s resolved=%s",
                    dependency.name, dependency.requirement, dependency.resolved,
                )
        config.tempdir.mkdir(parents=True, exist_ok=True)
        tempdir_active = True
        with stage("reference"):
            fasta_records = read_fasta(config.allele_target)
            records = {record.identifier: record for record in fasta_records}
            species = allele_species(
                config.seqdef, records, config.distinguish_subspecies
            )
        with stage("kallisto"):
            quant, index, index_kmer_size = run_quant(
                config, config.allele_target
            )
            calls = derive_candidates(config, quant)
            calls_path = config.output / "kallisto_candidates.csv"
            calls.to_csv(calls_path, index=False)
        with stage("candidate_pseudobam"):
            pseudobam = run_candidate_pseudobam(
                config, calls, records, index_kmer_size
            )
        with stage("runtime_clustering"):
            assignments = cluster_candidates(calls, records, config)
            assignments_path = config.output / "runtime_clusters.csv"
            assignments.to_csv(assignments_path, index=False)
            clusters = cluster_abundance(calls, assignments)
            cluster_path = config.output / "runtime_cluster_abundance.csv"
            clusters.to_csv(cluster_path, index=False)
        with stage("cluster_read_assignment"):
            read_paths, read_stats = assign_and_extract_reads(pseudobam, assignments, config)
            read_stats_path = config.output / "runtime_cluster_read_counts.csv"
            read_stats.to_csv(read_stats_path, index=False)
        with stage("cluster_consensus"):
            consensus_calls = run_cluster_consensus(
                clusters, assignments, read_paths, read_stats, records, species,
                config.allele_target, config,
            )
            consensus_path = config.output / "runtime_cluster_consensus.csv"
            _consensus_frame(
                consensus_calls, config.output, persist_artefacts=config.keep_temp,
            ).to_csv(consensus_path, index=False)
        with stage("consensus_merging"):
            groups = merge_consensus_calls(consensus_calls)
            groups_path = config.output / "consensus_groups.csv"
            groups.to_csv(groups_path, index=False)
        with stage("species_calling"):
            final, em = final_species_calls(
                config.sample_id, groups, config.min_profile_loci,
                min_unambiguous_loci=config.min_unambiguous_loci,
                require_accepted_consensus=config.require_accepted_consensus,
            )
            final = add_rst_matches(
                final, groups, config.seqdef, limit=10,
                distinguish_subspecies=config.distinguish_subspecies,
            )
            em_path = config.output / "em_fallback.csv"
            em.insert(0, "sample_id", config.sample_id)
            em.to_csv(em_path, index=False)
        with stage("species_allele_fastas"):
            final, allele_fasta_manifest = write_direct_species_alleles(
                final, groups, config.output
            )
            final_path = config.output / "species_calls.csv"
            allele_fasta_manifest_path = config.output / "species_allele_fastas.csv"
            final.to_csv(final_path, index=False)
            allele_fasta_manifest.to_csv(allele_fasta_manifest_path, index=False)
        with stage("called_species_plot"):
            # Imported lazily: matplotlib and seaborn cost seconds of startup,
            # which samples that fail before this stage should not pay.
            from .plotting import write_called_species_pi_pdf

            called_species_plot_path = config.output / "called_species_pi.pdf"
            write_called_species_pi_pdf(
                final, called_species_plot_path, config.sample_id,
            )
        with stage("manifest"):
            manifest_path = config.output / "rmlst_manifest.json"
            manifest = {
                "schema": "bridge.rmlst-run",
                "schema_version": 11,
                "package": "bridge",
                "package_version": __version__,
                "created_utc": datetime.now(timezone.utc).isoformat(),
                "sample_id": config.sample_id,
                "reference": str(config.allele_target),
                "reference_sha256": sha256(config.allele_target),
                "seqdef": str(config.seqdef),
                "reads": [str(path) for path in config.reads],
                "kallisto_index": str(index),
                "kallisto_index_source": (
                    "supplied" if config.kallisto_index is not None else "automatic"
                ),
                "kallisto_index_kmer_size": index_kmer_size,
                "candidate_rule": (
                    f"per locus: positive in at least {config.support_threshold} of "
                    f"{config.bootstraps} bootstraps; if empty, relax to that "
                    "locus's highest observed support"
                ),
                "runtime_clusters": int(len(clusters)),
                "consensus_groups": int(len(groups)),
                "reported_species": final.loc[final["reported"], "species"].tolist(),
                "tools": {
                    name: executable(value) for name, value in ({
                        "kallisto": config.kallisto, "mash": config.mash,
                        "samtools": config.samtools,
                        "viral_consensus": config.viral_consensus,
                        "blastn": config.blastn,
                    } | (
                        {"minimap2": config.minimap2}
                        if config.short_read_aligner == "minimap2"
                        else {"bwa_mem2": config.bwa_mem2}
                    ) | (
                        {} if config.viralconsensus_only
                        else {"cliquesnv": config.cliquesnv}
                    )).items()
                } | {
                    "makeblastdb": companion_executable(
                        config.makeblastdb, config.blastn
                    )
                },
                "parameters": _json_config(config),
                "runtime_seconds": perf_counter() - started,
                # Per-stage and per-tool wall clock. The cluster_consensus stage
                # alone covers the aligner, CliqueSNV, ViralConsensus, BLAST and
                # a lot of in-process Python, so a single stage total cannot say
                # where a sample's time went; these can.
                "runtime_breakdown": runtime_accounting(),
            }
            manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        outputs = {
            "calls": calls_path,
            "clusters": assignments_path,
            "cluster_abundance": cluster_path,
            "read_counts": read_stats_path,
            "consensus": consensus_path,
            "consensus_groups": groups_path,
            "species_calls": final_path,
            "called_species_plot": called_species_plot_path,
            "species_allele_fastas": config.output / "species_allele_fastas",
            "species_allele_fastas_manifest": allele_fasta_manifest_path,
            "em": em_path,
            "manifest": manifest_path,
            "log": log_path,
        }
        if config.keep_temp:
            retained_outputs = outputs
        else:
            _remove_intermediate_outputs(config.output)
            logger.info(
                "intermediate_outputs_status=removed output=%s", config.output
            )
            retained_outputs = {
                "species_calls": final_path,
                "called_species_plot": called_species_plot_path,
                "species_allele_fastas": config.output / "species_allele_fastas",
                "manifest": manifest_path,
                "log": log_path,
            }
        logger.info(
            "pipeline_status=completed sample=%s runtime_seconds=%.3f reported_species=%d",
            config.sample_id, perf_counter() - started, int(final["reported"].sum()),
            extra={"console": True},
        )
        return retained_outputs
    except InsufficientMaterial as error:
        logger.warning(
            "pipeline_status=insufficient_material sample=%s stage=%s reason=%s "
            "runtime_seconds=%.3f",
            config.sample_id, error.stage or "unknown", error,
            perf_counter() - started, extra={"console": True},
        )
        return _write_insufficient_material(config, error, started)
    except Exception as error:
        logger.exception(
            "pipeline_status=failed sample=%s error=%s",
            config.sample_id, error,
        )
        raise
    finally:
        try:
            if tempdir_active and config.keep_temp:
                logger.info("temporary_work_status=retained path=%s", config.tempdir)
            elif tempdir_active and config.tempdir.is_dir():
                shutil.rmtree(config.tempdir)
                logger.info("temporary_work_status=removed path=%s", config.tempdir)
        finally:
            close_logging(log_handlers)


def _write_insufficient_material(
    config: Config, error: InsufficientMaterial, started: float
) -> dict[str, Path]:
    """Record a genuine negative as a result instead of crashing.

    Writes the same two deliverables a successful run leaves behind - an empty
    ``species_calls.csv`` and a manifest - so the sample reads as finished. The
    batch driver uses ``species_calls.csv`` as its completion sentinel, so
    without this every negative control is resubmitted, and repays the full
    quant, on every re-run.
    """
    final_path = config.output / "species_calls.csv"
    pd.DataFrame(columns=FINAL_COLUMNS).to_csv(final_path, index=False)
    manifest_path = config.output / "rmlst_manifest.json"
    manifest_path.write_text(
        json.dumps({
            "schema": "bridge.rmlst-run",
            "schema_version": 11,
            "package": "bridge",
            "package_version": __version__,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "sample_id": config.sample_id,
            "outcome": "insufficient_material",
            "outcome_stage": error.stage,
            "outcome_reason": str(error),
            "reference": str(config.allele_target),
            "seqdef": str(config.seqdef),
            "reads": [str(path) for path in config.reads],
            "reported_species": [],
            "parameters": _json_config(config),
            "runtime_seconds": perf_counter() - started,
            "runtime_breakdown": runtime_accounting(),
        }, indent=2) + "\n",
        encoding="utf-8",
    )
    return {"species_calls": final_path, "manifest": manifest_path,
            "log": config.output / "rmlst.log"}


def _remove_intermediate_outputs(output: Path) -> None:
    """Remove pipeline-owned diagnostics while leaving final deliverables intact."""
    for name in INTERMEDIATE_OUTPUT_NAMES:
        path = output / name
        if path.is_symlink() or path.is_file():
            path.unlink()
        elif path.is_dir():
            shutil.rmtree(path)


def _consensus_frame(
    calls: list[ConsensusCall], output: Path, persist_artefacts: bool = True
) -> pd.DataFrame:
    """Flatten consensus calls, optionally persisting their FASTA/BLAST files.

    ``persist_artefacts`` is false on a normal run: ``cluster_consensus`` is in
    ``INTERMEDIATE_OUTPUT_NAMES``, so without ``--keeptmp`` any copies would be
    deleted a few stages later anyway.
    """
    persistent = output / "cluster_consensus"
    rows = []
    for call in calls:
        data = asdict(call)
        data["source_members"] = ";".join(call.source_members)
        data["compatible_alleles"] = ";".join(call.compatible_alleles)
        data["compatible_species"] = ";".join(call.compatible_species)
        for field in ("consensus_path", "blast_path"):
            source = data[field]
            if not source:
                data[field] = ""
                continue
            source = Path(source)
            if not source.is_file():
                data[field] = ""
                continue
            if not persist_artefacts:
                # The path is recorded as-is; it lives under the temporary
                # directory and is removed with it.
                data[field] = str(source)
                continue
            destination = persistent / call.locus_id / call.runtime_cluster_id / source.name
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
            data[field] = str(destination.resolve())
        rows.append(data)
    return pd.DataFrame(rows)


def _json_config(config: Config) -> dict[str, object]:
    ignored = {"allele_target", "seqdef", "reads", "sample_id", "output", "tempdir"}
    return {
        key: str(value) if isinstance(value, Path) else value
        for key, value in asdict(config).items() if key not in ignored
    }


def _validate(config: Config) -> None:
    if config.short_read_aligner not in {"bwa-mem2", "minimap2"}:
        raise ValueError("short-read aligner must be bwa-mem2 or minimap2")
    if len(config.reads) != 2:
        raise ValueError("Exactly two paired-end FASTQ files are currently required")
    for path in (config.allele_target, config.seqdef, *config.reads):
        if not path.is_file():
            raise FileNotFoundError(path)
    if config.kallisto_index is not None and not config.kallisto_index.is_file():
        raise FileNotFoundError(config.kallisto_index)
    if config.threads < 1 or config.bootstraps < 1:
        raise ValueError("threads and bootstraps must be positive")
    if config.memory_limit_gb is not None and (
        not isfinite(config.memory_limit_gb) or config.memory_limit_gb <= 0
    ):
        raise ValueError("memory limit must be positive")
    if (
        config.cliquesnv_time_limit_seconds is not None
        and config.cliquesnv_time_limit_seconds <= 0
    ):
        raise ValueError("CliqueSNV runtime limit must be positive")
    if not 0 <= config.support_threshold <= config.bootstraps:
        raise ValueError("support threshold must be between 0 and bootstraps")
    if not 0 <= config.subconsensus_min_proportion < 1:
        raise ValueError("subconsensus minimum proportion must be in [0, 1)")
    if config.subconsensus_min_fragments < 1:
        raise ValueError("subconsensus minimum fragments must be positive")
    if config.subconsensus_min_depth <= 0:
        raise ValueError("subconsensus minimum depth must be positive")
    if not 0 < config.subconsensus_min_breadth <= 1:
        raise ValueError("subconsensus minimum breadth must be in (0, 1]")
    if config.subconsensus_min_end_distance < 0:
        raise ValueError("subconsensus minimum end distance cannot be negative")
    if config.minor_allele_min_depth < 1:
        raise ValueError("minor allele minimum depth must be positive")
    if config.read_assignment not in {"compatible", "exclusive"}:
        raise ValueError("read assignment must be compatible or exclusive")
