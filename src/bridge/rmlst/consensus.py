from __future__ import annotations

from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
import hashlib
import logging
from math import ceil
import os
from pathlib import Path

import pandas as pd

from .fasta import FastaRecord, write_fasta
from .kallisto import sha256
from .haplotypes import (
    HaplotypeCandidate,
    alignment_evidence_from_sam,
    evaluate_haplotype,
    maximum_minor_allele_depth,
    read_cliquesnv_haplotypes,
)
from .models import Config, ConsensusCall
from .seqdef import remove_mixed_placeholders
from bridge.tools import companion_executable, executable, run


GIB = 1024 ** 3
MIB = 1024 ** 2
CLIQUESNV_MEMORY_FRACTION = 0.80
MINIMUM_JAVA_HEAP_MB = 64
MAX_TERMINAL_SEQUENCE_DIFFERENCE = 15
MINIMAP2_SHORT_READ_PARAMETERS = (
    "-ax", "sr", "-k", "17", "-w", "7", "-m", "20", "-s", "30",
)
BWA_MEM2_MIN_MAPPING_QUALITY = 20
MINIMAP2_MIN_MAPPING_QUALITY = 10


def detected_memory_limit_bytes() -> tuple[int, str]:
    """Return the smallest usable cgroup or physical-memory limit."""
    candidates: list[tuple[int, str]] = []
    for path in (
        Path("/sys/fs/cgroup/memory.max"),
        Path("/sys/fs/cgroup/memory/memory.limit_in_bytes"),
    ):
        try:
            value = path.read_text(encoding="utf-8").strip()
            if value != "max":
                limit = int(value)
                if limit > 0:
                    candidates.append((limit, f"cgroup:{path}"))
        except (FileNotFoundError, OSError, ValueError):
            pass

    try:
        physical = os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")
        if physical > 0:
            candidates.append((physical, "physical_memory"))
    except (AttributeError, OSError, TypeError, ValueError):
        pass

    if not candidates:
        raise RuntimeError(
            "Unable to detect available memory; supply --memory-limit-gb"
        )
    return min(candidates, key=lambda value: value[0])


def cliquesnv_heap_megabytes(
    workers: int,
    memory_limit_gb: float | None = None,
) -> tuple[int, float, str]:
    """Allocate 80% of the overall budget across concurrent Java tasks."""
    if workers < 1:
        raise ValueError("CliqueSNV worker count must be positive")
    if memory_limit_gb is not None:
        total_bytes = int(memory_limit_gb * GIB)
        source = "configured"
    else:
        total_bytes, source = detected_memory_limit_bytes()
    heap_mb = int(total_bytes * CLIQUESNV_MEMORY_FRACTION) // workers // MIB
    if heap_mb < MINIMUM_JAVA_HEAP_MB:
        raise ValueError(
            "Memory budget is too small for the configured CliqueSNV "
            f"parallelism: {workers} workers would receive {heap_mb} MB each"
        )
    return heap_mb, total_bytes / GIB, source


def cliquesnv_launcher(path: str, heap_mb: int) -> tuple[str, ...]:
    """Build a launcher for either the Conda wrapper or a direct JAR."""
    initial_heap_mb = min(512, heap_mb)
    java_options = (f"-Xms{initial_heap_mb}m", f"-Xmx{heap_mb}m")
    if Path(path).suffix.lower() == ".jar":
        return (executable("java"), *java_options, "-jar", path)
    return (path, *java_options)


def compatible_alleles(
    consensus: str,
    alleles: dict[str, str],
    max_terminal_difference: int = MAX_TERMINAL_SEQUENCE_DIFFERENCE,
) -> tuple[str, ...]:
    """Return alleles agreeing across a bounded, ungapped terminal alignment.

    An allele may contain additional sequence at either end. A shorter allele
    is also compatible when only non-callable consensus bases are discarded
    from the ends. Internal gaps are never permitted, and every overlapping
    callable consensus base must agree.
    """
    if max_terminal_difference < 0:
        raise ValueError("max_terminal_difference must be non-negative")
    query = consensus.upper()
    callable_positions = [
        (index, base) for index, base in enumerate(query) if base in "ACGT"
    ]
    matches = []
    for allele, sequence in alleles.items():
        sequence = sequence.upper()
        difference = len(sequence) - len(query)
        if abs(difference) > max_terminal_difference:
            continue

        compatible = False
        if difference >= 0:
            for left_extension in range(difference + 1):
                if all(
                    sequence[left_extension + index] == base
                    for index, base in callable_positions
                ):
                    compatible = True
                    break
        else:
            trimming = -difference
            for left_trimming in range(trimming + 1):
                right_trimming = trimming - left_trimming
                right_boundary = len(query) - right_trimming
                discarded = query[:left_trimming]
                if right_trimming:
                    discarded += query[right_boundary:]
                if any(base in "ACGT" for base in discarded):
                    continue
                if all(
                    sequence[index - left_trimming] == base
                    for index, base in callable_positions
                    if left_trimming <= index < right_boundary
                ):
                    compatible = True
                    break

        if compatible:
            matches.append(allele)
    return tuple(sorted(matches))


def pad_terminal_consensus(
    consensus: str,
    alleles: dict[str, str],
    max_padding: int = MAX_TERMINAL_SEQUENCE_DIFFERENCE,
) -> str:
    """Pad a short terminal truncation with N when it fits a known allele."""
    query = consensus.upper()
    best = query
    best_score = (0, 0, 0)
    lengths = sorted({len(sequence) for sequence in alleles.values()})
    # Fast path for the common case: a consensus already as long as the
    # longest allele has nothing to pad to.
    if not lengths or lengths[-1] <= len(query):
        return query
    for target_length in lengths:
        padding = target_length - len(query)
        if padding <= 0 or padding > max_padding:
            continue
        for left_padding in range(padding + 1):
            padded = (
                "N" * left_padding
                + query
                + "N" * (padding - left_padding)
            )
            matches = compatible_alleles(
                padded, alleles, max_terminal_difference=0
            )
            score = (len(matches), -padding, -left_padding)
            if matches and score > best_score:
                best = padded
                best_score = score
    return best


def mask_weak_variants(sequence: str, reasons: tuple[str, ...]) -> str:
    """Replace unsupported variant calls with N instead of losing the locus."""
    masked = list(sequence.upper())
    prefix = "weak_variant_evidence_at_"
    for reason in reasons:
        if not reason.startswith(prefix):
            continue
        position = int(reason[len(prefix):]) - 1
        if 0 <= position < len(masked):
            masked[position] = "N"
    return "".join(masked)


def composite_consensus_possible(
    consensus: str,
    member_sequences: dict[str, str] | list[str] | tuple[str, ...],
) -> bool:
    """Return whether cluster members can compose an absent consensus.

    This deliberately conservative, alignment-free check is restricted to a
    fully callable consensus and same-length cluster members. A consensus is a
    possible composite when it is not itself a member sequence and, at every
    position, at least one member supplies the consensus base.
    """
    query = consensus.upper()
    if not query or any(base not in "ACGT" for base in query):
        return False
    values = (
        member_sequences.values()
        if isinstance(member_sequences, dict)
        else member_sequences
    )
    members = [
        str(sequence).upper()
        for sequence in values
        if len(str(sequence)) == len(query)
    ]
    if not members or query in members:
        return False
    return all(
        any(member[position] == base for member in members)
        for position, base in enumerate(query)
    )


def alignment_min_mapping_quality(config: Config) -> int:
    return (
        MINIMAP2_MIN_MAPPING_QUALITY
        if config.short_read_aligner == "minimap2"
        else BWA_MEM2_MIN_MAPPING_QUALITY
    )


def _align_cluster_reads(
    config: Config,
    aligner: str,
    reference: Path,
    read1: Path,
    read2: Path,
    sam: Path,
    threads: int,
) -> None:
    """Write one cluster SAM with the configured short-read backend."""
    if config.short_read_aligner == "minimap2":
        run([
            aligner, *MINIMAP2_SHORT_READ_PARAMETERS,
            "-t", str(threads), str(reference), str(read1), str(read2),
        ], "minimap2", stdout_path=sam)
        return
    run([aligner, "index", str(reference)], "bwa-mem2")
    run([
        aligner, "mem", "-t", str(threads),
        str(reference), str(read1), str(read2),
    ], "bwa-mem2", stdout_path=sam)


def run_cluster_consensus(
    clusters: pd.DataFrame,
    assignments: pd.DataFrame,
    read_paths: dict[str, tuple[Path, Path]],
    read_stats: pd.DataFrame,
    records: dict[str, FastaRecord],
    species_by_allele: dict[str, set[str]],
    reference: Path,
    config: Config,
) -> list[ConsensusCall]:
    aligner = executable(
        config.minimap2
        if config.short_read_aligner == "minimap2"
        else config.bwa_mem2
    )
    viral_consensus = executable(config.viral_consensus)
    stats = read_stats.set_index("runtime_cluster_id").to_dict("index") if not read_stats.empty else {}
    members = assignments.groupby("runtime_cluster_id")["Member_Template"].apply(
        lambda values: tuple(sorted(set(map(str, values))))
    ).to_dict()
    sequences_by_locus: dict[str, dict[str, str]] = {}
    for allele, record in records.items():
        sequences_by_locus.setdefault(allele.split("_", 1)[0], {})[allele] = record.sequence
    root = config.tempdir / "consensus"
    rows = list(clusters.itertuples(index=False))
    per_cluster_threads = min(2, config.threads)
    workers = max(
        1, min(len(rows), max(1, config.threads // per_cluster_threads))
    )
    logging.getLogger("bridge.rmlst").info(
        "consensus_parallelism_status=configured workers=%d "
        "threads_per_cluster=%d total_thread_budget=%d aligner=%s "
        "minimum_mapping_quality=%d",
        workers, per_cluster_threads, config.threads,
        config.short_read_aligner, alignment_min_mapping_quality(config),
    )
    if config.viralconsensus_only:
        cliquesnv: tuple[str, ...] = ()
        logging.getLogger("bridge.rmlst").info(
            "cliquesnv_status=disabled reason=viralconsensus_only"
        )
    else:
        cliquesnv_path = executable(config.cliquesnv)
        heap_mb, total_memory_gb, memory_source = cliquesnv_heap_megabytes(
            workers, config.memory_limit_gb
        )
        cliquesnv = cliquesnv_launcher(cliquesnv_path, heap_mb)
        logging.getLogger("bridge.rmlst").info(
            "cliquesnv_memory_status=configured source=%s total_memory_gb=%.3f "
            "reserved_fraction=%.2f workers=%d heap_per_task_mb=%d "
            "maximum_aggregate_heap_gb=%.3f",
            memory_source, total_memory_gb, 1.0 - CLIQUESNV_MEMORY_FRACTION,
            workers, heap_mb, workers * heap_mb / 1024,
        )
        logging.getLogger("bridge.rmlst").info(
            "cliquesnv_time_limit_status=configured seconds=%s scope=per_invocation "
            "implementation=native_tl",
            (
                config.cliquesnv_time_limit_seconds
                if config.cliquesnv_time_limit_seconds is not None
                else "cliquesnv_default"
            ),
        )

    def generate(row) -> list[dict[str, object]]:
        return _generate_cluster_consensus(
            row, stats, members, read_paths, records, sequences_by_locus,
            root, config, aligner, cliquesnv, viral_consensus,
            per_cluster_threads,
        )

    with ThreadPoolExecutor(max_workers=workers) as executor:
        drafts = [
            draft
            for cluster_drafts in executor.map(generate, rows)
            for draft in cluster_drafts
        ]

    pending = [draft for draft in drafts if draft["status"] == "pending_full_reference"]
    if pending:
        try:
            _batch_full_reference_search(
                pending, reference, sequences_by_locus, root, config
            )
        except Exception:
            logging.getLogger("bridge.rmlst").exception(
                "allele_search_status=batch_failed queries=%d", len(pending)
            )
            for draft in pending:
                draft["status"] = "alignment_failed"
                draft["allele_search_scope"] = "full_reference_batch_failed"

    results: list[ConsensusCall] = []
    for draft in drafts:
        compatible = tuple(draft["compatible_alleles"])
        cluster_members = tuple(draft["source_members"])
        if compatible and not draft["composite_from_members"]:
            species = remove_mixed_placeholders({
                value for allele in compatible
                for value in species_by_allele.get(allele, set())
            })
        else:
            species = remove_mixed_placeholders({
                value for allele in cluster_members
                for value in species_by_allele.get(allele, set())
            })
        draft["compatible_species"] = tuple(sorted(species))
        draft.pop("directory")
        results.append(ConsensusCall(**draft))
    return results


def _viralconsensus_candidate(
    cluster_id: str,
    sam: Path,
    representative_path: Path,
    directory: Path,
    viral_consensus: str,
    config: Config,
) -> HaplotypeCandidate:
    viral_path = directory / "viralconsensus.fasta"
    run([
        viral_consensus, "-i", str(sam), "-r", str(representative_path),
        "-o", str(viral_path), "-q", "20", "-d",
        str(max(1, ceil(config.subconsensus_min_depth))),
    ], "viralconsensus")
    return HaplotypeCandidate(
        f"{cluster_id}_viralconsensus",
        _read_single_consensus(viral_path), 1.0, viral_path,
    )


def _generate_cluster_consensus(
    row,
    stats: dict[str, dict[str, object]],
    members: dict[str, tuple[str, ...]],
    read_paths: dict[str, tuple[Path, Path]],
    records: dict[str, FastaRecord],
    sequences_by_locus: dict[str, dict[str, str]],
    root: Path,
    config: Config,
    aligner: str,
    cliquesnv: tuple[str, ...],
    viral_consensus: str,
    per_cluster_threads: int,
) -> list[dict[str, object]]:
    cluster_id = str(row.runtime_cluster_id)
    locus = str(row.locus_id)
    representative = str(row.representative)
    cluster_stats = stats.get(cluster_id, {})
    assigned = int(cluster_stats.get("assigned_fragments", 0))
    exclusive = int(cluster_stats.get("exclusive_fragments", 0))
    ambiguous = int(cluster_stats.get("ambiguous_fragments", 0))
    directory = root / locus / cluster_id
    directory.mkdir(parents=True, exist_ok=True)
    representative_path = write_fasta(
        [FastaRecord(cluster_id, cluster_id, records[representative].sequence)],
        directory / "representative.fasta",
    )
    if assigned < config.min_fragments or cluster_id not in read_paths:
        return [_empty_cluster_draft(
            row, cluster_id, locus, representative, members[cluster_id], assigned,
            exclusive, ambiguous, len(records[representative].sequence), directory,
            config, "insufficient_reads",
        )]

    try:
        sam = directory / "alignment.sam"
        read1, read2 = read_paths[cluster_id]
        _align_cluster_reads(
            config, aligner, representative_path, read1, read2, sam,
            per_cluster_threads,
        )
        reference_sequence = records[representative].sequence.upper()
        evidence = alignment_evidence_from_sam(
            sam, len(reference_sequence),
            min_mapping_quality=alignment_min_mapping_quality(config),
            min_end_distance=config.subconsensus_min_end_distance,
        )
        maximum_minor_depth = maximum_minor_allele_depth(evidence)
        if (
            not config.viralconsensus_only
            and maximum_minor_depth >= config.minor_allele_min_depth
        ):
            consensus_method = "cliquesnv"
            logging.getLogger("bridge.rmlst").info(
                "consensus_route=cliquesnv cluster=%s maximum_minor_allele_depth=%d "
                "threshold=%d strand_filter=false frequency_filter=false",
                cluster_id, maximum_minor_depth, config.minor_allele_min_depth,
            )
            haplotype_dir = directory / "cliquesnv"
            haplotype_dir.mkdir(parents=True, exist_ok=True)
            command = [
                *cliquesnv, "-m", "snv-illumina", "-in", str(sam),
                "-outDir", str(haplotype_dir), "-threads", str(per_cluster_threads),
                "-fdf", "extended9", "-tf", str(config.subconsensus_min_proportion),
                "-cm", "fast", "-el", "200", "-ignoreDeletion"
            ]
            if config.cliquesnv_time_limit_seconds is not None:
                command.extend(["-tl", str(config.cliquesnv_time_limit_seconds)])
            try:
                run(command, "cliquesnv")
                haplotype_paths = sorted(
                    set(haplotype_dir.rglob("*.fasta"))
                    | set(haplotype_dir.rglob("*.fa"))
                    | set(haplotype_dir.rglob("*.fna"))
                )
                if not haplotype_paths:
                    raise ValueError(
                        "CliqueSNV did not produce a haplotype FASTA in "
                        f"{haplotype_dir}"
                    )
                candidates = read_cliquesnv_haplotypes(haplotype_paths)
            except Exception:
                logging.getLogger("bridge.rmlst").exception(
                    "consensus_fallback=viralconsensus cluster=%s "
                    "reason=cliquesnv_failed",
                    cluster_id,
                )
                consensus_method = "viralconsensus"
                candidates = [_viralconsensus_candidate(
                    cluster_id, sam, representative_path, directory,
                    viral_consensus, config,
                )]
            if (
                consensus_method == "cliquesnv"
                and len(candidates) == 1
                and abs(candidates[0].proportion - 1.0) <= 1e-9
            ):
                logging.getLogger("bridge.rmlst").info(
                    "consensus_fallback=viralconsensus cluster=%s "
                    "reason=cliquesnv_single_full_proportion_haplotype",
                    cluster_id,
                )
                consensus_method = "viralconsensus"
                candidates = [_viralconsensus_candidate(
                    cluster_id, sam, representative_path, directory,
                    viral_consensus, config,
                )]
        else:
            consensus_method = "viralconsensus"
            logging.getLogger("bridge.rmlst").info(
                "consensus_route=viralconsensus cluster=%s "
                "maximum_minor_allele_depth=%d threshold=%d "
                "strand_filter=false frequency_filter=false forced=%s",
                cluster_id, maximum_minor_depth, config.minor_allele_min_depth,
                str(config.viralconsensus_only).lower(),
            )
            candidates = [_viralconsensus_candidate(
                cluster_id, sam, representative_path, directory,
                viral_consensus, config,
            )]
        cluster_sequences = {
            allele: sequences_by_locus[locus][allele]
            for allele in members[cluster_id]
        }
        drafts = _retained_candidate_drafts(
            row, candidates, consensus_method, cluster_id, locus, representative,
            members[cluster_id], cluster_sequences, reference_sequence, evidence,
            assigned, exclusive, ambiguous, directory, config,
        )
        if not drafts and consensus_method == "cliquesnv":
            logging.getLogger("bridge.rmlst").info(
                "consensus_fallback=viralconsensus cluster=%s "
                "reason=no_cliquesnv_haplotype_passed_filters",
                cluster_id,
            )
            candidates = [_viralconsensus_candidate(
                cluster_id, sam, representative_path, directory,
                viral_consensus, config,
            )]
            consensus_method = "viralconsensus"
            drafts = _retained_candidate_drafts(
                row, candidates, consensus_method, cluster_id, locus,
                representative, members[cluster_id], cluster_sequences,
                reference_sequence, evidence, assigned, exclusive, ambiguous,
                directory, config,
            )
        if drafts:
            return drafts
        return [_empty_cluster_draft(
            row, cluster_id, locus, representative, members[cluster_id], assigned,
            exclusive, ambiguous, len(reference_sequence), directory, config,
            "no_haplotype_passed_filters",
        )]
    except Exception:
        logging.getLogger("bridge.rmlst").exception(
            "subconsensus_status=failed cluster=%s", cluster_id
        )
        return [_empty_cluster_draft(
            row, cluster_id, locus, representative, members[cluster_id], assigned,
            exclusive, ambiguous, len(records[representative].sequence), directory,
            config, "haplotype_reconstruction_failed",
        )]


def _retained_candidate_drafts(
    row,
    candidates: list[HaplotypeCandidate],
    consensus_method: str,
    cluster_id: str,
    locus: str,
    representative: str,
    source_members: tuple[str, ...],
    cluster_sequences: dict[str, str],
    reference_sequence: str,
    evidence,
    assigned: int,
    exclusive: int,
    ambiguous: int,
    directory: Path,
    config: Config,
) -> list[dict[str, object]]:
    drafts: list[dict[str, object]] = []
    for ordinal, candidate in enumerate(
        sorted(candidates, key=lambda value: (-value.proportion, value.identifier)),
        1,
    ):
        padded_sequence = pad_terminal_consensus(
            candidate.sequence, cluster_sequences
        )
        if padded_sequence != candidate.sequence:
            logging.getLogger("bridge.rmlst").info(
                "subconsensus_status=terminally_padded cluster=%s candidate=%s "
                "original_length=%d padded_length=%d",
                cluster_id, candidate.identifier, len(candidate.sequence),
                len(padded_sequence),
            )
            candidate = HaplotypeCandidate(
                candidate.identifier, padded_sequence, candidate.proportion,
                candidate.source_path,
            )
        cluster_matches = compatible_alleles(
            candidate.sequence, cluster_sequences
        )
        if all(base in "ACGT" for base in candidate.sequence.upper()):
            exact_cluster_matches = tuple(
                allele
                for allele in cluster_matches
                if len(cluster_sequences[allele]) == len(candidate.sequence)
            )
        else:
            exact_cluster_matches = ()
        evaluation = evaluate_haplotype(
            candidate, reference_sequence, evidence, assigned,
            min_proportion=config.subconsensus_min_proportion,
            min_supporting_fragments=config.subconsensus_min_fragments,
            min_depth=config.subconsensus_min_depth,
            min_breadth=config.subconsensus_min_breadth,
            validate_variants=not bool(exact_cluster_matches),
        )
        masked_sequence = mask_weak_variants(
            candidate.sequence, evaluation.reasons
        )
        if masked_sequence != candidate.sequence:
            logging.getLogger("bridge.rmlst").info(
                "subconsensus_status=weak_variants_masked cluster=%s candidate=%s "
                "masked_positions=%d",
                cluster_id, candidate.identifier,
                sum(
                    reason.startswith("weak_variant_evidence_at_")
                    for reason in evaluation.reasons
                ),
            )
            candidate = HaplotypeCandidate(
                candidate.identifier, masked_sequence, candidate.proportion,
                candidate.source_path,
            )
            cluster_matches = compatible_alleles(
                candidate.sequence, cluster_sequences
            )
            exact_cluster_matches = ()
            evaluation = evaluate_haplotype(
                candidate, reference_sequence, evidence, assigned,
                min_proportion=config.subconsensus_min_proportion,
                min_supporting_fragments=config.subconsensus_min_fragments,
                min_depth=config.subconsensus_min_depth,
                min_breadth=config.subconsensus_min_breadth,
            )
        subconsensus_id = f"{cluster_id}_SC{ordinal:04d}"
        if not evaluation.retained:
            logging.getLogger("bridge.rmlst").info(
                "subconsensus_status=filtered cluster=%s candidate=%s "
                "proportion=%.6f supporting_fragments=%.3f "
                "projected_median_depth=%.3f callable_breadth=%.6f reasons=%s",
                cluster_id, candidate.identifier, candidate.proportion,
                evaluation.supporting_fragments,
                evaluation.projected_median_depth, evaluation.callable_breadth,
                ";".join(evaluation.reasons),
            )
            continue
        consensus_path = write_fasta(
            [FastaRecord(
                subconsensus_id,
                f"{subconsensus_id} Freq={candidate.proportion:.6f}",
                candidate.sequence,
            )],
            directory / f"{subconsensus_id}.fasta",
        )
        draft = _candidate_draft(
            row, cluster_id, subconsensus_id, locus, representative,
            source_members, assigned, exclusive, ambiguous, candidate.sequence,
            consensus_path, directory, config, candidate.proportion,
            evaluation.supporting_fragments, evaluation.projected_median_depth,
            evaluation.callable_breadth, evaluation.variant_evidence_passed,
            consensus_method,
        )
        draft["composite_from_members"] = (
            config.filter_composite_consensus
            and not bool(exact_cluster_matches)
            and composite_consensus_possible(
                candidate.sequence, cluster_sequences
            )
        )
        if draft["composite_from_members"]:
            logging.getLogger("bridge.rmlst").info(
                "consensus_status=composite_from_members cluster=%s "
                "subconsensus=%s members=%d",
                cluster_id, subconsensus_id, len(cluster_sequences),
            )
        if exact_cluster_matches:
            draft["compatible_alleles"] = exact_cluster_matches
            draft["allele_search_scope"] = "cluster_members_exact"
            draft["status"] = "complete_exact_match"
            logging.getLogger("bridge.rmlst").info(
                "allele_search_status=cluster_exact cluster=%s subconsensus=%s "
                "hits=%d variant_validation=exact_match full_blast=skipped",
                cluster_id, subconsensus_id, len(exact_cluster_matches),
            )
        drafts.append(draft)
    return drafts


def _candidate_draft(
    row,
    cluster_id: str,
    subconsensus_id: str,
    locus: str,
    representative: str,
    source_members: tuple[str, ...],
    assigned: int,
    exclusive: int,
    ambiguous: int,
    sequence: str,
    consensus_path: Path,
    directory: Path,
    config: Config,
    proportion: float,
    supporting_fragments: float,
    projected_median_depth: float,
    projected_breadth: float,
    variant_evidence_passed: bool,
    consensus_method: str,
) -> dict[str, object]:
    callable_bases = sum(base in "ACGT" for base in sequence)
    callable_breadth = callable_bases / len(sequence) if sequence else 0.0
    return {
        "sample_id": config.sample_id,
        "locus_id": locus,
        "runtime_cluster_id": cluster_id,
        "representative": representative,
        "source_members": source_members,
        "est_counts": float(row.est_counts) * proportion,
        "assigned_fragments": round(assigned * proportion),
        "exclusive_fragments": round(exclusive * proportion),
        "ambiguous_fragments": round(ambiguous * proportion),
        "status": "pending_full_reference",
        "consensus": sequence,
        "callable_bases": callable_bases,
        "allele_length": len(sequence),
        "callable_breadth": callable_breadth,
        "compatible_alleles": (),
        "allele_search_scope": "full_reference_batched",
        "consensus_path": consensus_path,
        "blast_path": None,
        "subconsensus_id": subconsensus_id,
        "haplotype_proportion": proportion,
        "supporting_fragments": supporting_fragments,
        "projected_median_depth": projected_median_depth,
        "projected_breadth": projected_breadth,
        "variant_evidence_passed": variant_evidence_passed,
        "consensus_method": consensus_method,
        "composite_from_members": False,
        "directory": directory,
    }


def _empty_cluster_draft(
    row,
    cluster_id: str,
    locus: str,
    representative: str,
    source_members: tuple[str, ...],
    assigned: int,
    exclusive: int,
    ambiguous: int,
    allele_length: int,
    directory: Path,
    config: Config,
    status: str,
) -> dict[str, object]:
    return {
        "sample_id": config.sample_id,
        "locus_id": locus,
        "runtime_cluster_id": cluster_id,
        "representative": representative,
        "source_members": source_members,
        "est_counts": 0.0,
        "assigned_fragments": assigned,
        "exclusive_fragments": exclusive,
        "ambiguous_fragments": ambiguous,
        "status": status,
        "consensus": "",
        "callable_bases": 0,
        "allele_length": allele_length,
        "callable_breadth": 0.0,
        "compatible_alleles": (),
        "allele_search_scope": "not_searched",
        "consensus_path": None,
        "blast_path": None,
        "subconsensus_id": f"{cluster_id}_SC0000",
        "haplotype_proportion": 0.0,
        "supporting_fragments": 0.0,
        "projected_median_depth": 0.0,
        "projected_breadth": 0.0,
        "variant_evidence_passed": False,
        "consensus_method": "none",
        "composite_from_members": False,
        "directory": directory,
    }


def _read_single_consensus(path: Path) -> str:
    """Read exactly one plain FASTA record without allele-ID normalization."""
    records: list[str] = []
    sequence: list[str] = []
    saw_header = False
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if saw_header:
                    records.append("".join(sequence).upper())
                saw_header = True
                sequence = []
            elif not saw_header:
                raise ValueError(f"Sequence before FASTA header in {path}")
            else:
                sequence.append(line)
    if saw_header:
        records.append("".join(sequence).upper())
    if len(records) != 1:
        raise ValueError(
            f"ViralConsensus must produce exactly one FASTA record in {path}"
        )
    invalid = sorted(set(records[0]) - set("ACGTN"))
    if not records[0] or invalid:
        raise ValueError(
            f"Invalid ViralConsensus sequence in {path}: "
            f"{''.join(invalid) or 'empty'}"
        )
    return records[0]


def _batch_full_reference_search(
    pending: list[dict[str, object]],
    reference: Path,
    sequences_by_locus: dict[str, dict[str, str]],
    root: Path,
    config: Config,
) -> None:
    logger = logging.getLogger("bridge.rmlst")
    blastn = executable(config.blastn)
    makeblastdb = companion_executable(config.makeblastdb, config.blastn)
    batch = root / "full_reference_blast"
    batch.mkdir(parents=True, exist_ok=True)
    query = write_fasta([
        FastaRecord(
            str(draft["subconsensus_id"]),
            str(draft["subconsensus_id"]),
            str(draft["consensus"]),
        )
        for draft in pending
    ], batch / "queries.fasta")
    cache_root = config.blastdb_cache_dir or config.output.parent
    database = _ensure_blast_database(reference, cache_root, makeblastdb)
    combined = batch / "blast.tsv"
    run([
        blastn, "-task", "blastn", "-query", str(query),
        "-db", str(database), "-dust", "no", "-word_size", "7",
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qlen slen bitscore",
        "-out", str(combined), "-num_threads", str(config.threads),
    ], "blastn")

    lines_by_query: dict[str, list[str]] = defaultdict(list)
    with combined.open(encoding="utf-8") as handle:
        for line in handle:
            query_id = line.split("\t", 1)[0]
            lines_by_query[query_id].append(line)
    for draft in pending:
        cluster_id = str(draft["runtime_cluster_id"])
        subconsensus_id = str(draft["subconsensus_id"])
        locus = str(draft["locus_id"])
        blast_path = Path(draft["directory"]) / f"blast_{subconsensus_id}.tsv"
        blast_path.write_text("".join(lines_by_query[subconsensus_id]), encoding="utf-8")
        compatible = compatible_alleles(
            str(draft["consensus"]), sequences_by_locus[locus]
        )
        draft["blast_path"] = blast_path
        draft["compatible_alleles"] = compatible
        breadth = float(draft["callable_breadth"])
        consensus = str(draft["consensus"])
        all_matches_same_length = all(
            len(sequences_by_locus[locus][allele]) == len(consensus)
            for allele in compatible
        )
        if not compatible:
            draft["status"] = "no_compatible_match"
        elif (
            breadth == 1.0
            and all(base in "ACGT" for base in consensus.upper())
            and all_matches_same_length
        ):
            draft["status"] = "complete_exact_match"
        elif len(compatible) == 1:
            draft["status"] = "unique_partial_match"
        else:
            draft["status"] = "ambiguous_partial_match"
    logger.info(
        "allele_search_status=batch_completed queries=%d database=%s",
        len(pending), database,
    )


def _ensure_blast_database(reference: Path, cache_root: Path, makeblastdb: str) -> Path:
    # kallisto.sha256 memoises on (path, mtime, size), so the manifest's digest
    # of the same immutable reference is reused here instead of reading the
    # whole file a second time.
    directory = cache_root / ".bridge_blastdb"
    directory.mkdir(parents=True, exist_ok=True)
    prefix = directory / f"reference_{sha256(reference)[:16]}"
    required = [Path(str(prefix) + suffix) for suffix in (".nhr", ".nin", ".nsq")]
    if not all(path.is_file() for path in required):
        run([
            makeblastdb, "-in", str(reference), "-dbtype", "nucl",
            "-parse_seqids", "-out", str(prefix),
        ], "makeblastdb")
    else:
        logging.getLogger("makeblastdb").info("database_status=reused path=%s", prefix)
    return prefix
