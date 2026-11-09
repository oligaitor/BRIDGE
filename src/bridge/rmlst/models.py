from __future__ import annotations

from dataclasses import dataclass, field
import os
from pathlib import Path


@dataclass(frozen=True)
class Config:
    allele_target: Path
    seqdef: Path
    reads: tuple[Path, ...]
    sample_id: str
    output: Path
    tempdir: Path
    threads: int = 4
    memory_limit_gb: float | None = None
    keep_temp: bool = False
    verbose: bool = False
    kallisto: str = "kallisto"
    kallisto_index: Path | None = None
    mash: str = "mash"
    short_read_aligner: str = "minimap2"
    bwa_mem2: str = "bwa-mem2"
    minimap2: str = "minimap2"
    samtools: str = "samtools"
    cliquesnv: str = field(
        default_factory=lambda: os.environ.get(
            "BRIDGE_CLIQUESNV", "cliquesnv"
        )
    )
    cliquesnv_time_limit_seconds: int | None = None
    viral_consensus: str = "viral_consensus"
    viralconsensus_only: bool = False
    blastn: str = "blastn"
    makeblastdb: str = "makeblastdb"
    # Where the reusable full-reference BLAST database is cached. None means
    # `output.parent`, correct when -o is <run_root>/<SAMPLE>. Callers that nest
    # the output deeper (bridge writes <run_root>/<SAMPLE>/rmlst) must name the
    # run root, or every sample builds its own copy.
    blastdb_cache_dir: Path | None = None
    bootstraps: int = 100
    support_threshold: int = 100
    seed: int = 1000
    use_bootstrap_mean: bool = True
    mash_k: int = 21
    mash_threshold: float = 0.01
    min_fragments: int = 1
    subconsensus_min_proportion: float = 0.01
    # Keep these defaults in step with cli.py.
    subconsensus_min_fragments: int = 20
    subconsensus_min_depth: float = 10.0
    subconsensus_min_breadth: float = 0.80
    subconsensus_min_end_distance: int = 5
    minor_allele_min_depth: int = 5
    min_profile_loci: int = 3
    # Unambiguous loci a species needs before it may be reported. See
    # species.direct_profile_evidence for why this is the bar that separates a
    # real organism from a near-neighbour riding on shared evidence.
    min_unambiguous_loci: int = 3
    # Require an exact allele match before a locus counts as profile evidence.
    # False counts every locus that resolves to one species, which is what a
    # real low-abundance organism looks like: broad locus support, imperfect
    # assembly. See species.DEFAULT_REQUIRE_ACCEPTED_CONSENSUS.
    require_accepted_consensus: bool = False
    read_assignment: str = "compatible"
    filter_composite_consensus: bool = True
    distinguish_subspecies: bool = False


@dataclass(frozen=True)
class FastaRecord:
    identifier: str
    description: str
    sequence: str


@dataclass(frozen=True)
class ConsensusCall:
    sample_id: str
    locus_id: str
    runtime_cluster_id: str
    representative: str
    source_members: tuple[str, ...]
    est_counts: float
    assigned_fragments: int
    exclusive_fragments: int
    ambiguous_fragments: int
    status: str
    consensus: str
    callable_bases: int
    allele_length: int
    callable_breadth: float
    compatible_alleles: tuple[str, ...]
    compatible_species: tuple[str, ...]
    allele_search_scope: str
    consensus_path: Path | None
    blast_path: Path | None
    subconsensus_id: str = ""
    haplotype_proportion: float = 1.0
    supporting_fragments: float = 0.0
    projected_median_depth: float = 0.0
    projected_breadth: float = 0.0
    variant_evidence_passed: bool = False
    consensus_method: str = ""
    composite_from_members: bool = False
