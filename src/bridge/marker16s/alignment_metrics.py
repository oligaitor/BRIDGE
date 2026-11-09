"""Identity and coverage for the 16S arm, from a real alignment.

kallisto pseudoalignment decides *which* organism the reads belong to but
produces no alignment, so on its own it cannot say how well the reads matched
the reference it named - the percent identity and query coverage a Sanger 16S
report would give. This stage adds them, after the calls are made and without
changing them:

1. the MIMt sequences of every called genus are written to a subset FASTA
   (one streaming pass over the reference; ``load_reference`` keeps headers
   only);
2. the sample's reads are aligned to that subset once with minimap2, using the
   same short-read parameters the rMLST arm uses, then sorted and indexed;
3. per genus, from the primary alignments: how many reads landed, their mean
   identity (``1 - NM / aligned columns``), and the reference most of them
   chose; for that reference, breadth (fraction of positions covered at least
   once) and mean depth;
4. per genus, a consensus of the reads on that top reference (ViralConsensus,
   as the rMLST arm uses) is BLASTed against the genus's own references, and
   the best hit's percent identity and query coverage are reported.

Near-identical 16S references give MAPQ of ~0, so nothing is filtered on
mapping quality; ``-F 0x904`` (unmapped, secondary, supplementary) is what
makes each read count once, against the reference minimap2 chose for it.

Everything here is optional at runtime: a missing tool skips the stage with
one warning and the calls file carries blanks; a failure on one genus is
recorded in that genus's ``alignment_status`` and the others proceed.
"""

from __future__ import annotations

import logging
import re
import subprocess
from collections import Counter
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Iterable

from bridge.tools import companion_executable, executable, run
from .reference import ReferenceTaxon, parse_lineage_header

# Kept equal to consensus.MINIMAP2_SHORT_READ_PARAMETERS (a test pins this);
# repeated here so the 16S arm does not import the rMLST arm's heavier modules.
MINIMAP2_SHORT_READ_PARAMETERS = (
    "-ax", "sr", "-k", "17", "-w", "7", "-m", "20", "-s", "30",
)
# Same field list as the rMLST arm's allele search, so the two BLAST outputs
# read alike.
BLAST_OUTFMT = ("6 qseqid sseqid pident length mismatch gapopen qstart qend "
                "sstart send qlen slen bitscore")
PRIMARY_ONLY_FLAGS = "0x904"     # drop unmapped, secondary, supplementary
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")

METRIC_COLUMNS = (
    "n_aligned_reads", "mean_read_identity", "ref_breadth", "mean_depth",
    "top_reference", "consensus_length", "consensus_pident", "consensus_qcov",
    "consensus_top_hit", "alignment_status",
)


@dataclass(frozen=True)
class AlignmentTools:
    minimap2: str
    samtools: str
    viral_consensus: str
    blastn: str
    makeblastdb: str


@dataclass(frozen=True)
class GenusAlignmentMetrics:
    genus: str
    n_references: int = 0
    n_aligned_reads: int = 0
    mean_read_identity: float | None = None   # 0-1
    top_reference: str | None = None          # MIMt target_id
    ref_breadth: float | None = None          # 0-1, of the top reference
    mean_depth: float | None = None
    consensus_length: int | None = None
    consensus_pident: float | None = None     # percent, BLAST pident
    consensus_qcov: float | None = None       # percent of the consensus aligned
    consensus_top_hit: str | None = None
    status: str = "ok"                        # ok | no_reads | failed:<stage>


def resolve_tools(*, minimap2: str = "minimap2", samtools: str = "samtools",
                  viral_consensus: str = "viral_consensus", blastn: str = "blastn",
                  makeblastdb: str = "makeblastdb") -> AlignmentTools | None:
    """Every binary the stage needs, or None (with one warning) if any is missing."""
    try:
        return AlignmentTools(
            minimap2=executable(minimap2), samtools=executable(samtools),
            viral_consensus=executable(viral_consensus), blastn=executable(blastn),
            makeblastdb=companion_executable(makeblastdb, blastn),
        )
    except (FileNotFoundError, OSError) as error:
        logging.getLogger("marker16s").warning(
            "alignment_metrics_status=skipped reason=missing_tool detail=%s", error)
        return None


# --- pure helpers (unit-tested without any binary) ---------------------------

def read_identity(cigar: str, nm: int) -> tuple[int, float]:
    """Aligned columns and identity for one read: ``1 - NM / columns``.

    Columns are the M/=/X/I/D operations - the length of the alignment as
    BLAST would count it - so NM (mismatches plus inserted and deleted bases)
    over them is the fraction of the alignment that disagrees.
    """
    columns = sum(int(count) for count, op in CIGAR_RE.findall(cigar) if op in "M=XID")
    if columns <= 0:
        return 0, 0.0
    return columns, max(0.0, 1.0 - nm / columns)


def _nm_tag(fields: list[str]) -> int | None:
    for field in fields[11:]:
        if field.startswith("NM:i:"):
            try:
                return int(field[5:])
            except ValueError:
                return None
    return None


def parse_primary_alignments(lines: Iterable[str], genus_of: dict[str, str]
                             ) -> dict[str, tuple[int, float, Counter]]:
    """SAM lines (primary alignments only) -> per genus: (reads, mean identity, refs).

    ``genus_of`` maps a reference id to its genus; reads on any other reference
    are ignored.
    """
    n: Counter = Counter()
    identity_sum: dict[str, float] = {}
    refs: dict[str, Counter] = {}
    for line in lines:
        if not line or line.startswith("@"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 11:
            continue
        rname, cigar = fields[2], fields[5]
        genus = genus_of.get(rname)
        if genus is None or cigar == "*":
            continue
        nm = _nm_tag(fields)
        if nm is None:
            continue
        columns, identity = read_identity(cigar, nm)
        if columns == 0:
            continue
        n[genus] += 1
        identity_sum[genus] = identity_sum.get(genus, 0.0) + identity
        refs.setdefault(genus, Counter())[rname] += 1
    return {genus: (n[genus], identity_sum[genus] / n[genus], refs[genus]) for genus in n}


def parse_depth(lines: Iterable[str], length: int) -> tuple[float, float]:
    """``samtools depth -a`` output -> (breadth at >=1x, mean depth) over ``length``."""
    covered = 0
    total = 0
    for line in lines:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        try:
            depth = int(parts[2])
        except ValueError:
            continue
        total += depth
        if depth > 0:
            covered += 1
    if length <= 0:
        return 0.0, 0.0
    return min(1.0, covered / length), total / length


def parse_blast_top_hit(lines: Iterable[str]) -> tuple[str, float, float] | None:
    """Best-bitscore row of an outfmt-6 table -> (subject id, pident, query coverage %)."""
    best: tuple[float, str, float, float] | None = None
    for line in lines:
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 13:
            continue
        try:
            pident = float(cols[2])
            qstart, qend, qlen = int(cols[6]), int(cols[7]), int(cols[10])
            bitscore = float(cols[12])
        except ValueError:
            continue
        qcov = (abs(qend - qstart) + 1) / qlen * 100 if qlen else 0.0
        if best is None or bitscore > best[0]:
            best = (bitscore, cols[1], pident, qcov)
    if best is None:
        return None
    return best[1], best[2], best[3]


def read_fasta_records(path: Path) -> Iterable[tuple[str, str]]:
    """(header without '>', sequence) for each record, streaming."""
    header: str | None = None
    chunks: list[str] = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:], []
            elif header is not None:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def safe_name(genus: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", genus) or "genus"


def write_genus_subsets(fasta: Path, reference: dict[str, ReferenceTaxon],
                        genera: Iterable[str], out_dir: Path
                        ) -> tuple[Path, dict[str, Path], dict[str, int]]:
    """One pass over the MIMt FASTA: the called genera's sequences, combined and per genus.

    Returns the combined FASTA, ``{genus: its FASTA}`` and ``{target_id: length}``.
    """
    wanted = {genus for genus in genera if genus}
    by_lower = {genus.lower(): genus for genus in wanted}
    out_dir.mkdir(parents=True, exist_ok=True)
    combined = out_dir / "references.fasta"
    per_genus: dict[str, Path] = {}
    handles: dict[str, object] = {}
    lengths: dict[str, int] = {}
    try:
        with open(combined, "w", encoding="utf-8") as combined_handle:
            for header, sequence in read_fasta_records(fasta):
                target_id = header.split()[0] if header.strip() else ""
                taxon = reference.get(target_id) or parse_lineage_header(header)
                genus = by_lower.get((taxon.genus or "").lower())
                if genus is None:
                    continue
                record = f">{target_id}\n{sequence}\n"
                combined_handle.write(record)
                if genus not in handles:
                    path = out_dir / safe_name(genus) / "references.fasta"
                    path.parent.mkdir(parents=True, exist_ok=True)
                    handles[genus] = open(path, "w", encoding="utf-8")
                    per_genus[genus] = path
                handles[genus].write(record)  # type: ignore[attr-defined]
                lengths[target_id] = len(sequence)
    finally:
        for handle in handles.values():
            handle.close()  # type: ignore[attr-defined]
    return combined, per_genus, lengths


def write_single_record(source: Path, target_id: str, out: Path) -> Path:
    for header, sequence in read_fasta_records(source):
        if header.split()[0] == target_id:
            out.write_text(f">{target_id}\n{sequence}\n", encoding="utf-8")
            return out
    raise ValueError(f"{target_id} not found in {source}")


# --- steps that call the binaries --------------------------------------------

def align_reads(tools: AlignmentTools, subset_fasta: Path, reads: list[Path],
                out_dir: Path, threads: int) -> Path:
    sam = out_dir / "alignment.sam"
    run([tools.minimap2, *MINIMAP2_SHORT_READ_PARAMETERS, "-t", str(threads),
         str(subset_fasta), *map(str, reads)], "minimap2", stdout_path=sam)
    bam = out_dir / "alignment.bam"
    run([tools.samtools, "sort", "-@", str(threads), "-o", str(bam), str(sam)], "samtools")
    run([tools.samtools, "index", str(bam)], "samtools")
    sam.unlink(missing_ok=True)
    return bam


def primary_alignments(tools: AlignmentTools, bam: Path, out: Path) -> Path:
    run([tools.samtools, "view", "-F", PRIMARY_ONLY_FLAGS, str(bam)], "samtools", stdout_path=out)
    return out


def reference_depth(tools: AlignmentTools, bam: Path, target_id: str, length: int,
                    out: Path) -> tuple[float, float]:
    run([tools.samtools, "depth", "-a", "-r", target_id, str(bam)], "samtools", stdout_path=out)
    with open(out, encoding="utf-8") as handle:
        return parse_depth(handle, length)


def single_reference_sam(source: Path, target_id: str, dest: Path) -> int:
    """Copy ``source`` to ``dest`` keeping only ``target_id`` in the header.

    ``samtools view <bam> <region>`` keeps the whole header, so the SAM for one
    reference still declares every reference in the genus subset, and
    ViralConsensus refuses a file with more than one reference. Records are
    all on ``target_id`` already (the region query); a mate placed on another
    reference is cleared to ``*``/0 so the file names no reference the header
    no longer declares.
    Returns the number of alignment records written.
    """
    records = 0
    with open(source, encoding="utf-8") as src, open(dest, "w", encoding="utf-8") as out:
        for line in src:
            if line.startswith("@"):
                if line.startswith("@SQ") and f"\tSN:{target_id}\t" not in line.rstrip("\n") + "\t":
                    continue
                out.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) > 7 and fields[6] not in ("=", "*", target_id):
                fields[6], fields[7] = "*", "0"
            out.write("\t".join(fields) + "\n")
            records += 1
    return records


def consensus_for_reference(tools: AlignmentTools, bam: Path, target_id: str,
                            single_reference: Path, out_dir: Path) -> Path:
    all_refs_sam = out_dir / "top_reference.all_refs.sam"
    run([tools.samtools, "view", "-h", "-F", PRIMARY_ONLY_FLAGS, "-o", str(all_refs_sam),
         str(bam), target_id], "samtools")
    top_sam = out_dir / "top_reference.sam"
    single_reference_sam(all_refs_sam, target_id, top_sam)
    all_refs_sam.unlink()
    consensus = out_dir / "consensus.fasta"
    run([tools.viral_consensus, "-i", str(top_sam), "-r", str(single_reference),
         "-o", str(consensus), "-q", "20", "-d", "1"], "viralconsensus")
    return consensus


def blast_consensus(tools: AlignmentTools, consensus: Path, genus_fasta: Path,
                    out_dir: Path, threads: int) -> tuple[str, float, float] | None:
    db = out_dir / "blastdb" / "references"
    db.parent.mkdir(parents=True, exist_ok=True)
    run([tools.makeblastdb, "-in", str(genus_fasta), "-dbtype", "nucl", "-out", str(db)],
        "makeblastdb")
    table = out_dir / "blast.tsv"
    run([tools.blastn, "-task", "blastn", "-query", str(consensus), "-db", str(db),
         "-dust", "no", "-word_size", "7", "-max_target_seqs", "5",
         "-outfmt", BLAST_OUTFMT, "-out", str(table), "-num_threads", str(threads)],
        "blastn")
    with open(table, encoding="utf-8") as handle:
        return parse_blast_top_hit(handle)


def _consensus_length(path: Path) -> int:
    for _, sequence in read_fasta_records(path):
        return len(sequence)
    return 0


def alignment_metrics_for_calls(calls, reference: dict[str, ReferenceTaxon], fasta: Path,
                                reads: list[Path], out_dir: Path, *,
                                tools: AlignmentTools, threads: int = 4,
                                strengths: tuple[str, ...] = ("species_strong", "genus_only"),
                                ) -> dict[str, GenusAlignmentMetrics]:
    """Run the stage for every called genus; return metrics keyed by genus."""
    logger = logging.getLogger("marker16s")
    genera = []
    for call in calls:
        if call.genus and call.strength in strengths and call.genus not in genera:
            genera.append(call.genus)
    if not genera:
        logger.info("alignment_metrics_status=nothing_to_align")
        return {}

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    combined, per_genus_fasta, lengths = write_genus_subsets(fasta, reference, genera, out_dir)
    genus_of = {target_id: reference[target_id].genus or ""
                for target_id in lengths if target_id in reference}
    for target_id in lengths:
        genus_of.setdefault(target_id, parse_lineage_header(target_id).genus or "")
    n_references = Counter(genus_of.values())
    results: dict[str, GenusAlignmentMetrics] = {}

    try:
        bam = align_reads(tools, combined, reads, out_dir, threads)
        primary = primary_alignments(tools, bam, out_dir / "primary.sam")
    except (subprocess.CalledProcessError, OSError) as error:
        logger.error("alignment_metrics_status=failed stage=align detail=%s", error)
        return {genus: GenusAlignmentMetrics(genus, n_references[genus], status="failed:align")
                for genus in genera}

    with open(primary, encoding="utf-8") as handle:
        stats = parse_primary_alignments(handle, genus_of)

    for genus in genera:
        found = stats.get(genus)
        if not found:
            results[genus] = GenusAlignmentMetrics(genus, n_references[genus], status="no_reads")
            continue
        n_reads, mean_identity, refs = found
        top_reference = refs.most_common(1)[0][0]
        metrics = GenusAlignmentMetrics(genus, n_references[genus], n_reads, mean_identity,
                                        top_reference)
        genus_dir = out_dir / safe_name(genus)
        genus_dir.mkdir(parents=True, exist_ok=True)
        try:
            breadth, depth = reference_depth(tools, bam, top_reference,
                                             lengths.get(top_reference, 0), genus_dir / "depth.tsv")
            metrics = replace(metrics, ref_breadth=breadth, mean_depth=depth)
        except (subprocess.CalledProcessError, OSError) as error:
            logger.warning("alignment_metrics_status=failed genus=%s stage=depth detail=%s",
                           genus, error)
            results[genus] = replace(metrics, status="failed:depth")
            continue
        try:
            single = write_single_record(per_genus_fasta[genus], top_reference,
                                         genus_dir / "top_reference.fasta")
            consensus = consensus_for_reference(tools, bam, top_reference, single, genus_dir)
            hit = blast_consensus(tools, consensus, per_genus_fasta[genus], genus_dir, threads)
            metrics = replace(metrics, consensus_length=_consensus_length(consensus))
            if hit is not None:
                subject, pident, qcov = hit
                metrics = replace(metrics, consensus_pident=pident, consensus_qcov=qcov,
                                  consensus_top_hit=subject)
        except (subprocess.CalledProcessError, OSError, ValueError) as error:
            logger.warning("alignment_metrics_status=failed genus=%s stage=consensus detail=%s",
                           genus, error)
            metrics = replace(metrics, status="failed:consensus")
        results[genus] = metrics
        logger.info(
            "alignment_metrics genus=%s reads=%d identity=%.4f breadth=%s depth=%s pident=%s qcov=%s status=%s",
            genus, metrics.n_aligned_reads, metrics.mean_read_identity or 0.0,
            f"{metrics.ref_breadth:.3f}" if metrics.ref_breadth is not None else "-",
            f"{metrics.mean_depth:.1f}" if metrics.mean_depth is not None else "-",
            f"{metrics.consensus_pident:.2f}" if metrics.consensus_pident is not None else "-",
            f"{metrics.consensus_qcov:.1f}" if metrics.consensus_qcov is not None else "-",
            metrics.status)
    return results


def metric_cells(metrics: GenusAlignmentMetrics | None) -> dict[str, str]:
    """The CSV cells for one genus; every cell blank where nothing was measured."""
    if metrics is None:
        return {column: "" for column in METRIC_COLUMNS}

    def fmt(value, spec: str) -> str:
        return "" if value is None else format(value, spec)

    return {
        "n_aligned_reads": str(metrics.n_aligned_reads) if metrics.n_aligned_reads else "",
        "mean_read_identity": fmt(metrics.mean_read_identity, ".4f"),
        "ref_breadth": fmt(metrics.ref_breadth, ".4f"),
        "mean_depth": fmt(metrics.mean_depth, ".1f"),
        "top_reference": metrics.top_reference or "",
        "consensus_length": "" if metrics.consensus_length is None else str(metrics.consensus_length),
        "consensus_pident": fmt(metrics.consensus_pident, ".2f"),
        "consensus_qcov": fmt(metrics.consensus_qcov, ".1f"),
        "consensus_top_hit": metrics.consensus_top_hit or "",
        "alignment_status": metrics.status,
    }
