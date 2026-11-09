from __future__ import annotations

import csv
import hashlib
import logging
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from .errors import InsufficientMaterial
from .fasta import normalize_allele, write_fasta
from .models import Config, FastaRecord
from bridge.tools import executable, log_tool_output, run


CALL_COLUMNS = [
    "sample_id", "locus_id", "Template", "length", "eff_length",
    "est_counts", "initial_est_counts", "bootstrap_mean", "tpm",
    "bootstrap_support", "bootstrap_reps", "requested_bootstrap_threshold",
    "effective_bootstrap_threshold", "bootstrap_threshold_relaxed",
    "min_bootstrap_count", "max_bootstrap_count",
]


# The allele reference is large and immutable for the life of a run, but three
# callers want its digest: the manifest, the automatic index name, and the BLAST
# database cache key. Memoising on identity + mtime + size means it is read once
# and still re-hashed if the file is genuinely replaced.
_SHA256_CACHE: dict[tuple[str, int, int], str] = {}


def sha256(path: Path) -> str:
    try:
        stat = path.stat()
        key = (str(path.resolve()), stat.st_mtime_ns, stat.st_size)
    except OSError:
        key = None
    if key is not None and key in _SHA256_CACHE:
        return _SHA256_CACHE[key]
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    result = digest.hexdigest()
    if key is not None:
        _SHA256_CACHE[key] = result
    return result


def adjacent_index(reference: Path) -> Path | None:
    candidates = (
        reference.with_suffix(".kallisto.idx"),
        reference.with_suffix(".idx"),
        Path(str(reference) + ".kallisto.idx"),
        Path(str(reference) + ".idx"),
    )
    return next((path for path in candidates if path.is_file()), None)


def ensure_index(reference: Path, work: Path, kallisto: str) -> Path:
    found = adjacent_index(reference)
    if found is not None:
        logging.getLogger("kallisto").info("index_status=reused path=%s", found)
        return found
    index = work / f"reference_{sha256(reference)[:12]}.kallisto.idx"
    if not index.is_file():
        run([kallisto, "index", "-i", str(index), str(reference)], "kallisto")
    return index


def index_kmer_size(index: Path, kallisto: str) -> int:
    """Read the k-mer size embedded in a Kallisto index."""
    command = [kallisto, "inspect", str(index)]
    logging.getLogger("kallisto").info("command=%s", " ".join(command))
    completed = subprocess.run(command, capture_output=True, text=True)
    if completed.returncode:
        log_tool_output(
            "kallisto", stdout=completed.stdout, stderr=completed.stderr,
        )
        raise subprocess.CalledProcessError(
            completed.returncode, command,
            output=completed.stdout, stderr=completed.stderr,
        )
    output = completed.stdout + "\n" + completed.stderr
    match = re.search(r"(?:\[index\]\s*)?k-mer length:\s*(\d+)", output)
    if match is None:
        raise ValueError(f"Unable to determine k-mer size from Kallisto index: {index}")
    return int(match.group(1))


def _supplied_quant_directory(root: Path, index: Path) -> Path:
    """Keep cached quantification separate for each supplied index revision."""
    stat = index.stat()
    signature = f"{index.resolve()}\0{stat.st_size}\0{stat.st_mtime_ns}"
    digest = hashlib.sha256(signature.encode()).hexdigest()[:12]
    return root / f"quant_index_{digest}"


def run_quant(config: Config, reference: Path) -> tuple[Path, Path, int | None]:
    kallisto = executable(config.kallisto)
    root = config.tempdir / "kallisto"
    root.mkdir(parents=True, exist_ok=True)
    if config.kallisto_index is not None:
        index = config.kallisto_index
        kmer_size = index_kmer_size(index, kallisto)
        logging.getLogger("kallisto").info(
            "index_status=supplied path=%s kmer_size=%d", index, kmer_size,
        )
        quant = _supplied_quant_directory(root, index)
    else:
        index = ensure_index(reference, root, kallisto)
        kmer_size = None
        quant = root / "quant"
    abundance = quant / "abundance.tsv"
    bootstrap_files = list(quant.glob("bs_abundance_*.tsv")) if quant.is_dir() else []
    complete = abundance.is_file() and len(bootstrap_files) == config.bootstraps
    if not complete:
        if quant.exists():
            import shutil
            shutil.rmtree(quant)
        quant.mkdir(parents=True)
        command = [
            kallisto, "quant", "-i", str(index), "--plaintext",
            "-b", str(config.bootstraps), f"--seed={config.seed}",
            "-t", str(config.threads), "-o", str(quant),
        ]
        if len(config.reads) == 1:
            raise ValueError("Single-end reads are not yet supported by bridge.rmlst")
        command.extend(map(str, config.reads))
        run(command, "kallisto")
    else:
        logging.getLogger("kallisto").info("quant_status=reused path=%s", quant)
    return quant, index, kmer_size


def run_candidate_pseudobam(
    config: Config,
    candidates: pd.DataFrame,
    records: dict[str, FastaRecord],
    index_kmer_size: int | None = None,
) -> Path:
    """Generate pseudobAM against candidates, avoiding a huge full-index BAM header."""
    kallisto = executable(config.kallisto)
    root = config.tempdir / "kallisto" / "candidate_pseudobam"
    quant = root / "quant"
    pseudobam = quant / "pseudoalignments.bam"
    logger = logging.getLogger("kallisto")
    selected = sorted(set(candidates["Template"].astype(str)))
    missing = set(selected) - set(records)
    if missing:
        raise ValueError(f"Candidate alleles absent from FASTA: {sorted(missing)[:5]}")

    # kallisto cannot build a de Bruijn graph from a sequence shorter than k,
    # and an empty sequence is written as a bare header - malformed FASTA that
    # `kallisto index` can die on with SIGSEGV rather than an error. Drop those
    # here, loudly.
    minimum_length = index_kmer_size or 1
    usable = [a for a in selected if len(records[a].sequence) >= minimum_length]
    unusable = [a for a in selected if len(records[a].sequence) < minimum_length]
    if unusable:
        logger.warning(
            "candidate_alleles_status=dropped count=%d minimum_length=%d alleles=%s",
            len(unusable), minimum_length, ",".join(unusable[:10]),
            extra={"console": True},
        )
    if not usable:
        raise InsufficientMaterial(
            f"No candidate allele reached the {minimum_length} bp kallisto needs",
            stage="candidate_pseudobam",
        )
    lengths = [len(records[a].sequence) for a in usable]
    logger.info(
        "candidate_alleles_status=written count=%d shortest=%d longest=%d",
        len(usable), min(lengths), max(lengths),
    )

    candidate_fasta = write_fasta(
        [records[allele] for allele in usable], root / "candidate_alleles.fasta"
    )
    selected = usable
    kmer_suffix = f".k{index_kmer_size}" if index_kmer_size is not None else ""
    index = root / f"candidate_alleles{kmer_suffix}.kallisto.idx"
    if not index.is_file():
        command = [kallisto, "index", "-i", str(index)]
        if index_kmer_size is not None:
            command.extend(["-k", str(index_kmer_size)])
        command.append(str(candidate_fasta))
        # One retry: this step has died with a signal rather than an error, and
        # a crash that does not repeat costs a second here instead of the whole
        # sample. A crash that does repeat keeps its input for inspection -
        # tempdir is usually the job's scratch space and is deleted with it.
        for attempt in (1, 2):
            try:
                run(command, "kallisto")
                break
            except subprocess.CalledProcessError as error:
                index.unlink(missing_ok=True)
                if attempt == 1:
                    logger.warning(
                        "candidate_index_status=retrying targets=%d error=%s",
                        len(usable), error, extra={"console": True},
                    )
                    continue
                kept = _preserve_failed_input(candidate_fasta, config)
                raise RuntimeError(
                    f"kallisto index failed twice on {len(usable)} candidate alleles "
                    f"({min(lengths)}-{max(lengths)} bp): {error}. The input it "
                    f"choked on has been kept at {kept}."
                ) from error

    complete = pseudobam.is_file() and _bam_has_eof(pseudobam)
    if not complete:
        if quant.exists():
            import shutil
            shutil.rmtree(quant)
        quant.mkdir(parents=True)
        command = [
            kallisto, "quant", "-i", str(index), "--plaintext", "--pseudobam",
            "-t", str(config.threads), "-o", str(quant),
        ]
        command.extend(map(str, config.reads))
        run(command, "kallisto", stdout_path=pseudobam)
    else:
        logging.getLogger("kallisto").info(
            "candidate_pseudobam_status=reused targets=%d path=%s",
            len(selected), pseudobam,
        )
    logging.getLogger("kallisto").info(
        "candidate_pseudobam_status=completed targets=%d path=%s",
        len(selected), pseudobam,
    )
    return pseudobam


def _preserve_failed_input(candidate_fasta: Path, config: Config) -> Path:
    """Copy a crashing input somewhere it will survive the job.

    ``tempdir`` is normally the scheduler's per-job scratch, which is deleted
    the moment the job ends, so the file that caused the failure is gone before
    anyone can look at it. The output folder outlives the job.
    """
    import shutil

    destination = config.output / "failed_candidate_alleles.fasta"
    try:
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(candidate_fasta, destination)
        return destination
    except OSError as error:  # pragma: no cover - a full or read-only filesystem
        logging.getLogger("kallisto").warning(
            "could not preserve %s: %s", candidate_fasta, error)
        return candidate_fasta


def _bam_has_eof(path: Path) -> bool:
    # Standard BGZF EOF marker required for a cleanly closed BAM.
    marker = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")
    if path.stat().st_size < len(marker):
        return False
    with path.open("rb") as handle:
        handle.seek(-len(marker), 2)
        return handle.read() == marker


def _read_abundance(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    targets: list[str] = []
    rows: list[dict[str, str]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            targets.append(normalize_allele(row["target_id"]))
            rows.append(row)
    if len(targets) != len(set(targets)):
        raise ValueError(f"Duplicate normalized targets in {path}")
    return targets, rows


def derive_candidates(config: Config, quant: Path) -> pd.DataFrame:
    targets, abundance = _read_abundance(quant / "abundance.tsv")
    target_count = len(targets)
    support = np.zeros(target_count, dtype=np.uint16)
    sums = np.zeros(target_count, dtype=np.float64)
    minimum = np.full(target_count, np.inf, dtype=np.float64)
    maximum = np.zeros(target_count, dtype=np.float64)
    logger = logging.getLogger("kallisto")
    for replicate in range(config.bootstraps):
        path = quant / f"bs_abundance_{replicate}.tsv"
        # Bootstrap files share the index-defined target order. Reading only
        # the numeric column with Pandas' C parser avoids normalizing millions
        # of repeated allele identifiers in Python.
        values = pd.read_csv(
            path,
            sep="\t",
            usecols=["est_counts"],
            dtype={"est_counts": "float64"},
            memory_map=True,
        )["est_counts"].to_numpy(copy=False)
        if len(values) != target_count:
            raise ValueError(
                f"Bootstrap target count differs in {path}: "
                f"expected {target_count}, found {len(values)}"
            )
        sums += values
        support += values > 0
        np.minimum(minimum, values, out=minimum)
        np.maximum(maximum, values, out=maximum)
        if (replicate + 1) % 10 == 0 or replicate + 1 == config.bootstraps:
            logger.info(
                "bootstrap_parse_status=running completed=%d total=%d",
                replicate + 1, config.bootstraps,
            )
    records: list[dict[str, object]] = []
    loci = [target.split("_", 1)[0] for target in targets]
    maximum_support_by_locus: dict[str, int] = {}
    for locus, n in zip(loci, support):
        maximum_support_by_locus[locus] = max(
            maximum_support_by_locus.get(locus, 0), int(n)
        )
    effective_threshold_by_locus = {
        locus: min(config.support_threshold, maximum_support)
        for locus, maximum_support in maximum_support_by_locus.items()
    }
    for target, row, n, total, low, high in zip(
        targets, abundance, support, sums, minimum, maximum
    ):
        locus = target.split("_", 1)[0]
        effective_threshold = effective_threshold_by_locus[locus]
        initial = float(row["est_counts"])
        mean = total / config.bootstraps
        # A zero effective threshold rescues a locus only from a positive
        # initial estimate; it must not turn every zero-abundance reference
        # allele into a runtime candidate.
        if n < effective_threshold or (effective_threshold == 0 and initial <= 0):
            continue
        records.append({
            "sample_id": config.sample_id,
            "locus_id": locus,
            "Template": target,
            "length": int(row["length"]),
            "eff_length": float(row["eff_length"]),
            "est_counts": mean if config.use_bootstrap_mean else initial,
            "initial_est_counts": initial,
            "bootstrap_mean": mean,
            "tpm": float(row["tpm"]),
            "bootstrap_support": n,
            "bootstrap_reps": config.bootstraps,
            "requested_bootstrap_threshold": config.support_threshold,
            "effective_bootstrap_threshold": effective_threshold,
            "bootstrap_threshold_relaxed": effective_threshold < config.support_threshold,
            "min_bootstrap_count": low,
            "max_bootstrap_count": high,
        })
    result = pd.DataFrame.from_records(records, columns=CALL_COLUMNS)
    if result.empty:
        raise InsufficientMaterial(
            "No alleles passed the Kallisto candidate rule",
            stage="kallisto",
        )
    relaxed_loci = sum(
        threshold < config.support_threshold
        for threshold in effective_threshold_by_locus.values()
    )
    logging.getLogger("kallisto").info(
        "candidate_status=completed requested_threshold=%d retained=%d loci=%d "
        "relaxed_loci=%d minimum_effective_threshold=%d",
        config.support_threshold, len(result), result["locus_id"].nunique(),
        relaxed_loci, min(effective_threshold_by_locus.values(), default=0),
    )
    return result
