"""The arm itself: KMA against ResFinder, for one sample.

    <out>/<SAMPLE>/amr/
        kma.res  kma.mapstat  kma.mat.gz  kma.fsa  kma.aln   KMA's own output
        amr_hits.csv                                         every template KMA reported
        amr_manifest.json                                    what ran, on what, how many reads
        amr.log

Why KMA
-------
ResFinder is highly redundant - hundreds of blaTEM, blaCTX-M and qnrB alleles
differing by a few bases - and the question per gene is "which allele, at what
identity and coverage". KMA was built for exactly that database and exactly
that question (it is what ResFinder itself runs on reads): with ``-1t1`` each
read goes to the one template it fits best, and ConClave resolves the ties a
plain aligner would split. kallisto, which the rMLST arm uses, answers a
different question - how much of each of many near-identical references is
present - and reports no identity.

The read denominator
--------------------
``amr_manifest.json`` records ``total_read_pairs``, counted from R1 here. It
is the RPM denominator for the resistance gate and must mean the same thing in
the sample and the control, which counting it the same way in both guarantees.
The per-gene numerator is ``fragments`` - KMA's mapstat ``fragmentCount``,
read pairs for paired input, so a gene's ``fragments`` should never exceed
``total_read_pairs``.
"""

from __future__ import annotations

import csv
import gzip
import json
import logging
import shlex
import subprocess
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from time import perf_counter

from bridge import __version__
from .calls import classify
from .database import family_of, load_family_overrides, load_phenotypes, parse_template, read_manifest
from .phenotype import annotate
from .settings import AmrSettings

LOGGER = logging.getLogger("bridge.amr.kma")

PREFIX = "kma"
HITS_NAME = "amr_hits.csv"
MANIFEST_NAME = "amr_manifest.json"

# KMA .res header -> our column. Percentages stay percentages (0-100).
RES_COLUMNS = {
    "template": "template",
    "score": "score",
    "expected": "expected",
    "template_length": "template_length",
    "template_identity": "template_identity",
    "template_coverage": "template_coverage",
    "query_identity": "query_identity",
    "query_coverage": "query_coverage",
    "depth": "depth",
    "q_value": "q_value",
    "p_value": "p_value",
}
MAPSTAT_COLUMNS = {
    "readcount": "reads",
    "fragmentcount": "fragments",
    "refcoveredpositions": "ref_covered_positions",
    "snpsum": "snp_sum",
    "insertsum": "insert_sum",
    "deletionsum": "deletion_sum",
}

HIT_COLUMNS = [
    "sample", "template", "gene", "variant", "accession", "family",
    # KMA, verbatim
    "template_length", "template_identity", "template_coverage",
    "query_identity", "query_coverage", "depth", "fragments", "reads",
    "score", "expected", "q_value", "p_value",
    "ref_covered_positions", "snp_sum", "insert_sum", "deletion_sum",
    # our reading of it (recomputed by the report from the columns above)
    "tier", "allele_level", "est_differences",
    # the database's, fixed per DB version
    "drug_class", "antibiotics", "mechanism", "phenotype_notes",
    "family_spectrum_varies", "family_possible_antibiotics",
    "flags",
]


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def _key(name: str) -> str:
    return name.strip().lstrip("#").strip().lower()


def parse_res(path: Path) -> dict[str, dict]:
    """KMA ``.res``: one row per template, tab-separated, header starts ``#Template``."""
    rows: dict[str, dict] = {}
    if not Path(path).is_file():
        return rows
    with open(path, encoding="utf-8", errors="replace") as handle:
        header = None
        for raw in handle:
            line = raw.rstrip("\n")
            if not line.strip():
                continue
            cells = line.split("\t")
            if header is None:
                if line.startswith("#"):
                    header = [_key(cell) for cell in cells]
                continue
            record = {}
            for name, value in zip(header, cells):
                if name in RES_COLUMNS:
                    record[RES_COLUMNS[name]] = value.strip()
            if record.get("template"):
                rows[record["template"]] = record
    return rows


def parse_mapstat(path: Path) -> tuple[dict[str, str], dict[str, dict]]:
    """KMA ``.mapstat`` (``-ef``): ``##`` metadata lines, a ``# refSequence`` header, rows."""
    meta: dict[str, str] = {}
    rows: dict[str, dict] = {}
    if not Path(path).is_file():
        return meta, rows
    with open(path, encoding="utf-8", errors="replace") as handle:
        header = None
        for raw in handle:
            line = raw.rstrip("\n")
            if not line.strip():
                continue
            if line.startswith("##"):
                parts = line[2:].strip().split("\t", 1)
                if len(parts) == 2:
                    meta[parts[0].strip()] = parts[1].strip()
                continue
            if line.startswith("#"):
                header = [_key(cell) for cell in line.split("\t")]
                continue
            if header is None:
                continue
            cells = line.split("\t")
            record = {"template": cells[0].strip()}
            for name, value in zip(header[1:], cells[1:]):
                if name in MAPSTAT_COLUMNS:
                    record[MAPSTAT_COLUMNS[name]] = value.strip()
            rows[record["template"]] = record
    return meta, rows


def merge_hits(sample: str, res: dict[str, dict], mapstat: dict[str, dict],
               overrides=None) -> list[dict]:
    """One row per template in ``.res``, with its mapstat counts and parsed name."""
    hits = []
    for template, record in res.items():
        name = parse_template(template)
        hit = {column: "" for column in HIT_COLUMNS}
        hit.update(record)
        hit.update({key: value for key, value in mapstat.get(template, {}).items()
                    if key != "template"})
        hit.update(sample=sample, template=template, gene=name.gene, variant=name.variant,
                   accession=name.accession, family=family_of(name.gene, overrides))
        hits.append(hit)
    hits.sort(key=lambda hit: (hit["family"], hit["gene"], hit["template"]))
    return hits


def write_hits(hits: list[dict], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=HIT_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(hits)
    return path


def read_hits(path: Path, overrides: list | None = None) -> list[dict]:
    """The arm's ``amr_hits.csv``.

    ``family`` was written when the arm ran. Given the current
    ``family_overrides``, it is recomputed from the gene name, so an override
    added later reaches an existing run through ``bridge amr-report`` without
    re-running KMA, as thresholds do.
    """
    if not Path(path).is_file():
        return []
    with open(path, newline="", encoding="utf-8") as handle:
        hits = list(csv.DictReader(handle))
    if overrides is not None:
        for hit in hits:
            if hit.get("gene"):
                hit["family"] = family_of(hit["gene"], overrides)
    return hits


# ---------------------------------------------------------------------------
# Running
# ---------------------------------------------------------------------------

def count_read_pairs(r1: Path) -> int:
    """FASTQ records in R1 (= read pairs), gzip or plain."""
    opener = gzip.open if str(r1).endswith(".gz") else open
    newlines = 0
    with opener(r1, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 22), b""):
            newlines += block.count(b"\n")
    return newlines // 4


def kma_command(kma: str, r1: Path, r2: Path, db: Path, prefix: Path, threads: int,
                extra_args: str) -> list[str]:
    return [kma, "-ipe", str(r1), str(r2), "-o", str(prefix), "-t_db", str(db),
            "-t", str(max(1, int(threads))), *shlex.split(extra_args)]


def kma_version(kma: str) -> str | None:
    try:
        result = subprocess.run([kma, "-v"], capture_output=True, text=True, check=False)
    except OSError:
        return None
    text = (result.stdout or result.stderr or "").strip()
    return text.splitlines()[0] if text else None


def _check_db(settings: AmrSettings) -> None:
    if settings.db is None:
        raise ValueError("references.amr_db (the KMA index prefix) must be set")
    if not Path(f"{settings.db}.name").is_file():
        raise FileNotFoundError(
            f"no KMA index at {settings.db} ({settings.db}.name missing); "
            f"build it with: bridge build-amr-db --db-dir <resfinder_db>")
    if settings.phenotypes is None or not settings.phenotypes.is_file():
        raise FileNotFoundError(f"references.amr_phenotypes not found: {settings.phenotypes}")


def run_arm(settings: AmrSettings, sample: str, r1: Path, r2: Path, out_dir: Path, *,
            threads: int = 4) -> Path:
    """Run KMA for one sample and write ``amr_hits.csv``. Returns its path."""
    from bridge.tools import configure_logging, close_logging, executable, run

    _check_db(settings)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    handlers = configure_logging(out_dir / "amr.log")
    started = perf_counter()
    try:
        kma = executable(settings.kma)
        overrides = load_family_overrides(settings.family_overrides)
        db = load_phenotypes(settings.phenotypes, overrides)
        LOGGER.info("sample=%s phenotype_entries=%d", sample, len(db), extra={"console": True})

        total = count_read_pairs(Path(r1))
        LOGGER.info("sample=%s total_read_pairs=%d", sample, total, extra={"console": True})

        prefix = out_dir / PREFIX
        command = kma_command(kma, r1, r2, settings.db, prefix, threads, settings.kma_args)
        run(command, "kma")

        res = parse_res(prefix.with_suffix(".res"))
        meta, mapstat = parse_mapstat(prefix.with_suffix(".mapstat"))
        hits = merge_hits(sample, res, mapstat, overrides)
        classify(hits, settings.thresholds)
        for hit in hits:
            annotate(hit, db)
        hits_path = write_hits(hits, out_dir / HITS_NAME)

        manifest = {
            "sample": sample,
            "bridge": __version__,
            "finished": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "runtime_seconds": round(perf_counter() - started, 3),
            "reads": [str(Path(r1).resolve()), str(Path(r2).resolve())],
            "total_read_pairs": total,
            "kma": kma,
            "kma_version": kma_version(kma),
            "command": command,
            "kma_mapstat_meta": meta,
            "db_prefix": str(settings.db),
            "db_manifest": read_manifest(settings.db_manifest),
            "phenotypes": str(settings.phenotypes),
            "thresholds": asdict(settings.thresholds),
            "templates_reported": len(hits),
            "templates_above_floor": sum(1 for hit in hits if hit["tier"] != "below"),
        }
        (out_dir / MANIFEST_NAME).write_text(json.dumps(manifest, indent=2) + "\n",
                                             encoding="utf-8")
        LOGGER.info("sample=%s templates=%d above_floor=%d -> %s", sample, len(hits),
                    manifest["templates_above_floor"], hits_path, extra={"console": True})
        return hits_path
    finally:
        close_logging(handlers)
