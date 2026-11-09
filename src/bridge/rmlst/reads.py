from __future__ import annotations

import gzip
import logging
import subprocess
from collections import defaultdict
from contextlib import ExitStack
from itertools import zip_longest
from pathlib import Path
from typing import Iterator, TextIO

import pandas as pd

from .fasta import normalize_allele
from .models import Config
from bridge.tools import executable


def _read_name(value: str) -> str:
    name = value.strip().lstrip("@").split()[0]
    return name[:-2] if name.endswith(("/1", "/2")) else name


def _fastq(path: Path) -> Iterator[tuple[str, tuple[str, str, str, str]]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        while header := handle.readline():
            sequence, plus, quality = handle.readline(), handle.readline(), handle.readline()
            if not quality:
                raise ValueError(f"Truncated FASTQ: {path}")
            yield _read_name(header), (header, sequence, plus, quality)


def assign_and_extract_reads(
    pseudobam: Path,
    assignments: pd.DataFrame,
    config: Config,
) -> tuple[dict[str, tuple[Path, Path]], pd.DataFrame]:
    samtools = executable(config.samtools)
    target_cluster = assignments.set_index("Member_Template")["runtime_cluster_id"].to_dict()
    command = [samtools, "view", str(pseudobam)]
    process = subprocess.Popen(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert process.stdout is not None and process.stderr is not None
    qname_clusters: dict[str, set[str]] = defaultdict(set)
    qname_unknown: set[str] = set()
    for line in process.stdout:
        fields = line.split("\t", 4)
        if len(fields) < 4:
            continue
        qname, flag_text, reference = fields[:3]
        if int(flag_text) & 4 or reference == "*":
            continue
        try:
            allele = normalize_allele(reference)
        except ValueError:
            qname_unknown.add(_read_name(qname))
            continue
        name = _read_name(qname)
        cluster = target_cluster.get(allele)
        if cluster is None:
            qname_unknown.add(name)
        else:
            qname_clusters[name].add(cluster)
    stderr = process.stderr.read()
    returncode = process.wait()
    for line in stderr.replace("\r", "\n").splitlines():
        if line.strip():
            logging.getLogger("samtools").info("%s", line)
    if returncode:
        raise subprocess.CalledProcessError(returncode, command)

    selected: dict[str, set[str]] = {}
    stats: dict[str, dict[str, int]] = defaultdict(
        lambda: {"assigned_fragments": 0, "exclusive_fragments": 0, "ambiguous_fragments": 0}
    )
    for name, clusters in qname_clusters.items():
        exclusive = len(clusters) == 1 and name not in qname_unknown
        if config.read_assignment == "exclusive" and not exclusive:
            continue
        selected[name] = clusters
        for cluster in clusters:
            stats[cluster]["assigned_fragments"] += 1
            stats[cluster]["exclusive_fragments" if exclusive else "ambiguous_fragments"] += 1

    root = config.tempdir / "cluster_reads"
    root.mkdir(parents=True, exist_ok=True)
    paths = {
        cluster: (root / f"{cluster}_1.fastq", root / f"{cluster}_2.fastq")
        for cluster in sorted(stats)
    }
    with ExitStack() as stack:
        handles: dict[str, tuple[TextIO, TextIO]] = {
            cluster: (
                stack.enter_context(pair[0].open("w", encoding="utf-8")),
                stack.enter_context(pair[1].open("w", encoding="utf-8")),
            )
            for cluster, pair in paths.items()
        }
        first, second = _fastq(config.reads[0]), _fastq(config.reads[1])
        for left, right in zip_longest(first, second):
            if left is None or right is None:
                raise ValueError("Paired FASTQ files contain different numbers of records")
            if left[0] != right[0]:
                raise ValueError(f"FASTQ mates are not synchronized: {left[0]} != {right[0]}")
            for cluster in selected.get(left[0], set()):
                out1, out2 = handles[cluster]
                out1.writelines(left[1])
                out2.writelines(right[1])

    rows = [{"runtime_cluster_id": cluster, **values} for cluster, values in sorted(stats.items())]
    logging.getLogger("bridge.rmlst").info(
        "read_assignment_status=completed fragments=%d clusters=%d",
        len(selected), len(rows),
    )
    return paths, pd.DataFrame(rows)
