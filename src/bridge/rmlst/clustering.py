from __future__ import annotations

import logging
import subprocess
from collections import defaultdict
from pathlib import Path

import pandas as pd

from .errors import InsufficientMaterial
from .fasta import FastaRecord, write_fasta
from .models import Config
from bridge.tools import executable, log_tool_output, run


class UnionFind:
    def __init__(self, values: list[str]):
        self.parent = {value: value for value in values}

    def find(self, value: str) -> str:
        while self.parent[value] != value:
            self.parent[value] = self.parent[self.parent[value]]
            value = self.parent[value]
        return value

    def union(self, left: str, right: str) -> None:
        left, right = self.find(left), self.find(right)
        if left != right:
            self.parent[max(left, right)] = min(left, right)


def _retryable_mash_failure(
    returncode: int,
    stdout: str | None = None,
    stderr: str | None = None,
) -> bool:
    output = "\n".join(text for text in (stdout, stderr) if text).lower()
    singularity_mount_collision = (
        "container creation failed" in output
        and "failed to mount squashfs filesystem: file exists" in output
    )
    # Mash normally uses small positive exit codes; an otherwise-silent 255
    # comes from a Singularity wrapper whose container failed to start before
    # Mash itself ran.
    return singularity_mount_collision or returncode == 255


def cluster_candidates(
    calls: pd.DataFrame,
    records: dict[str, FastaRecord],
    config: Config,
) -> pd.DataFrame:
    mash = executable(config.mash)
    root = config.tempdir / "mash"
    root.mkdir(parents=True, exist_ok=True)
    output: list[dict[str, object]] = []
    for locus, frame in calls.groupby("locus_id", sort=True):
        names = sorted(set(frame["Template"].astype(str)))
        if missing := set(names) - set(records):
            raise ValueError(f"Candidate alleles absent from FASTA: {sorted(missing)[:5]}")
        distances: dict[tuple[str, str], float] = {}
        union = UnionFind(names)
        if len(names) > 1:
            fasta = write_fasta(
                [FastaRecord(name, name, records[name].sequence) for name in names],
                root / f"{locus}_candidates.fasta",
            )
            prefix = root / f"{locus}_candidates"
            sketch = Path(str(prefix) + ".msh")
            sketch.unlink(missing_ok=True)
            sketch_command = [
                mash, "sketch", "-i", "-p", str(config.threads),
                "-k", str(config.mash_k), "-s", "1000", "-o", str(prefix),
                str(fasta),
            ]
            try:
                run(sketch_command, "mash")
            except subprocess.CalledProcessError as error:
                if not _retryable_mash_failure(
                    error.returncode, error.output, error.stderr,
                ):
                    raise
                logging.getLogger("mash").warning(
                    "transient Mash/container failure exit_code=%d; removing any "
                    "partial sketch and retrying once",
                    error.returncode,
                )
                sketch.unlink(missing_ok=True)
                run(sketch_command, "mash")
            command = [
                mash, "triangle", "-E", "-p", str(config.threads),
                "-d", str(config.mash_threshold), str(sketch),
            ]
            logging.getLogger("mash").info(
                "command=%s", " ".join(map(str, command)),
            )
            for attempt in (1, 2):
                completed = subprocess.run(command, text=True, capture_output=True)
                if not completed.returncode:
                    break
                log_tool_output(
                    "mash", stdout=completed.stdout, stderr=completed.stderr,
                )
                if attempt == 1 and _retryable_mash_failure(
                    completed.returncode, completed.stdout, completed.stderr,
                ):
                    logging.getLogger("mash").warning(
                        "transient Mash/container failure exit_code=%d; "
                        "retrying triangle once",
                        completed.returncode,
                    )
                    continue
                raise subprocess.CalledProcessError(
                    completed.returncode, command,
                    output=completed.stdout, stderr=completed.stderr,
                )
            for line in completed.stderr.replace("\r", "\n").splitlines():
                if line.strip():
                    logging.getLogger("mash").info("%s", line)
            for line in completed.stdout.splitlines():
                fields = line.split()
                if len(fields) < 3:
                    continue
                left, right, text = fields[:3]
                distance = float(text)
                if left == right:
                    continue
                distances[tuple(sorted((left, right)))] = distance
                if distance <= config.mash_threshold:
                    union.union(left, right)
        groups: dict[str, list[str]] = defaultdict(list)
        for name in names:
            groups[union.find(name)].append(name)
        ordered = sorted(tuple(sorted(members)) for members in groups.values())
        for ordinal, members in enumerate(ordered, 1):
            representative = _medoid(members, distances, config.mash_threshold)
            cluster_id = f"{locus}_RC{ordinal:04d}"
            for member in members:
                output.append({
                    "locus_id": locus,
                    "runtime_cluster_id": cluster_id,
                    "representative": representative,
                    "Member_Template": member,
                    "cluster_member_count": len(members),
                })
    result = pd.DataFrame(output)
    if result.empty:
        raise InsufficientMaterial(
            "Runtime clustering produced no clusters",
            stage="runtime_clustering",
        )
    return result


def _medoid(
    members: tuple[str, ...],
    distances: dict[tuple[str, str], float],
    missing: float,
) -> str:
    if len(members) == 1:
        return members[0]
    return min(
        members,
        key=lambda member: (
            sum(
                distances.get(tuple(sorted((member, other))), missing)
                for other in members if other != member
            ),
            member,
        ),
    )


def cluster_abundance(calls: pd.DataFrame, assignments: pd.DataFrame) -> pd.DataFrame:
    mapped = calls.merge(
        assignments,
        left_on=["locus_id", "Template"],
        right_on=["locus_id", "Member_Template"],
        validate="one_to_one",
    )
    rows = []
    for keys, frame in mapped.groupby(
        ["sample_id", "locus_id", "runtime_cluster_id", "representative"],
        sort=True,
    ):
        sample, locus, cluster, representative = keys
        rows.append({
            "sample_id": sample,
            "locus_id": locus,
            "runtime_cluster_id": cluster,
            "representative": representative,
            "est_counts": float(frame["est_counts"].sum()),
            "source_members": ";".join(sorted(frame["Template"].astype(str))),
        })
    return pd.DataFrame(rows)
