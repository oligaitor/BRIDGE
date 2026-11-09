"""Read lists and read-pair discovery.

The rules:

* the sample name is the name of the directory holding R1;
* R2 is R1 with the first ``_R1_`` replaced by ``_R2_`` (or ``_R1.`` by
  ``_R2.``) in the BASENAME only, so no parent directory name is load-bearing;
* a list may carry ``#`` comments, blank lines and Windows line endings;
* a sample listed twice is taken once.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable


@dataclass(frozen=True)
class ReadPair:
    sample: str
    r1: Path
    r2: Path

    def size_bytes(self) -> int:
        return self.r1.stat().st_size + self.r2.stat().st_size


@dataclass(frozen=True)
class ReadProblem:
    sample: str
    reason: str
    source: str


def mate_basename(r1_basename: str) -> str | None:
    """The R2 file name for an R1 name, or None when there is no R1 token."""
    if "_R1_" in r1_basename:
        return r1_basename.replace("_R1_", "_R2_", 1)
    if "_R1." in r1_basename:
        return r1_basename.replace("_R1.", "_R2.", 1)
    return None


def pair_from_r1(r1: Path, sample: str | None = None) -> ReadPair | ReadProblem:
    r1 = Path(r1)
    name = sample or r1.parent.name
    mate = mate_basename(r1.name)
    if mate is None:
        return ReadProblem(name, f"no '_R1_' or '_R1.' token in {r1.name}", str(r1))
    r2 = r1.with_name(mate)
    if not r1.is_file():
        return ReadProblem(name, f"R1 not found: {r1}", str(r1))
    if not r2.is_file():
        return ReadProblem(name, f"R2 mate not found: {r2}", str(r1))
    return ReadPair(name, r1, r2)


def _collect(candidates: Iterable[Path],
             name_of: Callable[[Path], str]) -> tuple[list[ReadPair], list[ReadProblem]]:
    pairs: list[ReadPair] = []
    problems: list[ReadProblem] = []
    seen: set[str] = set()
    for r1 in candidates:
        result = pair_from_r1(r1, name_of(r1))
        if isinstance(result, ReadProblem):
            problems.append(result)
            continue
        if result.sample in seen:
            problems.append(ReadProblem(
                result.sample, "listed more than once; the first entry is used",
                str(r1)))
            continue
        seen.add(result.sample)
        pairs.append(result)
    return pairs, problems


def parse_read_list(path: Path) -> tuple[list[ReadPair], list[ReadProblem]]:
    """One absolute R1 path per line; the sample name is the parent directory of R1."""
    lines = []
    for raw in Path(path).read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if line and not line.startswith("#"):
            lines.append(Path(line))
    return _collect(lines, lambda r1: r1.parent.name)


_FLAT_NAME = re.compile(r"^(?P<sample>.+?)(?:_S\d+)?(?:_L\d+)?_R1[._]")


def flat_sample_name(r1_basename: str) -> str:
    """Sample name from a flat-layout R1 file name.

    Strips the Illumina sample-sheet suffix as well as the read token, so
    ``SAMPLE-01_S59_L001_R1_001.fastq.gz`` names ``SAMPLE-01`` and
    ``MIX03_80-20_R1.fastq.gz`` names ``MIX03_80-20``.
    """
    match = _FLAT_NAME.match(r1_basename)
    return match.group("sample") if match else r1_basename.split("_R1", 1)[0]


def discover_reads(reads_dir: Path) -> tuple[list[ReadPair], list[ReadProblem]]:
    """Find pairs under a directory.

    ``<dir>/<SAMPLE>/*_R1*.f*q.gz`` is tried first (sample = directory name);
    otherwise ``<dir>/*_R1*.f*q.gz`` (sample from the file name, see
    :func:`flat_sample_name`).
    """
    reads_dir = Path(reads_dir)
    nested = sorted(reads_dir.glob("*/*_R1*.f*q.gz"))
    if nested:
        return _collect(nested, lambda r1: r1.parent.name)
    flat = sorted(reads_dir.glob("*_R1*.f*q.gz"))
    return _collect(flat, lambda r1: flat_sample_name(r1.name))
