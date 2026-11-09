from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable

from .models import FastaRecord


ALLELE_RE = re.compile(r"^(BACT\d{6}_\d+)")


def normalize_allele(value: str) -> str:
    match = ALLELE_RE.match(str(value).strip().lstrip(">"))
    if match is None:
        raise ValueError(f"Cannot normalize allele identifier: {value!r}")
    return match.group(1)


def read_fasta(path: str | Path) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    header: str | None = None
    sequence: list[str] = []
    with Path(path).expanduser().resolve().open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append(_record(header, sequence))
                header, sequence = line[1:].strip(), []
            elif header is None:
                raise ValueError(f"Sequence before FASTA header in {path}")
            else:
                sequence.append(line.upper())
    if header is not None:
        records.append(_record(header, sequence))
    if not records:
        raise ValueError(f"No FASTA records in {path}")
    identifiers = [record.identifier for record in records]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError(f"Duplicate normalized allele identifiers in {path}")
    return records


def _record(header: str, chunks: list[str]) -> FastaRecord:
    return FastaRecord(normalize_allele(header), header, "".join(chunks))


def write_fasta(records: Iterable[FastaRecord], destination: str | Path) -> Path:
    path = Path(destination)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(f">{record.description}\n")
            for start in range(0, len(record.sequence), 80):
                handle.write(record.sequence[start : start + 80] + "\n")
    return path


def write_sequence(identifier: str, sequence: str, destination: str | Path) -> Path:
    return write_fasta([FastaRecord(identifier, identifier, sequence)], destination)
