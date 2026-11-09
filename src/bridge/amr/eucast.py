"""EUCAST Expected Resistant Phenotypes: what an organism is resistant to by nature.

The table (``data/eucast_expected_resistant.tsv``) lists, per organism group,
the agents EUCAST says the organism is intrinsically resistant to. Two uses:

* a predicted-phenotype row for each such agent, sourced to EUCAST rather than
  to a gene - so a report never presents Klebsiella's ampicillin resistance
  as though a detected gene explained it;
* context for the host table: blaSHV in K. pneumoniae IS that ampicillin
  resistance.

Every row carries ``verified``. The shipped rows were transcribed for this
pipeline and have all been checked against EUCAST ERP v1.2 (``yes``); a row
added later stays ``no`` until it has been. :func:`expected_for` returns the
flag with each match and the report carries it into every row it produces.
Rows are matched by regular expression on the organism name,
case-insensitively, so a group ("Enterobacter cloacae complex") can list its
member species.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

FIELDS = ("organism_group", "taxon_regex", "agents", "source", "verified", "note")


@dataclass(frozen=True)
class ExpectedRow:
    organism_group: str
    pattern: re.Pattern
    agents: tuple[str, ...]
    source: str
    verified: bool
    note: str = ""


def load_expected(path: Path) -> list[ExpectedRow]:
    rows: list[ExpectedRow] = []
    header = None
    path = Path(path)
    with open(path, encoding="utf-8") as handle:
        for lineno, raw in enumerate(handle, 1):
            line = raw.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            cells = [cell.strip() for cell in line.split("\t")]
            if header is None:
                header = cells
                if tuple(header) != FIELDS:
                    raise ValueError(f"{path.name}:{lineno}: header is {header}, "
                                     f"expected {list(FIELDS)}")
                continue
            cells += [""] * (len(FIELDS) - len(cells))
            record = dict(zip(FIELDS, cells))
            if not record["source"]:
                raise ValueError(f"{path.name}:{lineno}: every row needs a source "
                                 f"(EUCAST document, version and table)")
            agents = tuple(agent.strip() for agent in record["agents"].split(";") if agent.strip())
            if not agents:
                raise ValueError(f"{path.name}:{lineno}: no agents listed")
            rows.append(ExpectedRow(
                organism_group=record["organism_group"],
                pattern=re.compile(record["taxon_regex"], re.IGNORECASE),
                agents=agents,
                source=record["source"],
                verified=record["verified"].lower() in ("yes", "true", "y"),
                note=record["note"],
            ))
    if header is None:
        raise ValueError(f"{path.name}: no header row found")
    return rows


def expected_for(organism: str, rows: list[ExpectedRow]) -> list[tuple[str, ExpectedRow]]:
    """Every (agent, row) the organism is expected to resist. Empty for no organism."""
    name = str(organism or "").strip()
    if not name:
        return []
    return [(agent, row) for row in rows if row.pattern.search(name) for agent in row.agents]
