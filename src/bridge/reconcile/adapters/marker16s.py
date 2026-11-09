"""``marker16s_calls.csv`` -> ``list[BlastStrengthResult]``.

Thin by design. The strength decision was made in
``bridge.marker16s.calls``, where the abundance data is; this reads the row
and hands it to :func:`..rank_marker16s.rank_marker16s` to translate. Nothing
re-derives the verdict, so there is only ever one copy of the rule.
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path

from ..config import Config
from ..rank_marker16s import Marker16SRow, rank_marker16s
from ..strength import BlastStrength, BlastStrengthResult

REQUIRED_COLUMNS = {"sample_id", "taxon", "genus", "strength"}


def _number(row: dict, key: str, cast, default):
    try:
        value = row.get(key, "")
        return cast(value) if str(value).strip() != "" else default
    except (TypeError, ValueError):
        return default


def load_marker16s_calls(path: Path) -> list[Marker16SRow]:
    """Parse ``marker16s_calls.csv``."""
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = REQUIRED_COLUMNS - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"{path} is not a marker16s_calls.csv: missing {sorted(missing)}"
            )
        return [
            Marker16SRow(
                sample_id=(row.get("sample_id") or "").strip(),
                taxon=(row.get("taxon") or "").strip(),
                genus=(row.get("genus") or "").strip(),
                species=(row.get("species") or "").strip(),
                rank=(row.get("rank") or "").strip(),
                strength=(row.get("strength") or "").strip(),
                abundance_frac=_number(row, "abundance_frac", float, 0.0),
                genus_frac=_number(row, "genus_frac", float, 0.0),
                reads_assigned=_number(row, "reads_assigned", float, 0.0),
                raw_reads=_number(row, "raw_reads", float, 0.0),
                copy_number=_number(row, "copy_number", float, 0.0),
                within_genus_share=_number(row, "within_genus_share", float, 0.0),
                margin=_number(row, "margin", float, 0.0),
                n_references=_number(row, "n_references", int, 0),
                competing_species=(row.get("competing_species") or "").strip(),
                reason=(row.get("reason") or "").strip(),
                # Alignment metrics; absent columns and blank cells are
                # None, not 0 - see Marker16SRow.
                n_aligned_reads=_number(row, "n_aligned_reads", int, None),
                mean_read_identity=_number(row, "mean_read_identity", float, None),
                ref_breadth=_number(row, "ref_breadth", float, None),
                mean_depth=_number(row, "mean_depth", float, None),
                consensus_length=_number(row, "consensus_length", int, None),
                consensus_pident=_number(row, "consensus_pident", float, None),
                consensus_qcov=_number(row, "consensus_qcov", float, None),
                consensus_top_hit=(row.get("consensus_top_hit") or "").strip(),
                alignment_status=(row.get("alignment_status") or "").strip(),
            )
            for row in reader
        ]


def marker16s_results(
    rows: list[Marker16SRow],
    cfg: Config,
    *,
    include_silent: bool = False,
) -> list[BlastStrengthResult]:
    """Build the reconciler's 16S side.

    ``silent_weak`` genera are dropped by default. A genus below the reporting
    floor is the arm declining to comment, and the reconciler already treats an
    absent arm as "no opinion" - passing thousands of trace genera through would
    add noise without adding information. Set ``include_silent`` when
    diagnosing why an expected organism was not reported.
    """
    logger = logging.getLogger("bridge.reconcile.adapters.marker16s")
    results = []
    for index, row in enumerate(rows):
        # The index is this row's position in the CSV, carried so the report can
        # fetch this organism's own read counts back out of the arm file.
        result = rank_marker16s(row, cfg, row_index=index)
        if result.strength is BlastStrength.SILENT_WEAK and not include_silent:
            continue
        results.append(result)

    logger.info(
        "adapter_status=complete rows=%d results=%d species_strong=%d genus_only=%d",
        len(rows), len(results),
        sum(1 for r in results if r.strength is BlastStrength.SPECIES_STRONG),
        sum(1 for r in results if r.strength is BlastStrength.GENUS_ONLY),
    )
    return results
