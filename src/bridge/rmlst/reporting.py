from __future__ import annotations

import re
import shutil
from pathlib import Path

import pandas as pd

from .fasta import write_fasta
from .models import FastaRecord


ACCEPTED_CONSENSUS = {
    "complete_exact_match", "unique_partial_match", "ambiguous_partial_match",
}
FASTA_MANIFEST_COLUMNS = [
    "sample_id", "species", "call_method", "allele_count", "locus_count",
    "alleles", "fasta_path",
]


def write_direct_species_alleles(
    species_calls: pd.DataFrame,
    consensus_groups: pd.DataFrame,
    output: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Write supporting read-derived consensuses for each reported direct call."""
    calls = species_calls.copy()
    calls["allele_fasta"] = ""
    calls["allele_fasta_allele_count"] = 0
    root = output / "species_allele_fastas"
    if root.exists():
        shutil.rmtree(root)
    root.mkdir(parents=True, exist_ok=True)

    direct = calls[
        calls["reported"].map(_truthy)
        & calls["call_method"].ne("em_fallback")
    ]
    accepted = consensus_groups[
        consensus_groups["consensus_status"].isin(ACCEPTED_CONSENSUS)
    ]
    manifest: list[dict[str, object]] = []
    used_names: set[str] = set()
    for index, call in direct.iterrows():
        species = str(call["species"])
        consensus_records: list[tuple[str, str, str, str, str]] = []
        for group in accepted.itertuples(index=False):
            supported_species = {
                value for value in str(group.compatible_species).split(";") if value
            }
            if species not in supported_species:
                continue
            sequence = str(group.consensus)
            if not sequence or sequence == "nan":
                continue
            consensus_records.append((
                str(group.locus_id),
                str(group.consensus_group_id),
                str(group.consensus_status),
                str(group.compatible_alleles),
                sequence,
            ))
        consensus_records.sort(key=lambda value: (value[0], value[1]))
        slug = _unique_slug(species, used_names)
        fasta = root / f"{slug}.fasta"
        write_fasta([
            FastaRecord(
                group_id,
                (
                    f"{group_id} locus={locus} species={species} status={status} "
                    f"{_database_match_description(status, compatible_alleles)}"
                ),
                sequence,
            )
            for locus, group_id, status, compatible_alleles, sequence
            in consensus_records
        ], fasta)
        resolved = str(fasta.resolve())
        calls.at[index, "allele_fasta"] = resolved
        calls.at[index, "allele_fasta_allele_count"] = len(consensus_records)
        manifest.append({
            "sample_id": str(call["sample_id"]),
            "species": species,
            "call_method": str(call["call_method"]),
            "allele_count": len(consensus_records),
            "locus_count": len({record[0] for record in consensus_records}),
            "alleles": ";".join(record[1] for record in consensus_records),
            "fasta_path": resolved,
        })
    return calls, pd.DataFrame.from_records(
        manifest, columns=FASTA_MANIFEST_COLUMNS
    )


def _truthy(value: object) -> bool:
    return value is True or str(value).strip().casefold() == "true"


def _unique_slug(species: str, used: set[str]) -> str:
    base = re.sub(r"[^A-Za-z0-9._-]+", "_", species).strip("._-") or "species"
    candidate = base
    ordinal = 2
    while candidate in used:
        candidate = f"{base}_{ordinal}"
        ordinal += 1
    used.add(candidate)
    return candidate


def _database_match_description(status: str, value: str) -> str:
    alleles = [allele for allele in str(value).split(";") if allele and allele != "nan"]
    match_type = "exact" if status == "complete_exact_match" else "partial"
    if len(alleles) == 1:
        return f"database_match={match_type} allele={alleles[0]}"
    return f"database_match={match_type} allele_count={len(alleles)}"
