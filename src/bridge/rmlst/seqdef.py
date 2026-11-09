from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import pandas as pd


LOCUS_RE = re.compile(r"^BACT\d{6}$")
NO_SUBSPECIES_LABEL = "[no subspecies annotation]"
SUBSPECIES_SUFFIX_RE = re.compile(
    rf"(?:\s+subsp\.\s+.+|\s+{re.escape(NO_SUBSPECIES_LABEL)})$"
)


def taxon_label(species: str, subspecies: str, distinguish_subspecies: bool) -> str:
    species = str(species).strip()
    if not distinguish_subspecies:
        return species
    subspecies = str(subspecies).strip()
    if not subspecies:
        return f"{species} {NO_SUBSPECIES_LABEL}"
    return f"{species} subsp. {subspecies}"


def base_species_label(value: str) -> str:
    """Remove the display-only subspecies qualifier from a taxon label."""
    return SUBSPECIES_SUFFIX_RE.sub("", str(value).strip())


def allele_species(
    seqdef: str | Path,
    wanted: Iterable[str] | None = None,
    distinguish_subspecies: bool = False,
) -> dict[str, set[str]]:
    selected = None if wanted is None else set(wanted)
    result: dict[str, set[str]] = defaultdict(set)
    with Path(seqdef).expanduser().resolve().open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "species" not in reader.fieldnames:
            raise ValueError("Sequence definition requires a species column")
        if distinguish_subspecies and "subspecies" not in reader.fieldnames:
            raise ValueError(
                "Sequence definition requires a subspecies column when "
                "--distinguish-subspecies is used"
            )
        loci = [name for name in reader.fieldnames if LOCUS_RE.fullmatch(name)]
        for row in reader:
            species = taxon_label(
                row["species"], row.get("subspecies", ""), distinguish_subspecies
            )
            for locus in loci:
                number = row[locus].strip()
                if not number or number == "N":
                    continue
                allele = f"{locus}_{number}"
                if selected is None or allele in selected:
                    result[allele].add(species)
    return dict(result)


def remove_mixed_placeholders(species: set[str]) -> set[str]:
    definitive = {
        value for value in species
        if not base_species_label(value).lower().endswith(" sp.")
    }
    return definitive or set(species)


def add_rst_matches(
    species_calls: pd.DataFrame,
    consensus_groups: pd.DataFrame,
    seqdef: str | Path,
    limit: int = 10,
    distinguish_subspecies: bool = False,
) -> pd.DataFrame:
    """Append rSTs matching any species-assigned allele option at each locus."""
    result = species_calls.copy()
    result["matching_rSTs"] = ""
    result["matching_rST_count"] = 0
    if result.empty or consensus_groups.empty:
        return result

    accepted = consensus_groups[consensus_groups["consensus_status"].isin({
        "complete_exact_match", "unique_partial_match", "ambiguous_partial_match",
    })]
    constraints: dict[str, dict[str, set[str]]] = defaultdict(
        lambda: defaultdict(set)
    )
    for row in accepted.itertuples(index=False):
        locus = str(row.locus_id)
        allele_numbers = {
            allele.split("_", 1)[1]
            for allele in str(row.compatible_alleles).split(";")
            if allele.startswith(f"{locus}_")
        }
        if not allele_numbers:
            continue
        for species in str(row.compatible_species).split(";"):
            if species:
                constraints[species][locus].update(allele_numbers)
    if not constraints:
        return result

    loci = sorted({locus for values in constraints.values() for locus in values})
    path = Path(seqdef).expanduser().resolve()
    header = pd.read_csv(path, sep="\t", nrows=0).columns
    required = ["rST", "species", *[locus for locus in loci if locus in header]]
    if "rST" not in header or "species" not in header:
        raise ValueError("Sequence definition requires rST and species columns")
    if distinguish_subspecies:
        if "subspecies" not in header:
            raise ValueError(
                "Sequence definition requires a subspecies column when "
                "--distinguish-subspecies is used"
            )
        required.append("subspecies")
    profiles = pd.read_csv(
        path, sep="\t", usecols=required, dtype=str, keep_default_na=False,
    )
    if distinguish_subspecies:
        profiles["taxon"] = [
            taxon_label(species, subspecies, distinguish_subspecies)
            for species, subspecies in zip(
                profiles["species"],
                profiles.get("subspecies", pd.Series("", index=profiles.index)),
            )
        ]
    else:
        # Without subspecies mode taxon_label is str(species).strip(); the
        # vectorised form avoids a per-row Python pass over the whole table.
        profiles["taxon"] = profiles["species"].astype(str).str.strip()
    wanted_species = set(result["species"].astype(str)) & set(constraints)
    profiles = profiles[profiles["taxon"].isin(wanted_species)]

    for index, call in result.iterrows():
        species = str(call["species"])
        evidence = constraints.get(species, {})
        species_profiles = profiles[profiles["taxon"].eq(species)]
        candidates = species_profiles
        for locus, allowed in evidence.items():
            if locus in candidates.columns:
                # A consensus group can contain alleles from several species.
                # Constrain this species only with allele options that occur in
                # at least one of its seqdef profiles at this locus.
                species_options = set(species_profiles[locus]) - {"", "N"}
                assigned_options = allowed & species_options
                if not assigned_options:
                    continue
                candidates = candidates[candidates[locus].isin(assigned_options)]
            else:
                candidates = candidates.iloc[0:0]
            if candidates.empty:
                break
        identifiers = sorted(
            set(candidates["rST"].astype(str)) - {"", "N"}, key=_rst_sort_key
        )
        result.at[index, "matching_rSTs"] = ";".join(identifiers[:limit])
        result.at[index, "matching_rST_count"] = len(identifiers)
    return result


def _rst_sort_key(value: str) -> tuple[int, int | str]:
    return (0, int(value)) if value.isdigit() else (1, value)
