"""What a clinician needs from a run: ``clinical_summary_all.csv``.

One row per reported organism (``status`` report); a sample with none gets a
single row carrying its outcome statement (see ``outcome.py``). Deliberately
short - the full evidence stays in the other tables:

``organism``             as reported ("Genus sp." for a genus-rank call)
``resistance_genes``     reported resistance genes attributed to this organism
                         that are NOT species-typical, grouped by drug class:
                         ``Glycopeptide: vanA | Aminoglycoside: aac(6')-aph(2'')``;
                         "No acquired resistance genes detected" when there are
                         none. A gene whose allele was not resolved enough to
                         name its spectrum is marked "(possible)".
``species_typical_genes`` genes the host table calls intrinsic to this species
                         or genus (blaPDC, fosA and catB7 in P. aeruginosa). They
                         confirm the identification and add no resistance beyond
                         ``expected_resistance``, so they are kept apart from it
``expected_resistance``  EUCAST expected resistant phenotype of the organism
                         (intrinsic, no gene needed), as agents
``unlinked_resistance``  reported genes that could not be tied to one organism
                         (several reported and none chosen, or none reported),
                         in the same format - shown, never silently dropped or
                         pinned on the wrong organism
``statement``            the sample's outcome statement

Grouped by class rather than by agent: ResFinder lists every agent a gene
touches (butirosin, lividomycin...), which buries the one line that matters.
There is no "susceptible" anywhere: the panel captures acquired genes only, so
absence of a gene is not evidence of susceptibility (see bridge.amr.report).
"""

from __future__ import annotations

from .outcome import by_sample, display_name

CLINICAL_NAME = "clinical_summary_all.csv"
CLINICAL_COLUMNS = [
    "sample", "run_id", "outcome", "organism", "resistance_genes",
    "species_typical_genes", "expected_resistance", "unlinked_resistance", "statement",
]
UNLINKED = ("unassigned", "no_organism", "no_id")
NO_ACQUIRED = "No acquired resistance genes detected"
# ResFinder spells a few classes two ways (tet(39) is "Tetracyclines"); one
# class must not be split into two groups on the same line.
CLASS_ALIASES = {"Tetracyclines": "Tetracycline"}


def gene_label(call: dict) -> str:
    gene = (call.get("gene") or "").strip()
    if not gene:
        return ""
    possible = call.get("phenotype_uncertain") == "True" and call.get("host_verdict") != "intrinsic"
    return f"{gene} (possible)" if possible else gene


def by_class(calls: list[dict]) -> str:
    """``Class: gene, gene | Class: gene``, classes in order of first appearance."""
    grouped: dict[str, list[str]] = {}
    for call in calls:
        drug_class = (call.get("drug_class") or "").strip() or "Unclassified"
        drug_class = CLASS_ALIASES.get(drug_class, drug_class)
        label = gene_label(call)
        if label and label not in grouped.setdefault(drug_class, []):
            grouped[drug_class].append(label)
    return " | ".join(f"{name}: {', '.join(genes)}" for name, genes in grouped.items())


def expected_agents(phenotype_rows: list[dict], organism: str) -> str:
    agents = sorted({row["agent"] for row in phenotype_rows
                     if row.get("organism") == organism and row.get("eucast_expected") == "True"
                     and row.get("agent")})
    return ", ".join(agents)


def build_clinical(outcomes: list[dict], id_rows: list[dict], amr_rows: list[dict],
                   phenotype_rows: list[dict]) -> list[dict]:
    ids, amr, pheno = by_sample(id_rows), by_sample(amr_rows), by_sample(phenotype_rows)
    out = []
    for outcome in outcomes:
        sample = outcome["sample"]
        if outcome["outcome"] in ("negative_control", "not_run"):
            continue
        reported_calls = [call for call in amr.get(sample, []) if call.get("status") == "report"]
        shared = {
            "sample": sample, "run_id": outcome["run_id"], "outcome": outcome["outcome"],
            "unlinked_resistance": by_class([call for call in reported_calls
                                             if call.get("attribution") in UNLINKED]),
            "statement": outcome["statement"],
        }
        organisms = ([row for row in ids.get(sample, []) if row.get("status") == "report"]
                     if outcome["outcome"] == "detected" else [])
        if not organisms:
            out.append({**shared, "organism": "", "resistance_genes": "",
                        "species_typical_genes": "", "expected_resistance": ""})
            continue
        for row in organisms:
            name = (row.get("organism") or "").strip()
            own = [call for call in reported_calls
                   if (call.get("host_organism") or "").strip() == name]
            typical = [call for call in own if call.get("host_verdict") == "intrinsic"]
            acquired = [call for call in own if call.get("host_verdict") != "intrinsic"]
            out.append({**shared, "organism": display_name(row),
                        "resistance_genes": by_class(acquired) or NO_ACQUIRED,
                        "species_typical_genes": by_class(typical),
                        "expected_resistance": expected_agents(pheno.get(sample, []), name)})
    return out
