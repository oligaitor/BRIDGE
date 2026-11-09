"""Is a resistance gene above what the run's negative control carried?

The identification arms' gate (``bridge.reconcile.rpm``), reused unchanged: RPM-r
for reagent background, abs-r for cross-contamination from a strong positive
on the same run, the same 10x convention and the same guards, with the gene
(``blaTEM-1B``) in the place of the organism.

Gene, not family
----------------
The identification side falls back from species to genus, because a shallow
blank resolves an organism less finely than a deep sample. The same fallback
is wrong for resistance genes. A blank can carry vector-type sequence
(aph(3')-Ia, blaTEM-1A) at a large share of its reads, and a family fallback
turns that into a bar no sample can clear - an exact blaTEM-1B, or
P. aeruginosa's own aph(3')-IIb at tens of thousands of fragments, would be
called background against it. Those are different
sequences, and KMA (``-1t1``) puts a read on the template it fits best in the
sample exactly as in the blank, so contaminating vector reads in a sample land
on the vector allele and are gated there. A different allele in the sample is
a different molecule. So the gate is per gene. That a relative of the gene was
in the blank is still worth a reader's eye: it is reported as the flag
``family_in_control``, which changes no arithmetic (the same move as
``arms_with_background`` on the identification side).

The denominator is ``total_read_pairs`` from ``amr/amr_manifest.json``,
counted by the arm itself in the sample and the control alike; the numerator
is KMA's ``fragments`` (read pairs for paired input).
"""

from __future__ import annotations

import json
from pathlib import Path

from bridge.reconcile.rpm import (
    CohortEvidence, RpmResult, baselines_from_rows, build_baselines, evaluate,
    lookup_baseline, ungated,
)

from .calls import _num, collapse_genes
from .kma import HITS_NAME, MANIFEST_NAME, read_hits
from .settings import GateSettings, Thresholds

ARM = "amr"


def hits_file(sample_dir: Path) -> Path:
    return Path(sample_dir) / "amr" / HITS_NAME


def total_read_pairs(sample_dir: Path) -> int | None:
    path = Path(sample_dir) / "amr" / MANIFEST_NAME
    if not path.is_file():
        return None
    try:
        value = json.loads(path.read_text(encoding="utf-8")).get("total_read_pairs")
        return int(value) if value is not None else None
    except (ValueError, TypeError, OSError):
        return None


def control_profile(sample_dir: Path, overrides: list | None = None) -> dict[str, float]:
    """Fragments per gene, and per family, in a control. Every row counts.

    Like ``baseline.rmlst_profile``: a baseline is a picture of what the
    reagents carry, so a gene below the detection floor in the blank still
    counts. Family totals are written to the baseline file for the reader and
    for the ``family_in_control`` flag; the gate itself looks genes up only.
    """
    genes: dict[str, float] = {}
    families: dict[str, float] = {}
    for hit in read_hits(hits_file(sample_dir), overrides):
        fragments = _num(hit.get("fragments"))
        if fragments <= 0:
            continue
        gene, family = hit.get("gene", ""), hit.get("family", "")
        if gene:
            genes[gene] = genes.get(gene, 0.0) + fragments
        if family:
            families[family] = families.get(family, 0.0) + fragments
    profile = dict(genes)
    for family, fragments in families.items():
        profile.setdefault(family, fragments)
    return profile


def _evaluate(gene: str, reads: float, total: int, baselines, run_id: str,
              gate: GateSettings) -> RpmResult:
    return evaluate(gene, ARM, reads, total,
                    lookup_baseline(baselines, run_id, ARM, gene),
                    rpm_r_min=gate.rpm_r_min, rpm_min=gate.rpm_min,
                    min_ntc_rpm=gate.min_ntc_rpm, abs_r_min=gate.abs_r_min,
                    min_sample_reads=gate.min_sample_reads)


def family_in_control(hit: dict, baselines, run_id: str) -> bool:
    """The blank held this gene's family but not this gene."""
    if baselines is None:
        return False
    gene, family = hit.get("gene", ""), hit.get("family", "")
    if not family or family == gene:
        return False
    return (lookup_baseline(baselines, run_id, ARM, gene) is None
            and lookup_baseline(baselines, run_id, ARM, family) is not None)


def build_amr_baselines(out_root: Path, controls: list[str], samples: list[Path],
                        run_id: str, gate: GateSettings, thresholds: Thresholds, *,
                        cohort: bool = True,
                        overrides: list | None = None) -> list[dict[str, object]]:
    """Baseline rows for the resistance genes, same two passes as the ID side.

    First reagent-only rows; then each gene in the blank that a sample on the
    run credibly carries (confident, and clear of the reagent-only baseline)
    is attributed to cross-contamination and gated on absolute reads.
    """
    profiles: dict[str, dict[str, tuple[float, int]]] = {}
    for name in controls:
        folder = Path(out_root) / name
        total = total_read_pairs(folder)
        if total is None:
            continue
        profiles[name] = {gene: (reads, total)
                          for gene, reads in control_profile(folder, overrides).items()}
    if not profiles:
        return []

    def rows_with(evidence: CohortEvidence | None):
        return build_baselines(profiles, run_id, ARM, min_control_reads=gate.min_control_reads,
                               cohort=evidence, max_ntc_source_frac=gate.max_ntc_source_frac)

    reagent_rows = rows_with(None)
    if not cohort or not reagent_rows:
        return reagent_rows
    reagent = baselines_from_rows(reagent_rows)
    genes: dict[str, tuple[str, float]] = {}
    for sample_dir in samples:
        total = total_read_pairs(sample_dir)
        if not total:
            continue
        for hit in collapse_genes(read_hits(hits_file(sample_dir), overrides), thresholds):
            if hit.get("tier") != "confident":
                continue
            gene = hit.get("gene", "")
            reads = _num(hit.get("fragments"))
            if not gene or reads <= 0:
                continue
            if not _evaluate(gene, reads, total, reagent, run_id, gate).passed:
                continue
            if gene not in genes or reads > genes[gene][1]:
                genes[gene] = (Path(sample_dir).name, reads)
    return rows_with(CohortEvidence(genes, {}))


def gate_hit(hit: dict, total: int | None, baselines, run_id: str,
             gate: GateSettings) -> RpmResult | None:
    """One hit's verdict; ``ungated`` without a baseline; None without a denominator."""
    if total is None or total <= 0:
        return None
    gene = hit.get("gene", "")
    reads = _num(hit.get("fragments"))
    if baselines is None or not gate.enabled:
        return ungated(gene, ARM, reads, total)
    return _evaluate(gene, reads, total, baselines, run_id, gate)
