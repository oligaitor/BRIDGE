"""From ``amr_hits.csv`` to a report: gated, attributed, interpreted.

    bridge amr-report --out <results> --run-id RUN1     (also run by bridge finalise)

Per sample, under ``<SAMPLE>/amr_report/``:

``amr_calls.csv``
    One row per gene above the detection floor: the allele and how firmly it
    was resolved, the contamination verdict, which organism it is attributed
    to and on what grounds, and the antibiotics it affects. ``status`` is
    ``report`` or ``review``, with ``flags`` saying why. A sample with nothing
    above the floor gets one ``no_call`` row, so "absent from the table" never
    has to mean both "nothing found" and "never ran".

``predicted_phenotype.csv``
    One row per (organism, agent) with evidence of NON-susceptibility: an
    acquired determinant, or EUCAST's expected resistant phenotype.
    ``attribution`` lists every basis the row rests on: ``eucast_expected``
    (no gene - the reported organism resists the agent by nature) and/or the
    gene attribution tiers below, e.g. ``eucast_expected; sole_organism``.
    There is no "S" anywhere, by design. The current panel captures acquired
    genes only, so absence of a gene is not evidence of susceptibility - and
    EUCAST's own position is that genotype does not yet predict
    susceptibility for most organism-agent pairs.

Run level: ``amr_summary_all.csv``, ``predicted_phenotype_all.csv`` and
``amr_baselines.csv`` beside the identification tables.

Attributing a gene to an organism
---------------------------------
Targeted short reads cannot link a plasmid gene to its host, so attribution
is done in stated tiers and every row says which one it rests on:

Only organisms the identification reports (``status`` report: a reported
basis, clear of the reagent background) are candidates, and only they get a
predicted phenotype.

``sole_organism``   the ID arms reported exactly one organism;
``intrinsic_host``  several organisms, and the curated host table says the
                    gene is intrinsic to exactly one of them (blaACT and
                    Enterobacter);
``unassigned``      several organisms and nothing to choose between them.
                    The candidates are listed with the ratio of gene depth to
                    each organism's rMLST depth - as evidence, not a decision:
                    the resistance and rMLST probes capture at different
                    efficiencies, so the ratio is not calibrated;
``no_organism``     the ID arms ran and reported nothing with ``status`` report
                    - nothing, or only organisms held as review /
                    review_override (``host_note`` names them);
``no_id``           no identification results for this sample (AMR run alone).
"""

from __future__ import annotations

import csv
import logging
from dataclasses import dataclass, field
from pathlib import Path

from bridge.reconcile.reconcile import CALL_BASES
from bridge.reconcile.rpm import baselines_from_rows, rpm_detail, write_baselines

from .calls import _num, collapse_genes
from .database import load_family_overrides, normalise_agent
from .eucast import expected_for, load_expected
from .gate import (
    build_amr_baselines, family_in_control, gate_hit, hits_file, total_read_pairs,
)
from .host import load_host_context, resolve
from .kma import read_hits
from .phenotype import is_uncertain
from .settings import AmrSettings

LOGGER = logging.getLogger("bridge.amr.report")

REPORT_DIR = "amr_report"
CALLS_NAME = "amr_calls.csv"
PHENOTYPE_NAME = "predicted_phenotype.csv"
BASELINE_NAME = "amr_baselines.csv"
SUMMARY_ALL = "amr_summary_all.csv"
PHENOTYPE_ALL = "predicted_phenotype_all.csv"

CALL_COLUMNS = [
    # --- what -----------------------------------------------------------------
    "sample", "run_id", "gene", "allele_level", "family", "drug_class", "status",
    # --- how firmly -----------------------------------------------------------
    "tier", "template", "query_identity", "template_identity", "template_coverage", "depth",
    "fragments", "est_differences",
    # --- contamination gate ---------------------------------------------------
    "rpm", "rpm_r", "abs_r", "gate_band", "gate_mechanism", "gate_reason",
    # --- whose ----------------------------------------------------------------
    "attribution", "host_organism", "host_verdict", "host_citation", "host_note",
    "candidates", "depth_ratio_to_candidates",
    # --- what it does ---------------------------------------------------------
    "antibiotics", "phenotype_uncertain", "possible_antibiotics", "mechanism",
    "phenotype_notes",
    # --- provenance -----------------------------------------------------------
    "flags",
]

PHENOTYPE_COLUMNS = [
    "sample", "run_id", "organism", "attribution", "agent", "prediction",
    "acquired_evidence", "intrinsic_gene_evidence", "possible_evidence",
    "eucast_expected", "eucast_group", "eucast_verified", "notes",
]

# Strongest first; one row per (organism, agent) takes the first that applies.
# intrinsic_gene_evidence is not here: an intrinsic gene only annotates a row
# EUCAST already made, so it never decides a prediction.
PREDICTIONS = (
    ("eucast_expected", "resistant (EUCAST expected phenotype)"),
    ("acquired_evidence", "non-susceptible predicted (acquired determinant)"),
    ("possible_evidence", "possible non-susceptibility (allele not resolved)"),
)


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

def id_summary(sample_dir: Path) -> list[dict] | None:
    """The identification rows for this sample, or None if the ID arms never ran."""
    path = Path(sample_dir) / "bacterial_ID" / "identification_summary.csv"
    if not path.is_file():
        return None
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def reported_organisms(rows: list[dict]) -> list[dict]:
    """The organisms the identification side REPORTS: ``status`` report.

    That is a reported basis (rmlst / 16s_fallback) AND clear of the run's
    reagent background AND not one of several 16S fallback calls. Rows with
    ``status`` review - an organism the gate put at or below the blank's level,
    a demoted near neighbour, a 16S-only genus - or review_override are not
    candidates for a gene and get no predicted phenotype: listing one would
    assert an organism the identification itself did not. A summary without a
    ``status`` column is read on basis alone.
    """
    out = []
    for row in rows:
        if row.get("basis") not in CALL_BASES or not (row.get("organism") or "").strip():
            continue
        status = (row.get("status") or "").strip()
        if status and status != "report":
            continue
        out.append(row)
    return out


# ---------------------------------------------------------------------------
# Attribution
# ---------------------------------------------------------------------------

def attribute(hit: dict, id_rows: list[dict] | None, host_table) -> dict:
    """Which organism a gene belongs to, and the host verdict for it."""
    family = hit.get("family", "")
    out = {"attribution": "", "host_organism": "", "host_verdict": "", "host_citation": "",
           "host_note": "", "candidates": "", "depth_ratio_to_candidates": ""}
    if id_rows is None:
        out.update(attribution="no_id", host_verdict="undetermined",
                   host_note="no identification results for this sample")
        return out
    organisms = reported_organisms(id_rows)
    if not organisms:
        held = [row["organism"].strip() for row in id_rows
                if row.get("basis") in CALL_BASES and (row.get("organism") or "").strip()]
        out.update(attribution="no_organism", host_verdict="undetermined",
                   host_note=("the ID arms reported no organism" if not held else
                              "no organism reported; all held for review (background "
                              "or several 16S fallback calls): "
                              + "; ".join(held) + ")"))
        return out

    verdicts = []
    for row in organisms:
        name = row["organism"].strip()
        verdict, citation, note, _ = resolve(family, name, host_table)
        verdicts.append((name, verdict, citation, note, row))

    depth = _num(hit.get("depth"))
    ratios = []
    for name, _, _, _, row in verdicts:
        organism_depth = _num(row.get("median_projected_depth"))
        ratios.append(f"{name}: {depth / organism_depth:.3g}" if organism_depth > 0
                      else f"{name}: -")
    out["depth_ratio_to_candidates"] = "; ".join(ratios)

    chosen = None
    if len(verdicts) == 1:
        chosen, out["attribution"] = verdicts[0], "sole_organism"
    else:
        intrinsic = [item for item in verdicts if item[1] == "intrinsic"]
        if len(intrinsic) == 1:
            chosen, out["attribution"] = intrinsic[0], "intrinsic_host"
        else:
            out["attribution"] = "unassigned"
            out["candidates"] = "; ".join(f"{name} ({verdict})"
                                          for name, verdict, *_ in verdicts)
            out["host_verdict"] = "undetermined"
            out["host_note"] = (f"{len(verdicts)} organisms reported; "
                                + ("intrinsic to more than one" if intrinsic
                                   else "not intrinsic to any"))
    if chosen is not None:
        name, verdict, citation, note, _ = chosen
        out.update(host_organism=name, host_verdict=verdict, host_citation=citation,
                   host_note=note)
    return out


# ---------------------------------------------------------------------------
# One sample
# ---------------------------------------------------------------------------

def _flags(*groups) -> list[str]:
    flags: list[str] = []
    for group in groups:
        for flag in str(group or "").split(";"):
            flag = flag.strip()
            if flag and flag not in flags:
                flags.append(flag)
    return flags


def call_rows(sample: str, run_id: str, hits: list[dict], total: int | None, baselines,
              id_rows: list[dict] | None, host_table, settings: AmrSettings) -> list[dict]:
    # Re-read under today's thresholds, one call per gene name.
    hits = collapse_genes(hits, settings.thresholds)
    rows = []
    for hit in hits:
        if hit.get("tier") == "below":
            continue
        row = {column: "" for column in CALL_COLUMNS}
        row.update({key: hit.get(key, "") for key in (
            "gene", "allele_level", "family", "drug_class", "tier", "template",
            "query_identity", "template_identity", "template_coverage", "depth", "fragments",
            "est_differences", "antibiotics", "mechanism", "phenotype_notes")})
        row.update(sample=sample, run_id=run_id)
        flags = _flags(hit.get("flags"))

        verdict = gate_hit(hit, total, baselines, run_id, settings.gate)
        if verdict is None:
            flags.append("amr_gate_no_read_total")
        else:
            detail = rpm_detail(verdict)
            row.update(rpm=f"{verdict.rpm:.4g}",
                       rpm_r="" if verdict.rpm_r is None else f"{verdict.rpm_r:.4g}",
                       abs_r="" if verdict.abs_r is None else f"{verdict.abs_r:.4g}",
                       gate_band=detail["band"], gate_mechanism=detail["mechanism"] or "",
                       gate_reason=verdict.reason)
            if verdict.rule == "not_applied":
                flags.append("amr_gate_not_applied")
            elif not verdict.passed:
                flags.append("amr_below_background")
            if family_in_control(hit, baselines, run_id):
                # A relative of this gene was in the blank; this gene was not.
                # Informational: the gate is per gene (see gate.py).
                flags.append("family_in_control")

        row.update(attribute(hit, id_rows, host_table))
        if row["host_verdict"] == "not_expected_in_taxon":
            flags.append("not_expected_in_taxon")

        if is_uncertain(hit):
            row["phenotype_uncertain"] = "True"
            row["possible_antibiotics"] = hit.get("family_possible_antibiotics", "")
            flags.append("phenotype_uncertain")
        else:
            row["phenotype_uncertain"] = "False"

        reasons = []
        if row["tier"] != "confident":
            reasons.append("partial")
        if "amr_below_background" in flags:
            reasons.append("background")
        if "not_expected_in_taxon" in flags:
            reasons.append("not_expected_in_taxon")
        # A nearest allele can't be gated on its own: its name is just KMA's
        # closest template, so reagent reads of one allele in the blank can be
        # named as a different allele in a sample, which then has no blank level
        # to be compared with. An exact allele keeps its own per-allele gate.
        if row["allele_level"] != "exact" and "family_in_control" in flags:
            flags.append("nearest_family_in_control")
            reasons.append("nearest_family_in_control")
        row["status"] = "review" if reasons else "report"
        row["flags"] = "; ".join(flags)
        rows.append(row)

    if not rows:
        row = {column: "" for column in CALL_COLUMNS}
        row.update(sample=sample, run_id=run_id, status="no_call",
                   flags="no resistance gene above the detection floor")
        rows.append(row)
    return rows


def phenotype_rows(sample: str, run_id: str, calls: list[dict], id_rows: list[dict] | None,
                   eucast_rows) -> list[dict]:
    """Aggregate to one row per (organism, agent). Only ``report`` calls contribute."""
    table: dict[tuple[str, str], dict] = {}

    def add(record: dict, column: str, value: str) -> None:
        present = [item for item in record[column].split("; ") if item]
        if value and value not in present:
            present.append(value)
        record[column] = "; ".join(present)

    def entry(organism: str, attribution: str, agent: str) -> dict:
        """The (organism, agent) row; ``attribution`` collects every basis it rests on."""
        key = (organism, normalise_agent(agent))
        if key not in table:
            table[key] = {column: "" for column in PHENOTYPE_COLUMNS}
            table[key].update(sample=sample, run_id=run_id, organism=organism,
                              agent=agent.strip().lower())
        add(table[key], "attribution", attribution)
        return table[key]

    for row in reported_organisms(id_rows or []):
        organism = row["organism"].strip()
        for agent, expected in expected_for(organism, eucast_rows):
            # No gene here: the organism itself was reported, and EUCAST says
            # it resists this agent by nature.
            record = entry(organism, "eucast_expected", agent)
            record["eucast_expected"] = "True"
            add(record, "eucast_group", expected.organism_group)
            record["eucast_verified"] = str(expected.verified)
            if expected.note:
                add(record, "notes", expected.note)

    for call in calls:
        if call.get("status") != "report":
            continue
        organism = call.get("host_organism") or ""
        attribution = call.get("attribution") or ""
        label = call["gene"] + ("" if call.get("allele_level") == "exact" else " (nearest)")
        if call.get("host_verdict") == "intrinsic":
            # A gene the species always carries says nothing the species name
            # does not: its clinical meaning is EUCAST's expected phenotype, not
            # ResFinder's agent list. blaPDC is in every P. aeruginosa, yet wild
            # type is susceptible to ceftazidime and cefepime; fosA is in every
            # P. aeruginosa, yet EUCAST lists no expected fosfomycin resistance.
            # So an intrinsic gene only annotates a row EUCAST already made for
            # this organism (catB7 beside chloramphenicol) and never makes one.
            for agent in _flags(call.get("antibiotics"), call.get("possible_antibiotics")):
                record = table.get((organism, normalise_agent(agent)))
                if record is not None and record["eucast_expected"]:
                    add(record, "intrinsic_gene_evidence", label)
                    add(record, "attribution", attribution)
            continue
        if call.get("phenotype_uncertain") == "True":
            # The named antibiotics are the nearest allele's, and the family
            # disagrees with itself about spectrum: nothing here is a prediction.
            column = "possible_evidence"
            agents = _flags(call.get("antibiotics"), call.get("possible_antibiotics"))
        else:
            column = "acquired_evidence"
            agents = _flags(call.get("antibiotics"))
        note = f"candidates: {call['candidates']}" if call.get("candidates") else ""
        for agent in agents:
            record = entry(organism, attribution, agent)
            add(record, column, label)
            add(record, "notes", note)

    rows = []
    for record in table.values():
        for column, prediction in PREDICTIONS:
            if record[column]:
                record["prediction"] = prediction
                break
        rows.append(record)
    rows.sort(key=lambda r: (r["organism"] == "", r["organism"], r["agent"]))
    return rows


def _write(rows: list[dict], columns: list[str], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return path


# ---------------------------------------------------------------------------
# A run
# ---------------------------------------------------------------------------

@dataclass
class AmrResult:
    baseline_path: Path | None = None
    baseline_rows: int = 0
    reported: list[str] = field(default_factory=list)
    skipped: list[tuple[str, str]] = field(default_factory=list)
    failed: list[tuple[str, str]] = field(default_factory=list)
    summary_path: Path | None = None
    phenotype_path: Path | None = None
    unverified_eucast_rows_used: int = 0


def amr_sample_dirs(out_root: Path) -> list[Path]:
    return sorted(path.parent.parent for path in Path(out_root).glob("*/amr/amr_hits.csv"))


def finalise_amr(settings: AmrSettings, out_root: Path, run_id: str, *,
                 controls: list[str], cohort: bool = True,
                 pooled_rows: list[dict] | None = None) -> AmrResult:
    """Gate, attribute and interpret every sample with an AMR arm result.

    ``pooled_rows`` are ``amr`` rows of a cross-run baseline (run_id ``*``,
    from ``bridge pool-baseline``). They are written beside the run's own
    rows and used where this run's blank never saw the gene - or as the only
    background, where the run has no usable blank.
    """
    out_root = Path(out_root).resolve()
    result = AmrResult()
    host_table = load_host_context(settings.host_context)
    eucast_rows = load_expected(settings.eucast_table)
    # Families from today's overrides, not the ones written when the arm ran.
    overrides = load_family_overrides(settings.family_overrides)
    control_set = set(controls)
    folders = amr_sample_dirs(out_root)
    samples = [folder for folder in folders if folder.name not in control_set]

    baselines = None
    present = [name for name in controls if hits_file(out_root / name).is_file()]
    pooled = list(pooled_rows or [])
    baseline_path = out_root / BASELINE_NAME
    if settings.gate.enabled and (present or pooled):
        rows = (build_amr_baselines(out_root, present, samples, run_id, settings.gate,
                                    settings.thresholds, cohort=cohort, overrides=overrides)
                if present else [])
        rows += pooled
        write_baselines(rows, baseline_path)
        result.baseline_path, result.baseline_rows = baseline_path, len(rows)
        baselines = baselines_from_rows(rows)
    elif settings.gate.enabled:
        print("AMR WARNING: no negative control has an AMR arm result; the resistance "
              "gate is skipped and every call is flagged amr_gate_not_applied.")

    all_calls: list[dict] = []
    all_pheno: list[dict] = []
    for folder in samples:
        name = folder.name
        try:
            hits = read_hits(hits_file(folder), overrides)
            id_rows = id_summary(folder)
            calls = call_rows(name, run_id, hits, total_read_pairs(folder), baselines,
                              id_rows, host_table, settings)
            pheno = phenotype_rows(name, run_id, calls, id_rows, eucast_rows)
        except Exception as error:  # noqa: BLE001 - one sample must not stop the run
            result.failed.append((name, f"{type(error).__name__}: {error}"))
            continue
        report_dir = folder / REPORT_DIR
        _write(calls, CALL_COLUMNS, report_dir / CALLS_NAME)
        _write(pheno, PHENOTYPE_COLUMNS, report_dir / PHENOTYPE_NAME)
        result.unverified_eucast_rows_used += sum(
            1 for row in pheno if row["eucast_expected"] and row["eucast_verified"] != "True")
        all_calls += calls
        all_pheno += pheno
        result.reported.append(name)
    for name in controls:
        if name in {folder.name for folder in folders}:
            result.skipped.append((name, "negative control"))

    result.summary_path = _write(all_calls, CALL_COLUMNS, out_root / SUMMARY_ALL)
    result.phenotype_path = _write(all_pheno, PHENOTYPE_COLUMNS, out_root / PHENOTYPE_ALL)
    return result


def describe(result: AmrResult, out_root: Path) -> str:
    lines = ["", f"AMR samples: {len(result.reported)}"]
    if result.baseline_path:
        lines.append(f"AMR baseline: {result.baseline_path}  ({result.baseline_rows} rows)")
    if result.failed:
        lines.append(f"AMR FAILED : {len(result.failed)}")
        lines += [f"             {name}: {why}" for name, why in result.failed]
    if result.unverified_eucast_rows_used:
        lines.append(f"AMR NOTE   : {result.unverified_eucast_rows_used} predicted-phenotype "
                     f"row(s) rest on EUCAST rows not yet verified (eucast_verified=False)")
    lines += [
        f"AMR results: {out_root}/<SAMPLE>/{REPORT_DIR}/{CALLS_NAME}",
        f"             {out_root}/{SUMMARY_ALL}         <- every gene call in the run",
        f"             {out_root}/{PHENOTYPE_ALL}  <- organism x agent",
    ]
    return "\n".join(lines)
