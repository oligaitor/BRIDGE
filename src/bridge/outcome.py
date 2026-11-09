"""One outcome per sequenced sample: ``sample_outcomes_all.csv``.

    bridge finalise --out <dir> --run-id RUN          (writes it as step 5)
    bridge outcomes --out <dir> --run-id RUN          (rebuild it alone)
    bridge outcomes --tables <run folder>             (from the *_all.csv tables alone)

The identification table only says ``no_call`` for a sample whose arms called
nothing at all. A sample whose calls were all held, a sample with no reads, a
sample whose arm crashed and a sample that never ran were all left to be read
from absence - and absence cannot tell "negative" from "failed". This table
gives every sample on the run's sample list (``sample_bacteria.list``, else
``sample.list``, else the result folders) exactly one row, decided in this
order:

``negative_control``  the run's blank(s); never reported
``not_run``           on the sample list, but no results folder
``detected``          at least one organism with ``status`` report
``held_for_review``   none reported, but a ``review_override`` call: several
                      16S fallback calls that may be unmodelled background, or
                      may be real. A scientist decides. Split further
                      (``held_split.py``) into
``dna_detected_unconfirmed`` one dominant, high-identity call not common on
                      its run: bacterial DNA, most consistent with that genus
``background_profile`` every call common on the run, recurring across runs, or
                      weak: reported as a negative, flagged
``insufficient_reads`` nothing reported and fewer than ``outcome.min_read_pairs``
                      read pairs (zero-read libraries included): too shallow to
                      call a negative
``analysis_failed``   the identification arms did not finish, or reconcile failed
``not_detected``      everything else: nothing called, or only calls held as
                      background / 16S-only / rMLST alternatives

A ``detected`` or ``held_for_review`` sample below the depth floor keeps its
outcome and is flagged ``low_depth``: a positive on few reads is still a
positive. The statements are in ``config.yaml`` (``outcome.statements``) so
the wording can change without code; ``{organisms}`` is filled in for
``detected``.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path

from bridge.reconcile.baseline import is_control_name

from .held_split import BACKGROUND, HELD, UNCONFIRMED, RunContext, SplitSettings, split_held, \
    split_settings

OUTCOMES_NAME = "sample_outcomes_all.csv"
# The list this pipeline writes. A run folder can also hold another tool's
# ``sample.list``: that is read when ours is absent, but never written to.
SAMPLE_LIST_NAME = "sample_bacteria.list"
OTHER_SAMPLE_LIST = "sample.list"

OUTCOME_COLUMNS = [
    "sample", "run_id", "outcome", "statement", "reported_organisms",
    "n_report", "n_review", "n_review_override", "total_read_pairs",
    "amr_outcome", "reason", "flags",
]

DEFAULT_MIN_READ_PAIRS = 10_000
DEFAULT_STATEMENTS = {
    "negative_control": "Negative control (not reported)",
    "not_run": "Not analysed",
    "detected": "Bacteria detected: {organisms}",
    "dna_detected_unconfirmed": "Bacterial DNA detected, most consistent with {genus}; not "
                                "identified to organism level (16S only, not confirmed by "
                                "rMLST). May represent contamination; interpret with clinical "
                                "context.",
    "held_for_review": "Result pending scientist review",
    "background_profile": "No bacteria present above reporting thresholds",
    "insufficient_reads": "Insufficient sequencing data for a reliable result",
    "analysis_failed": "Result unavailable: analysis did not complete",
    "not_detected": "No bacteria present above reporting thresholds",
}
OUTCOME_ORDER = tuple(DEFAULT_STATEMENTS)

# amr_outcome: what the resistance arm says about the same sample.
AMR_REPORTED = "genes_reported"
AMR_HELD = "genes_held_for_review"
AMR_NONE = "no_gene_above_floor"
AMR_NOT_RUN = "not_run"


@dataclass(frozen=True)
class OutcomeSettings:
    min_read_pairs: int = DEFAULT_MIN_READ_PAIRS
    statements: dict = field(default_factory=lambda: dict(DEFAULT_STATEMENTS))
    split: SplitSettings = field(default_factory=lambda: SplitSettings(enabled=False))


def outcome_settings(cfg, recurrence: Path | None = None) -> OutcomeSettings:
    """The ``outcome:`` block of the config, falling back to the defaults key by key.

    ``recurrence`` overrides ``outcome.held_split.recurrence_table``.
    """
    statements = dict(DEFAULT_STATEMENTS)
    configured = cfg.get("outcome.statements") or {}
    if isinstance(configured, dict):
        statements.update({key: str(value) for key, value in configured.items()
                           if key in DEFAULT_STATEMENTS and value not in (None, "")})
    return OutcomeSettings(
        min_read_pairs=int(cfg.get("outcome.min_read_pairs", DEFAULT_MIN_READ_PAIRS)),
        statements=statements, split=split_settings(cfg, recurrence))


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

def read_sample_list(path: Path) -> list[str]:
    """Sample names from a read list (sample = R1's parent directory), in order.

    Unlike ``reads.parse_read_list`` this never touches the FASTQs, so it works
    on a copy of the list far from the reads.
    """
    names: list[str] = []
    for raw in Path(path).read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        name = Path(line.replace("\\", "/")).parent.name
        if name and name not in names:
            names.append(name)
    return names


def record_sample_list(out_root: Path, pairs) -> tuple[Path, int]:
    """Keep ``<out>/sample_bacteria.list`` as the run's full set of samples (``bridge run``).

    Only ever this file: a ``sample.list`` in the same folder may be someone
    else's and is left alone.

    Written when absent; otherwise only samples not already on it are
    appended, so re-running a subset (``--arms 16s`` on a few samples) never
    shrinks the list finalise uses to spot a sample that never ran.
    """
    path = Path(out_root) / SAMPLE_LIST_NAME
    existing = read_sample_list(path) if path.is_file() else []
    new = [pair for pair in pairs if pair.sample not in existing]
    if new:
        text = path.read_text(encoding="utf-8") if path.is_file() else ""
        if text and not text.endswith("\n"):
            text += "\n"
        path.write_text(text + "".join(f"{pair.r1}\n" for pair in new), encoding="utf-8")
    return path, len(new)


def read_csv(path: Path | None) -> list[dict[str, str]]:
    if path is None or not Path(path).is_file():
        return []
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def read_depth_csv(path: Path) -> dict[str, int]:
    """``sample,total_read_pairs`` (extra columns ignored; blanks skipped)."""
    depth: dict[str, int] = {}
    for row in read_csv(path):
        value = (row.get("total_read_pairs") or "").strip()
        if row.get("sample") and value:
            depth[row["sample"].strip()] = int(float(value))
    return depth


def by_sample(rows: list[dict]) -> dict[str, list[dict]]:
    grouped: dict[str, list[dict]] = {}
    for row in rows:
        grouped.setdefault(row.get("sample", ""), []).append(row)
    return grouped


def display_name(row: dict) -> str:
    """How an organism is named for a reader: a genus-rank call as "Genus sp."."""
    organism = (row.get("organism") or "").strip()
    if row.get("reported_rank") == "genus" and organism and " " not in organism:
        return f"{organism} sp."
    return organism


# ---------------------------------------------------------------------------
# The rule
# ---------------------------------------------------------------------------

def amr_outcome(rows: list[dict]) -> str:
    statuses = {row.get("status") for row in rows}
    if "report" in statuses:
        return AMR_REPORTED
    if "review" in statuses:
        return AMR_HELD
    if "no_call" in statuses:
        return AMR_NONE
    return AMR_NOT_RUN


def _statement(settings: OutcomeSettings, outcome: str, organisms: list[str],
               genus: str = "") -> str:
    template = settings.statements.get(outcome, DEFAULT_STATEMENTS[outcome])
    try:
        return template.format(organisms=", ".join(organisms), genus=genus)
    except (KeyError, IndexError, ValueError):
        return template  # a stray brace in edited wording must not stop a run


def build_outcomes(run_id: str, universe: list[str], *,
                   id_rows: list[dict], amr_rows: list[dict],
                   depth: dict[str, int | None],
                   settings: OutcomeSettings,
                   controls: set[str] | None = None,
                   control_prefixes: tuple[str, ...] = (),
                   present: set[str] | None = None,
                   failures: dict[str, str] | None = None,
                   listed: set[str] | None = None) -> list[dict]:
    """One row per sample in ``universe``.

    ``present`` is the set of samples with a results folder; ``None`` means
    unknown (tables mode), and then nothing can be ``not_run``. ``failures``
    maps a sample to why finalise could not reconcile it; its identification
    rows are then ignored, since any that exist are left over from an earlier
    finalise. ``listed`` is the sample list, to flag samples found in the
    results but not on it.
    """
    controls = set(controls or ())
    failures = dict(failures or {})
    ids = by_sample(id_rows)
    amr = by_sample(amr_rows)
    context = (RunContext.build(run_id, [row for row in id_rows
                                         if row.get("sample") not in failures
                                         and row.get("sample") not in controls],
                                settings.split, control_prefixes)
               if settings.split.enabled else None)
    out = []
    for sample in universe:
        rows = [] if sample in failures else ids.get(sample, [])
        statuses = [row.get("status", "") for row in rows]
        reported = [display_name(row) for row in rows if row.get("status") == "report"]
        reads = depth.get(sample)
        low = reads is not None and reads < settings.min_read_pairs
        flags: list[str] = []
        reason = ""
        genus = ""

        if sample in controls or is_control_name(sample, control_prefixes):
            outcome = "negative_control"
        elif present is not None and sample not in present:
            outcome, reason = "not_run", "no results folder"
        elif reported:
            outcome = "detected"
        elif "review_override" in statuses:
            outcome = "held_for_review"
            reason = "only several 16S fallback calls (review_override)"
            if context is not None:
                outcome, genus, reason = split_held(
                    sample, [row for row in rows if row.get("status") == "review_override"],
                    context, settings.split)
                if outcome == BACKGROUND:
                    flags.append("background_profile")
                    if low:     # a negative in all but name: too shallow to call one
                        outcome = "insufficient_reads"
                        reason = f"{reads} read pairs < {settings.min_read_pairs}; {reason}"
        elif low:
            outcome = "insufficient_reads"
            reason = f"{reads} read pairs < {settings.min_read_pairs}"
        elif sample in failures or not rows:
            outcome = "analysis_failed"
            reason = failures.get(sample) or "no identification rows"
        else:
            outcome = "not_detected"
            if "review" in statuses:
                flags.append("held_calls_not_reported")

        if low and outcome in ("detected", HELD, UNCONFIRMED):
            flags.append("low_depth")
        if reads is None and outcome not in ("negative_control", "not_run"):
            flags.append("depth_unknown")
            if outcome == "analysis_failed":
                reason += "; depth unknown"
        if listed is not None and sample not in listed:
            flags.append("not_on_sample_list")

        out.append({
            "sample": sample, "run_id": run_id, "outcome": outcome,
            "statement": _statement(settings, outcome, reported, genus),
            "reported_organisms": "; ".join(reported),
            "n_report": str(statuses.count("report")),
            "n_review": str(statuses.count("review")),
            "n_review_override": str(statuses.count("review_override")),
            "total_read_pairs": "" if reads is None else str(reads),
            "amr_outcome": amr_outcome(amr.get(sample, [])),
            "reason": reason, "flags": "; ".join(flags),
        })
    return out


def universe_of(listed: list[str], *extra: list[str]) -> list[str]:
    """The sample list first, in its order; then anything else seen, sorted."""
    seen = list(listed)
    others = sorted({name for group in extra for name in group if name} - set(seen))
    return seen + others


def count_outcomes(rows: list[dict]) -> dict[str, int]:
    counts = {name: 0 for name in OUTCOME_ORDER}
    for row in rows:
        counts[row["outcome"]] = counts.get(row["outcome"], 0) + 1
    return {name: n for name, n in counts.items() if n}


def write_rows(rows: list[dict], columns: list[str], path: Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return path


# ---------------------------------------------------------------------------
# A run: from a results tree, or from the run-level tables alone
# ---------------------------------------------------------------------------

@dataclass
class OutcomeResult:
    outcomes_path: Path | None = None
    clinical_path: Path | None = None
    counts: dict[str, int] = field(default_factory=dict)
    sample_list: Path | None = None


def write_run_tables(dest: Path, run_id: str, universe: list[str], *,
                     id_rows: list[dict], amr_rows: list[dict], phenotype_rows: list[dict],
                     depth: dict, settings: OutcomeSettings, **kwargs) -> OutcomeResult:
    """Build and write ``sample_outcomes_all.csv`` and ``clinical_summary_all.csv``."""
    from .clinical import CLINICAL_COLUMNS, CLINICAL_NAME, build_clinical
    outcomes = build_outcomes(run_id, universe, id_rows=id_rows, amr_rows=amr_rows,
                              depth=depth, settings=settings, **kwargs)
    failures = kwargs.get("failures") or {}
    clinical = build_clinical(outcomes,
                              [row for row in id_rows if row.get("sample") not in failures],
                              amr_rows, phenotype_rows)
    dest = Path(dest)
    return OutcomeResult(
        outcomes_path=write_rows(outcomes, OUTCOME_COLUMNS, dest / OUTCOMES_NAME),
        clinical_path=write_rows(clinical, CLINICAL_COLUMNS, dest / CLINICAL_NAME),
        counts=count_outcomes(outcomes))


def sample_list_path(out_root: Path, explicit: Path | None = None) -> Path | None:
    """``--sample-list``, else ``sample_bacteria.list``, else ``sample.list`` (read only)."""
    if explicit is not None:
        return Path(explicit)
    for name in (SAMPLE_LIST_NAME, OTHER_SAMPLE_LIST):
        candidate = Path(out_root) / name
        if candidate.is_file():
            return candidate
    return None


def outcomes_from_tables(folder: Path, settings: OutcomeSettings, *,
                         run_id: str | None = None, sample_list: Path | None = None,
                         depth_csv: Path | None = None, dest: Path | None = None,
                         controls: list[str] | None = None,
                         control_prefixes: tuple[str, ...] = ()) -> OutcomeResult:
    """Tables mode: a run folder holding only the ``*_all.csv`` tables (+ sample.list).

    There are no sample folders to look at, so nothing can be ``not_run`` and
    a missing identification row cannot be told apart from a zero-read
    library unless ``depth_csv`` says how deep each sample was.
    """
    folder = Path(folder)
    id_rows = read_csv(folder / "identification_summary_all.csv")
    amr_rows = read_csv(folder / "amr_summary_all.csv")
    phenotype_rows = read_csv(folder / "predicted_phenotype_all.csv")
    list_path = sample_list_path(folder, sample_list)
    listed = read_sample_list(list_path) if list_path else []
    run_id = run_id or next((row["run_id"] for row in id_rows if row.get("run_id")),
                            folder.name)
    universe = universe_of(listed, [row.get("sample", "") for row in id_rows],
                           [row.get("sample", "") for row in amr_rows])
    result = write_run_tables(
        dest or folder, run_id, universe, id_rows=id_rows, amr_rows=amr_rows,
        phenotype_rows=phenotype_rows,
        depth=read_depth_csv(depth_csv) if depth_csv else {},
        settings=settings, controls=set(controls or ()), control_prefixes=control_prefixes,
        listed=set(listed) if list_path else None)
    result.sample_list = list_path
    return result


ARM_MARKERS = ("rmlst", "16s", "amr", "marker16s", "bacterial_ID", "quant", "species_calls.csv")


def looks_like_sample(folder: Path) -> bool:
    """Whether a results-tree folder holds any arm output.

    A run folder also collects things that are not samples - backups of the
    run tables made by hand - and they must not be scored as failed samples.
    """
    return any((Path(folder) / name).exists() for name in ARM_MARKERS)


def tree_inputs(out_root: Path, controls: set[str]) -> tuple[set[str], dict, dict, list[dict]]:
    """For a results tree: folders present, depth, failures, identification rows.

    Mirrors finalise step 2's skip rules without reconciling anything, so
    ``bridge outcomes --out`` can rebuild the table from what is on disk.
    """
    from bridge.reconcile import baseline
    from bridge.amr.gate import total_read_pairs

    from .finalise import find_arm_outputs, has_amr, sample_dirs

    present, depth, failures, id_rows = set(), {}, {}, []
    for sample_dir in sample_dirs(out_root):
        name = sample_dir.name
        present.add(name)
        depth[name] = sample_depth(sample_dir, baseline.total_reads, total_read_pairs)
        if name in controls:
            continue
        rmlst, s16 = find_arm_outputs(sample_dir)
        summary = sample_dir / "bacterial_ID" / "identification_summary.csv"
        if rmlst is None:
            failures[name] = ("no rmlst/species_calls.csv (arm incomplete)"
                              if s16 is not None or not has_amr(sample_dir)
                              else "AMR arm only (no identification)")
        elif not summary.is_file():
            failures[name] = "not reconciled (no bacterial_ID/identification_summary.csv)"
        else:
            id_rows += read_csv(summary)
    return present, depth, failures, id_rows


def sample_depth(sample_dir: Path, id_total, amr_total) -> int | None:
    """Read pairs: kallisto's ``n_processed``, else the AMR manifest (0 for an empty library)."""
    value = id_total(sample_dir)
    return value if value is not None else amr_total(sample_dir)


def outcomes_from_tree(out_root: Path, run_id: str, settings: OutcomeSettings, *,
                       controls: list[str], control_prefixes: tuple[str, ...] = (),
                       sample_list: Path | None = None,
                       failures: dict[str, str] | None = None,
                       id_rows: list[dict] | None = None,
                       use_amr: bool = True) -> OutcomeResult:
    """Tree mode (finalise step 5, or ``bridge outcomes --out``).

    finalise passes its own ``failures`` and combined ``id_rows``; on its own
    they are read from the tree. ``use_amr`` False ignores an
    ``amr_summary_all.csv`` this finalise did not write (it may be stale).
    """
    out_root = Path(out_root)
    control_set = set(controls)
    present, depth, found_failures, found_rows = tree_inputs(out_root, control_set)
    if failures is None:
        failures = found_failures
    if id_rows is None:
        id_rows = found_rows
    list_path = sample_list_path(out_root, sample_list)
    listed = read_sample_list(list_path) if list_path else []
    # A folder off the list with no arm output is not a sample (a backup folder).
    present = {name for name in present
               if name in listed or looks_like_sample(out_root / name)}
    failures = {name: why for name, why in failures.items() if name in present}
    amr_rows = read_csv(out_root / "amr_summary_all.csv") if use_amr else []
    phenotype_rows = read_csv(out_root / "predicted_phenotype_all.csv") if use_amr else []
    universe = universe_of(listed, sorted(present))
    result = write_run_tables(
        out_root, run_id, universe, id_rows=id_rows, amr_rows=amr_rows,
        phenotype_rows=phenotype_rows, depth=depth, settings=settings,
        controls=control_set, control_prefixes=control_prefixes, present=present,
        failures=failures, listed=set(listed) if list_path else None)
    result.sample_list = list_path
    return result


def describe(result: OutcomeResult) -> list[str]:
    lines = [f"Outcomes   : {result.outcomes_path}"
             + (f"  (samples from {result.sample_list.name})" if result.sample_list
                else "  (no sample list: samples = result folders)")]
    lines += [f"             {name:<18} {n}" for name, n in result.counts.items()]
    lines.append(f"Clinical   : {result.clinical_path}")
    return lines
