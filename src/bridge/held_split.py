"""Splitting ``held_for_review``: background, unconfirmed bacterial DNA, or truly mixed.

A held sample has no organism confirmed by rMLST and two or more 16S fallback
calls above its run's blank (``review_override``). On clinical runs this can
be a large share of samples, too many to review by hand, and most carry the
same environmental / reagent genera again and again. Each held sample is
re-read from its override calls, in this order:

``dna_detected_unconfirmed``
    one call dominates (>= ``dominant_share_min`` of the override 16S reads,
    >= ``dominant_reads_min`` reads, read identity >= ``min_read_identity``)
    and its genus is not common on this run. Reported as bacterial DNA, most
    consistent with that genus, not identified to organism level.
``background_profile``
    every call is background-compatible: its genus is common on this run
    (>= ``run_prevalence_min`` of the run's samples), or recurs across other
    runs (>= ``cross_prevalence_min`` of the samples on >= ``cross_runs_min``
    other runs, from the recurrence table), or the call itself is weak
    (read identity < ``min_read_identity``). Reported as a negative, flagged.
``held_for_review``
    everything else, e.g. one organism's 16S spread over a genus and its
    family: several strong calls that recur nowhere.

Why each piece is there: a kit's background genera can be rare on one run yet
recur across many, so only cross-run evidence sees them; spurious extra
fallback calls tend to have low read identity, hence the identity floor.
Recurrence cannot be the whole story - common pathogens recur too - so a
dominant call escapes cross-run recurrence; only being common on its own run
(a batch contaminant) makes it background. No genus is named in code: all of
it is evidence from the data, not a hard-coded contaminant list.

The recurrence table (``bridge recurrence``) holds each run's genus
prevalence; a run is never compared with itself.
"""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from bridge.reconcile.baseline import is_control_name

RECURRENCE_COLUMNS = ["run_id", "genus", "n_samples_with", "n_samples", "prevalence"]
UNCONFIRMED = "dna_detected_unconfirmed"
BACKGROUND = "background_profile"
HELD = "held_for_review"


@dataclass(frozen=True)
class SplitSettings:
    enabled: bool = True
    run_prevalence_min: float = 0.30
    cross_prevalence_min: float = 0.10
    cross_runs_min: int = 3
    min_read_identity: float = 0.99
    dominant_share_min: float = 0.80
    dominant_reads_min: float = 50
    recurrence: tuple = field(default_factory=tuple)   # rows of the recurrence table


def _num(value) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


# ---------------------------------------------------------------------------
# Prevalence
# ---------------------------------------------------------------------------

def genera_by_sample(id_rows: list[dict], control_prefixes: tuple[str, ...] = ()) -> dict[str, set]:
    """Genera seen in each non-control sample, in any row (any status)."""
    by: dict[str, set] = defaultdict(set)
    for row in id_rows:
        sample, genus = row.get("sample", ""), (row.get("genus") or "").strip()
        if genus and not is_control_name(sample, control_prefixes):
            by[sample].add(genus)
    return by


def prevalence_rows(run_id: str, id_rows: list[dict],
                    control_prefixes: tuple[str, ...] = ()) -> list[dict]:
    """One run's rows of the recurrence table.

    The denominator is the run's non-control samples with at least one genus
    row - the samples that could have carried it.
    """
    by = genera_by_sample(id_rows, control_prefixes)
    counts = Counter(genus for genera in by.values() for genus in genera)
    n = len(by)
    return [{"run_id": run_id, "genus": genus, "n_samples_with": str(k), "n_samples": str(n),
             "prevalence": f"{k / n:.4f}"} for genus, k in sorted(counts.items())]


def read_recurrence(path: Path | None) -> tuple:
    if path is None or not Path(path).is_file():
        return ()
    with open(path, newline="", encoding="utf-8") as handle:
        return tuple(csv.DictReader(handle))


def recurrent_elsewhere(run_id: str, settings: SplitSettings) -> set[str]:
    """Genera at >= cross_prevalence_min on >= cross_runs_min runs other than this one."""
    runs: dict[str, set] = defaultdict(set)
    for row in settings.recurrence:
        if row.get("run_id") == run_id:
            continue
        if (_num(row.get("prevalence")) or 0) >= settings.cross_prevalence_min:
            runs[row.get("genus", "")].add(row.get("run_id"))
    return {genus for genus, seen in runs.items() if len(seen) >= settings.cross_runs_min}


# ---------------------------------------------------------------------------
# The split
# ---------------------------------------------------------------------------

@dataclass
class RunContext:
    """What a run tells the split: who carries which genus, and what recurs elsewhere."""
    genera: dict[str, set]
    n_samples: int
    counts: Counter
    elsewhere: set[str]

    @classmethod
    def build(cls, run_id: str, id_rows: list[dict], settings: SplitSettings,
              control_prefixes: tuple[str, ...] = ()) -> "RunContext":
        genera = genera_by_sample(id_rows, control_prefixes)
        return cls(genera=genera, n_samples=len(genera),
                   counts=Counter(g for gs in genera.values() for g in gs),
                   elsewhere=recurrent_elsewhere(run_id, settings))

    def run_prevalence(self, sample: str, genus: str) -> float:
        """Share of the run's OTHER samples carrying the genus."""
        own = 1 if genus in self.genera.get(sample, set()) else 0
        others = self.n_samples - (1 if sample in self.genera else 0)
        return (self.counts[genus] - own) / others if others > 0 else 0.0


def split_held(sample: str, override_rows: list[dict], context: RunContext,
               settings: SplitSettings) -> tuple[str, str, str]:
    """(outcome, genus for the statement, reason) for one held sample."""
    if not override_rows:
        return HELD, "", "no review_override rows"
    calls = sorted(override_rows, key=lambda row: -(_num(row.get("16s_raw_reads")) or 0))
    total = sum(_num(row.get("16s_raw_reads")) or 0 for row in calls)
    top = calls[0]
    top_reads = _num(top.get("16s_raw_reads")) or 0
    top_genus = (top.get("genus") or "").strip()

    def identity(row):
        return _num(row.get("16s_read_identity"))

    def weak(row):
        value = identity(row)                    # missing is unknown, not weak
        return value is not None and value < settings.min_read_identity

    def run_common(genus):
        return context.run_prevalence(sample, genus) >= settings.run_prevalence_min

    share = top_reads / total if total else 0.0
    dominant = (share >= settings.dominant_share_min
                and top_reads >= settings.dominant_reads_min
                and identity(top) is not None and identity(top) >= settings.min_read_identity)
    if dominant and top_genus and not run_common(top_genus):
        return (UNCONFIRMED, top_genus,
                f"{top_genus} {share:.0%} of {total:.0f} 16S fallback reads "
                f"(read identity {identity(top):.3f}); not common on this run")

    evidence = []
    for row in calls:
        genus = (row.get("genus") or "").strip()
        if genus and run_common(genus):
            evidence.append(f"{genus}: common on run")
        elif genus in context.elsewhere:
            evidence.append(f"{genus}: recurs on other runs")
        elif weak(row):
            evidence.append(f"{genus or row.get('organism', '')}: weak (read identity "
                            f"{identity(row):.3f})")
        else:
            return HELD, "", (f"{genus} is not background-compatible "
                              f"({_num(row.get('16s_raw_reads')) or 0:.0f} reads)")
    return BACKGROUND, "", "; ".join(evidence)


def split_settings(cfg, recurrence: Path | None = None) -> SplitSettings:
    """``outcome.held_split`` from the config; the table from ``recurrence`` or the config."""
    def get(key, default):
        return cfg.get(f"outcome.held_split.{key}", default)

    path = recurrence
    if path is None:
        value = get("recurrence_table", "")
        if value:
            path = Path(str(value)).expanduser()
            source = getattr(cfg, "source_path", None)
            if not path.is_absolute() and source is not None:
                path = Path(source).parent / path
    return SplitSettings(
        enabled=bool(get("enabled", True)),
        run_prevalence_min=float(get("run_prevalence_min", 0.30)),
        cross_prevalence_min=float(get("cross_prevalence_min", 0.10)),
        cross_runs_min=int(get("cross_runs_min", 3)),
        min_read_identity=float(get("min_read_identity", 0.99)),
        dominant_share_min=float(get("dominant_share_min", 0.80)),
        dominant_reads_min=float(get("dominant_reads_min", 50)),
        recurrence=read_recurrence(path))


# ---------------------------------------------------------------------------
# bridge recurrence
# ---------------------------------------------------------------------------

def build_recurrence(folders: list[Path], out: Path,
                     control_prefixes: tuple[str, ...] = ()) -> tuple[Path, int]:
    """The recurrence table from runs' ``identification_summary_all.csv``.

    Each folder is a run (a results tree or a copy of its tables); its run_id
    is read from the table, else the folder name. Returns (path, runs).
    """
    rows: list[dict] = []
    runs = 0
    for folder in folders:
        table = Path(folder) / "identification_summary_all.csv"
        if not table.is_file():
            continue
        with open(table, newline="", encoding="utf-8") as handle:
            id_rows = list(csv.DictReader(handle))
        run_id = next((row["run_id"] for row in id_rows if row.get("run_id")), Path(folder).name)
        rows += prevalence_rows(run_id, id_rows, control_prefixes)
        runs += 1
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=RECURRENCE_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    return out, runs
