"""RPM-r: is this organism above its own reagent background?

A coverage gate ("is the reference evenly covered") has no equivalent in a
kallisto-based pipeline. This asks a different question - "is there more of
this than in a blank" - but catches the same failure: a confident-looking call
resting on almost no material.

Why a ratio rather than a raw count
-----------------------------------
Reads-per-million normalises for how deeply a sample was sequenced, but not for
reagent background. Extraction kits and PCR reagents carry a reproducible,
organism-specific bacterial signature (the "kitome"), so a fixed read threshold
is simultaneously too strict for organisms that never appear in reagents and too
lenient for those that always do. The clinical-metagenomics literature settles
this by dividing by the organism's own level in a negative control:

    RPM-r = RPM_sample / RPM_negative_control

with a reporting threshold around 10 for organisms that have background, and a
small absolute RPM floor for those that have none.

Two different control problems, two different guards
----------------------------------------------------
These are easy to conflate and need opposite fixes.

*A very clean control* gives an organism a near-zero RPM_NTC. Dividing by it
makes the ratio enormous and everything passes - false positives. The guard is a
**floor on the denominator**, :data:`MIN_NTC_RPM`; published assays use 127.

*A very shallow control* does the reverse. A negative with 79 total reads makes
a single stray read worth 12,658 RPM, an enormous denominator, so genuine
organisms fail - false negatives. A floor cannot help here, because the value is
already far above it. The guard is to **refuse to build a baseline from a
control that was not sequenced deeply enough to measure anything**, which is
:data:`MIN_CONTROL_READS` in :func:`build_baselines`. A control below it is
excluded and the pooled baseline is used instead.

Copy number must not reach the ratio
------------------------------------
The 16S arm divides each taxon's reads by its estimated rrn copy number, which
is right for **composition** - it turns reads into something closer to genome
equivalents, so two different organisms can be compared. RPM-r is not
composition. It compares **one organism to itself** across two samples, so copy
number is a constant on both sides and cancels. Feeding it the corrected figure
therefore changes nothing... until the two sides resolve to different ranks,
and then it changes everything.

A negative control is shallow by nature, so it resolves to genus where the
sample resolves to species, and the genus fallback in :func:`lookup_baseline`
is built to exploit exactly that. A species copy number of 3.00 against a genus
estimate of 4.08 would then shift RPM-r by 36% on identical reads - the sample
being better characterised than the blank, which is not evidence that it is
cleaner. Worse, it pushes the wrong way: the genus fallback exists because a genus total
"is at least as high as any species-level value, making the gate stricter
rather than looser", and dividing that higher total by a higher divisor undoes
precisely that.

So the gate is fed **raw** reads on both sides - here and in
``baseline.py``, which must agree or the ratio is meaningless. The
corrected figures remain in ``marker16s_calls.csv`` for the questions they do
answer.

Scope - read this before using the numbers for anything else
------------------------------------------------------------
RPM-r is valid for **one organism compared across samples**, because capture
efficiency is constant for a given organism and so a change is real. That covers
the gate below, and longitudinal trending.

It is **not** valid for comparing **different organisms within one sample**.
RPM-r rescales every organism by its own private constant, so

    RPM-r(A)/RPM-r(B) = RPM(A)/RPM(B) x RPM_NTC(B)/RPM_NTC(A)

and that second term is a contamination artefact with no abundance meaning. A
heavy kitome organism and a never-seen-in-reagents organism at identical true
abundance get RPM-r values differing by whatever their background levels differ
by. Within-sample composition belongs to the rMLST arm's ``pi``, used ordinally.
"""

from __future__ import annotations

import csv
import logging
import math
from dataclasses import dataclass
from pathlib import Path

# --- Defaults, from the clinical metagenomics literature ------------------

# Report an organism that HAS negative-control background only if it is this
# many times above it.
DEFAULT_RPM_R_MIN = 10.0

# For an organism with NO background at all, the ratio is undefined, so an
# absolute floor applies instead.
DEFAULT_RPM_MIN = 0.05

# Denominator floor: stops a very CLEAN control letting everything through.
MIN_NTC_RPM = 127.0

# Minimum depth for a control to be usable at all: stops a very SHALLOW control
# blocking everything. Below this, one stray read dominates the RPM and the
# control measures nothing. At 10,000 reads a single read is 100 RPM, which is
# already below MIN_NTC_RPM and therefore cannot distort the ratio.
MIN_CONTROL_READS = 10_000

# --- Two mechanisms, two comparators ----------------------------------------
#
# What a negative control holds gets there in two ways, and they scale
# differently with sequencing depth.
#
# *Reagent background* (the kitome) is a roughly fixed number of molecules per
# library. In a blank those molecules ARE the library, so they are a large
# fraction of its reads; in a deep, low-biomass sample they are still a large
# fraction. A fraction is what RPM measures, so RPM-r is the right comparator.
#
# *Cross-contamination* - carry-over or index hopping from a strong positive on
# the same run - arrives as an absolute number of reads per library, not a
# fraction. A blank of tens of thousands of reads beside samples of millions
# turns a large absolute excess into a small RPM-r: scaling the blank's count
# up by the depth ratio inflates it ~30-100x and marks genuine organisms as
# background. For these the *absolute* count is the comparator.
#
# Which mechanism applies is a run-level fact: an organism in the blank that an
# rMLST-confirmed sample on the same run carries in bulk is cross-contamination;
# one no sample credibly carries is reagent. ``build_baselines`` records the
# mechanism per organism from that evidence (see ``baseline.scan_cohort``), and
# ``evaluate`` picks the rule from it. A file without the column is read as
# reagent throughout.
MECHANISM_REAGENT = "reagent"
MECHANISM_CROSS = "cross_contamination"
MECHANISMS = (MECHANISM_REAGENT, MECHANISM_CROSS)

# Cross-contamination rule: sample reads must be this many times the blank's
# ABSOLUTE reads - the same 10x convention as RPM-r, on the right quantity...
DEFAULT_ABS_R_MIN = 10.0
# ...and at least this many reads outright, so a blank that saw 0.5 of a read
# cannot be cleared by 5.
DEFAULT_MIN_SAMPLE_READS = 10.0
# A blank holding more than this fraction of its source sample's reads is
# unusually high for cross-contamination alone; the row is noted so a reviewer
# can weigh a reagent component too.
MAX_NTC_SOURCE_FRAC = 0.05

BASELINE_COLUMNS = [
    "run_id", "arm", "organism", "rpm_ntc", "reads_ntc", "total_reads_ntc",
    "n_controls", "control_ids", "source", "mechanism", "source_sample",
    "source_reads", "ntc_over_source_frac", "note",
]


@dataclass(frozen=True)
class Baseline:
    """One organism's level in the run's negative control(s), per arm."""

    rpm_ntc: float
    reads_ntc: float | None = None        # absolute reads in the control that set rpm_ntc
    total_reads_ntc: int | None = None    # that control's depth
    mechanism: str = MECHANISM_REAGENT
    source_sample: str | None = None      # strongest credible source on the run, if cross
    source_reads: float | None = None
    n_controls: int = 1
    pooled: bool = False                  # a run_id `*` row: another run's blank, not this one's


def as_baseline(value: "Baseline | float | int | None") -> Baseline | None:
    """Accept a bare RPM where a :class:`Baseline` is wanted (older callers, tests)."""
    if value is None:
        return None
    if isinstance(value, Baseline):
        return value
    return Baseline(rpm_ntc=float(value))


class CohortEvidence:
    """Per arm: organism -> (sample, reads) of its strongest credible source on the run.

    Keyed by lowercased species binomial and by lowercased genus. Lookup falls
    from species to genus - the same direction :func:`lookup_baseline` falls -
    so a blank that resolved only ``Enterobacter`` is attributed from the
    sample that carried *Enterobacter hormaechei*.
    """

    def __init__(self, species: dict[str, tuple[str, float]] | None = None,
                 genera: dict[str, tuple[str, float]] | None = None) -> None:
        self.species = {k.strip().lower(): v for k, v in (species or {}).items()}
        self.genera = {k.strip().lower(): v for k, v in (genera or {}).items()}

    def source_for(self, organism: str) -> tuple[str, float] | None:
        name = organism.strip().lower()
        if not name:
            return None
        if name in self.species:
            return self.species[name]
        if name in self.genera:
            return self.genera[name]
        genus = name.split()[0]
        return self.genera.get(genus)

    def __len__(self) -> int:
        return len(self.species) + len(self.genera)


def _optional_float(value) -> float | None:
    try:
        text = str(value).strip()
        return float(text) if text else None
    except (TypeError, ValueError):
        return None


def _optional_int(value) -> int | None:
    number = _optional_float(value)
    return int(number) if number is not None else None


def baseline_from_row(row: dict) -> Baseline | None:
    """One CSV row (new or legacy schema) -> :class:`Baseline`, or None if unusable."""
    rpm_ntc = _optional_float(row.get("rpm_ntc"))
    if rpm_ntc is None:
        return None
    mechanism = (str(row.get("mechanism") or "").strip().lower()) or MECHANISM_REAGENT
    if mechanism not in MECHANISMS:
        logging.getLogger("bridge.reconcile.rpm").warning(
            "baseline row for %s has unknown mechanism %r; treated as reagent",
            row.get("organism"), mechanism)
        mechanism = MECHANISM_REAGENT
    return Baseline(
        rpm_ntc=rpm_ntc,
        reads_ntc=_optional_float(row.get("reads_ntc")),
        total_reads_ntc=_optional_int(row.get("total_reads_ntc")),
        mechanism=mechanism,
        source_sample=(str(row.get("source_sample") or "").strip() or None),
        source_reads=_optional_float(row.get("source_reads")),
        n_controls=_optional_int(row.get("n_controls")) or 1,
        pooled=str(row.get("run_id", "")).strip() == "*",
    )


def baselines_from_rows(rows: list[dict]) -> dict[tuple[str, str, str], Baseline]:
    """Rows as :func:`build_baselines` returns them -> the lookup :func:`evaluate` uses."""
    baselines: dict[tuple[str, str, str], Baseline] = {}
    for row in rows:
        baseline = baseline_from_row(row)
        if baseline is None:
            continue
        key = (str(row.get("run_id", "*")).strip(),
               str(row.get("arm", "")).strip().lower(),
               str(row.get("organism", "")).strip().lower())
        baselines[key] = baseline
    return baselines


@dataclass(frozen=True)
class RpmResult:
    """One organism's standing against its reagent background."""

    organism: str
    arm: str
    reads: float
    total_reads: int
    rpm: float
    rpm_ntc: float | None      # None where the organism is absent from controls
    rpm_r: float | None        # None where there is no background to divide by
    passed: bool
    reason: str
    # Both comparators are always computed where they can be; `rule`
    # says which one decided. "rpm_r" | "abs_r" | "floor" | "not_applied".
    reads_ntc: float | None = None
    total_reads_ntc: int | None = None
    abs_r: float | None = None            # sample reads / control reads (absolute)
    mechanism: str | None = None
    rule: str = "rpm_r"
    source_sample: str | None = None
    pooled: bool = False                  # the baseline came from the cross-run pool

    @property
    def log10_rpm_r(self) -> float | None:
        if self.rpm_r is None or self.rpm_r <= 0:
            return None
        return math.log10(self.rpm_r)

    @property
    def governing_ratio(self) -> float | None:
        """The ratio the rule actually tested."""
        return self.abs_r if self.rule == "abs_r" else self.rpm_r

    @property
    def band(self) -> str:
        """Ordinal band, which is how this should be read.

        Run-to-run, kallisto's ``est_counts`` move by up to ~27% (about 0.1
        log10) while the reporting threshold spans a full log10, so the decision
        sits an order of magnitude above the noise. Banding keeps it that way;
        quoting a continuous ratio would not. The band is on whichever ratio
        the rule tested - RPM-r for reagent background, abs-r for
        cross-contamination.
        """
        if self.rule == "not_applied":
            return "ungated"
        if not self.passed:
            return "background"
        ratio = self.governing_ratio
        if ratio is None:
            return "no_control_background"
        if ratio >= 100:
            return "high"
        if ratio >= 10:
            return "moderate"
        return "low"


def compute_rpm(reads: float, total_reads: int) -> float:
    """Reads per million of the sample's total raw reads."""
    if total_reads <= 0:
        raise ValueError(f"total_reads must be positive, got {total_reads}")
    return (reads / total_reads) * 1_000_000


def load_baselines(path: Path) -> dict[tuple[str, str, str], Baseline]:
    """Read per-run, per-arm, per-organism negative-control levels.

    Keyed ``(run_id, arm, organism_lowercased)``. Rows with ``run_id`` of
    ``"*"`` are a pooled fallback used where a given run's control never saw
    that organism.

    Per-run granularity is deliberate: reagent lots change between runs and the
    kitome changes with them, so a control is only properly comparable to the
    samples it was extracted alongside.

    A file with only ``rpm_ntc`` loads with ``reads_ntc=None`` and every
    mechanism ``reagent``, so the gate falls back to RPM-r alone.
    """
    with open(path, newline="", encoding="utf-8") as handle:
        return baselines_from_rows(list(csv.DictReader(handle)))


def lookup_baseline(
    baselines: dict[tuple[str, str, str], "Baseline | float"],
    run_id: str,
    arm: str,
    organism: str,
    genus: str | None = None,
) -> "Baseline | float | None":
    """This run's control level for an organism.

    Tried in order: this run's exact name, the pooled exact name, then the same
    two at **genus** level.

    The genus fallback matters more than it looks. A negative control is shallow
    by nature, so its organisms are often only resolvable to genus - the
    control reports ``Corynebacterium`` (genus_only) while the sample calls
    ``Corynebacterium resistens``. Matching on the full name alone, that call
    reads as having no reagent background at all, and the gate silently stops
    protecting against exactly the contamination it was built to catch.

    Falling back to the genus errs the safe way: a genus baseline aggregates
    every species in it, so it is at least as high as any species-level value,
    making the gate stricter rather than looser.
    """
    arm = arm.strip().lower()
    names = [organism.strip().lower()]
    if genus:
        genus = genus.strip().lower()
        if genus not in names:
            names.append(genus)
    for name in names:
        for key in ((run_id.strip(), arm, name), ("*", arm, name)):
            if key in baselines:
                return baselines[key]
    return None


def arms_with_background(
    baselines: dict[tuple[str, str, str], float],
    run_id: str,
    organism: str,
    genus: str | None = None,
) -> list[str]:
    """Which arms' controls saw this organism at all - regardless of arm asked.

    Presence in a negative control is a fact about the tube, not about the
    marker used to detect it. The *ratio* is arm-specific and must stay that way
    - rMLST counts fragments over ~53 ribosomal protein loci, 16S counts
    copy-number-corrected reads over one gene, and dividing one by the other
    means nothing. But when an arm reports ``no_control_background`` for an
    organism the *other* arm plainly saw in the blank, that phrase is false as
    English even though it is true as arithmetic.

    For example, a species with no rMLST baseline passes that arm on the 0.05
    RPM floor, which nothing can fail, while the 16S arm scores the same
    organism barely clear of a genus baseline in the blank. The rMLST row
    alone would say there was no background.

    Matching is deliberately wider here than in :func:`lookup_baseline`, in the
    one direction that function cannot safely go. ``lookup_baseline`` will read
    a *genus* baseline for a *species* call, because a genus total is at least
    as large as any species within it and so only makes the gate stricter. It
    will not do the reverse. But for a flag - which changes no arithmetic - the
    reverse is exactly the case that matters: a call of *Staphylococcus* should
    not read "no control background" when the blank's rMLST arm held
    *Staphylococcus aureus*.

    So an entry counts if it names the organism, names the genus, or sits
    **within** the genus.

    Returns the arms sorted, so the flag text is stable between runs.
    """
    names = {organism.strip().lower()}
    within = f"{genus.strip().lower()} " if genus else None
    if genus:
        names.add(genus.strip().lower())
    wanted = {run_id.strip(), "*"}
    return sorted({
        arm for (run, arm, name) in baselines
        if run in wanted
        and (name in names or (within and name.startswith(within)))
    })


def rpm_detail(
    result: "RpmResult",
    other_arm_background: list[str] | None = None,
) -> dict[str, object]:
    """One organism's gate verdict, shaped for the JSON record and the CSV.

    ``other_arm_background`` names arms whose control saw this organism when
    this arm's did not - see :func:`arms_with_background`. Empty for the normal
    case, which keeps the record unchanged where there is nothing to say.
    """
    detail: dict[str, object] = {
        "rpm": result.rpm,
        "rpm_r": result.rpm_r,
        "band": result.band,
        "passed": result.passed,
        "reason": result.reason,
        # The control's own figures and the second comparator, so a
        # reader can see "15,488 reads against 78 in the blank" next to an
        # RPM-r of 6.2 and know which one the rule used.
        "rpm_ntc": result.rpm_ntc,
        "reads_ntc": result.reads_ntc,
        "total_reads_ntc": result.total_reads_ntc,
        "abs_r": result.abs_r,
        "mechanism": result.mechanism,
        "rule": result.rule,
    }
    if result.source_sample:
        detail["source_sample"] = result.source_sample
    if result.pooled:
        detail["pooled_baseline"] = True
    if other_arm_background:
        detail["other_arm_background"] = list(other_arm_background)
    return detail


def ungated(organism: str, arm: str, reads: float, total_reads: int) -> RpmResult:
    """The verdict when there is no baseline at all: RPM computed, nothing tested.

    RPM itself never needed the control - only the ratio does - so a run with
    no blank still gets every ``rpm_*`` column filled and only the ratios
    blank. ``passed`` is True because
    nothing was failed; the band says ``ungated`` so it cannot be read as clean.
    """
    return RpmResult(
        organism=organism, arm=arm, reads=reads, total_reads=total_reads,
        rpm=compute_rpm(reads, total_reads), rpm_ntc=None, rpm_r=None,
        passed=True, reason="no baseline supplied; RPM computed, gate not applied",
        rule="not_applied",
    )


def evaluate(
    organism: str,
    arm: str,
    reads: float,
    total_reads: int,
    rpm_ntc: "Baseline | float | None",
    *,
    rpm_r_min: float = DEFAULT_RPM_R_MIN,
    rpm_min: float = DEFAULT_RPM_MIN,
    min_ntc_rpm: float = MIN_NTC_RPM,
    abs_r_min: float = DEFAULT_ABS_R_MIN,
    min_sample_reads: float = DEFAULT_MIN_SAMPLE_READS,
) -> RpmResult:
    """Decide whether one organism clears its background in the control.

    ``rpm_ntc`` is the organism's :class:`Baseline` (preferred) or, for older
    callers, its bare control RPM - which is read as a reagent baseline with no
    read count. Three outcomes:

    * ``floor`` - absent from the control: RPM against the absolute floor;
    * ``abs_r`` - cross-contamination with a known control read count:
      ``reads / reads_ntc`` against ``abs_r_min``, and ``reads`` against
      ``min_sample_reads``;
    * ``rpm_r`` - everything else: ``RPM / max(rpm_ntc, min_ntc_rpm)``
      against ``rpm_r_min``, the original rule.

    Both ratios are computed whenever they can be, whichever one decides.
    """
    rpm = compute_rpm(reads, total_reads)
    baseline = as_baseline(rpm_ntc)

    if baseline is None or baseline.rpm_ntc <= 0:
        # No background for this organism: the ratio is undefined, so fall back
        # to the absolute floor. Anything seen at all is notable here.
        passed = rpm >= rpm_min
        return RpmResult(
            organism=organism, arm=arm, reads=reads, total_reads=total_reads,
            rpm=rpm, rpm_ntc=None, rpm_r=None, passed=passed,
            reason=(f"absent from controls; RPM {rpm:.3g} "
                    f"{'>=' if passed else '<'} floor {rpm_min}"),
            rule="floor",
        )

    denominator = max(float(baseline.rpm_ntc), min_ntc_rpm)
    rpm_r = rpm / denominator
    floored = denominator > baseline.rpm_ntc
    abs_r = (reads / max(float(baseline.reads_ntc), 1.0)
             if baseline.reads_ntc is not None else None)
    common = dict(
        organism=organism, arm=arm, reads=reads, total_reads=total_reads,
        rpm=rpm, rpm_ntc=baseline.rpm_ntc, rpm_r=rpm_r,
        reads_ntc=baseline.reads_ntc, total_reads_ntc=baseline.total_reads_ntc,
        abs_r=abs_r, mechanism=baseline.mechanism,
        source_sample=baseline.source_sample, pooled=baseline.pooled,
    )
    pooled_note = " [pooled baseline]" if baseline.pooled else ""

    if baseline.mechanism == MECHANISM_CROSS and abs_r is not None:
        enough = reads >= min_sample_reads
        clear = abs_r >= abs_r_min
        reason = (f"abs-r {abs_r:.3g} {'>=' if clear else '<'} {abs_r_min:g} "
                  f"({reads:.0f} vs {baseline.reads_ntc:.0f} NTC reads) "
                  f"[{MECHANISM_CROSS}]")
        if not enough:
            reason += f"; sample reads {reads:.0f} < floor {min_sample_reads:g}"
        return RpmResult(**common, passed=clear and enough, reason=reason + pooled_note,
                         rule="abs_r")

    passed = rpm_r >= rpm_r_min
    reason = (f"RPM-r {rpm_r:.3g} {'>=' if passed else '<'} {rpm_r_min}"
              + (f" (NTC RPM {baseline.rpm_ntc:.3g} floored to {denominator:.0f})"
                 if floored else "")
              + f" [{baseline.mechanism}]")
    if baseline.mechanism == MECHANISM_CROSS:
        reason += " (reads_ntc unknown; RPM-r used)"
    return RpmResult(**common, passed=passed, reason=reason + pooled_note, rule="rpm_r")


def build_baselines(
    control_profiles: dict[str, dict[str, tuple[float, int]]],
    run_id: str,
    arm: str,
    *,
    min_control_reads: int = MIN_CONTROL_READS,
    cohort: CohortEvidence | None = None,
    max_ntc_source_frac: float = MAX_NTC_SOURCE_FRAC,
) -> list[dict[str, object]]:
    """Summarise negative controls into baseline rows ready to be written.

    ``control_profiles`` maps control sample id -> organism -> (reads, total).
    Where a run has more than one control, the **maximum** RPM across them is
    taken, not the mean: the question a baseline answers is "how high does this
    organism get in a blank", and one heavily contaminated control is exactly
    the case a mean would hide. The read count and depth recorded beside it
    are that same control's.

    Controls shallower than ``min_control_reads`` are excluded entirely. Their
    RPM values are dominated by single stray reads and would block genuine
    organisms - see the module docstring.

    ``cohort`` - what the run's samples credibly carry - decides each row's
    ``mechanism``: an organism with a credible source is ``cross_contamination``
    and will be gated on absolute reads; the rest is ``reagent``. Without it
    every row is reagent.
    """
    logger = logging.getLogger("bridge.reconcile.rpm")
    per_organism: dict[str, list[tuple[float, float, int, str]]] = {}
    for sample, profile in control_profiles.items():
        totals = {total for _, total in profile.values()}
        depth = max(totals) if totals else 0
        if depth < min_control_reads:
            logger.warning(
                "control %s has only %d reads (< %d); excluded from the baseline "
                "- too shallow for its RPM to mean anything",
                sample, depth, min_control_reads,
            )
            continue
        for organism, (reads, total) in profile.items():
            if total <= 0:
                logger.warning("control %s has no reads; skipped", sample)
                continue
            per_organism.setdefault(organism, []).append(
                (compute_rpm(reads, total), reads, total, sample))

    rows = []
    for organism, values in sorted(per_organism.items()):
        rpm, reads, total, control = max(values, key=lambda item: item[0])
        row: dict[str, object] = {
            "run_id": run_id,
            "arm": arm,
            "organism": organism,
            "rpm_ntc": f"{rpm:.6g}",
            "reads_ntc": f"{reads:.6g}",
            "total_reads_ntc": total,
            "n_controls": len(values),
            "control_ids": "; ".join(sorted({item[3] for item in values})),
            "source": "max_across_controls",
            "mechanism": MECHANISM_REAGENT,
            "source_sample": "",
            "source_reads": "",
            "ntc_over_source_frac": "",
            "note": "",
        }
        source = cohort.source_for(organism) if cohort is not None else None
        if source is not None:
            source_sample, source_reads = source
            row["mechanism"] = MECHANISM_CROSS
            row["source_sample"] = source_sample
            row["source_reads"] = f"{source_reads:.6g}"
            if source_reads > 0:
                frac = reads / source_reads
                row["ntc_over_source_frac"] = f"{frac:.4g}"
                if frac > max_ntc_source_frac:
                    row["note"] = (
                        f"blank holds {frac:.0%} of its source's reads, high for "
                        f"cross-contamination alone (> {max_ntc_source_frac:.0%}); "
                        f"a reagent component is possible")
        rows.append(row)
    logger.info("baseline_status=built run=%s arm=%s organisms=%d cross_contamination=%d",
                run_id, arm, len(rows),
                sum(1 for row in rows if row["mechanism"] == MECHANISM_CROSS))
    return rows


def write_baselines(rows: list[dict[str, object]], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=BASELINE_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    return path
