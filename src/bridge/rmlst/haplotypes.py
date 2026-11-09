from __future__ import annotations

from dataclasses import dataclass
import logging
from pathlib import Path
import re
from statistics import median


FREQUENCY_RE = re.compile(
    r"(?:^|\s)(?:freq(?:uency)?)[=:](?P<frequency>(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)",
    re.IGNORECASE,
)
EXTENDED_FREQUENCY_RE = re.compile(
    r"_(?P<frequency>(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)$"
)
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


@dataclass(frozen=True)
class HaplotypeCandidate:
    identifier: str
    sequence: str
    proportion: float
    source_path: Path


@dataclass
class BaseEvidence:
    total: int = 0
    forward: int = 0
    reverse: int = 0
    interior: int = 0


@dataclass(frozen=True)
class AlignmentEvidence:
    depth: tuple[int, ...]
    bases: tuple[dict[str, BaseEvidence], ...]


@dataclass(frozen=True)
class HaplotypeEvaluation:
    retained: bool
    reasons: tuple[str, ...]
    supporting_fragments: float
    projected_median_depth: float
    callable_breadth: float
    variant_evidence_passed: bool


def maximum_minor_allele_depth(evidence: AlignmentEvidence) -> int:
    """Return the largest second-highest nucleotide count at any position."""
    maximum = 0
    for observations in evidence.bases:
        depths = sorted(
            (observation.total for observation in observations.values()),
            reverse=True,
        )
        if len(depths) >= 2:
            maximum = max(maximum, depths[1])
    return maximum


def read_cliquesnv_haplotypes(paths: list[Path]) -> list[HaplotypeCandidate]:
    """Read valid CliqueSNV FASTA records and skip malformed alternatives."""
    candidates: list[HaplotypeCandidate] = []
    seen: set[tuple[str, float]] = set()
    logger = logging.getLogger("cliquesnv")
    for path in sorted(paths):
        identifier: str | None = None
        description = ""
        sequence: list[str] = []

        def append_record() -> None:
            if identifier is None:
                return
            match = FREQUENCY_RE.search(description) or EXTENDED_FREQUENCY_RE.search(
                description
            )
            if match is None:
                logger.warning(
                    "haplotype_status=skipped identifier=%s path=%s "
                    "reason=missing_frequency",
                    identifier, path,
                )
                return
            proportion = float(match.group("frequency"))
            assembled = "".join(sequence).upper()
            if not 0.0 <= proportion <= 1.0:
                logger.warning(
                    "haplotype_status=skipped identifier=%s path=%s "
                    "reason=frequency_out_of_range value=%s",
                    identifier, path, proportion,
                )
                return
            invalid = sorted(set(assembled) - set("ACGTN"))
            if not assembled or invalid:
                logger.warning(
                    "haplotype_status=skipped identifier=%s path=%s "
                    "reason=%s invalid_characters=%s",
                    identifier, path,
                    "empty_sequence" if not assembled else "invalid_sequence",
                    "".join(invalid) or "none",
                )
                return
            key = (assembled, proportion)
            if key not in seen:
                candidates.append(
                    HaplotypeCandidate(identifier, assembled, proportion, path)
                )
                seen.add(key)

        with path.open(encoding="utf-8") as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    append_record()
                    description = line[1:]
                    identifier = description.split(None, 1)[0]
                    sequence = []
                elif identifier is None:
                    raise ValueError(f"Sequence before FASTA header in {path}")
                else:
                    sequence.append(line)
        append_record()
    return candidates


def alignment_evidence_from_sam(
    path: Path,
    reference_length: int,
    *,
    min_mapping_quality: int = 20,
    min_base_quality: int = 20,
    min_end_distance: int = 5,
) -> AlignmentEvidence:
    """Collect usable depth and strand/read-position support from a SAM file."""
    depth = [0] * reference_length
    bases: list[dict[str, BaseEvidence]] = [dict() for _ in range(reference_length)]
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            flag = int(fields[1])
            mapping_quality = int(fields[4])
            if flag & (0x4 | 0x100 | 0x800) or mapping_quality < min_mapping_quality:
                continue
            reference_position = int(fields[3]) - 1
            query_position = 0
            sequence = fields[9].upper()
            qualities = fields[10]
            reverse = bool(flag & 0x10)
            operations = CIGAR_RE.findall(fields[5])
            if not operations or "".join(length + code for length, code in operations) != fields[5]:
                continue
            for length_text, operation in operations:
                length = int(length_text)
                if operation in "M=X":
                    for offset in range(length):
                        ref_index = reference_position + offset
                        query_index = query_position + offset
                        if not 0 <= ref_index < reference_length or query_index >= len(sequence):
                            continue
                        quality = ord(qualities[query_index]) - 33 if qualities != "*" else 0
                        base = sequence[query_index]
                        if quality < min_base_quality or base not in "ACGT":
                            continue
                        depth[ref_index] += 1
                        observation = bases[ref_index].setdefault(base, BaseEvidence())
                        observation.total += 1
                        if reverse:
                            observation.reverse += 1
                        else:
                            observation.forward += 1
                        distance = min(query_index, len(sequence) - query_index - 1)
                        if distance >= min_end_distance:
                            observation.interior += 1
                    reference_position += length
                    query_position += length
                elif operation in "IS":
                    query_position += length
                elif operation in "DN":
                    reference_position += length
                elif operation in "HP":
                    continue
    return AlignmentEvidence(tuple(depth), tuple(bases))


def evaluate_haplotype(
    candidate: HaplotypeCandidate,
    reference: str,
    evidence: AlignmentEvidence,
    assigned_fragments: int,
    *,
    min_proportion: float,
    min_supporting_fragments: int,
    min_depth: float,
    min_breadth: float,
    min_interior_fraction: float = 0.5,
    validate_variants: bool = True,
) -> HaplotypeEvaluation:
    """Apply abundance and coverage filters, plus optional variant validation."""
    proportion = candidate.proportion
    supporting_fragments = assigned_fragments * proportion
    projected_depth = [value * proportion for value in evidence.depth]
    projected_median_depth = float(median(projected_depth)) if projected_depth else 0.0
    callable_breadth = (
        sum(
            value >= min_depth and index < len(candidate.sequence)
            and candidate.sequence[index] in "ACGT"
            for index, value in enumerate(projected_depth)
        ) / len(projected_depth)
        if projected_depth else 0.0
    )
    reasons: list[str] = []
    # Deliberately strict: exactly 1% is not "over 1%".
    if proportion <= min_proportion:
        reasons.append("proportion_not_over_threshold")
    if supporting_fragments < min_supporting_fragments:
        reasons.append("insufficient_supporting_fragments")
    if projected_median_depth < min_depth:
        reasons.append("insufficient_projected_depth")
    if callable_breadth < min_breadth:
        reasons.append("insufficient_projected_breadth")

    variant_evidence_passed = True
    if validate_variants and (
        len(candidate.sequence) != len(reference)
        or len(evidence.depth) != len(reference)
    ):
        variant_evidence_passed = False
        reasons.append("unphaseable_length_change")
    elif validate_variants:
        for index, (reference_base, haplotype_base) in enumerate(
            zip(reference.upper(), candidate.sequence)
        ):
            if haplotype_base not in "ACGT" or haplotype_base == reference_base:
                continue
            observation = evidence.bases[index].get(haplotype_base, BaseEvidence())
            interior_fraction = (
                observation.interior / observation.total if observation.total else 0.0
            )
            if (
                observation.forward < 1
                or observation.reverse < 1
                or interior_fraction < min_interior_fraction
            ):
                variant_evidence_passed = False
                reasons.append(f"weak_variant_evidence_at_{index + 1}")
    return HaplotypeEvaluation(
        retained=not reasons,
        reasons=tuple(reasons),
        supporting_fragments=supporting_fragments,
        projected_median_depth=projected_median_depth,
        callable_breadth=callable_breadth,
        variant_evidence_passed=variant_evidence_passed,
    )
