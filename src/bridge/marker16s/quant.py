"""One kallisto pass against the 16S reference.

Deliberately much leaner than the rMLST arm's quantification, and the omissions
are the point:

No bootstraps
    The rMLST arm runs 100 bootstrap replicates because bootstrap support is
    what selects its *candidate alleles* - the set of sequences worth
    assembling. This arm selects nothing; it reports the abundances kallisto
    produces. With no candidate rule to drive, bootstraps would cost roughly a
    hundredfold for a number nothing reads. (They are also the source of the
    rMLST arm's run-to-run instability, so not running them is a second win.)

No --pseudobam
    That exists to attribute individual fragments to runtime clusters. There is
    no clustering here.

No clustering or consensus
    A single marker gene cannot separate close relatives however carefully it is
    assembled, so reconstructing haplotypes would add cost and an illusion of
    resolution. Where 16S genuinely cannot resolve, the right answer is to say
    so - which is what the GENUS_ONLY rank in :mod:`.calls` is for.

The result is a stage measured in seconds rather than the rMLST arm's minutes.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import pandas as pd

from bridge.tools import executable, run, stage

# Columns kallisto writes to abundance.tsv with --plaintext.
ABUNDANCE_COLUMNS = ["target_id", "length", "eff_length", "est_counts", "tpm"]

# The read counts kallisto records in run_info.json, copied into the manifest.
# n_processed is the RPM denominator the reconciler uses.
RUN_INFO_FIELDS = ("n_processed", "n_pseudoaligned", "n_unique",
                   "p_pseudoaligned", "p_unique")


def read_run_info(output: Path) -> dict[str, float | int]:
    """The read counts from ``<output>/quant/run_info.json``; empty if absent."""
    path = Path(output) / "quant" / "run_info.json"
    if not path.is_file():
        return {}
    try:
        info = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return {key: info[key] for key in RUN_INFO_FIELDS if key in info}


def quant_16s(
    reads: list[Path],
    index: Path,
    output: Path,
    *,
    threads: int = 8,
    kallisto: str = "kallisto",
    reuse: bool = True,
) -> pd.DataFrame:
    """Quantify one sample against the 16S index; return the abundance table.

    ``reuse`` skips the quant when a complete ``abundance.tsv`` is already
    present, which makes re-running the downstream steps cheap during
    development and after a partial failure.
    """
    logger = logging.getLogger("marker16s")
    quant_dir = output / "quant"
    abundance = quant_dir / "abundance.tsv"

    if reuse and abundance.is_file() and abundance.stat().st_size > 0:
        logger.info("quant_status=reused path=%s", abundance)
    else:
        if len(reads) != 2:
            raise ValueError(
                "The 16S arm expects paired reads; got "
                f"{len(reads)} file(s). Single-end support is not implemented."
            )
        for path in reads:
            if not Path(path).is_file():
                raise FileNotFoundError(f"Read file not found: {path}")

        quant_dir.mkdir(parents=True, exist_ok=True)
        command = [
            executable(kallisto), "quant",
            "-i", str(index),
            "--plaintext",
            "-t", str(threads),
            "-o", str(quant_dir),
            *map(str, reads),
        ]
        with stage("marker16s_quant"):
            run(command, "kallisto")

    if not abundance.is_file():
        raise FileNotFoundError(
            f"kallisto produced no abundance.tsv in {quant_dir}. "
            "Check the run log for the kallisto error."
        )

    frame = pd.read_csv(abundance, sep="\t")
    missing = set(ABUNDANCE_COLUMNS) - set(frame.columns)
    if missing:
        raise ValueError(
            f"abundance.tsv is missing expected columns {sorted(missing)}; "
            f"found {list(frame.columns)}"
        )
    logger.info(
        "quant_status=parsed targets=%d assigned_reads=%.0f",
        len(frame), float(frame["est_counts"].sum()),
    )
    return frame
