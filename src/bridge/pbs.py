"""PBS submission: one job per sample, then one finalise job that waits for them.

* resource tiers by combined FASTQ size;
* a sample whose job is still queued or running is left alone;
* a sample with every requested arm complete is skipped (so ``--arms amr``
  on an existing results tree submits only the samples still missing it,
  and clears only ``amr/`` before resubmitting);
* a refusal to touch a directory that looks like a Castanet output tree;
* ``--max-submit`` to stay under a per-user queue limit.

Every job script is rendered in full and kept under ``<out>/oe/scripts/``, so
what ran is on disk. The finalise job is submitted with
``-W depend=afterany:<all sample jobs>``, so one ``bridge run`` ends with the
run-level table. ``afterany`` rather than ``afterok``: one failed sample must
not hold up the whole run, and finalise already skips incomplete samples.
"""

from __future__ import annotations

import getpass
import os
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass
from importlib import resources
from pathlib import Path

from bridge.tools import companion_executable

from . import __version__
from .reads import ReadPair
from .settings import (
    ARMS, REPO_ROOT, Settings, Tier, arms_arg, tier_for,
)

# Directories that exist only in a Castanet results folder (bacterial_ID/ is
# not one: `bridge finalise` writes it into our own tree).
CASTANET_MARKERS = ("consensus_sequences", "graphcalls", "hashes", "qc")
ARM_RESULT = {"rmlst": Path("rmlst") / "species_calls.csv",
              "16s": Path("16s") / "marker16s_calls.csv",
              "amr": Path("amr") / "amr_hits.csv"}

FINALISE_TEMPLATE = """#!/bin/bash
#PBS -N {job_name}
#PBS -P {project}
#PBS -q {queue}
#PBS -l walltime={walltime}
#PBS -l ncpus={ncpus}
#PBS -l mem={mem_gb}GB
{storage_line}
#PBS -j oe
{depend_line}
#
# Rendered by bridge {version}: baseline + reconcile + run table for {run_id}.
module purge >/dev/null 2>&1 || true
{python_module_line}
{pythonpath_line}
set -euo pipefail

{bridge_command} finalise \\
    --config "{config}" \\
    --out "{out_root}" \\
    --run-id "{run_id}"{control_args}
"""


@dataclass(frozen=True)
class Submission:
    sample: str
    action: str          # submitted | dry_run | skipped_complete | skipped_inflight | deferred | failed | refused
    detail: str = ""
    job_id: str | None = None


def bridge_command() -> str:
    """The environment's own ``bridge``, by absolute path.

    The job then does not care what ``python`` the tool modules put on PATH.
    """
    bindir = Path(sys.executable).resolve().parent
    for name in ("bridge", "bridge.exe"):
        candidate = bindir / name
        if candidate.is_file():
            return shlex.quote(str(candidate))
    return f"{shlex.quote(sys.executable)} -m bridge.cli"


def arms_complete(sample_dir: Path, arms: tuple[str, ...]) -> bool:
    return all((Path(sample_dir) / ARM_RESULT[arm]).is_file() for arm in arms)


def castanet_marker(sample_dir: Path) -> str | None:
    for marker in CASTANET_MARKERS:
        if (Path(sample_dir) / marker).is_dir():
            return marker
    return None


def _user() -> str:
    return os.environ.get("USER") or os.environ.get("USERNAME") or getpass.getuser()


def active_jobs() -> set[str]:
    """This user's queued and running job ids; empty when PBS is not here."""
    if shutil.which("qselect") is None:
        return set()
    try:
        result = subprocess.run(["qselect", "-u", _user()], capture_output=True,
                                text=True, check=False)
    except OSError:
        return set()
    return {line.strip() for line in result.stdout.splitlines() if line.strip()}


def qsub(script: Path, cwd: Path) -> str:
    result = subprocess.run(["qsub", str(script)], capture_output=True, text=True,
                            check=True, cwd=str(cwd))
    return result.stdout.strip()


def _storage_line(settings: Settings) -> str:
    storage = settings.hpc.storage
    return f"#PBS -l storage={storage}" if storage else "# (no storage directive configured)"


def _module_line(settings: Settings) -> str:
    module = settings.hpc.python_module
    return f"module load {module}" if module else "# (no python module configured)"


def _tools_line(settings: Settings) -> str:
    env = settings.tools.tools_env
    if env is None:
        return "# (tools.tools_env is null: binaries expected on PATH)"
    return f"source {shlex.quote(str(env))}"


def _pythonpath_line() -> str:
    """Pin PYTHONPATH to this repository's code (or clear it).

    The tool modules can put their own Python packages on PYTHONPATH, and a
    PYTHONPATH entry outranks an editable install - so without this a job
    could import an older copy of the code while reporting the venv's version.
    """
    src = REPO_ROOT / "src"
    return f'export PYTHONPATH={shlex.quote(str(src)) if src.is_dir() else ""!s}'


def prepare_blast_cache(settings: Settings, out_root: Path) -> str:
    """Build the rMLST arm's shared BLAST database once, before any job needs it.

    The rMLST arm caches the full-reference database under
    ``<out>/.bridge_blastdb`` and builds it unlocked on first use. With
    hundreds of concurrent jobs several can build the same files at once and
    corrupt them, so it is built here up front. Returns a one-line status for
    the console.
    """
    reference = settings.references.rmlst_alleles
    if reference is None or not reference.is_file():
        return f"skipped: rMLST allele reference not found ({reference})"
    try:
        makeblastdb = companion_executable(settings.tools.makeblastdb, settings.tools.blastn)
    except (FileNotFoundError, OSError):
        return ("WARNING: makeblastdb not on PATH (source tools.env first). Each job "
                "will try to build the BLAST cache itself, and concurrent builds can "
                "corrupt it.")
    from bridge.rmlst.consensus import _ensure_blast_database
    try:
        prefix = _ensure_blast_database(reference, Path(out_root).resolve(), makeblastdb)
    except (subprocess.CalledProcessError, OSError) as error:
        return f"WARNING: makeblastdb failed ({error}); each job will try to build it"
    return f"ready: {prefix}"


def render_sample_script(settings: Settings, pair: ReadPair, out_root: Path,
                         run_id: str, tier: Tier, *, arms: tuple[str, ...],
                         size_bytes: int, keep_temp: bool = False) -> str:
    template = resources.files("bridge").joinpath("templates/sample.pbs").read_text(encoding="utf-8")
    return template.format(
        job_name=f"{pair.sample}_{run_id}"[:230],
        project=settings.hpc.project, queue=settings.hpc.queue,
        walltime=tier.walltime, ncpus=tier.ncpus, mem_gb=tier.mem_gb,
        jobfs_gb=tier.jobfs_gb, storage_line=_storage_line(settings),
        version=__version__, sample=pair.sample, tier=tier.name,
        size_mb=f"{size_bytes / 1e6:.0f}",
        python_module_line=_module_line(settings),
        tools_env_line=_tools_line(settings),
        pythonpath_line=_pythonpath_line(),
        bridge_command=bridge_command(),
        config=settings.path, r1=pair.r1, r2=pair.r2, out_root=out_root,
        arms=arms_arg(tuple(arms)),
        keep_temp_arg=" \\\n    --keep-temp" if keep_temp else "",
    )


def render_finalise_script(settings: Settings, out_root: Path, run_id: str, *,
                           depends_on: list[str], controls: list[str]) -> str:
    control_args = "".join(f' \\\n    --control "{name}"' for name in controls)
    depend = (f"#PBS -W depend=afterany:{':'.join(depends_on)}" if depends_on
              else "# (no sample jobs to wait for)")
    return FINALISE_TEMPLATE.format(
        job_name=f"{run_id}_finalise"[:230],
        project=settings.hpc.project, queue=settings.hpc.queue,
        walltime=settings.hpc.finalise_walltime, ncpus=settings.hpc.finalise_ncpus,
        mem_gb=settings.hpc.finalise_mem_gb, storage_line=_storage_line(settings),
        depend_line=depend, version=__version__, run_id=run_id,
        python_module_line=_module_line(settings),
        pythonpath_line=_pythonpath_line(),
        bridge_command=bridge_command(), config=settings.path,
        out_root=out_root, control_args=control_args,
    )


def submit_samples(settings: Settings, pairs: list[ReadPair], out_root: Path,
                   run_id: str, *, arms: tuple[str, ...] = ARMS,
                   max_submit: int = 0, dry_run: bool = False,
                   force: bool = False, keep_temp: bool = False) -> list[Submission]:
    """Submit one job per sample that needs one.

    ``force`` resubmits samples whose requested arm output already exists
    (their arm folders are removed first); a sample whose job is still queued
    or running is left alone either way.
    """
    out_root = Path(out_root).resolve()
    oe_dir = out_root / "oe"
    scripts_dir = oe_dir / "scripts"
    jobids_dir = oe_dir / ".jobids"
    for directory in (scripts_dir, jobids_dir):
        directory.mkdir(parents=True, exist_ok=True)

    running = active_jobs() if not dry_run else set()
    submissions: list[Submission] = []
    submitted = 0

    for pair in pairs:
        sample_dir = out_root / pair.sample
        jid_file = jobids_dir / pair.sample

        if jid_file.is_file():
            previous = jid_file.read_text().strip()
            if previous and previous in running:
                submissions.append(Submission(pair.sample, "skipped_inflight",
                                              f"job {previous} still queued or running", previous))
                continue

        if not force and arms_complete(sample_dir, arms):
            submissions.append(Submission(pair.sample, "skipped_complete",
                                          "requested arm output(s) already present"))
            continue

        if max_submit and submitted >= max_submit:
            submissions.append(Submission(pair.sample, "deferred",
                                          f"--max-submit {max_submit} reached; re-run to submit the rest"))
            continue

        marker = castanet_marker(sample_dir)
        if marker:
            submissions.append(Submission(pair.sample, "refused",
                                          f"{sample_dir} contains {marker}/ and looks like a Castanet tree"))
            continue

        size = pair.size_bytes()
        tier = tier_for(settings.hpc, size)
        script_text = render_sample_script(settings, pair, out_root, run_id, tier,
                                           arms=arms, size_bytes=size, keep_temp=keep_temp)
        script = scripts_dir / f"{pair.sample}.pbs"
        script.write_text(script_text, encoding="utf-8")

        if dry_run:
            submissions.append(Submission(pair.sample, "dry_run",
                                          f"tier {tier.name}: {script}"))
            continue

        # A partial previous attempt is removed so the arms start clean; the
        # Castanet guard above is what makes this safe.
        for arm in arms:
            arm_dir = sample_dir / ARM_RESULT[arm].parent
            if arm_dir.exists():
                shutil.rmtree(arm_dir)
        try:
            job_id = qsub(script, cwd=oe_dir)
        except (subprocess.CalledProcessError, OSError) as error:
            detail = getattr(error, "stderr", "") or str(error)
            submissions.append(Submission(pair.sample, "failed", f"qsub failed: {detail.strip()}"))
            jid_file.unlink(missing_ok=True)
            continue
        jid_file.write_text(job_id + "\n")
        submitted += 1
        submissions.append(Submission(
            pair.sample, "submitted",
            f"tier {tier.name} ({size / 1e6:.0f} MB, {tier.ncpus} cpu, {tier.walltime})",
            job_id))
    return submissions


def submit_finalise(settings: Settings, out_root: Path, run_id: str, *,
                    depends_on: list[str], controls: list[str],
                    dry_run: bool = False) -> tuple[Path, str | None]:
    out_root = Path(out_root).resolve()
    oe_dir = out_root / "oe"
    scripts_dir = oe_dir / "scripts"
    scripts_dir.mkdir(parents=True, exist_ok=True)
    script = scripts_dir / f"finalise_{run_id}.pbs"
    script.write_text(render_finalise_script(settings, out_root, run_id,
                                             depends_on=depends_on, controls=controls),
                      encoding="utf-8")
    if dry_run:
        return script, None
    return script, qsub(script, cwd=oe_dir)
