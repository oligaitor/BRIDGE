"""``bridge check`` - is this environment able to run the pipeline?

Python version, the venv, the Python packages, every binary the requested arms
call, every reference file, and (for PBS) the scheduler. Prints one line per
item and exits non-zero if anything required is missing.
"""

from __future__ import annotations

import os
import shutil
import sys
from importlib import metadata
from pathlib import Path

from bridge.tools import companion_executable, executable

from .settings import ARMS, REPO_ROOT, Settings

PYTHON_PACKAGES = ("pandas", "numpy", "matplotlib", "seaborn", "PyYAML")

CLIQUESNV_HINT = ("if it is a JAR rather than a program on PATH: set tools.cliquesnv "
                  "in config.yaml to its path, or export BRIDGE_CLIQUESNV in tools.env")


def _line(status: str, name: str, detail: str) -> str:
    return f"{status:<8} {name:<22} {detail}"


def _resolve(requested: str, companion: str | None = None) -> tuple[bool, str]:
    """Locate one binary. ``companion`` looks beside another tool, as BLAST's
    makeblastdb sits beside blastn."""
    try:
        if companion is not None:
            return True, companion_executable(requested, companion)
        return True, executable(requested)
    except (FileNotFoundError, OSError) as error:
        return False, str(error)


def run_check(settings: Settings, *, arms: tuple[str, ...] = ARMS,
              scheduler: str = "pbs", out: Path | None = None) -> int:
    lines: list[str] = []
    problems = 0

    version = ".".join(map(str, sys.version_info[:3]))
    if sys.version_info < (3, 10):
        lines.append(_line("MISSING", "python",
                           f"{version} - need 3.10+"))
        problems += 1
    else:
        lines.append(_line("ok", "python", f"{version}  {sys.executable}"))
    in_venv = sys.prefix != getattr(sys, "base_prefix", sys.prefix)
    lines.append(_line(
        "ok" if in_venv else "note", "environment",
        sys.prefix if in_venv else
        "not a virtual environment (setup_env.sh creates one; fine if this is a conda env)"))

    for package in PYTHON_PACKAGES:
        try:
            lines.append(_line("ok", package, metadata.version(package)))
        except metadata.PackageNotFoundError:
            lines.append(_line("MISSING", package, "pip install -e . (or setup_env.sh)"))
            problems += 1

    # Which copy of the pipeline code is actually imported. A PYTHONPATH left
    # behind by a tool module can point at an older copy; the job scripts pin
    # PYTHONPATH to this repository, and so should the shell.
    import bridge
    code_dir = Path(bridge.__file__).resolve().parent.parent
    lines.append(_line("ok" if code_dir == REPO_ROOT / "src" else "note", "bridge",
                       f"{bridge.__version__}  {code_dir}"
                       + ("" if code_dir == REPO_ROOT / "src" else
                          f"  (expected {REPO_ROOT / 'src'}; check PYTHONPATH)")))
    pythonpath = os.environ.get("PYTHONPATH", "")
    if pythonpath and str(REPO_ROOT / "src") not in pythonpath:
        lines.append(_line("note", "PYTHONPATH", f"{pythonpath}  (source tools.env to pin it to this repository)"))

    # Every tool gets its own line, so the report shows what does work too.
    tools = settings.tools
    wanted: list[tuple[str, str, str | None, bool, str]] = []
    if "amr" in arms:
        wanted += [
            ("kma", tools.kma, None, True, ""),
            # Only `bridge build-amr-db` needs it, and it falls back to
            # `kma index`, so its absence is a note, not a failure.
            ("kma_index", "kma_index", tools.kma, False,
             "build-amr-db will use 'kma index' instead"),
        ]
    if "rmlst" in arms:
        wanted += [
            ("kallisto", tools.kallisto, None, True, ""),
            ("mash", tools.mash, None, True, ""),
            ("minimap2", tools.minimap2, None, True, ""),
            ("samtools", tools.samtools, None, True, ""),
            ("viral_consensus", tools.viral_consensus, None, True, ""),
            ("blastn", tools.blastn, None, True, ""),
            ("makeblastdb", tools.makeblastdb, tools.blastn, True, ""),
            ("cliquesnv", tools.cliquesnv, None, True, CLIQUESNV_HINT),
        ]
    elif "16s" in arms:
        # The 16S arm needs only kallisto; the rest drive the optional
        # identity/coverage stage, which skips itself if they are absent.
        wanted.append(("kallisto", tools.kallisto, None, True, ""))
        wanted += [(name, requested, None, False,
                    "16S identity and coverage columns will be blank")
                   for name, requested in (("minimap2", tools.minimap2),
                                           ("samtools", tools.samtools),
                                           ("viral_consensus", tools.viral_consensus),
                                           ("blastn", tools.blastn))]

    cliquesnv_is_jar = False
    missing_tools: list[str] = []
    for name, requested, companion, required, hint in wanted:
        found, detail = _resolve(requested, companion)
        if found:
            lines.append(_line("ok", name, detail))
            if name == "cliquesnv" and Path(detail).suffix.lower() == ".jar":
                cliquesnv_is_jar = True
        elif required:
            lines.append(_line("MISSING", name, detail + (f"  ({hint})" if hint else "")))
            missing_tools.append(name)
            problems += 1
        else:
            lines.append(_line("note", name, f"not found: {hint}"))

    if cliquesnv_is_jar:
        found, detail = _resolve("java")
        lines.append(_line("ok" if found else "MISSING", "java",
                           detail if found else f"{detail}  (needed to launch the CliqueSNV JAR)"))
        problems += 0 if found else 1

    # Several tools missing at once is one fault, not several: the modules that
    # provide them are not loaded. Say that, rather than listing symptoms.
    if len(missing_tools) > 1:
        lines.append(_line("hint", "tools",
                           f"{len(missing_tools)} tools missing together - the modules that "
                           f"provide them are probably not loaded in this shell"))
        if tools.tools_env is not None:
            lines.append(_line("hint", "", f"source {tools.tools_env}"))
        lines.append(_line("hint", "", "then: module list   (to see what is loaded)"))
    elif problems and tools.tools_env is not None:
        lines.append(_line("hint", "tools",
                           f"source {tools.tools_env} (then re-run this check)"))

    required = set(settings.references.required_for(arms))
    for name, path in settings.references.items():
        if path is None:
            if name == "s16_copy_numbers":
                lines.append(_line("note", name,
                                   "not set: the 16S arm runs uncorrected (every taxon divided by 4.2)"))
            elif name in required:
                lines.append(_line("MISSING", name, f"not set in {settings.path}"))
                problems += 1
            continue
        if name == "amr_db":
            # A KMA index PREFIX, not a file: the index is <prefix>.name,
            # <prefix>.comp.b, ... beside the manifest build-amr-db wrote.
            present = Path(f"{path}.name").is_file()
            manifest = path.parent / "db_manifest.json"
            if present:
                lines.append(_line("ok", name, f"{path}  (manifest: "
                                   f"{'present' if manifest.is_file() else 'MISSING - rebuild with bridge build-amr-db'})"))
            elif name in required:
                lines.append(_line("MISSING", name, f"no KMA index at {path} - build it with: "
                                   f"bridge build-amr-db --db-dir <resfinder_db>"))
                problems += 1
            else:
                lines.append(_line("note", name, f"no KMA index yet (not needed for {'/'.join(arms)}): {path}"))
            continue
        if path.exists():
            lines.append(_line("ok", name, str(path)))
        elif name in required or name == "s16_copy_numbers":
            lines.append(_line("MISSING", name, f"not found: {path}"))
            problems += 1
        else:
            lines.append(_line("note", name, f"not found (not needed for {'/'.join(arms)}): {path}"))

    if tools.tools_env is not None:
        present = tools.tools_env.is_file()
        lines.append(_line("ok" if present else "MISSING", "tools_env", str(tools.tools_env)))
        if not present:
            problems += 1

    if scheduler == "pbs":
        qsub = shutil.which("qsub")
        if qsub:
            lines.append(_line("ok", "qsub", qsub))
        else:
            lines.append(_line("note", "qsub",
                               "not on PATH - PBS submission unavailable here; --scheduler local still works"))

    if out is not None:
        try:
            out.mkdir(parents=True, exist_ok=True)
            probe = out / ".bridge_write_test"
            probe.write_text("ok")
            probe.unlink()
            lines.append(_line("ok", "output", str(out)))
        except OSError as error:
            lines.append(_line("MISSING", "output", f"not writable: {out} ({error})"))
            problems += 1

    print(f"bridge check  (config: {settings.path})")
    for line in lines:
        print("  " + line)
    print()
    print("Everything needed is present." if not problems
          else f"{problems} problem(s) - fix the MISSING lines above.")
    return 1 if problems else 0
