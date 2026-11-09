from __future__ import annotations

from dataclasses import dataclass
from importlib import metadata
from pathlib import Path
import re
import sys

from .models import Config
from bridge.tools import companion_executable, executable


PYTHON_REQUIREMENTS = {
    "pandas": (2, 0),
    "numpy": (1, 24),
    "matplotlib": (3, 7),
    "seaborn": (0, 13),
}


@dataclass(frozen=True)
class DependencyStatus:
    name: str
    requirement: str
    resolved: str


class DependencyError(RuntimeError):
    """Raised when one or more required dependencies are unavailable."""

    def __init__(self, failures: tuple[str, ...]):
        self.failures = failures
        super().__init__(
            "Dependency check failed:\n  - " + "\n  - ".join(failures)
        )


def _numeric_version(version: str) -> tuple[int, ...]:
    parts: list[int] = []
    for component in version.split("."):
        match = re.match(r"\d+", component)
        if match is None:
            break
        parts.append(int(match.group()))
    return tuple(parts)


def check_dependencies(
    *,
    kallisto: str = "kallisto",
    mash: str = "mash",
    short_read_aligner: str = "minimap2",
    bwa_mem2: str = "bwa-mem2",
    minimap2: str = "minimap2",
    samtools: str = "samtools",
    cliquesnv: str = "cliquesnv",
    viral_consensus: str = "viral_consensus",
    viralconsensus_only: bool = False,
    blastn: str = "blastn",
    makeblastdb: str = "makeblastdb",
) -> tuple[DependencyStatus, ...]:
    """Resolve every dependency needed by the selected pipeline route."""
    if short_read_aligner not in {"bwa-mem2", "minimap2"}:
        raise ValueError("short-read aligner must be bwa-mem2 or minimap2")

    statuses: list[DependencyStatus] = []
    failures: list[str] = []

    python_requirement = ">=3.10"
    python_version = ".".join(map(str, sys.version_info[:3]))
    if sys.version_info < (3, 10):
        failures.append(
            f"python {python_requirement}: found {python_version}"
        )
    else:
        statuses.append(DependencyStatus("python", python_requirement, python_version))

    for package, minimum in PYTHON_REQUIREMENTS.items():
        requirement = f">={'.'.join(map(str, minimum))}"
        try:
            version = metadata.version(package)
        except metadata.PackageNotFoundError:
            failures.append(f"{package} {requirement}: package is not installed")
            continue
        if _numeric_version(version) < minimum:
            failures.append(f"{package} {requirement}: found {version}")
            continue
        statuses.append(DependencyStatus(package, requirement, version))

    commands = [
        ("kallisto", kallisto),
        ("mash", mash),
        ("samtools", samtools),
        ("viral_consensus", viral_consensus),
        ("blastn", blastn),
        (
            short_read_aligner,
            minimap2 if short_read_aligner == "minimap2" else bwa_mem2,
        ),
    ]
    if not viralconsensus_only:
        commands.extend((("cliquesnv", cliquesnv), ("java", "java")))

    for name, requested in commands:
        try:
            resolved = executable(requested)
        except (FileNotFoundError, OSError) as error:
            failures.append(f"{name} ({requested}): {error}")
        else:
            statuses.append(DependencyStatus(name, requested, resolved))

    try:
        resolved_makeblastdb = companion_executable(makeblastdb, blastn)
    except (FileNotFoundError, OSError) as error:
        failures.append(f"makeblastdb ({makeblastdb}): {error}")
    else:
        statuses.append(
            DependencyStatus("makeblastdb", makeblastdb, resolved_makeblastdb)
        )

    if not viralconsensus_only:
        cliquesnv_status = next(
            (status for status in statuses if status.name == "cliquesnv"), None
        )
        if cliquesnv_status is not None and Path(cliquesnv_status.resolved).suffix.lower() == ".jar":
            java_available = any(status.name == "java" for status in statuses)
            if not java_available and not any(
                failure.startswith("java ") for failure in failures
            ):
                failures.append("java (java): required to launch the CliqueSNV JAR")

    if failures:
        raise DependencyError(tuple(failures))
    return tuple(statuses)


def check_config_dependencies(config: Config) -> tuple[DependencyStatus, ...]:
    return check_dependencies(
        kallisto=config.kallisto,
        mash=config.mash,
        short_read_aligner=config.short_read_aligner,
        bwa_mem2=config.bwa_mem2,
        minimap2=config.minimap2,
        samtools=config.samtools,
        cliquesnv=config.cliquesnv,
        viral_consensus=config.viral_consensus,
        viralconsensus_only=config.viralconsensus_only,
        blastn=config.blastn,
        makeblastdb=config.makeblastdb,
    )
