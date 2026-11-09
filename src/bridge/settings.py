"""Everything a run needs to know, read from one ``config.yaml``.

Besides the arms' thresholds it holds ``references`` (the databases),
``tools`` (binary names) and ``hpc`` (the scheduler block: project, storage, queue,
size tiers, control-name prefixes).
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

from bridge.reconcile.baseline import DEFAULT_CONTROL_PREFIXES, is_control_name
from bridge.reconcile.config import Config, ConfigError, load_config

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = REPO_ROOT / "config.yaml"
ARMS = ("rmlst", "16s", "amr")
# The identification arms: what `--arms both` means.
ID_ARMS = ("rmlst", "16s")


@dataclass(frozen=True)
class References:
    rmlst_alleles: Path | None
    rmlst_index: Path | None
    rmlst_seqdef: Path | None
    s16_reference: Path | None
    s16_index: Path | None
    s16_copy_numbers: Path | None
    # AMR arm: the KMA index PREFIX (not a file) and ResFinder's phenotypes.txt
    amr_db: Path | None = None
    amr_phenotypes: Path | None = None

    def items(self) -> list[tuple[str, Path | None]]:
        return [(name, getattr(self, name)) for name in (
            "rmlst_alleles", "rmlst_index", "rmlst_seqdef",
            "s16_reference", "s16_index", "s16_copy_numbers",
            "amr_db", "amr_phenotypes")]

    def required_for(self, arms: tuple[str, ...]) -> list[str]:
        names = []
        if "rmlst" in arms:
            names += ["rmlst_alleles", "rmlst_index", "rmlst_seqdef"]
        if "16s" in arms:
            names += ["s16_reference", "s16_index"]
        if "amr" in arms:
            names += ["amr_db", "amr_phenotypes"]
        return names


@dataclass(frozen=True)
class Tools:
    kallisto: str = "kallisto"
    mash: str = "mash"
    samtools: str = "samtools"
    minimap2: str = "minimap2"
    bwa_mem2: str = "bwa-mem2"
    blastn: str = "blastn"
    makeblastdb: str = "makeblastdb"
    viral_consensus: str = "viral_consensus"
    # CliqueSNV may be a JAR rather than a program on PATH, located through
    # BRIDGE_CLIQUESNV; the rMLST arm (bridge.rmlst.models.Config) honours that
    # variable, so this layer must too.
    cliquesnv: str = "cliquesnv"
    kma: str = "kma"
    tools_env: Path | None = None


@dataclass(frozen=True)
class Tier:
    name: str
    max_bytes: int | None      # None = no upper bound
    walltime: str
    mem_gb: int
    ncpus: int
    jobfs_gb: int


DEFAULT_TIERS = (
    Tier("small", 10_000_000, "04:00:00", 16, 4, 20),
    Tier("medium", 100_000_000, "06:00:00", 32, 8, 20),
    Tier("large", 1_000_000_000, "12:00:00", 64, 16, 40),
    Tier("xlarge", None, "12:00:00", 100, 25, 100),
)


@dataclass(frozen=True)
class Hpc:
    project: str
    queue: str
    storage: str
    python_module: str | None
    controls_prefixes: tuple[str, ...]
    cliquesnv_time_limit_seconds: int | None
    tiers: tuple[Tier, ...]
    finalise_walltime: str
    finalise_mem_gb: int
    finalise_ncpus: int


@dataclass(frozen=True)
class Settings:
    cfg: Config
    path: Path
    references: References
    tools: Tools
    hpc: Hpc


def _path(value, base: Path | None = None) -> Path | None:
    if value in (None, "", False):
        return None
    path = Path(str(value)).expanduser()
    if base is not None and not path.is_absolute():
        path = base / path
    return path


def _prefixes(value) -> tuple[str, ...]:
    if value in (None, ""):
        return DEFAULT_CONTROL_PREFIXES
    if isinstance(value, str):
        return tuple(value.split())
    return tuple(str(item) for item in value)


def _tiers(cfg: Config) -> tuple[Tier, ...]:
    raw = cfg.get("hpc.tiers")
    if not raw:
        return DEFAULT_TIERS
    if not isinstance(raw, dict):
        raise ConfigError("hpc.tiers must be a mapping of tier name -> settings")
    tiers = []
    for name, spec in raw.items():
        if not isinstance(spec, dict):
            raise ConfigError(f"hpc.tiers.{name} must be a mapping")
        cap = spec.get("max_bytes")
        tiers.append(Tier(
            name=str(name),
            max_bytes=None if cap in (None, 0) else int(cap),
            walltime=str(spec.get("walltime", "12:00:00")),
            mem_gb=int(spec.get("mem_gb", 16)),
            ncpus=int(spec.get("ncpus", 4)),
            jobfs_gb=int(spec.get("jobfs_gb", 20)),
        ))
    tiers.sort(key=lambda tier: (tier.max_bytes is None, tier.max_bytes or 0))
    return tuple(tiers)


def _hpc(cfg: Config) -> Hpc:
    limit = cfg.get("hpc.cliquesnv_time_limit_seconds")
    module = cfg.get("hpc.python_module")
    return Hpc(
        project=str(cfg.get("hpc.project", "PROJECT")),
        queue=str(cfg.get("hpc.queue", "normal")),
        storage=str(cfg.get("hpc.storage", "")),
        python_module=str(module) if module else None,
        controls_prefixes=_prefixes(cfg.get("hpc.controls_prefixes")),
        cliquesnv_time_limit_seconds=int(limit) if limit else None,
        tiers=_tiers(cfg),
        finalise_walltime=str(cfg.get("hpc.finalise.walltime", "02:00:00")),
        finalise_mem_gb=int(cfg.get("hpc.finalise.mem_gb", 8)),
        finalise_ncpus=int(cfg.get("hpc.finalise.ncpus", 1)),
    )


def load_settings(path: str | Path | None = None) -> Settings:
    config_path = Path(path).expanduser().resolve() if path else DEFAULT_CONFIG
    cfg = load_config(config_path)
    base = config_path.parent
    references = References(
        rmlst_alleles=_path(cfg.get("references.rmlst_alleles")),
        rmlst_index=_path(cfg.get("references.rmlst_index")),
        rmlst_seqdef=_path(cfg.get("references.rmlst_seqdef")),
        s16_reference=_path(cfg.get("references.s16_reference")),
        s16_index=_path(cfg.get("references.s16_index")),
        s16_copy_numbers=_path(cfg.get("references.s16_copy_numbers")),
        amr_db=_path(cfg.get("references.amr_db")),
        amr_phenotypes=_path(cfg.get("references.amr_phenotypes")),
    )
    tools = Tools(
        kallisto=str(cfg.get("tools.kallisto", "kallisto")),
        mash=str(cfg.get("tools.mash", "mash")),
        samtools=str(cfg.get("tools.samtools", "samtools")),
        minimap2=str(cfg.get("tools.minimap2", "minimap2")),
        bwa_mem2=str(cfg.get("tools.bwa_mem2", "bwa-mem2")),
        blastn=str(cfg.get("tools.blastn", "blastn")),
        makeblastdb=str(cfg.get("tools.makeblastdb", "makeblastdb")),
        viral_consensus=str(cfg.get("tools.viral_consensus", "viral_consensus")),
        # config.yaml wins; then BRIDGE_CLIQUESNV; then the bare name.
        cliquesnv=str(cfg.get("tools.cliquesnv")
                      or os.environ.get("BRIDGE_CLIQUESNV")
                      or "cliquesnv"),
        kma=str(cfg.get("tools.kma") or "kma"),
        tools_env=_path(cfg.get("tools.tools_env"), base),
    )
    return Settings(cfg=cfg, path=config_path, references=references,
                    tools=tools, hpc=_hpc(cfg))


def tier_for(hpc: Hpc, size_bytes: int) -> Tier:
    for tier in hpc.tiers:
        if tier.max_bytes is None or size_bytes < tier.max_bytes:
            return tier
    return hpc.tiers[-1]


def is_control(name: str, prefixes: tuple[str, ...]) -> bool:
    """Case-insensitive prefix match (``NEG`` catches Neg_, NEG_ and NEGMATRIX_)."""
    return is_control_name(name, prefixes)


def arms_from(value: str) -> tuple[str, ...]:
    """``all`` | ``both`` (rMLST + 16S) | one arm | a comma list (``rmlst,amr``).

    Always returned in pipeline order, so the same set of arms renders the
    same job script however it was typed.
    """
    value = (value or "").strip()
    if value in ("all", ""):
        return ARMS
    if value == "both":
        return ID_ARMS
    wanted = {part.strip() for part in value.split(",") if part.strip()}
    unknown = sorted(wanted - set(ARMS))
    if unknown or not wanted:
        raise ValueError(f"--arms must be all, both, or any of {', '.join(ARMS)} "
                         f"(comma-separated), got {value!r}")
    return tuple(arm for arm in ARMS if arm in wanted)


def arms_arg(arms: tuple[str, ...]) -> str:
    """The inverse of :func:`arms_from`, for a job script's ``--arms``."""
    arms = tuple(arm for arm in ARMS if arm in arms)
    if arms == ARMS:
        return "all"
    if arms == ID_ARMS:
        return "both"
    return ",".join(arms)
