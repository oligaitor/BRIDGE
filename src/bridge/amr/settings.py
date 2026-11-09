"""The AMR arm's slice of ``config.yaml``.

Read from the same file as the identification arms, so one config describes
a whole run, whether the arm is run from ``bridge amr`` or from
``bridge sample``.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bridge.reconcile.baseline import DEFAULT_CONTROL_PREFIXES
from bridge.reconcile.config import Config, load_config
from bridge.reconcile.rpm import (
    DEFAULT_ABS_R_MIN, DEFAULT_MIN_SAMPLE_READS, DEFAULT_RPM_MIN, DEFAULT_RPM_R_MIN,
    MAX_NTC_SOURCE_FRAC, MIN_CONTROL_READS, MIN_NTC_RPM,
)

DATA_DIR = Path(__file__).resolve().parent / "data"
REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_CONFIG = REPO_ROOT / "config.yaml"

# What `bridge amr build-db` writes and what the arm then reads. The index
# prefix is `references.amr_db`; the manifest sits beside it.
DB_MANIFEST_NAME = "db_manifest.json"

# What config.yaml's amr.kma_args documents; used when that key is absent.
DEFAULT_KMA_ARGS = "-1t1 -mem_mode -ef -and -matrix -nf"


@dataclass(frozen=True)
class Thresholds:
    """When a KMA template counts as detected, and how firmly.

    The floor is ResFinder's own default for reads (identity 0.90, coverage
    0.60); everything is provisional until the blood-culture AST validation.
    """

    min_identity: float = 90.0        # Query_Identity (identity where covered), percent
    min_coverage: float = 60.0        # Template_Coverage, percent
    confident_coverage: float = 80.0
    min_depth: float = 2.0            # KMA Depth, for `confident`
    min_fragments: float = 10.0       # mapstat fragmentCount, for `confident`
    exact_identity: float = 100.0     # identity AND coverage for an `exact` allele


@dataclass(frozen=True)
class GateSettings:
    enabled: bool = True
    rpm_r_min: float = DEFAULT_RPM_R_MIN
    rpm_min: float = DEFAULT_RPM_MIN
    min_ntc_rpm: float = MIN_NTC_RPM
    abs_r_min: float = DEFAULT_ABS_R_MIN
    min_sample_reads: float = DEFAULT_MIN_SAMPLE_READS
    min_control_reads: int = MIN_CONTROL_READS
    max_ntc_source_frac: float = MAX_NTC_SOURCE_FRAC


@dataclass(frozen=True)
class AmrSettings:
    kma: str
    db: Path | None                 # KMA index prefix
    phenotypes: Path | None         # ResFinder phenotypes.txt
    kma_args: str
    thresholds: Thresholds
    gate: GateSettings
    host_context: Path
    eucast_table: Path
    family_overrides: Path
    point_mutations: bool = False
    control_prefixes: tuple[str, ...] = DEFAULT_CONTROL_PREFIXES
    config_path: Path | None = None

    @property
    def db_manifest(self) -> Path | None:
        return self.db.parent / DB_MANIFEST_NAME if self.db is not None else None


def _path(value, base: Path | None = None) -> Path | None:
    if value in (None, "", False):
        return None
    path = Path(str(value)).expanduser()
    if base is not None and not path.is_absolute():
        path = base / path
    return path


def _prefixes(value) -> tuple[str, ...]:
    if value in (None, ""):
        return AmrSettings.control_prefixes
    if isinstance(value, str):
        return tuple(value.split())
    return tuple(str(item) for item in value)


def _data(cfg: Config, key: str, default_name: str, base: Path | None) -> Path:
    """A curated table: the config's path if set, else the one shipped here."""
    return _path(cfg.get(key), base) or DATA_DIR / default_name


def from_config(cfg: Config, path: Path | None = None) -> AmrSettings:
    base = Path(path).parent if path else None
    t = Thresholds()
    thresholds = Thresholds(
        min_identity=float(cfg.get("amr.min_identity", t.min_identity)),
        min_coverage=float(cfg.get("amr.min_coverage", t.min_coverage)),
        confident_coverage=float(cfg.get("amr.confident_coverage", t.confident_coverage)),
        min_depth=float(cfg.get("amr.min_depth", t.min_depth)),
        min_fragments=float(cfg.get("amr.min_fragments", t.min_fragments)),
        exact_identity=float(cfg.get("amr.exact_identity", t.exact_identity)),
    )
    g = GateSettings()
    # The gate reuses the identification side's rpm.* thresholds - the same
    # 10x convention, the same guards - unless amr.gate overrides one.
    def gate_value(key: str, default):
        return cfg.get(f"amr.gate.{key}", cfg.get(f"rpm.{key}", default))
    gate = GateSettings(
        enabled=bool(cfg.get("amr.gate.enabled", True)),
        rpm_r_min=float(gate_value("rpm_r_min", g.rpm_r_min)),
        rpm_min=float(gate_value("rpm_min", g.rpm_min)),
        min_ntc_rpm=float(gate_value("min_ntc_rpm", g.min_ntc_rpm)),
        abs_r_min=float(gate_value("abs_r_min", g.abs_r_min)),
        min_sample_reads=float(gate_value("min_sample_reads", g.min_sample_reads)),
        min_control_reads=int(cfg.get("amr.gate.min_control_reads", g.min_control_reads)),
        max_ntc_source_frac=float(gate_value("max_ntc_source_frac", g.max_ntc_source_frac)),
    )
    return AmrSettings(
        kma=str(cfg.get("tools.kma") or "kma"),
        db=_path(cfg.get("references.amr_db")),
        phenotypes=_path(cfg.get("references.amr_phenotypes")),
        kma_args=str(cfg.get("amr.kma_args") or DEFAULT_KMA_ARGS),
        thresholds=thresholds,
        gate=gate,
        host_context=_data(cfg, "amr.host_context", "gene_host_context.tsv", base),
        eucast_table=_data(cfg, "amr.eucast_table", "eucast_expected_resistant.tsv", base),
        family_overrides=_data(cfg, "amr.family_overrides", "family_overrides.tsv", base),
        point_mutations=bool(cfg.get("amr.point_mutations", False)),
        control_prefixes=_prefixes(cfg.get("hpc.controls_prefixes")),
        config_path=Path(path) if path else None,
    )


def load_amr_settings(path: str | Path | None = None) -> AmrSettings:
    config_path = Path(path).expanduser().resolve() if path else DEFAULT_CONFIG
    return from_config(load_config(config_path), config_path)
