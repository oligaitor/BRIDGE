"""Load, validate, and expose the pipeline configuration.

``config.yaml`` is the single source of truth for every threshold and the real
column/prefix names. This module loads it, checks that required keys exist (so a
typo fails fast rather than silently mis-thresholding a clinical call), and
exposes values by dotted path. The raw dict is preserved for the run manifest so
the exact effective configuration of every run is auditable.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

try:
    import yaml
except ModuleNotFoundError:  # pragma: no cover - exercised where PyYAML is absent
    yaml = None


class ConfigError(ValueError):
    """Raised when the configuration is missing required keys or is malformed."""


# ---------------------------------------------------------------------------
# Minimal YAML fallback
# ---------------------------------------------------------------------------
# The arms never need PyYAML, only this loader does, so where it is not
# installed this parses the small YAML subset the config actually uses: nested
# mappings by indentation, scalar values, empty inline mappings, and comments.
#
# It is deliberately STRICT. Anything outside that subset - a list, a block
# scalar, an anchor - raises instead of being guessed at. These values are
# clinical thresholds, and a parser that quietly mis-reads one is far worse than
# a parser that refuses to start. If the config ever grows a construct this
# cannot handle, the error names the line and the fix is to install PyYAML.

_SCALARS = {"true": True, "false": False, "yes": True, "no": False,
            "null": None, "~": None, "": None}
_INT_RE = re.compile(r"^[+-]?\d+$")
_FLOAT_RE = re.compile(r"^[+-]?(\d+\.\d*|\.\d+|\d+)([eE][+-]?\d+)?$")


def _strip_comment(line: str) -> str:
    """Remove a trailing ``#`` comment that is not inside quotes."""
    out, quote = [], None
    for index, char in enumerate(line):
        if quote:
            out.append(char)
            if char == quote:
                quote = None
        elif char in "\"'":
            quote = char
            out.append(char)
        elif char == "#" and (index == 0 or line[index - 1] in " \t"):
            break
        else:
            out.append(char)
    return "".join(out).rstrip()


def _scalar(text: str, lineno: int) -> Any:
    text = text.strip()
    if len(text) >= 2 and text[0] == text[-1] and text[0] in "\"'":
        return text[1:-1]
    if text == "{}":
        return {}
    lowered = text.lower()
    if lowered in _SCALARS:
        return _SCALARS[lowered]
    if _INT_RE.match(text):
        return int(text)
    if _FLOAT_RE.match(text):
        return float(text)
    if text.startswith(("[", "{", "&", "*", "|", ">")):
        raise ConfigError(
            f"line {lineno}: the built-in YAML reader supports only mappings "
            f"and scalars, got {text!r}. Install PyYAML to use this config."
        )
    return text


def _mini_yaml(text: str) -> dict:
    """Parse the mapping-and-scalars subset of YAML. Raises on anything else."""
    root: dict[str, Any] = {}
    # Each entry is (indent, mapping-at-that-indent).
    stack: list[tuple[int, dict]] = [(-1, root)]

    for lineno, raw in enumerate(text.splitlines(), start=1):
        line = _strip_comment(raw.replace("\t", "    "))
        if not line.strip():
            continue
        if line.lstrip().startswith("- "):
            raise ConfigError(
                f"line {lineno}: lists are not supported by the built-in YAML "
                f"reader. Install PyYAML to use this config."
            )

        indent = len(line) - len(line.lstrip())
        body = line.strip()
        if ":" not in body:
            raise ConfigError(f"line {lineno}: expected 'key: value', got {body!r}")

        key, _, value = body.partition(":")
        key, value = key.strip(), value.strip()
        if not key:
            raise ConfigError(f"line {lineno}: empty key")

        while stack and indent <= stack[-1][0]:
            stack.pop()
        if not stack:
            raise ConfigError(f"line {lineno}: inconsistent indentation")
        parent = stack[-1][1]

        if value == "":
            child: dict[str, Any] = {}
            parent[key] = child
            stack.append((indent, child))
        else:
            parent[key] = _scalar(value, lineno)
    return root


# Dotted keys that must be present for the pipeline to run safely. Optional
# tuning knobs (e.g. rmlst.loci_floor) are deliberately not listed.
REQUIRED_KEYS: tuple[str, ...] = (
    "taxonomy.genus_species_separator",
    "marker16s.genus_min_frac",
    "marker16s.species_dominance",
    "marker16s.species_margin",
    "rpm.enabled",
)

_MISSING = object()


class Config:
    """Read-only view over the loaded configuration dictionary.

    Access values by dotted path::

        cfg.get("rpm.rpm_r_min")                        # -> 10.0
        cfg.get("rmlst.loci_floor")                     # -> None
        cfg.get("does.not.exist", default=42)           # -> 42 (no error)

    ``cfg.require(path)`` is like ``get`` but raises ``ConfigError`` if the key
    is absent, for call sites that must not silently fall back to a default.
    """

    def __init__(self, data: dict[str, Any], source_path: Path | None = None) -> None:
        self._data = data
        self.source_path = source_path

    def get(self, dotted_key: str, default: Any = None) -> Any:
        node: Any = self._data
        for part in dotted_key.split("."):
            if not isinstance(node, dict) or part not in node:
                return default
            node = node[part]
        return node

    def require(self, dotted_key: str) -> Any:
        value = self.get(dotted_key, _MISSING)
        if value is _MISSING:
            raise ConfigError(f"Required config key missing: {dotted_key!r}")
        return value

    def override(self, dotted_key: str, value: Any) -> None:
        """Set a value by dotted path (for CLI overrides). Creates intermediate maps.

        The change is recorded in :meth:`as_dict`, so a CLI override is captured
        verbatim in each run's manifest.
        """
        parts = dotted_key.split(".")
        node = self._data
        for part in parts[:-1]:
            child = node.setdefault(part, {})
            if not isinstance(child, dict):
                raise ConfigError(f"Cannot override {dotted_key!r}: {part!r} is not a mapping")
            node = child
        node[parts[-1]] = value

    def as_dict(self) -> dict[str, Any]:
        """Deep-ish copy of the raw config, for writing into a run manifest."""
        import copy

        return copy.deepcopy(self._data)


def load_config(path: str | Path) -> Config:
    """Load and validate ``config.yaml`` from ``path``."""
    path = Path(path)
    if not path.is_file():
        raise ConfigError(f"Config file not found: {path}")

    text = path.read_text(encoding="utf-8")
    # PyYAML where available; the strict built-in reader otherwise, so the
    # reconciliation step needs nothing installed beyond the standard library.
    data = yaml.safe_load(text) if yaml is not None else _mini_yaml(text)

    if not isinstance(data, dict):
        raise ConfigError(f"Config root must be a mapping, got {type(data).__name__}")

    cfg = Config(data, source_path=path)
    _validate(cfg)
    return cfg


def _validate(cfg: Config) -> None:
    missing = [key for key in REQUIRED_KEYS if cfg.get(key, _MISSING) is _MISSING]
    if missing:
        raise ConfigError(
            "Config is missing required keys: " + ", ".join(sorted(missing))
        )

    # A handful of sanity checks that would otherwise surface as confusing
    # downstream errors or, worse, silently wrong calls.
    for numeric_key in (
        "marker16s.genus_min_frac",
        "marker16s.species_dominance",
        "marker16s.species_margin",
        "rpm.rpm_r_min",
        "rpm.rpm_min",
    ):
        value = cfg.get(numeric_key)
        if not isinstance(value, (int, float)):
            raise ConfigError(f"Config key {numeric_key!r} must be a number, got {value!r}")
    # Optional knobs: checked only when present, so a config written before
    # they existed still loads with the module defaults.
    for numeric_key in (
        "rpm.min_ntc_rpm",
        "rpm.abs_r_min",
        "rpm.min_sample_reads",
        "rpm.max_ntc_source_frac",
        "marker16s.min_genus_reads",
        "reconcile.multi_fallback_review_min",
    ):
        value = cfg.get(numeric_key, _MISSING)
        if value is _MISSING or value is None:
            continue
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ConfigError(f"Config key {numeric_key!r} must be a number, got {value!r}")
