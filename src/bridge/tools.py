from __future__ import annotations

import logging
import shutil
import subprocess
import threading
from collections import defaultdict
from contextlib import contextmanager
from pathlib import Path
from time import perf_counter
from typing import Iterator, Sequence


# ---------------------------------------------------------------------------
# Runtime accounting
# ---------------------------------------------------------------------------
# Stage timings alone cannot attribute cost (cluster_consensus wraps the
# aligner, CliqueSNV, ViralConsensus and BLAST), so every subprocess is timed
# here and the totals are written to the run manifest.

_ACCOUNTING_LOCK = threading.Lock()
_TOOL_SECONDS: dict[str, float] = defaultdict(float)
_TOOL_CALLS: dict[str, int] = defaultdict(int)
_STAGE_SECONDS: dict[str, float] = {}


def reset_runtime_accounting() -> None:
    with _ACCOUNTING_LOCK:
        _TOOL_SECONDS.clear()
        _TOOL_CALLS.clear()
        _STAGE_SECONDS.clear()


def _record_tool(tool: str, seconds: float) -> None:
    with _ACCOUNTING_LOCK:
        _TOOL_SECONDS[tool] += seconds
        _TOOL_CALLS[tool] += 1


def _record_stage(name: str, seconds: float) -> None:
    with _ACCOUNTING_LOCK:
        _STAGE_SECONDS[name] = round(seconds, 3)


def runtime_accounting() -> dict[str, dict[str, object]]:
    """Per-tool and per-stage totals for the run so far.

    Tool seconds are wall time summed across calls, so with concurrent clusters
    they can exceed the stage's own elapsed time. That is the intended reading:
    it is the work done, not the time spent waiting for it.
    """
    with _ACCOUNTING_LOCK:
        tools = {
            tool: {
                "calls": _TOOL_CALLS[tool],
                "seconds": round(_TOOL_SECONDS[tool], 3),
            }
            for tool in sorted(_TOOL_SECONDS)
        }
        return {"stages": dict(_STAGE_SECONDS), "tools": tools}


class _ConsoleFilter(logging.Filter):
    def __init__(self, verbose: bool) -> None:
        super().__init__()
        self.verbose = verbose

    def filter(self, record: logging.LogRecord) -> bool:
        return (
            self.verbose
            or record.levelno >= logging.WARNING
            or bool(getattr(record, "console", False))
        )


def configure_logging(log_path: Path, verbose: bool = False) -> tuple[logging.Handler, ...]:
    """Write all INFO records to a run log and selected records to stderr."""
    log_path.parent.mkdir(parents=True, exist_ok=True)
    root = logging.getLogger()
    root.setLevel(logging.INFO)
    for handler in tuple(root.handlers):
        if getattr(handler, "_rmlst_handler", False):
            root.removeHandler(handler)
            handler.close()

    formatter = logging.Formatter(
        "%(asctime)s [%(name)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    file_handler = logging.FileHandler(log_path, mode="w", encoding="utf-8")
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(formatter)
    file_handler._rmlst_handler = True  # type: ignore[attr-defined]

    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    console_handler.addFilter(_ConsoleFilter(verbose))
    console_handler._rmlst_handler = True  # type: ignore[attr-defined]

    handlers: tuple[logging.Handler, ...] = (file_handler, console_handler)
    for handler in handlers:
        root.addHandler(handler)
    return handlers


def close_logging(handlers: Sequence[logging.Handler]) -> None:
    root = logging.getLogger()
    for handler in handlers:
        root.removeHandler(handler)
        handler.flush()
        handler.close()


@contextmanager
def stage(name: str) -> Iterator[None]:
    logger = logging.getLogger("bridge.rmlst")
    started = perf_counter()
    logger.info("stage=%s status=started", name, extra={"console": True})
    try:
        yield
    except Exception:
        elapsed = perf_counter() - started
        _record_stage(name, elapsed)
        logger.error(
            "stage=%s status=failed runtime_seconds=%.3f",
            name, elapsed, extra={"console": True},
        )
        raise
    elapsed = perf_counter() - started
    _record_stage(name, elapsed)
    logger.info(
        "stage=%s status=completed runtime_seconds=%.3f",
        name, elapsed, extra={"console": True},
    )


def executable(value: str) -> str:
    resolved = shutil.which(value)
    if resolved is None:
        path = Path(value).expanduser()
        if not path.is_file():
            raise FileNotFoundError(f"Executable not found: {value}")
        return str(path.resolve())
    return resolved


def companion_executable(value: str, companion: str) -> str:
    """Resolve a tool, falling back to the installed directory of a companion."""
    try:
        return executable(value)
    except FileNotFoundError:
        if Path(value).name != value:
            raise
        sibling = Path(executable(companion)).with_name(value)
    return executable(str(sibling))


def log_tool_output(tool: str, *, stdout: str = "", stderr: str = "") -> None:
    """Write captured tool output to the run log, preserving its source stream."""
    logger = logging.getLogger(tool)
    for stream, output in (("stdout", stdout), ("stderr", stderr)):
        for line in output.replace("\r", "\n").splitlines():
            if line.strip():
                logger.info("%s: %s", stream, line)


def run(command: Sequence[str], tool: str, *, stdout_path: Path | None = None) -> None:
    logger = logging.getLogger(tool)
    logger.info("command=%s", " ".join(map(str, command)))
    started = perf_counter()
    try:
        _run(command, tool, logger, stdout_path)
    finally:
        elapsed = perf_counter() - started
        _record_tool(tool, elapsed)
        logger.info("tool=%s runtime_seconds=%.3f", tool, elapsed)


def _run(
    command: Sequence[str],
    tool: str,
    logger: logging.Logger,
    stdout_path: Path | None,
) -> None:
    if stdout_path is not None:
        stderr_lines: list[str] = []
        stdout_path.parent.mkdir(parents=True, exist_ok=True)
        with stdout_path.open("wb") as output:
            process = subprocess.Popen(
                list(command), stdout=output, stderr=subprocess.PIPE, text=True,
            )
            assert process.stderr is not None
            for line in process.stderr:
                stderr_lines.append(line)
                for message in line.replace("\r", "\n").splitlines():
                    if message.strip():
                        logger.info("%s", message)
            returncode = process.wait()
        if returncode:
            raise subprocess.CalledProcessError(
                returncode, list(command), stderr="".join(stderr_lines),
            )
    else:
        output_lines: list[str] = []
        process = subprocess.Popen(
            list(command), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
        )
        assert process.stdout is not None
        for line in process.stdout:
            output_lines.append(line)
            for message in line.replace("\r", "\n").splitlines():
                if message.strip():
                    logger.info("%s", message)
        returncode = process.wait()
        if returncode:
            raise subprocess.CalledProcessError(
                returncode, list(command), output="".join(output_lines),
            )
