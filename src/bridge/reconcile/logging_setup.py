"""Structured logging: file + console, ``--log-level`` configurable, sample-tagged.

Every record carries a ``sample`` field (``-`` when not in a per-sample
context) so the log can be filtered/grepped by sample ID.
"""

from __future__ import annotations

import logging
from pathlib import Path

LOGGER_NAME = "bridge.reconcile"


class _SampleFilter(logging.Filter):
    def filter(self, record: logging.LogRecord) -> bool:
        if not hasattr(record, "sample"):
            record.sample = "-"
        return True


def setup_logging(level: str = "INFO", log_file: str | Path | None = None) -> logging.Logger:
    """Configure and return the package logger. Idempotent (clears old handlers)."""
    logger = logging.getLogger(LOGGER_NAME)
    logger.setLevel(level.upper())
    logger.handlers.clear()
    logger.propagate = False

    fmt = logging.Formatter("%(asctime)s %(levelname)-7s [%(sample)s] %(message)s")
    sample_filter = _SampleFilter()

    console = logging.StreamHandler()
    console.setFormatter(fmt)
    console.addFilter(sample_filter)
    logger.addHandler(console)

    if log_file is not None:
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file, encoding="utf-8")
        file_handler.setFormatter(fmt)
        file_handler.addFilter(sample_filter)
        logger.addHandler(file_handler)

    return logger


def sample_logger(sample_id: str) -> logging.LoggerAdapter:
    """A logger adapter that tags every record with ``sample_id``."""
    return logging.LoggerAdapter(logging.getLogger(LOGGER_NAME), {"sample": sample_id})
