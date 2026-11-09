"""Shared exception taxonomy mapping to the agreed failure policy.

* :class:`SampleError` — a problem with one sample's data. The orchestrator
  flags the sample and continues the batch.
* :class:`InfraError` — an infrastructure failure (a reference database
  or a required binary missing). The orchestrator stops the run with a
  non-zero exit (fail-closed), because it will affect every sample.
"""

from __future__ import annotations


class PipelineError(Exception):
    """Base class for pipeline errors."""


class SampleError(PipelineError):
    """A per-sample data problem: flag and continue."""


class InfraError(PipelineError):
    """An infrastructure failure: fail-closed, stop the run."""
