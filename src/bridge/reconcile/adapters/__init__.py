"""Adapters: turn each arm's CSV into the typed lists the reconciler consumes.

The reconciler is source-agnostic: it decides on the strength enum, the genus
and the species name - never on a filename, a column or a numeric field. These
adapters build ``list[RmlstStrengthResult]`` and ``list[BlastStrengthResult]``
from the rMLST and 16S arms' call tables.
"""

from .rmlst import load_rmlst_calls, rmlst_results
from .marker16s import load_marker16s_calls, marker16s_results

__all__ = [
    "load_rmlst_calls",
    "rmlst_results",
    "load_marker16s_calls",
    "marker16s_results",
]
