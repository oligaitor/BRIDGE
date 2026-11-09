"""BRIDGE: dual-marker (rMLST + 16S) bacterial identification from targeted NGS.

Two arms (``bridge.rmlst`` for rMLST, ``bridge.marker16s`` for 16S), a
negative-control baseline, per-sample reconciliation (``bridge.reconcile``), an
optional AMR arm (``bridge.amr``) and a run-level table, all driven by the
``bridge`` command.
"""

__version__ = "1.0.0"
