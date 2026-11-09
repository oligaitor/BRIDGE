"""BRIDGE's resistance (AMR) arm: ResFinder alleles called from the same reads by KMA.

The arm (:mod:`bridge.amr.kma`) runs per sample on a compute node and writes
``<SAMPLE>/amr/amr_hits.csv``, carrying everything database-derived (family,
antibiotics, whether a family's alleles disagree about spectrum) so the report
never needs the ResFinder files; the report (:mod:`bridge.amr.report`) is the
cheap, re-runnable step - negative-control gate, organism attribution and
predicted phenotype. Only ResFinder's *acquired* genes are captured, so
point-mutation resistance (gyrA/parC, mgrB, porin loss) is invisible,
``amr.point_mutations`` stays false, and an absent gene is never read as
susceptibility.
"""
