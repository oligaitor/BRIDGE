"""BRIDGE's reconcile step: one call set per sample from the rMLST and 16S arms.

The rMLST arm (allele profiles across ~53 ribosomal protein loci, species
called locally) governs: if it resolved anything the call set is its own, and
the 16S marker arm (one kallisto pass against the curated MIMt database)
produces calls only for a sample where rMLST resolved nothing. The 16S arm
still runs on every sample and its numbers still reach the report; organisms
only 16S saw appear as labelled, non-reported rows, because the rMLST
reference covers a restricted set of target genera. Contamination is judged by
measurement rather than confidence tiers - the RPM-r and abs-r ratios against
the run's own negative control (:mod:`bridge.reconcile.rpm`), with
per-organism read counts and per-organism coverage of the loci behind each
call.
"""
