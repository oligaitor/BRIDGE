"""16S marker arm: abundance and taxonomic assignment against a curated 16S database.

An independent second opinion to the rMLST arm. Where rMLST resolves species by
assembling allele profiles across ~53 ribosomal protein loci, this arm does one
kallisto pass against MIMt and reads the taxonomy straight off the reference
header. No clustering, no consensus reconstruction, no haplotype calling - those
are what make the rMLST arm slow, and none of them help a single marker gene
that cannot separate close relatives however carefully it is assembled.

The two arms are deliberately independent: different reference database,
different algorithm, different failure modes. Agreement between them is
therefore evidence, which is the whole point of running both.
"""

from .reference import ReferenceTaxon, parse_lineage_header, load_reference
from .quant import quant_16s
from .aggregate import aggregate_taxa, TaxonAbundance
from .calls import marker16s_calls, write_marker16s_calls

__all__ = [
    "ReferenceTaxon",
    "parse_lineage_header",
    "load_reference",
    "quant_16s",
    "aggregate_taxa",
    "TaxonAbundance",
    "marker16s_calls",
    "write_marker16s_calls",
]
