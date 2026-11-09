from __future__ import annotations


class InsufficientMaterial(ValueError):
    """The sample carried too little rMLST evidence to call anything.

    This is a *result*, not a fault. A negative control, an extraction blank or
    a sample whose pathogen load sits below the assay's floor will legitimately
    produce no candidate alleles and no runtime clusters, and the pipeline
    should say so rather than crash. Raising this type lets the pipeline
    record the outcome and still write ``species_calls.csv``, which the batch
    driver treats as its completion sentinel.

    Subclasses ``ValueError`` so callers that catch ``ValueError`` still work.
    """

    def __init__(self, message: str, *, stage: str = "") -> None:
        super().__init__(message)
        self.stage = stage
