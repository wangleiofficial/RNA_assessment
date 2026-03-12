from __future__ import annotations


class RNAAssessmentError(Exception):
    """Base class for project-specific exceptions."""


class InvalidResidueRangeError(RNAAssessmentError):
    """Raised when a residue range specification is malformed."""


class NormalizationError(RNAAssessmentError):
    """Raised when a PDB file cannot be normalized safely."""


class SequenceMismatchError(RNAAssessmentError):
    """Raised when two indexed structures do not describe the same sequence."""


class ToolNotAvailableError(RNAAssessmentError):
    """Raised when an optional third-party tool cannot be resolved or executed."""


class MetricCalculationError(RNAAssessmentError):
    """Raised when a metric cannot be computed from the provided structures."""
