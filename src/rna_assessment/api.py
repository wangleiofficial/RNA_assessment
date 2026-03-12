from __future__ import annotations

from pathlib import Path

from .extraction import extract_pdb
from .metrics import (
    AssessmentResult,
    InteractionNetworkResult,
    LDDTResult,
    RMSDResult,
    calculate_assessment,
    calculate_interaction_network_fidelity,
    calculate_lddt,
    calculate_rmsd,
)
from .normalization import PDBNormalizer


def normalize_structure(
    input_file: str | Path,
    output_file: str | Path,
    normalizer: PDBNormalizer | None = None,
) -> Path:
    normalizer = normalizer or PDBNormalizer.from_defaults()
    return normalizer.normalize_or_raise(input_file, output_file)


def extract_structure(
    input_file: str | Path,
    residue_ranges: str,
    output_file: str | Path,
) -> Path:
    return extract_pdb(input_file, residue_ranges, output_file)


__all__ = [
    "AssessmentResult",
    "InteractionNetworkResult",
    "LDDTResult",
    "RMSDResult",
    "calculate_assessment",
    "calculate_interaction_network_fidelity",
    "calculate_lddt",
    "calculate_rmsd",
    "extract_structure",
    "normalize_structure",
]
