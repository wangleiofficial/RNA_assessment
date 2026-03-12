from __future__ import annotations

from pathlib import Path

from .alignment import ChainAlignment, StructureAlignment, StructureMatcher, infer_structure_alignment
from .benchmark import (
    BenchmarkEntry,
    BenchmarkJob,
    BenchmarkResult,
    ChainMappingResult,
    build_benchmark_jobs,
    describe_prepared_pair,
    load_benchmark_manifest,
    run_benchmark,
)
from .extraction import extract_pdb
from .metrics import (
    AssessmentResult,
    InteractionNetworkResult,
    LDDTResult,
    PreparedStructurePair,
    ResidueAssessment,
    RMSDResult,
    calculate_assessment,
    calculate_assessment_from_prepared,
    calculate_interaction_network_fidelity,
    calculate_lddt,
    calculate_rmsd,
    prepare_structure_pair,
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
    "BenchmarkEntry",
    "BenchmarkJob",
    "BenchmarkResult",
    "ChainAlignment",
    "ChainMappingResult",
    "InteractionNetworkResult",
    "LDDTResult",
    "PreparedStructurePair",
    "ResidueAssessment",
    "RMSDResult",
    "StructureAlignment",
    "StructureMatcher",
    "build_benchmark_jobs",
    "calculate_assessment",
    "calculate_assessment_from_prepared",
    "calculate_interaction_network_fidelity",
    "calculate_lddt",
    "calculate_rmsd",
    "describe_prepared_pair",
    "extract_structure",
    "infer_structure_alignment",
    "load_benchmark_manifest",
    "normalize_structure",
    "prepare_structure_pair",
    "run_benchmark",
]
