from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

from rna_assessment import MCAnnotateRunner, calculate_assessment, normalize_structure


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parent
    data_dir = root / "data"

    parser = argparse.ArgumentParser(
        description="Evaluate an RNA prediction against a reference structure.",
    )
    parser.add_argument(
        "reference",
        nargs="?",
        default=data_dir / "14_solution_0.pdb",
        type=Path,
        help="Reference PDB file.",
    )
    parser.add_argument(
        "prediction",
        nargs="?",
        default=data_dir / "14_ChenPostExp_2.pdb",
        type=Path,
        help="Predicted PDB file.",
    )
    parser.add_argument("--reference-index", type=Path, default=data_dir / "14_solution_0.index")
    parser.add_argument("--prediction-index", type=Path, default=data_dir / "14_ChenPostExp_2.index")
    parser.add_argument("--reference-annotation", type=Path)
    parser.add_argument("--prediction-annotation", type=Path)
    parser.add_argument("--mc-annotate", type=Path)
    parser.add_argument("--annotation-cache-dir", type=Path)
    parser.add_argument("--normalize-reference", action="store_true")
    parser.add_argument("--inclusion-radius", type=float, default=15.0)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    root = Path(__file__).resolve().parent
    output_dir = root / "output"
    output_dir.mkdir(exist_ok=True)

    original_reference_path = args.reference
    reference_path = original_reference_path
    if args.normalize_reference:
        reference_path = normalize_structure(
            original_reference_path,
            output_dir / f"{original_reference_path.stem}.normalized.pdb",
        )

    annotation_overrides = {}
    if args.reference_annotation:
        annotation_overrides[reference_path] = args.reference_annotation
    elif args.normalize_reference:
        original_annotation = original_reference_path.with_name(f"{original_reference_path.name}.mcout")
        if original_annotation.exists():
            annotation_overrides[reference_path] = original_annotation
    if args.prediction_annotation:
        annotation_overrides[args.prediction] = args.prediction_annotation

    annotator = MCAnnotateRunner(
        binary_path=args.mc_annotate,
        cache_dir=args.annotation_cache_dir,
        annotation_overrides=annotation_overrides or None,
    )
    result = calculate_assessment(
        reference_path,
        args.reference_index,
        args.prediction,
        args.prediction_index,
        annotator=annotator,
        inclusion_radius=args.inclusion_radius,
    )

    print(
        json.dumps(
            {
                "reference": str(reference_path),
                "prediction": str(args.prediction),
                "metrics": asdict(result),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
