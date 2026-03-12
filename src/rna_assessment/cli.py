from __future__ import annotations

import argparse
import json
from dataclasses import asdict

from .api import extract_structure, normalize_structure
from .exceptions import RNAAssessmentError
from .mc_annotate import MCAnnotateRunner
from .metrics import calculate_assessment, calculate_interaction_network_fidelity, calculate_lddt, calculate_rmsd


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="rna-assessment")
    subparsers = parser.add_subparsers(dest="command", required=True)

    normalize_parser = subparsers.add_parser("normalize", help="Normalize a PDB file.")
    normalize_parser.add_argument("input")
    normalize_parser.add_argument("output")

    extract_parser = subparsers.add_parser("extract", help="Extract indexed residues from a PDB file.")
    extract_parser.add_argument("input")
    extract_parser.add_argument("residues")
    extract_parser.add_argument("output")

    rmsd_parser = subparsers.add_parser("rmsd", help="Calculate RMSD and p-value.")
    rmsd_parser.add_argument("native")
    rmsd_parser.add_argument("native_index")
    rmsd_parser.add_argument("prediction")
    rmsd_parser.add_argument("prediction_index")
    rmsd_parser.add_argument("--pvalue-mode", default="-", choices=["+", "-"])

    inf_parser = subparsers.add_parser("inf", help="Calculate interaction network fidelity metrics.")
    inf_parser.add_argument("native")
    inf_parser.add_argument("native_index")
    inf_parser.add_argument("prediction")
    inf_parser.add_argument("prediction_index")
    inf_parser.add_argument("--mc-annotate")
    inf_parser.add_argument("--annotation-cache-dir")
    inf_parser.add_argument("--native-annotation")
    inf_parser.add_argument("--prediction-annotation")

    lddt_parser = subparsers.add_parser("lddt", help="Calculate all-atom lDDT.")
    _add_structure_pair_arguments(lddt_parser)
    lddt_parser.add_argument("--inclusion-radius", type=float, default=15.0)

    assess_parser = subparsers.add_parser(
        "assess",
        help="Calculate RMSD, P-value, INF and lDDT for a reference/prediction pair.",
    )
    _add_structure_pair_arguments(assess_parser)
    assess_parser.add_argument("--pvalue-mode", default="-", choices=["+", "-"])
    assess_parser.add_argument("--mc-annotate")
    assess_parser.add_argument("--annotation-cache-dir")
    assess_parser.add_argument("--native-annotation")
    assess_parser.add_argument("--prediction-annotation")
    assess_parser.add_argument("--inclusion-radius", type=float, default=15.0)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        if args.command == "normalize":
            path = normalize_structure(args.input, args.output)
            print(json.dumps({"output": str(path)}, indent=2))
            return 0

        if args.command == "extract":
            path = extract_structure(args.input, args.residues, args.output)
            print(json.dumps({"output": str(path)}, indent=2))
            return 0

        if args.command == "rmsd":
            result = calculate_rmsd(
                args.native,
                args.native_index,
                args.prediction,
                args.prediction_index,
                pvalue_mode=args.pvalue_mode,
            )
            print(json.dumps(asdict(result), indent=2))
            return 0

        if args.command == "lddt":
            result = calculate_lddt(
                args.native,
                args.native_index,
                args.prediction,
                args.prediction_index,
                inclusion_radius=args.inclusion_radius,
            )
            print(json.dumps(asdict(result), indent=2))
            return 0

        annotator = _build_annotator(args)
        if args.command == "inf":
            result = calculate_interaction_network_fidelity(
                args.native,
                args.native_index,
                args.prediction,
                args.prediction_index,
                annotator=annotator,
            )
            print(json.dumps(asdict(result), indent=2))
            return 0

        result = calculate_assessment(
            args.native,
            args.native_index,
            args.prediction,
            args.prediction_index,
            pvalue_mode=args.pvalue_mode,
            annotator=annotator,
            inclusion_radius=args.inclusion_radius,
        )
        print(json.dumps(asdict(result), indent=2))
        return 0
    except RNAAssessmentError as exc:
        parser.exit(1, f"{exc}\n")


def _add_structure_pair_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("native")
    parser.add_argument("prediction")
    parser.add_argument("--native-index")
    parser.add_argument("--prediction-index")


def _build_annotator(args: argparse.Namespace) -> MCAnnotateRunner:
    annotation_overrides = {}
    if getattr(args, "native_annotation", None):
        annotation_overrides[args.native] = args.native_annotation
    if getattr(args, "prediction_annotation", None):
        annotation_overrides[args.prediction] = args.prediction_annotation

    return MCAnnotateRunner(
        binary_path=getattr(args, "mc_annotate", None),
        cache_dir=getattr(args, "annotation_cache_dir", None),
        annotation_overrides=annotation_overrides or None,
    )


if __name__ == "__main__":
    raise SystemExit(main())
