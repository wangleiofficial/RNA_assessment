from __future__ import annotations

import argparse
import glob
import json
from dataclasses import asdict
from pathlib import Path

from .api import extract_structure, normalize_structure
from .benchmark import build_benchmark_jobs, describe_prepared_pair, load_benchmark_manifest, run_benchmark
from .exceptions import ManifestFormatError, RNAAssessmentError
from .mc_annotate import MCAnnotateRunner
from .metrics import (
    PreparedStructurePair,
    calculate_assessment,
    calculate_interaction_network_fidelity,
    calculate_lddt,
    calculate_rmsd,
    prepare_structure_pair,
)


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
    lddt_parser.add_argument("--per-residue", action="store_true")

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
    assess_parser.add_argument("--per-residue", action="store_true")

    map_parser = subparsers.add_parser(
        "map",
        help="Infer residue mapping and generated index specifications for a structure pair.",
    )
    _add_structure_pair_arguments(map_parser)

    benchmark_parser = subparsers.add_parser(
        "benchmark",
        help="Run batch assessment for one reference structure against many predictions.",
    )
    benchmark_parser.add_argument("native", nargs="?")
    benchmark_parser.add_argument("predictions", nargs="*")
    benchmark_parser.add_argument("--native-index")
    benchmark_parser.add_argument("--prediction-glob")
    benchmark_parser.add_argument("--manifest")
    benchmark_parser.add_argument("--pvalue-mode", default="-", choices=["+", "-"])
    benchmark_parser.add_argument("--mc-annotate")
    benchmark_parser.add_argument("--annotation-cache-dir")
    benchmark_parser.add_argument("--native-annotation")
    benchmark_parser.add_argument("--inclusion-radius", type=float, default=15.0)
    benchmark_parser.add_argument("--per-residue", action="store_true")
    benchmark_parser.add_argument(
        "--sort-by",
        choices=["input", "rmsd", "pvalue", "inf_all", "lddt"],
        default="input",
    )
    benchmark_parser.add_argument("--descending", action="store_true")

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
                include_per_residue=args.per_residue,
            )
            print(json.dumps(asdict(result), indent=2))
            return 0

        if args.command == "map":
            prepared = prepare_structure_pair(
                args.native,
                args.native_index,
                args.prediction,
                args.prediction_index,
            )
            entry = describe_prepared_pair(prepared, args.prediction, args.native)
            print(json.dumps(asdict(entry), indent=2))
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

        if args.command == "benchmark":
            jobs = _collect_benchmark_jobs(
                native=args.native,
                native_index=args.native_index,
                predictions=args.predictions,
                prediction_glob=args.prediction_glob,
                manifest=args.manifest,
                native_annotation=args.native_annotation,
            )
            if not jobs:
                parser.exit(2, "benchmark requires at least one prediction path, --prediction-glob, or --manifest.\n")
            result = run_benchmark(
                native_file=args.native,
                pvalue_mode=args.pvalue_mode,
                annotator=annotator,
                inclusion_radius=args.inclusion_radius,
                jobs=jobs,
                include_per_residue=args.per_residue,
            )
            print(json.dumps(asdict(_sort_benchmark_result(result, args.sort_by, args.descending)), indent=2))
            return 0

        result = calculate_assessment(
            args.native,
            args.native_index,
            args.prediction,
            args.prediction_index,
            pvalue_mode=args.pvalue_mode,
            annotator=annotator,
            inclusion_radius=args.inclusion_radius,
            include_per_residue=args.per_residue,
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


def _collect_prediction_paths(predictions: list[str], prediction_glob: str | None) -> list[str]:
    collected = list(predictions)
    if prediction_glob:
        collected.extend(sorted(glob.glob(prediction_glob)))
    deduplicated: list[str] = []
    seen: set[str] = set()
    for item in collected:
        if item in seen:
            continue
        seen.add(item)
        deduplicated.append(item)
    return deduplicated


def _collect_benchmark_jobs(
    native: str | None,
    native_index: str | None,
    predictions: list[str],
    prediction_glob: str | None,
    manifest: str | None,
    native_annotation: str | None,
):
    jobs = []
    if manifest:
        jobs.extend(load_benchmark_manifest(manifest))

    prediction_paths = _collect_prediction_paths(predictions, prediction_glob)
    if prediction_paths:
        if native is None:
            raise ManifestFormatError("benchmark path mode requires a reference structure.")
        jobs.extend(
            build_benchmark_jobs(
                native_file=native,
                predictions=prediction_paths,
                native_index=native_index,
            )
        )

    if native_annotation is not None:
        jobs = [
            job.__class__(
                prediction=job.prediction,
                native=job.native,
                label=job.label,
                native_index=job.native_index,
                prediction_index=job.prediction_index,
                native_annotation=native_annotation if job.native_annotation is None else job.native_annotation,
                prediction_annotation=job.prediction_annotation,
            )
            for job in jobs
        ]

    if native is not None:
        jobs = [
            job.__class__(
                prediction=job.prediction,
                native=native if job.native is None else job.native,
                label=job.label,
                native_index=native_index if job.native_index is None else job.native_index,
                prediction_index=job.prediction_index,
                native_annotation=job.native_annotation,
                prediction_annotation=job.prediction_annotation,
            )
            for job in jobs
        ]
    return jobs


def _sort_benchmark_result(result, sort_by: str, descending: bool):
    if sort_by == "input":
        return result

    succeeded = [entry for entry in result.entries if entry.metrics is not None]
    failed = [entry for entry in result.entries if entry.metrics is None]
    succeeded.sort(key=lambda entry: getattr(entry.metrics, sort_by), reverse=descending)
    return result.__class__(
        reference=result.reference,
        total_predictions=result.total_predictions,
        succeeded=result.succeeded,
        failed=result.failed,
        entries=tuple(succeeded + failed),
    )


if __name__ == "__main__":
    raise SystemExit(main())
