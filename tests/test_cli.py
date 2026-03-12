from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest

from .conftest import DATA_DIR, PROJECT_ROOT


def test_rmsd_cli_returns_json() -> None:
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "rna_assessment",
            "rmsd",
            str(DATA_DIR / "14_solution_0.pdb"),
            str(DATA_DIR / "14_solution_0.index"),
            str(DATA_DIR / "14_solution_0.pdb"),
            str(DATA_DIR / "14_solution_0.index"),
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    payload = json.loads(result.stdout)
    assert payload["rmsd"] == pytest.approx(0.0, abs=1e-8)


def test_normalize_cli_writes_file(tmp_path: Path) -> None:
    output_path = tmp_path / "normalized.pdb"
    subprocess.run(
        [
            sys.executable,
            "-m",
            "rna_assessment",
            "normalize",
            str(DATA_DIR / "14_solution_0.pdb"),
            str(output_path),
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    assert output_path.exists()


def test_assess_cli_returns_combined_metrics() -> None:
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "rna_assessment",
            "assess",
            str(DATA_DIR / "14_solution_0.pdb"),
            str(DATA_DIR / "14_ChenPostExp_2.pdb"),
            "--native-index",
            str(DATA_DIR / "14_solution_0.index"),
            "--prediction-index",
            str(DATA_DIR / "14_ChenPostExp_2.index"),
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    payload = json.loads(result.stdout)
    assert payload["rmsd"] == pytest.approx(7.751173243045826)
    assert payload["lddt"] == pytest.approx(0.6126129382795444)


def test_assess_cli_uses_sidecar_indices_when_flags_are_omitted() -> None:
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "rna_assessment",
            "assess",
            str(DATA_DIR / "14_solution_0.pdb"),
            str(DATA_DIR / "14_ChenPostExp_2.pdb"),
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    payload = json.loads(result.stdout)
    assert payload["rmsd"] == pytest.approx(7.751173243045826)
    assert payload["inf_all"] == pytest.approx(0.7282347248904991)


def test_lddt_cli_returns_json() -> None:
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "rna_assessment",
            "lddt",
            str(DATA_DIR / "14_solution_0.pdb"),
            str(DATA_DIR / "14_solution_0.pdb"),
            "--native-index",
            str(DATA_DIR / "14_solution_0.index"),
            "--prediction-index",
            str(DATA_DIR / "14_solution_0.index"),
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    payload = json.loads(result.stdout)
    assert payload["lddt"] == pytest.approx(1.0, abs=1e-8)


def test_assess_cli_can_emit_per_residue_report() -> None:
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "rna_assessment",
            "assess",
            str(DATA_DIR / "14_solution_0.pdb"),
            str(DATA_DIR / "14_ChenPostExp_2.pdb"),
            "--native-index",
            str(DATA_DIR / "14_solution_0.index"),
            "--prediction-index",
            str(DATA_DIR / "14_ChenPostExp_2.index"),
            "--per-residue",
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    payload = json.loads(result.stdout)
    assert len(payload["per_residue"]) == 60
    assert payload["per_residue"][0]["native_chain"] == "A"
    assert payload["per_residue"][0]["prediction_chain"] == "U"


def test_map_cli_reports_inferred_chain_mapping(tmp_path: Path) -> None:
    reference_path = tmp_path / "reference.pdb"
    prediction_path = tmp_path / "prediction.pdb"
    reference_path.write_text((DATA_DIR / "14_ChenPostExp_2.pdb").read_text(encoding="utf-8"), encoding="utf-8")
    prediction_path.write_text(
        _rewrite_chain_ids((DATA_DIR / "14_ChenPostExp_2.pdb").read_text(encoding="utf-8"), "Z"),
        encoding="utf-8",
    )

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "rna_assessment",
            "map",
            str(reference_path),
            str(prediction_path),
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    payload = json.loads(result.stdout)
    assert payload["matched_residues"] > 0
    assert payload["chain_mappings"][0]["native_chain"] == "U"
    assert payload["chain_mappings"][0]["prediction_chain"] == "Z"


def test_benchmark_cli_returns_multiple_results() -> None:
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "rna_assessment",
            "benchmark",
            str(DATA_DIR / "14_solution_0.pdb"),
            str(DATA_DIR / "14_solution_0.pdb"),
            str(DATA_DIR / "14_ChenPostExp_2.pdb"),
            "--native-index",
            str(DATA_DIR / "14_solution_0.index"),
            "--sort-by",
            "rmsd",
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    payload = json.loads(result.stdout)
    assert payload["total_predictions"] == 2
    assert payload["succeeded"] == 2
    assert payload["entries"][0]["metrics"]["rmsd"] == pytest.approx(0.0, abs=1e-8)


def test_benchmark_cli_supports_json_manifest_and_per_residue(tmp_path: Path) -> None:
    manifest_path = tmp_path / "benchmark.json"
    manifest_path.write_text(
        json.dumps(
            [
                {
                    "label": "self",
                    "native": str(DATA_DIR / "14_solution_0.pdb"),
                    "native_index": str(DATA_DIR / "14_solution_0.index"),
                    "prediction": str(DATA_DIR / "14_solution_0.pdb"),
                },
                {
                    "label": "chen",
                    "native": str(DATA_DIR / "14_solution_0.pdb"),
                    "native_index": str(DATA_DIR / "14_solution_0.index"),
                    "prediction": str(DATA_DIR / "14_ChenPostExp_2.pdb"),
                },
            ]
        ),
        encoding="utf-8",
    )

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "rna_assessment",
            "benchmark",
            "--manifest",
            str(manifest_path),
            "--per-residue",
            "--sort-by",
            "rmsd",
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    payload = json.loads(result.stdout)
    assert payload["total_predictions"] == 2
    assert payload["entries"][0]["label"] == "self"
    assert len(payload["entries"][1]["metrics"]["per_residue"]) == 60


def test_benchmark_cli_supports_csv_manifest(tmp_path: Path) -> None:
    manifest_path = tmp_path / "benchmark.csv"
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["label", "native", "native_index", "prediction", "prediction_index"],
        )
        writer.writeheader()
        writer.writerow(
            {
                "label": "self",
                "native": str(DATA_DIR / "14_solution_0.pdb"),
                "native_index": str(DATA_DIR / "14_solution_0.index"),
                "prediction": str(DATA_DIR / "14_solution_0.pdb"),
                "prediction_index": str(DATA_DIR / "14_solution_0.index"),
            }
        )

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "rna_assessment",
            "benchmark",
            "--manifest",
            str(manifest_path),
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    payload = json.loads(result.stdout)
    assert payload["total_predictions"] == 1
    assert payload["entries"][0]["label"] == "self"
    assert payload["entries"][0]["metrics"]["rmsd"] == pytest.approx(0.0, abs=1e-8)


def _rewrite_chain_ids(pdb_text: str, chain_id: str) -> str:
    rows = []
    for row in pdb_text.splitlines():
        if row.startswith(("ATOM  ", "HETATM")):
            rows.append(f"{row[:21]}{chain_id}{row[22:]}")
        else:
            rows.append(row)
    return "\n".join(rows) + "\n"
