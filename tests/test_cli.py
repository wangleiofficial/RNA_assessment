from __future__ import annotations

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
