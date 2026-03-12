from __future__ import annotations

import pytest

from rna_kit import MCAnnotateRunner, ToolNotAvailableError

from .conftest import DATA_DIR, write_mmcif_from_pdb


def test_precomputed_annotation_is_used_without_binary() -> None:
    runner = MCAnnotateRunner(binary_path="/definitely/missing")
    result = runner.load(DATA_DIR / "14_solution_0.pdb")

    assert len(result.residues) > 0
    assert len(result.interactions) > 0


def test_missing_annotation_raises_without_available_binary(tmp_path) -> None:
    pdb_path = tmp_path / "model.pdb"
    pdb_path.write_text((DATA_DIR / "14_solution_0.pdb").read_text(encoding="utf-8"), encoding="utf-8")

    runner = MCAnnotateRunner(binary_path=tmp_path / "missing-binary")
    with pytest.raises(ToolNotAvailableError):
        runner.load(pdb_path)


def test_annotation_override_path_is_supported(tmp_path) -> None:
    pdb_path = tmp_path / "model.pdb"
    annotation_path = tmp_path / "custom_output.mcout"
    pdb_path.write_text((DATA_DIR / "14_solution_0.pdb").read_text(encoding="utf-8"), encoding="utf-8")
    annotation_path.write_text((DATA_DIR / "14_solution_0.pdb.mcout").read_text(encoding="utf-8"), encoding="utf-8")

    runner = MCAnnotateRunner(
        binary_path=tmp_path / "missing-binary",
        annotation_overrides={pdb_path: annotation_path},
    )
    result = runner.load(pdb_path)

    assert len(result.interactions) > 0


def test_mc_annotate_runner_converts_mmcif_input_to_pdb(tmp_path) -> None:
    structure_path = write_mmcif_from_pdb(DATA_DIR / "14_solution_0.pdb", tmp_path / "model.cif")
    seen_path = tmp_path / "seen_path.txt"
    binary_path = tmp_path / "fake_mc_annotate.py"
    binary_path.write_text(
        "\n".join(
            [
                "#!/usr/bin/env python3",
                "import sys",
                "from pathlib import Path",
                f"Path({str(seen_path)!r}).write_text(sys.argv[1], encoding='utf-8')",
                "print('Residue conformations -------------------------------------------')",
                "print('A1 x A anti anti')",
                "print('A2 x U anti anti')",
                "print('---------------------------------------------------------------')",
                "print('Base-pairs ----------------------------------------------------')",
                "print('A1-A2 : A-U W/W pairing antiparallel cis')",
                "print('---------------------------------------------------------------')",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    binary_path.chmod(0o755)

    runner = MCAnnotateRunner(binary_path=binary_path)
    result = runner.load(structure_path)

    assert len(result.interactions) == 1
    assert seen_path.read_text(encoding="utf-8").strip().endswith(".pdb")
