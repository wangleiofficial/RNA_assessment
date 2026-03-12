from __future__ import annotations

from pathlib import Path

from rna_kit import MolProbityRunner, calculate_molprobity

from .conftest import DATA_DIR, write_mmcif_from_pdb


def test_molprobity_runner_parses_summary_output(tmp_path: Path) -> None:
    structure_path = tmp_path / "model.pdb"
    structure_path.write_text((DATA_DIR / "14_solution_0.pdb").read_text(encoding="utf-8"), encoding="utf-8")
    binary_path = _write_fake_molprobity_script(tmp_path)

    result = calculate_molprobity(
        structure_path,
        runner=MolProbityRunner(binary_path=binary_path),
    )

    assert result.clashscore == 5.42
    assert result.molprobity_score == 2.11
    assert result.bond_outliers == 1.0
    assert result.angle_outliers == 3.0
    assert result.pucker_outliers == 2.0
    assert result.suite_outliers == 4.0


def test_molprobity_runner_converts_mmcif_input_to_pdb(tmp_path: Path) -> None:
    mmcif_path = write_mmcif_from_pdb(DATA_DIR / "14_solution_0.pdb", tmp_path / "model.cif")
    binary_path = _write_fake_molprobity_script(tmp_path)

    calculate_molprobity(
        mmcif_path,
        runner=MolProbityRunner(binary_path=binary_path),
    )

    seen_path = (tmp_path / "seen_path.txt").read_text(encoding="utf-8").strip()
    assert seen_path.endswith(".pdb")


def _write_fake_molprobity_script(tmp_path: Path) -> Path:
    script_path = tmp_path / "fake_molprobity.py"
    seen_path = tmp_path / "seen_path.txt"
    script_path.write_text(
        "\n".join(
            [
                "#!/usr/bin/env python3",
                "import sys",
                "from pathlib import Path",
                "input_path = Path(sys.argv[1])",
                f"Path({str(seen_path)!r}).write_text(str(input_path), encoding='utf-8')",
                "print('All-atom clashscore = 5.42')",
                "print('MolProbity score = 2.11')",
                "print('Bad bonds = 1')",
                "print('Bad angles = 3')",
                "print('Pucker outliers = 2')",
                "print('Suite outliers = 4')",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    script_path.chmod(0o755)
    return script_path
