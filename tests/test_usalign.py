from __future__ import annotations

from pathlib import Path

from rna_kit import (
    USAlignResult,
    USAlignRunner,
    calculate_lddt,
    calculate_us_align,
    write_us_align_html,
)

from .conftest import DATA_DIR, write_mmcif_from_pdb

def test_usalign_runner_parses_metrics_and_extracts_superposed_model(tmp_path: Path) -> None:
    reference = tmp_path / "reference.pdb"
    prediction = tmp_path / "prediction.pdb"
    reference.write_text(_simple_pdb("U"), encoding="utf-8")
    prediction.write_text(_simple_pdb("A"), encoding="utf-8")
    fake_binary = _write_fake_usalign_script(tmp_path)
    output_dir = tmp_path / "alignment"

    result = calculate_us_align(
        reference,
        prediction,
        runner=USAlignRunner(binary_path=fake_binary),
        output_dir=output_dir,
    )

    assert result.aligned_length == 36
    assert result.rmsd == 3.72
    assert result.sequence_identity == 0.889
    assert result.tm_score_prediction == 0.29691
    assert result.tm_score_reference == 0.41234
    assert result.superposed_prediction_output is not None
    superposed_content = Path(result.superposed_prediction_output).read_text(encoding="utf-8")
    assert " A   1" in superposed_content
    assert " B   1" not in superposed_content
    assert Path(result.superposition_output).exists()


def test_write_us_align_html_embeds_viewer_and_metrics(tmp_path: Path) -> None:
    reference = tmp_path / "reference.pdb"
    superposed_prediction = tmp_path / "superposed_prediction.pdb"
    reference.write_text(_simple_pdb("U"), encoding="utf-8")
    superposed_prediction.write_text(_simple_pdb("A"), encoding="utf-8")
    output = tmp_path / "viewer.html"
    result = USAlignResult(
        reference=str(reference),
        prediction="prediction.pdb",
        reference_name="reference.pdb:U",
        prediction_name="prediction.pdb:A",
        reference_length=61,
        prediction_length=61,
        aligned_length=36,
        rmsd=3.72,
        sequence_identity=0.889,
        tm_score_reference=0.41234,
        tm_score_prediction=0.29691,
        alignment_prediction="AACGUU",
        alignment_markup="::::..",
        alignment_reference="AACGAA",
        binary_path="USalign",
        output_directory=str(tmp_path),
        reference_structure_output=str(reference),
        prediction_structure_output="prediction_input.pdb",
        superposition_output="sup.pdb",
        superposed_prediction_output=str(superposed_prediction),
    )

    html_path = write_us_align_html(
        reference,
        superposed_prediction,
        output,
        result=result,
        title="RNA Alignment Viewer",
    )

    content = html_path.read_text(encoding="utf-8")
    assert "3Dmol-min.js" in content
    assert "RNA Alignment Viewer" in content
    assert "TM-score (reference)" in content
    assert "AACGUU" in content
    assert ".viewer-panel" in content
    assert "box-sizing: border-box;" in content


def test_write_us_align_html_can_color_prediction_by_lddt(tmp_path: Path) -> None:
    reference = DATA_DIR / "14_solution_0.pdb"
    prediction = DATA_DIR / "14_ChenPostExp_2.pdb"
    output_dir = tmp_path / "alignment"
    runner = USAlignRunner(binary_path=_write_fake_usalign_script(tmp_path))
    result = calculate_us_align(reference, prediction, runner=runner, output_dir=output_dir)
    lddt_result = calculate_lddt(
        reference,
        DATA_DIR / "14_solution_0.index",
        prediction,
        DATA_DIR / "14_ChenPostExp_2.index",
        include_per_residue=True,
    )
    html_path = write_us_align_html(
        result.reference_structure_output or reference,
        result.superposed_prediction_output or prediction,
        tmp_path / "viewer_colored.html",
        result=result,
        per_residue=lddt_result.per_residue,
    )

    content = html_path.read_text(encoding="utf-8")
    assert "Prediction lDDT" in content
    assert "const residueScores = [" in content
    assert "viewer.setStyle(selection" in content


def test_usalign_runner_converts_mmcif_inputs_to_pdb(tmp_path: Path) -> None:
    reference = write_mmcif_from_pdb(DATA_DIR / "14_solution_0.pdb", tmp_path / "reference.cif")
    prediction = write_mmcif_from_pdb(DATA_DIR / "14_ChenPostExp_2.pdb", tmp_path / "prediction.cif")
    script_path = _write_fake_usalign_script(tmp_path)
    output_dir = tmp_path / "alignment"

    calculate_us_align(
        reference,
        prediction,
        runner=USAlignRunner(binary_path=script_path),
        output_dir=output_dir,
    )

    argv_log = (tmp_path / "usalign_seen_args.txt").read_text(encoding="utf-8")
    first_line, second_line = argv_log.splitlines()[:2]
    assert first_line.endswith(".pdb")
    assert second_line.endswith(".pdb")


def _write_fake_usalign_script(tmp_path: Path) -> Path:
    script_path = tmp_path / "fake_usalign.py"
    argv_log = tmp_path / "usalign_seen_args.txt"
    script_path.write_text(
        "\n".join(
            [
                "#!/usr/bin/env python3",
                "import sys",
                "from pathlib import Path",
                "args = sys.argv[1:]",
                "prediction = Path(args[0])",
                "reference = Path(args[1])",
                f"Path({str(argv_log)!r}).write_text(str(prediction) + '\\n' + str(reference) + '\\n', encoding='utf-8')",
                "prefix = None",
                "for index, arg in enumerate(args):",
                "    if arg == '-o':",
                "        prefix = Path(args[index + 1])",
                "        break",
                "if prefix is None:",
                "    raise SystemExit(2)",
                "prefix.parent.mkdir(parents=True, exist_ok=True)",
                f"PDB_TEXT = {repr(_combined_superposition_pdb())}",
                "prefix.with_suffix('.pdb').write_text(PDB_TEXT, encoding='utf-8')",
                "print(' ********************************************************************')",
                "print(' * US-align (Version 20241108)                                      *')",
                "print(' ********************************************************************')",
                "print()",
                "print(f'Name of Structure_1: {prediction}:A (to be superimposed onto Structure_2)')",
                "print(f'Name of Structure_2: {reference}:U')",
                "print('Length of Structure_1: 61 residues')",
                "print('Length of Structure_2: 61 residues')",
                "print()",
                "print('Aligned length= 36, RMSD=   3.72, Seq_ID=n_identical/n_aligned= 0.889')",
                "print('TM-score= 0.29691 (normalized by length of Structure_1: L=61, d0=2.17)')",
                "print('TM-score= 0.41234 (normalized by length of Structure_2: L=61, d0=2.17)')",
                "print('(You should use TM-score normalized by length of the reference structure)')",
                "print()",
                "print('(\":\" denotes residue pairs of d < 5.0 Angstrom, \".\" denotes other aligned residues)')",
                "print('AACGUU')",
                "print('::::..')",
                "print('AACGAA')",
                "print()",
                "print('#Total CPU time is  0.01 seconds')",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    script_path.chmod(0o755)
    return script_path


def _simple_pdb(chain_id: str) -> str:
    return (
        f"ATOM      1  P     G {chain_id}   1      10.000  10.000  10.000  1.00 20.00           P  \n"
        f"ATOM      2  O5'   G {chain_id}   1      11.000  10.000  10.000  1.00 20.00           O  \n"
        "TER\n"
        "END\n"
    )


def _combined_superposition_pdb() -> str:
    return (
        "ATOM      1  P     G A   1      10.000  10.000  10.000  1.00 20.00           P  \n"
        "ATOM      2  O5'   G A   1      11.000  10.000  10.000  1.00 20.00           O  \n"
        "TER\n"
        "ATOM      3  P     G B   1      20.000  20.000  20.000  1.00 20.00           P  \n"
        "ATOM      4  O5'   G B   1      21.000  20.000  20.000  1.00 20.00           O  \n"
        "TER\n"
        "END\n"
    )
