from __future__ import annotations

import copy
import math
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from Bio.PDB import PDBIO, Superimposer

from .exceptions import SequenceMismatchError, ToolNotAvailableError
from .mc_annotate import MCAnnotateRunner
from .structures import PDBStructure


@dataclass(frozen=True)
class RMSDResult:
    rmsd: float
    pvalue: float


@dataclass(frozen=True)
class InteractionNetworkResult:
    rmsd: float
    deformation_index: float
    inf_all: float
    inf_wc: float
    inf_nwc: float
    inf_stack: float


def erf(z: float) -> float:
    t = 1.0 / (1.0 + 0.5 * abs(z))
    ans = 1 - t * math.exp(
        -z * z
        - 1.26551223
        + t
        * (
            1.00002368
            + t
            * (
                0.37409196
                + t
                * (
                    0.09678418
                    + t
                    * (
                        -0.18628806
                        + t
                        * (
                            0.27886807
                            + t
                            * (
                                -1.13520398
                                + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277))
                            )
                        )
                    )
                )
            )
        )
    )
    return ans if z >= 0.0 else -ans


class PDBComparer:
    BACKBONE_ATOMS = ["C1'", "C2'", "C3'", "C4'", "C5'", "O2'", "O3'", "O4'", "O5'", "OP1", "OP2", "P"]
    HEAVY_ATOMS = ["C2", "C4", "C5", "C6", "C8", "N1", "N2", "N3", "N4", "N6", "N7", "N9", "O2", "O4", "O6"]
    ALL_ATOMS = BACKBONE_ATOMS + HEAVY_ATOMS

    def rmsd(self, src_struct: PDBStructure, trg_struct: PDBStructure, fit_pdb: str | Path | None = None) -> float:
        src_atoms, trg_atoms = self._get_atoms_struct(self.ALL_ATOMS, src_struct.res_sequence(), trg_struct.res_sequence())
        superimposer = Superimposer()
        superimposer.set_atoms(src_atoms, trg_atoms)

        fit_structure = copy.deepcopy(trg_struct.struct)
        superimposer.apply(fit_structure.get_atoms())

        if fit_pdb is not None:
            io = PDBIO()
            io.set_structure(fit_structure)
            io.save(str(fit_pdb))

        return float(superimposer.rms)

    def pvalue(self, rmsd_value: float, residue_count: int, param: str) -> float:
        if param == "+":
            a, b = 5.1, 15.8
        elif param == "-":
            a, b = 6.4, 12.7
        else:
            raise ValueError(f"Wrong p-value parameter '{param}'. Expected '+' or '-'.")
        expected_rmsd = a * (residue_count**0.41) - b
        z_score = (rmsd_value - expected_rmsd) / 1.8
        return (1.0 + erf(z_score / (2**0.5))) / 2.0

    def inf(
        self,
        src_struct: PDBStructure,
        trg_struct: PDBStructure,
        interaction_type: str = "ALL",
        annotator: MCAnnotateRunner | None = None,
    ) -> float:
        annotator = annotator or MCAnnotateRunner()
        src_interactions = self._select_interactions(annotator.indexed_interactions(src_struct), interaction_type)
        trg_interactions = self._select_interactions(annotator.indexed_interactions(trg_struct), interaction_type)

        true_positives = len(set(src_interactions) & set(trg_interactions))
        false_negatives = len(set(src_interactions) - set(trg_interactions))
        false_positives = len(set(trg_interactions) - set(src_interactions))

        if true_positives == 0 and (false_positives == 0 or false_negatives == 0):
            return -1.0

        precision = true_positives / float(true_positives + false_positives)
        sensitivity = true_positives / float(true_positives + false_negatives)
        return math.sqrt(precision * sensitivity)

    def mcq(self, model_file: str | Path, target_file: str | Path, jar_path: str | Path | None = None) -> float:
        jar = Path(jar_path) if jar_path else _default_jar_path("mcq.ws.client-0.0.1-SNAPSHOT-jar-with-dependencies.jar")
        return self._run_java_metric(
            [
                "-cp",
                str(jar),
                "pl.poznan.put.mcq.ws.client.Global",
                "-m",
                str(model_file),
                "-t",
                str(target_file),
            ],
            parser=lambda stdout: float(stdout.strip()),
            tool_name="MCQ",
        )

    def gdt(self, reference_file: str | Path, model_file: str | Path, jar_path: str | Path | None = None) -> float:
        jar = Path(jar_path) if jar_path else _default_jar_path("gdt.jar")
        return self._run_java_metric(
            ["-jar", str(jar), str(model_file), str(reference_file)],
            parser=_parse_gdt_output,
            tool_name="GDT",
        )

    def _run_java_metric(self, args: list[str], parser, tool_name: str) -> float:
        if shutil.which("java") is None:
            raise ToolNotAvailableError(f"{tool_name} requires a Java runtime.")
        try:
            result = subprocess.run(["java", *args], check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as exc:
            raise ToolNotAvailableError(
                f"{tool_name} failed: {exc.stderr.strip() or exc.stdout.strip()}"
            ) from exc
        return parser(result.stdout)

    def _get_atoms_residue(self, atom_list: list[str], src_residue, trg_residue) -> tuple[list[object], list[object]]:
        src_atoms: list[object] = []
        trg_atoms: list[object] = []
        src_candidates = [atom for atom in src_residue if atom.get_name() in atom_list]
        trg_candidates = [atom for atom in trg_residue if atom.get_name() in atom_list]

        for src_atom in src_candidates:
            src_name = src_atom.get_full_id()[4][0]
            for trg_atom in trg_candidates:
                trg_name = trg_atom.get_full_id()[4][0]
                if src_name == trg_name:
                    src_atoms.append(src_atom)
                    trg_atoms.append(trg_atom)
                    break
        return src_atoms, trg_atoms

    def _get_atoms_struct(
        self, atom_list: list[str], src_residues: list[object], trg_residues: list[object]
    ) -> tuple[list[object], list[object]]:
        if len(src_residues) != len(trg_residues):
            raise ValueError("Different number of residues.")
        src_atoms: list[object] = []
        trg_atoms: list[object] = []
        for src_residue, trg_residue in zip(src_residues, trg_residues):
            src_chunk, trg_chunk = self._get_atoms_residue(atom_list, src_residue, trg_residue)
            src_atoms.extend(src_chunk)
            trg_atoms.extend(trg_chunk)
        return src_atoms, trg_atoms

    def _select_interactions(
        self, interactions: list[tuple[str, int, int, str]], interaction_type: str
    ) -> list[tuple[str, int, int, str]]:
        if interaction_type == "ALL":
            return interactions
        if interaction_type == "PAIR":
            return [item for item in interactions if item[0] in {"PAIR_2D", "PAIR_3D"}]
        if interaction_type in {"PAIR_2D", "PAIR_3D", "STACK"}:
            return [item for item in interactions if item[0] == interaction_type]
        raise ValueError(
            f"Wrong interaction type '{interaction_type}'. Expected 'ALL', 'PAIR', 'PAIR_2D', 'PAIR_3D' or 'STACK'."
        )


def calculate_rmsd(
    native_file: str | Path,
    native_index: str | Path | None,
    prediction_file: str | Path,
    prediction_index: str | Path | None,
    pvalue_mode: str = "-",
) -> RMSDResult:
    native, prediction = _load_pair(native_file, native_index, prediction_file, prediction_index)
    comparer = PDBComparer()
    rmsd_value = comparer.rmsd(prediction, native)
    return RMSDResult(
        rmsd=rmsd_value,
        pvalue=comparer.pvalue(rmsd_value, len(prediction.raw_sequence()), pvalue_mode),
    )


def calculate_interaction_network_fidelity(
    native_file: str | Path,
    native_index: str | Path | None,
    prediction_file: str | Path,
    prediction_index: str | Path | None,
    annotator: MCAnnotateRunner | None = None,
) -> InteractionNetworkResult:
    native, prediction = _load_pair(native_file, native_index, prediction_file, prediction_index)
    comparer = PDBComparer()
    rmsd_value = comparer.rmsd(prediction, native)
    inf_all = comparer.inf(prediction, native, interaction_type="ALL", annotator=annotator)
    return InteractionNetworkResult(
        rmsd=rmsd_value,
        deformation_index=rmsd_value / inf_all,
        inf_all=inf_all,
        inf_wc=comparer.inf(prediction, native, interaction_type="PAIR_2D", annotator=annotator),
        inf_nwc=comparer.inf(prediction, native, interaction_type="PAIR_3D", annotator=annotator),
        inf_stack=comparer.inf(prediction, native, interaction_type="STACK", annotator=annotator),
    )


def _load_pair(
    native_file: str | Path,
    native_index: str | Path | None,
    prediction_file: str | Path,
    prediction_index: str | Path | None,
) -> tuple[PDBStructure, PDBStructure]:
    native = PDBStructure.from_file(native_file, index_name=native_index)
    prediction = PDBStructure.from_file(prediction_file, index_name=prediction_index)
    if prediction.raw_sequence() != native.raw_sequence():
        raise SequenceMismatchError(
            "Result sequence does not match the reference sequence: "
            f"reference='{native.raw_sequence()}', prediction='{prediction.raw_sequence()}'."
        )
    return native, prediction


def _default_jar_path(filename: str) -> Path:
    return Path(__file__).resolve().parents[2] / "third_party" / "lib" / filename


def _parse_gdt_output(stdout: str) -> float:
    lines = [line.strip() for line in stdout.splitlines() if line.strip()]
    if len(lines) < 2:
        raise ToolNotAvailableError("GDT output was shorter than expected.")
    value = lines[1].split(",")[-1]
    if value == "NaN":
        return 0.0
    return float(value)
