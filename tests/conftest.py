from __future__ import annotations

from pathlib import Path

from Bio.PDB import MMCIFIO, PDBParser


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "examples" / "data"


def write_mmcif_from_pdb(source_pdb: str | Path, output_cif: str | Path) -> Path:
    source_path = Path(source_pdb)
    output_path = Path(output_cif)
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("input", str(source_path))
    io = MMCIFIO()
    io.set_structure(structure)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    io.save(str(output_path))
    return output_path
