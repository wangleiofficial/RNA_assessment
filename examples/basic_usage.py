from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

from rna_assessment import calculate_interaction_network_fidelity, calculate_rmsd, normalize_structure


def main() -> None:
    root = Path(__file__).resolve().parent
    data_dir = root / "data"
    output_dir = root / "output"
    output_dir.mkdir(exist_ok=True)

    normalize_structure(
        data_dir / "14_solution_0.pdb",
        output_dir / "14_solution_0.normalized.pdb",
    )

    rmsd_result = calculate_rmsd(
        data_dir / "14_solution_0.pdb",
        data_dir / "14_solution_0.index",
        data_dir / "14_ChenPostExp_2.pdb",
        data_dir / "14_ChenPostExp_2.index",
    )
    inf_result = calculate_interaction_network_fidelity(
        data_dir / "14_solution_0.pdb",
        data_dir / "14_solution_0.index",
        data_dir / "14_ChenPostExp_2.pdb",
        data_dir / "14_ChenPostExp_2.index",
    )

    print(
        json.dumps(
            {
                "normalized_file": str(output_dir / "14_solution_0.normalized.pdb"),
                "rmsd": asdict(rmsd_result),
                "inf": asdict(inf_result),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
