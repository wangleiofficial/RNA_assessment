from __future__ import annotations

import pytest

from rna_assessment import calculate_interaction_network_fidelity, calculate_rmsd

from .conftest import DATA_DIR


def test_rmsd_self_comparison_is_zero() -> None:
    result = calculate_rmsd(
        DATA_DIR / "14_solution_0.pdb",
        DATA_DIR / "14_solution_0.index",
        DATA_DIR / "14_solution_0.pdb",
        DATA_DIR / "14_solution_0.index",
    )

    assert result.rmsd == pytest.approx(0.0, abs=1e-8)
    assert 0.0 <= result.pvalue <= 1.0


def test_inf_self_comparison_is_one() -> None:
    result = calculate_interaction_network_fidelity(
        DATA_DIR / "14_solution_0.pdb",
        DATA_DIR / "14_solution_0.index",
        DATA_DIR / "14_solution_0.pdb",
        DATA_DIR / "14_solution_0.index",
    )

    assert result.rmsd == pytest.approx(0.0, abs=1e-8)
    assert result.inf_all == pytest.approx(1.0, abs=1e-8)
    assert result.inf_wc == pytest.approx(1.0, abs=1e-8)
    assert result.inf_nwc == pytest.approx(1.0, abs=1e-8)
    assert result.inf_stack == pytest.approx(1.0, abs=1e-8)


def test_cross_structure_metrics_are_finite() -> None:
    result = calculate_interaction_network_fidelity(
        DATA_DIR / "14_solution_0.pdb",
        DATA_DIR / "14_solution_0.index",
        DATA_DIR / "14_ChenPostExp_2.pdb",
        DATA_DIR / "14_ChenPostExp_2.index",
    )

    assert result.rmsd == pytest.approx(7.751173243045826)
    assert result.deformation_index == pytest.approx(10.643784178530252)
    assert result.inf_all == pytest.approx(0.7282347248904991)
    assert result.inf_wc == pytest.approx(0.9375)
    assert result.inf_nwc == pytest.approx(0.25)
    assert result.inf_stack == pytest.approx(0.7082882469748285)
