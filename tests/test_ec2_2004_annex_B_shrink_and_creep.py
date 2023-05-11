""" Tests for EUROCODE 2-1-2:2004 Annex B and chapter 3"""

import pytest
from structuralcodes.codes.ec2_2004 import annex_b_shrink_and_creep


@pytest.mark.parametrize(
    'test_h_0, test_fck, test_RH, test_cement_class, test_t_dead_load, test_t_live_load, expected',
    [
        (138.5, 35, 50, 'R', 7, 28, [2.567, 2.130]),
        (135.5, 30, 55, 'N', 7, 28, [3.086, 2.375]),
    ],
)
def test_creep_returns_expected_values(
    test_h_0, test_fck, test_RH, test_cement_class, test_t_dead_load, test_t_live_load, expected
):
    """Test that annex_B_Creep returns expected values"""
    creep = annex_b_shrink_and_creep.annex_B_creep(
        test_h_0, test_fck, test_RH, test_cement_class, test_t_dead_load, test_t_live_load
    )
    assert pytest.approx(creep, 0.001) == expected


@pytest.mark.parametrize(
    'test_h_0, test_fck, test_RH, test_cement_class, test_t_S, expected',
    [
        (138.5, 35, 50, 'R', 28, 0.00065817),
        (135.5, 30, 55, 'N', 28, 0.00048494),
    ],
)
def test_shrinkage_returns_expected_values(
    test_h_0, test_fck, test_RH, test_cement_class, test_t_S, expected
):
    """Test that annex_B_Shrinkage returns expected values"""
    shrinkage = annex_b_shrink_and_creep.annex_B_shrinkage(
        test_h_0, test_fck, test_RH, test_cement_class, test_t_S
    )
    assert pytest.approx(shrinkage, 0.001) == expected
