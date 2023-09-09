""" Tests for EUROCODE 2-1-2:2004 Annex B and chapter 3"""
import pytest
from structuralcodes.codes.ec2_2004 import annex_b_shrink_and_creep


@pytest.mark.parametrize(
    'test_h_0, test_f_cm, test_RH, test_cement_class, test_t0, test_t, expected',
    [
        (138.5, 43, 50, 'R', 7, 18263, 2.567),
        (136.5, 38, 55, 'N', 7, 18263, 3.083),
    ],
)
def test_creep_returns_expected_values(
    test_h_0, test_f_cm, test_RH, test_cement_class, test_t0, test_t, expected
):
    """Test that annex_B_Creep returns expected values"""
    creep = annex_b_shrink_and_creep.annex_B_creep(
        test_h_0, test_f_cm, test_RH, test_cement_class, test_t0, test_t
    )
    assert pytest.approx(creep, 0.001) == expected


@pytest.mark.parametrize(
    'test_h_0, test_f_cm, test_cement_class, test_RH, test_t_S, test_t, expected',
    [
        (138.5, 43, 'R', 50, 28, 18263, 0.00065817),
        (136.8, 38, 'N', 55, 28, 18263, 0.00048427),
    ],
)
def test_shrinkage_returns_expected_values(
    test_h_0, test_f_cm, test_cement_class, test_RH, test_t_S, test_t, expected
):
    """Test that annex_B_Shrinkage returns expected values"""
    shrinkage = annex_b_shrink_and_creep.annex_B_shrinkage(
        test_h_0, test_f_cm, test_cement_class, test_RH, test_t_S, test_t
    )
    assert pytest.approx(shrinkage, 0.001) == expected


if __name__ == '__main__':
    test_shrinkage_returns_expected_values()
    test_creep_returns_expected_values()
