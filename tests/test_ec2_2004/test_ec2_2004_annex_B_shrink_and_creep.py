"""Tests for EUROCODE 2-1-2:2004 Annex B and chapter 3."""

import math

import pytest

from structuralcodes.codes.ec2_2004 import annex_b_shrink_and_creep


@pytest.mark.parametrize(
    'test_h_0, test_f_cm, test_RH, test_c_class, test_t0, test_t, expected',
    [
        (138.5, 43, 50, 'R', 7, 18263, 2.567),
        (136.5, 38, 55, 'N', 7, 18263, 3.083),
        (136.5, 28, 55, 'N', 7, 18263, 3.748),
    ],
)
def test_creep_returns_expected_values(
    test_h_0, test_f_cm, test_RH, test_c_class, test_t0, test_t, expected
):
    """Test that phi returns expected values."""
    creep = annex_b_shrink_and_creep.phi(
        test_h_0, test_f_cm, test_RH, test_c_class, test_t0, test_t
    )
    assert math.isclose(creep, expected, abs_tol=1e-3)


@pytest.mark.parametrize(
    'test_h_0, test_f_cm, test_c_class, test_RH, test_t_S, test_t, expected',
    [
        (138.5, 43, 'R ', 50, 28, 18263, 0.00065817),
        (138.5, 43, 'r', 50, 28, 18263, 0.00065817),
        (600.0, 43, 'N', 50, 28, 18263, 0.00038040),
        (600.0, 43, 'n', 50, 28, 18263, 0.00038040),
        (500.0, 43, 'N', 50, 28, 18263, 0.00038040),
        (136.8, 38, 'N ', 55, 28, 18263, 0.00048427),
        (100.0, 38, 'N', 55, 28, 18263, 0.00050943),
        (90.0, 38, 'N', 55, 28, 18263, 0.00050943),
        (136.8, 28, 'N', 55, 28, 18263, 0.00051442),
    ],
)
def test_shrinkage_returns_expected_values(
    test_h_0, test_f_cm, test_c_class, test_RH, test_t_S, test_t, expected
):
    """Test that eps_cs returns expected values."""
    shrinkage = annex_b_shrink_and_creep.eps_cs(
        test_h_0, test_f_cm, test_c_class, test_RH, test_t_S, test_t
    )
    assert math.isclose(shrinkage, expected, abs_tol=1e-6)


invalid_cement_classes = pytest.mark.parametrize(
    'cement_class',
    [
        'not a valid class',
        'rns',
    ],
)


@invalid_cement_classes
def test_creep_wrong_cement_class(cement_class):
    """Test that ValueError is raised if a wrong cement class is provided."""
    with pytest.raises(ValueError):
        annex_b_shrink_and_creep.phi(300, 38, 50, cement_class)


@invalid_cement_classes
def test_shrinkage_wrong_cement_class(cement_class):
    """Test that ValueError is raised if a wrong cement class is provided."""
    with pytest.raises(ValueError):
        annex_b_shrink_and_creep.eps_cd_0(cement_class, 38, 50)
