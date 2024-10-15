"""Tests for EUROCODE 2-1-2:2004 Annex B and chapter 3."""

import numpy as np
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
    assert np.isclose(creep, expected, atol=1e-3)


def test_creep_array_input():
    """Test calculating creep number with array input."""
    # Arrange
    h_0 = 352
    f_cm = 51
    RH = 67
    cement_class = 's'
    t0 = 8
    t_final = 23 * 365
    num_time = 13
    t = np.linspace(0, t_final, num_time)

    # Act
    creep = annex_b_shrink_and_creep.phi(
        h_0=h_0,
        f_cm=f_cm,
        RH=RH,
        cement_class=cement_class,
        t0=t0,
        t=t,
    )
    creep_final = annex_b_shrink_and_creep.phi(
        h_0=h_0,
        f_cm=f_cm,
        RH=RH,
        cement_class=cement_class,
        t0=t0,
        t=t[-1],
    )

    # Assert
    assert len(creep) == num_time
    assert np.isclose(creep[-1], creep_final)


@pytest.mark.parametrize(
    'test_h_0, test_f_cm, test_c_class, test_RH, test_t_S, test_t, expected',
    [
        (138.5, 43, 'R ', 50, 28, 18263, 0.0006560211272528516),
        (138.5, 43, 'r', 50, 28, 18263, 0.0006560211272528516),
        (600.0, 43, 'N', 50, 28, 18263, 0.0003704813729755042),
        (600.0, 43, 'n', 50, 28, 18263, 0.0003704813729755042),
        (500.0, 43, 'N', 50, 28, 18263, 0.00037280025633482355),
        (136.8, 38, 'N ', 55, 28, 18263, 0.0004825589169449677),
        (100.0, 38, 'N', 55, 28, 18263, 0.000508432488714812),
        (90.0, 38, 'N', 55, 28, 18263, 0.0005085792190408255),
        (136.8, 28, 'N', 55, 28, 18263, 0.0005127088169780931),
    ],
)
def test_shrinkage_returns_expected_values(
    test_h_0, test_f_cm, test_c_class, test_RH, test_t_S, test_t, expected
):
    """Test that eps_cs returns expected values."""
    shrinkage = annex_b_shrink_and_creep.eps_cs(
        test_h_0, test_f_cm, test_c_class, test_RH, test_t_S, test_t
    )
    assert np.isclose(shrinkage, expected, atol=1e-6)


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


def test_shrinkage_array_input():
    """Test calculating the shrinkage strain with an array as input."""
    # Arrange
    h_0 = 400
    f_cm = 53
    cement_class = 'n'
    RH = 70
    t_s = 14
    t_final = 35 * 365
    num_time = 10
    t = np.linspace(0, t_final, num_time)

    # Act
    shrinkage = annex_b_shrink_and_creep.eps_cs(
        h_0=h_0,
        f_cm=f_cm,
        cement_class=cement_class,
        RH=RH,
        t_S=t_s,
        t=t,
    )
    shrinkage_final = annex_b_shrink_and_creep.eps_cs(
        h_0=h_0,
        f_cm=f_cm,
        cement_class=cement_class,
        RH=RH,
        t_S=t_s,
        t=t[-1],
    )

    # Assert
    assert len(shrinkage) == num_time
    assert np.isclose(shrinkage[-1], shrinkage_final)
