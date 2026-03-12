"""Tests related to development of concrete property with time."""

import typing as t

import numpy as np
import pytest

from structuralcodes.codes.ec2_2004._concrete_material_properties import (
    Ecm,
    Ecm_time,
    beta_cc,
    beta_ct,
    beta_E,
    fcm_time,
    fctm,
    fctm_time,
    s_time_development,
)


@pytest.mark.parametrize('s', (0.2, 0.25, 0.38))
def test_unit_scale_at_28_days(s: float):
    """Test that all time development functions are unity at 28 days."""
    for fun in beta_cc, beta_ct, beta_E:
        assert np.isclose(fun(28, s), 1.0)


@pytest.mark.parametrize('s', (0.2, 0.25, 0.38))
def test_array_input(s: float):
    """Test that all time development functions can be called with ArrayLike
    input.
    """
    time = np.linspace(1e-3, 1000, 50)
    for fun in beta_cc, beta_ct, beta_E:
        beta = fun(time, s)
        assert not np.isscalar(beta)
        assert len(beta) == len(time)


@pytest.mark.parametrize('s', (0.2, 0.25, 0.38))
@pytest.mark.parametrize('t', (11, 28, 100))
def test_relation_tension_compression(t: float, s: float):
    """Test that the relation between development of compressive and tensile
    strength is as intended.
    """
    if t < 28:
        assert np.isclose(beta_ct(t, s), beta_cc(t, s))
    else:
        assert np.isclose(beta_ct(t, s), np.pow(beta_cc(t, s), 2 / 3))


@pytest.mark.parametrize('s', (0.2, 0.25, 0.38))
@pytest.mark.parametrize('t', (11, 28, 100))
def test_relation_compression_stiffness(t: float, s: float):
    """Test that the relation between development of compressive strength and
    Young's modulus is as intended.
    """
    assert np.isclose(beta_E(t, s), np.pow(beta_cc(t, s), 0.3))


@pytest.mark.parametrize(
    't, expected',
    (
        (9, 0.8261668305331398),
        (27, 0.995422968073841),
        (56, 1.075970779425791),
        (365, 1.198124627036262),
    ),
)
def test_compressive_strength_development(t, expected):
    """Test time development function at specific stages."""
    assert np.isclose(beta_cc(t, 0.25), expected)


@pytest.mark.parametrize(
    'cement_class, expected',
    (
        ('r', 0.20),
        ('N', 0.25),
        ('S', 0.38),
        (None, None),
    ),
)
def test_s_time_development(
    cement_class: t.Optional[str], expected: t.Optional[float]
):
    """Test the function for getting s based on the cement class."""
    if cement_class is not None:
        assert np.isclose(s_time_development(cement_class), expected)
    else:
        with pytest.raises(ValueError):
            s_time_development(cement_class)


@pytest.mark.parametrize('s', (0.2, 0.25, 0.38))
@pytest.mark.parametrize('fcm', (20, 35, 65))
def test_ecm_time(fcm: float, s: float):
    """Test the function for development of Young's modulus with time."""
    # Arrange
    t_max = 367
    E_reference = Ecm(fcm=fcm)
    _beta_E = beta_E(t_max, s)
    _beta_cc = beta_cc(t_max, s)
    fcm_t = _beta_cc * fcm
    Ecm_t_beta = E_reference * _beta_E

    # Act
    Ecm_t = Ecm_time(fcm, fcm_t, E_reference)

    # Assert
    assert np.isclose(Ecm_t, Ecm_t_beta)


@pytest.mark.parametrize('s', (0.2, 0.25, 0.38))
@pytest.mark.parametrize('fcm', (20, 35, 65))
def test_fcm_time(fcm: float, s: float):
    """Test the function for development of compressive strength with time."""
    # Arrange
    t_max = 367
    _beta_cc = beta_cc(t_max, s)
    fcm_t_beta = _beta_cc * fcm

    # Act
    fcm_t = fcm_time(fcm, _beta_cc)

    # Assert
    assert np.isclose(fcm_t, fcm_t_beta)


@pytest.mark.parametrize('t', (26, 35, 367))
@pytest.mark.parametrize('s', (0.2, 0.25, 0.38))
@pytest.mark.parametrize('fcm', (20, 35, 65))
def test_fctm_time(fcm: float, s: float, t: float):
    """Test the function for development of tensile strength with time."""
    # Arrange
    fctm_reference = fctm(fck=fcm - 8)
    _beta_cc = beta_cc(t, s)
    _beta_ct = beta_ct(t, s)
    alpha = 1 if t < 28 else 2 / 3
    fctm_t_beta = _beta_ct * fctm_reference

    # Act
    fctm_t = fctm_time(fctm_reference, _beta_cc, alpha)

    # Assert
    assert np.isclose(fctm_t, fctm_t_beta)
