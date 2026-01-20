"""Tests related to development of concrete property with time."""

import numpy as np
import pytest

from structuralcodes.codes.ec2_2004._concrete_material_properties import (
    beta_cc,
    beta_ct,
    beta_E,
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
