"""Test module for testing the AASHTO 2020 shear."""

import math

import pytest

from structuralcodes.codes.aashto_2020 import _section_5


@pytest.mark.parametrize(
    'sx, ag, expected',
    [
        (11.811, 0.62992, 12.93668),
        (11.811, 0.3937, 15.92183),
        (11.811, 0.82677, 11.18857),
    ],
)
def test_calc_s_xe(sx, ag, expected):
    """Test the calc_s_xe function."""
    assert math.isclose(_section_5.calc_s_xe(sx, ag), expected, rel_tol=0.005)


@pytest.mark.parametrize(
    'VkN, rho_l, bw, dv, expected',
    [
        (273.9709, 0.01, 1000, 300, 0.001522),
        (380.0253, 0.01, 1000, 300, 0.002111),
        (424.0828, 0.01, 1000, 300, 0.002356),
    ],
)
def test_calc_strain(VkN, rho_l, bw, dv, expected):
    """Test the calc_strain function."""
    assert math.isclose(
        _section_5.calc_strain(VkN, rho_l, bw, dv), expected, rel_tol=0.005
    )


@pytest.mark.parametrize(
    's_xe, strain, expected',
    [
        (12.93668, 0.001522, 2.2009),
        (12.93668, 0.002111, 1.8245),
        (12.93668, 0.002356, 1.7034),
    ],
)
def test_calc_beta(s_xe, strain, expected):
    """Test the calc_beta function."""
    assert math.isclose(
        _section_5.calc_beta(s_xe, strain), expected, rel_tol=0.005
    )


@pytest.mark.parametrize(
    'beta, fc_prime, expected',
    [
        (2.2009, 3.625, 0.9132),
        (1.8245, 10.15, 1.267),
        (1.7034, 14.50, 1.414),
    ],
)
def test_calc_tau(beta, fc_prime, expected):
    """Test the calc_tau function."""
    assert math.isclose(
        _section_5.calc_tau(beta, fc_prime), expected, rel_tol=0.005
    )


@pytest.mark.parametrize(
    'strain, expected',
    [
        (0.003894, 42.63),
        (0.004106, 43.37),
        (0.004312, 44.09),
    ],
)
def test_calc_theta(strain, expected):
    """Test the calc_theta function."""
    assert math.isclose(_section_5.calc_theta(strain), expected, rel_tol=0.005)


@pytest.mark.parametrize(
    'strain, expected',
    [
        (0.003894, 1.22432),
        (0.004106, 1.17651),
        (0.004312, 1.13375),
    ],
)
def test_calc_beta_with_reinforcement(strain, expected):
    """Test the calc_beta_with_reinforcement."""
    assert math.isclose(
        _section_5.calc_beta_with_reinforcement(strain),
        expected,
        rel_tol=0.005,
    )


@pytest.mark.parametrize(
    'Av, fy, dv, bw, cot_theta, s, expected',
    [
        (1.2985, 58, 11.811, 39.37, 1.086, 8.858, 1.618),
        (1.2985, 65.25, 11.811, 39.37, 1.058, 8.858, 1.7735),
        (1.2985, 72.5, 11.811, 39.37, 1.032, 8.858, 1.9217),
    ],
)
def test_calc_tau_s(Av, fy, dv, bw, cot_theta, s, expected):
    """Test the calc_tau_s function."""
    assert math.isclose(
        _section_5.calc_tau_s(Av, fy, dv, bw, cot_theta, s),
        expected,
        rel_tol=0.005,
    )


@pytest.mark.parametrize(
    'tau, tau_s, tau_p, expected',
    [
        (0.7184, 1.6190, 0, 2.336),
        (0.6904, 1.7735, 0, 2.464),
        (0.6653, 1.9217, 0, 2.587),
    ],
)
def test_calc_tau_nominal(tau, tau_s, tau_p, expected):
    """Test the calc_tau_nominal function."""
    assert math.isclose(
        _section_5.calc_tau_nominal(tau, tau_s, tau_p),
        expected,
        rel_tol=0.005,
    )
