"""Test module for testing the AASHTO 2024 shear."""

import math

import pytest

from structuralcodes.codes.aashto_2024 import _shear


@pytest.mark.parametrize(
    'sx, ag, expected',
    [
        (11.811, 0.62992, 12.93668),
        (11.811, 0.3937, 15.92183),
        (11.811, 0.82677, 11.18857),
    ],
)
def test_s_xe(sx, ag, expected):
    """Test the s_xe function."""
    assert math.isclose(_shear.s_xe(sx, ag), expected, rel_tol=0.005)


@pytest.mark.parametrize(
    'sx, ag',
    [
        (-0.5, 2),
        (2, -0.5),
        (-3, -4),
    ],
)
def test_s_xe_raises_errors(sx, ag):
    """Test s_xe errors."""
    with pytest.raises(ValueError):
        _shear.s_xe(sx, ag)


@pytest.mark.parametrize(
    'VkN, rho_l, bw, dv, expected',
    [
        (273.9709, 0.01, 1000, 300, 0.001522),
        (380.0253, 0.01, 1000, 300, 0.002111),
        (424.0828, 0.01, 1000, 300, 0.002356),
    ],
)
def test_eps(VkN, rho_l, bw, dv, expected):
    """Test the eps function."""
    assert math.isclose(
        _shear.eps(VkN, rho_l, bw, dv), expected, rel_tol=0.005
    )


@pytest.mark.parametrize(
    'VkN, rho_l, bw, dv',
    [(-50, -0.015, 2, 5), (50, 0.01, -20, -8), (-30, -0.02, -15, -20)],
)
def test_eps_raises_errors(VkN, rho_l, bw, dv):
    """Test eps errors."""
    with pytest.raises(ValueError):
        _shear.eps(VkN, rho_l, bw, dv)


@pytest.mark.parametrize(
    's_xe, strain, expected',
    [
        (12.93668, 0.001522, 2.2009),
        (12.93668, 0.002111, 1.8245),
        (12.93668, 0.002356, 1.7034),
    ],
)
def test_beta(s_xe, strain, expected):
    """Test the beta function."""
    assert math.isclose(
        _shear.beta_wo_rein(s_xe, strain),
        expected,
        rel_tol=0.005,
    )


@pytest.mark.parametrize(
    's_xe, strain',
    [
        (11, 0.0016),
        (85, 0.0017),
    ],
)
def test_beta_raise_errors(s_xe, strain):
    """Test beta_wo_rein errors."""
    with pytest.raises(ValueError):
        _shear.beta_wo_rein(s_xe, strain)


@pytest.mark.parametrize(
    'beta, fc_prime, bw, d, expected',
    [
        (2.2009, 3.625, 11.811, 39.37, 61.57),
        (1.8245, 10.15, 11.811, 39.37, 85.41),
        (1.7034, 14.50, 11.811, 39.37, 95.31),
    ],
)
def test_Vc(beta, fc_prime, bw, d, expected):
    """Test the Vc function."""
    assert math.isclose(
        _shear.Vc(beta, fc_prime, bw, d), expected, rel_tol=0.005
    )


@pytest.mark.parametrize(
    'VkN, bw, dv, rho_l, s_xe, strain, beta, fc_prime, tau_MPa, expected',
    [
        (
            296.9441,
            1000,
            300,
            0.01,
            12.9367,
            0.0016497,
            2.1068,
            3.625,
            0.8742,
            0.9132,
        ),
        (
            399.1398,
            1000,
            300,
            0.01,
            12.9367,
            0.0022174,
            1.7699,
            10.15,
            1.2289,
            1.2668,
        ),
        (
            414.5062,
            1000,
            300,
            0.01,
            12.9367,
            0.0023028,
            1.7284,
            14.50,
            1.4343,
            1.4136,
        ),
    ],
)
def test_converge(
    VkN, bw, dv, rho_l, s_xe, strain, beta, fc_prime, tau_MPa, expected
):
    """Test the convergence function."""
    assert math.isclose(
        _shear._converge(
            VkN, bw, dv, rho_l, s_xe, strain, beta, fc_prime, tau_MPa
        ),
        expected,
        rel_tol=0.005,
    )


@pytest.mark.parametrize(
    'VkN, bw, dv, rho_l, s_xe, strain, beta, fc_prime, tau_MPa',
    [(-2, -15, -30, -0.01, -10, 0.0016, -1.7, -15.5, 1.14)],
)
def test_converge_errors(
    VkN, bw, dv, rho_l, s_xe, strain, beta, fc_prime, tau_MPa
):
    """Test converge errors."""
    with pytest.raises(ValueError):
        _shear._converge(
            VkN, bw, dv, rho_l, s_xe, strain, beta, fc_prime, tau_MPa
        )


@pytest.mark.parametrize(
    'strain, expected',
    [
        (0.003894, 42.63),
        (0.004106, 43.37),
        (0.004312, 44.09),
    ],
)
def test_theta(strain, expected):
    """Test the theta function."""
    assert math.isclose(_shear.theta(strain), expected, rel_tol=0.005)


@pytest.mark.parametrize(
    'strain, expected',
    [
        (0.003894, 1.22432),
        (0.004106, 1.17651),
        (0.004312, 1.13375),
    ],
)
def test_beta_with_reinforcement(strain, expected):
    """Test the beta_with_reinforcement."""
    assert math.isclose(
        _shear.beta_with_reinforcement(strain),
        expected,
        rel_tol=0.005,
    )


@pytest.mark.parametrize(
    'Av, fy, dv, bw, cot_theta, s, expected',
    [
        (1.2985, 58.00, 11.811, 39.37, 1.086, 8.858, 109.09),
        (1.2985, 65.25, 11.811, 39.37, 1.058, 8.858, 119.58),
        (1.2985, 72.50, 11.811, 39.37, 1.032, 8.858, 129.57),
    ],
)
def test_Vs(Av, fy, dv, bw, cot_theta, s, expected):
    """Test the Vs_s function."""
    assert math.isclose(
        _shear.Vs(Av, fy, dv, bw, cot_theta, s),
        expected,
        rel_tol=0.005,
    )


@pytest.mark.parametrize(
    'Av, fy, dv, bw, cot_theta, s',
    [
        (-0.5, -30, -5, -3, 1.75, -2),
        (-0.5, 40, -5, -3, 1.75, 21),
        (8, -30, 21, -3, 1.75, -2),
    ],
)
def test_Vs_errors(Av, fy, dv, bw, cot_theta, s):
    """Test Vs errors."""
    with pytest.raises(ValueError):
        _shear.Vs(Av, fy, dv, bw, cot_theta, s)


@pytest.mark.parametrize(
    'V_c, V_s, V_p, expected',
    [
        (48.44, 109.09, 0, 157.53),
        (46.55, 119.58, 0, 166.13),
        (44.86, 129.57, 0, 174.43),
    ],
)
def test_Vn(V_c, V_s, V_p, expected):
    """Test the Vn function."""
    assert math.isclose(
        _shear.Vn(V_c, V_s, V_p),
        expected,
        rel_tol=0.005,
    )
