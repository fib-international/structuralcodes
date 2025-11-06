"""Tests for the section EC2 2023 8.2 Shear."""

import math

import pytest

from structuralcodes.codes.ec2_2023 import _section_8_2_shear


@pytest.mark.parametrize(
    'VEd, bw, d, expected',
    [
        (100.0, 0.3, 0.5, 100.0 * 1000.0 / (0.3 * 0.9 * 0.5)),
        (200.0, 0.5, 0.6, 200.0 * 1000.0 / (0.5 * 0.9 * 0.6)),
    ],
)
def test_tao_Ed(VEd, bw, d, expected):
    """Test shear_stress_linear_members."""
    assert _section_8_2_shear.tau_Ed(VEd, bw, d) == pytest.approx(expected)


@pytest.mark.parametrize(
    'VEd, bw, d',
    [
        (100.0, -0.3, 0.5),
        (100.0, 0.3, -0.5),
        (100.0, -0.3, -0.5),
    ],
)
def test_tau_Ed_value_errors(VEd, bw, d):
    """Test tao_Ed raises ValueError for negative bw or d."""
    with pytest.raises(ValueError):
        _section_8_2_shear.tau_Ed(VEd, bw, d)


@pytest.mark.parametrize(
    'vEd, d, expected',
    [
        (100.0, 500.0, 100.0 / (0.9 * 500.0)),
        (200.0, 400.0, 200.0 / (0.9 * 400.0)),
        (50.0, 250.0, 50.0 / (0.9 * 250.0)),
        (0.0, 300.0, 0.0),  # Zero shear force
    ],
)
def test_tao_Ed_planar_valid(vEd, d, expected):
    """Test tao_Ed_planar with valid inputs."""
    assert _section_8_2_shear.tau_Ed_planar(vEd, d) == pytest.approx(
        expected, rel=1e-9
    )


@pytest.mark.parametrize(
    'vEd, d',
    [
        (100.0, -500.0),
        (0.0, -1.0),
        (50.0, -0.01),
    ],
)
def test_tao_Ed_planar_invalid_depth(vEd, d):
    """Test tao_Ed_planar raises ValueError for negative d."""
    with pytest.raises(ValueError):
        _section_8_2_shear.tau_Ed_planar(vEd, d)


@pytest.mark.parametrize(
    'f_ck, d_lower, expected',
    [
        (30, 20, 36),  # f_ck <= 60: min(16 + 20, 40) = 36
        (40, 25, 40),  # f_ck <= 60: min(16 + 25, 40) = 40
        (70, 30, 38.0204),  # f_ck > 60: min(16 + 30 * (60/70)^2, 40) = 38.0204
        (80, 20, 27.25),  # f_ck > 60: min(16 + 20 * (60/80)^2, 40) = 27.25
    ],
)
def test_d_dg(f_ck, d_lower, expected):
    """Test the d_dg function with example values."""
    result = _section_8_2_shear.d_dg(f_ck, d_lower)
    assert result == pytest.approx(expected, rel=1e-3)


@pytest.mark.parametrize(
    'gamma_v, f_ck, f_yd, d, d_lower, expected',
    [
        (1.4, 30, 500, 500, 20, 0.51642),
        (1.5, 40, 600, 600, 25, 0.48888),
        (1.3, 70, 700, 700, 30, 0.62377),
    ],
)
def test_tau_rdc_min(gamma_v, f_ck, f_yd, d, d_lower, expected):
    """Test the calculate_tau_rdc_min function with example values."""
    d_dg = _section_8_2_shear.d_dg(f_ck, d_lower)
    result = _section_8_2_shear.tau_rdc_min(gamma_v, f_ck, f_yd, d, d_dg)
    assert result == pytest.approx(expected, rel=1e-4)


@pytest.mark.parametrize(
    'vEd_x, vEd_y, expected',
    [
        (3.0, 4.0, 5.0),  # Pythagorean triplet
        (0.0, 4.0, 4.0),  # One direction zero
        (3.0, 0.0, 3.0),  # One direction zero
    ],
)
def test_v_Ed(vEd_x, vEd_y, expected):
    """Test calculation of design shear force per unit width."""
    assert math.isclose(
        _section_8_2_shear.v_Ed(vEd_x, vEd_y),
        expected,
        rel_tol=1e-9,
    )


@pytest.mark.parametrize(
    'dx, dy, vEd_x, vEd_y, expected',
    [
        (300, 400, 3.0, 1.0, 300),  # Ratio <= 0.5
        (300, 400, 3.0, 3.0, 350),  # Ratio between 0.5 and 2
        (300, 400, 1.0, 3.0, 400),  # Ratio >= 2
    ],
)
def test_d_eff(dx, dy, vEd_x, vEd_y, expected):
    """Test calculation of effective depth based on shear force ratio."""
    assert math.isclose(
        _section_8_2_shear.d_eff(dx, dy, vEd_x, vEd_y), expected, rel_tol=1e-9
    )


@pytest.mark.parametrize(
    'dx, dy, vEd_x, vEd_y, expected',
    [
        (300, 400, 3.0, 0.0, 300),  # alpha_v = 0
        (300, 400, 0.0, 3.0, 400),  # alpha_v = 90 degrees
        (300, 400, 3.0, 3.0, 350),  # alpha_v = 45 degrees
    ],
)
def test_ed_eff_with_angle(dx, dy, vEd_x, vEd_y, expected):
    """Test calculation of effective depth based on angle alpha_v."""
    assert math.isclose(
        _section_8_2_shear.d_eff_angle(dx, dy, vEd_x, vEd_y),
        expected,
        rel_tol=1e-9,
    )


@pytest.mark.parametrize(
    'dx, dy, vEd_x, vEd_y',
    [
        (-300, 400, 3.0, 1.0),  # Negative dx
        (300, -400, 3.0, 1.0),  # Negative dy
        (-300, -400, 3.0, 1.0),  # Both dx and dy negative
    ],
)
def test_d_eff_with_angle_value_errors(dx, dy, vEd_x, vEd_y):
    """Test d_eff_with_angle raises ValueError for negative dx or dy."""
    with pytest.raises(ValueError):
        _section_8_2_shear.d_eff_angle(dx, dy, vEd_x, vEd_y)


@pytest.mark.parametrize(
    'gamma_v, rho_l, f_ck, d, d_dg, tau_rdc_min, expected',
    [
        (1.5, 0.02, 30, 500, 16, 0.3, 0.5468),
        (1.4, 0.03, 40, 450, 20, 0.5, 0.82366),
        (1.6, 0.025, 35, 600, 18, 0.1, 0.5690),
    ],
)
def test_calculate_tau_Rdc(
    gamma_v, rho_l, f_ck, d, d_dg, tau_rdc_min, expected
):
    """Test the calculation of the shear stress resistance."""
    result = _section_8_2_shear.tau_Rdc(
        gamma_v, rho_l, f_ck, d, d_dg, tau_rdc_min
    )
    assert math.isclose(
        result,
        expected,
        rel_tol=1e-3,
    )


@pytest.mark.parametrize(
    'f_ck, d_lower',
    [
        (-30, 20),
        (30, -20),
    ],
)
def test_d_dg_value_errors(f_ck, d_lower):
    """Test d_dg raises ValueError for negative arguments."""
    with pytest.raises(ValueError):
        _section_8_2_shear.d_dg(f_ck, d_lower)


@pytest.mark.parametrize(
    'gamma_v, f_ck, f_yd, d, d_dg',
    [
        (-1.0, 30, 500, 500, 36),
        (1.4, -30, 500, 500, 36),
        (1.4, 30, -500, 500, 36),
        (1.4, 30, 500, -500, 36),
        (1.4, 30, 500, 500, -36),
    ],
)
def test_tau_rdc_min_value_errors(gamma_v, f_ck, f_yd, d, d_dg):
    """Test tau_rdc_min raises ValueError for negative arguments."""
    with pytest.raises(ValueError):
        _section_8_2_shear.tau_rdc_min(gamma_v, f_ck, f_yd, d, d_dg)


@pytest.mark.parametrize(
    'A_sl, b_w, d, expected',
    [
        (300, 200, 500, 0.003),
        (400, 250, 600, 0.00266667),
        (500, 300, 700, 0.00238095),
    ],
)
def test_rho_l(A_sl, b_w, d, expected):
    """Test the calculation of the reinforcement ratio."""
    assert math.isclose(
        _section_8_2_shear.rho_l(A_sl, b_w, d), expected, rel_tol=1e-5
    )


@pytest.mark.parametrize(
    'A_sl, b_w, d',
    [
        (-1, 200, 500),  # Negative A_sl
        (300, -200, 500),  # Negative b_w
        (300, 200, -500),  # Negative d
        (0, 200, 500),  # Zero A_sl
        (300, 0, 500),  # Zero b_w
        (300, 200, 0),  # Zero d
    ],
)
def test_rho_l_value_errors(A_sl, b_w, d):
    """Test that rho_l raises ValueError for non-positive arguments."""
    with pytest.raises(ValueError):
        _section_8_2_shear.rho_l(A_sl, b_w, d)


@pytest.mark.parametrize(
    'gamma_v, rho_l, f_ck, d, d_dg, tau_rdc_min',
    [
        (-1.0, 0.02, 30, 500, 16, 0.3),  # Negative gamma_v
        (1.5, -0.02, 30, 500, 16, 0.3),  # Negative rho_l
        (1.5, 0.02, -30, 500, 16, 0.3),  # Negative f_ck
        (1.5, 0.02, 30, -500, 16, 0.3),  # Negative d
        (1.5, 0.02, 30, 500, -16, 0.3),  # Negative d_dg
        (0.0, 0.02, 30, 500, 16, 0.3),  # Zero gamma_v
        (1.5, 0.02, 30, 0.0, 16, 0.3),  # Zero d
    ],
)
def test_tau_Rdc_value_errors(gamma_v, rho_l, f_ck, d, d_dg, tau_rdc_min):
    """Test tau_Rdc raises ValueError for negative values or zero gamma_v/d."""
    with pytest.raises(ValueError):
        _section_8_2_shear.tau_Rdc(gamma_v, rho_l, f_ck, d, d_dg, tau_rdc_min)


@pytest.mark.parametrize(
    'gamma_v, rho_l, f_ck, d, d_dg, tau_rdc_min, expected',
    [
        (
            1.5,
            0.0,
            30,
            500,
            16,
            0.3,
            0.3,
        ),  # Zero rho_l - should return tau_rdc_min
        (
            1.5,
            0.02,
            0.0,
            500,
            16,
            0.3,
            0.3,
        ),  # Zero f_ck - should return tau_rdc_min
        (
            1.5,
            0.02,
            30,
            500,
            0.0,
            0.3,
            0.3,
        ),  # Zero d_dg - should return tau_rdc_min
    ],
)
def test_tau_Rdc_zero_values(
    gamma_v, rho_l, f_ck, d, d_dg, tau_rdc_min, expected
):
    """Test tau_Rdc handles zero values correctly (should not raise errors)."""
    result = _section_8_2_shear.tau_Rdc(
        gamma_v, rho_l, f_ck, d, d_dg, tau_rdc_min
    )
    assert result == pytest.approx(expected, rel=1e-3)


@pytest.mark.parametrize(
    'dx, dy, vEd_x, vEd_y',
    [
        (-300, 400, 3.0, 1.0),  # Negative dx
        (300, -400, 3.0, 1.0),  # Negative dy
        (-300, -400, 3.0, 1.0),  # Both dx and dy negative
    ],
)
def test_d_eff_value_errors(dx, dy, vEd_x, vEd_y):
    """Test d_eff raises ValueError for negative dx or dy."""
    with pytest.raises(ValueError):
        _section_8_2_shear.d_eff(dx, dy, vEd_x, vEd_y)


@pytest.mark.parametrize(
    'a_cs, d, expected',
    [
        (300, 200, 122.4744),
        (400, 250, 158.1138),
        (500, 300, 193.6491),
    ],
)
def test_a_v(a_cs, d, expected):
    """Test the a_v."""
    result = _section_8_2_shear.a_v(a_cs, d)
    assert math.isclose(result, expected, rel_tol=1e-5)


@pytest.mark.parametrize(
    'a_cs, d',
    [
        (-300, 200),  # Negative a_cs
        (300, -200),  # Negative d
        (0, 200),  # Zero a_cs
        (300, 0),  # Zero d
        (-300, -200),  # Both negative
        (0, 0),  # Both zero
    ],
)
def test_a_v_value_errors(a_cs, d):
    """Test that a_v raises ValueError for non-positive arguments."""
    with pytest.raises(ValueError):
        _section_8_2_shear.a_v(a_cs, d)


@pytest.mark.parametrize(
    'M_Ed, V_Ed, d, expected',
    [
        (1000, 100, 500, 10000),
        (2000, 200, 600, 10000),
        (1500, 150, 700, 10000),
    ],
)
def test_calculate_a_cs(M_Ed, V_Ed, d, expected):
    """Test the calculation of the effective shear span."""
    result = _section_8_2_shear.a_cs(M_Ed, V_Ed, d)
    assert math.isclose(result, expected, rel_tol=1e-5)


@pytest.mark.parametrize(
    'M_Ed, V_Ed, d',
    [
        (1000, 100, -500),  # Negative d
        (1000, 100, 0),  # Zero d
    ],
)
def test_a_cs_value_error_d(M_Ed, V_Ed, d):
    """Test a_cs raises ValueError for non-positive d."""
    with pytest.raises(ValueError):
        _section_8_2_shear.a_cs(M_Ed, V_Ed, d)


@pytest.mark.parametrize(
    'N_Ed, V_Ed, d, a_cs, expected',
    [
        (100, 50, 500, 1500, 1.2222),
        (200, 100, 600, 1600, 1.25),
        (-150, 75, 700, 1700, 0.7254),
    ],
)
def test_calculate_k_vp(N_Ed, V_Ed, d, a_cs, expected):
    """Test the calculation of the coefficient k_vp."""
    result = _section_8_2_shear.k_vp(N_Ed, V_Ed, d, a_cs)
    assert math.isclose(result, expected, rel_tol=1e-3)


@pytest.mark.parametrize(
    'N_Ed, V_Ed, d, a_cs',
    [
        (100, 50, -500, 1500),  # Negative d
        (100, 50, 500, -1500),  # Negative a_cs
        (100, 50, 0, 1500),  # Zero d
        (100, 50, 500, 0),  # Zero a_cs
        (100, 50, 0, 0),  # Both zero
        (100, 50, -500, -1500),  # Both negative
    ],
)
def test_k_vp_value_errors(N_Ed, V_Ed, d, a_cs):
    """Test k_vp raises ValueError for non-positive d or a_cs."""
    with pytest.raises(ValueError):
        _section_8_2_shear.k_vp(N_Ed, V_Ed, d, a_cs)


@pytest.mark.parametrize(
    'gamma_v, rho_l, f_ck, d, d_dg, expected',
    [
        (1.5, 0.02, 30, 500, 16, 0.5468),
        (1.4, 0.03, 40, 450, 20, 0.8236),
        (1.6, 0.025, 35, 600, 18, 0.5690),
    ],
)
def test_calculate_tau_Rdc_0(gamma_v, rho_l, f_ck, d, d_dg, expected):
    """Test the sh stress resistance wo/ axial force effects."""
    result = _section_8_2_shear.tau_Rdc_0(gamma_v, rho_l, f_ck, d, d_dg)
    assert math.isclose(
        result,
        expected,
        rel_tol=1e-3,
    )


@pytest.mark.parametrize(
    'gamma_v, rho_l, f_ck, d, d_dg',
    [
        (0, 0.02, 30, 500, 16),  # gamma_v zero
        (-1, 0.02, 30, 500, 16),  # gamma_v negative
        (1.5, 0, 30, 500, 16),  # rho_l zero
        (1.5, -0.01, 30, 500, 16),  # rho_l negative
        (1.5, 0.02, 0, 500, 16),  # f_ck zero
        (1.5, 0.02, -30, 500, 16),  # f_ck negative
        (1.5, 0.02, 30, 0, 16),  # d zero
        (1.5, 0.02, 30, -500, 16),  # d negative
        (1.5, 0.02, 30, 500, 0),  # d_dg zero
        (1.5, 0.02, 30, 500, -16),  # d_dg negative
    ],
)
def test_tau_Rdc_0_value_errors(gamma_v, rho_l, f_ck, d, d_dg):
    """Test tau_Rdc_0 raises ValueError for non-positive arguments."""
    with pytest.raises(ValueError):
        _section_8_2_shear.tau_Rdc_0(gamma_v, rho_l, f_ck, d, d_dg)


@pytest.mark.parametrize(
    'tau_Rdc_0, k1, sigma_cp, tau_Rdc_max, tau_rdc_min, expected',
    [
        (1, 0.5, 0.1, 2, 0.3, 0.95),
        (1, 0.6, 0.2, 2, 0.3, 0.88),
        (1, 0.4, 0.3, 2, 0.3, 0.88),
        (1, 0.5, 0.1, 0.8, 0.3, 0.8),  # Limited by tau_Rdc_max
        (1, 0.5, 0.5, 2, 0.8, 0.8),  # Limited by tau_rdc_min
    ],
)
def test_calculate_tau_Rdc_comp(
    tau_Rdc_0, k1, sigma_cp, tau_Rdc_max, tau_rdc_min, expected
):
    """Test the calculation of the shear considering comp normal forces."""
    assert math.isclose(
        _section_8_2_shear.tau_Rdc_comp(
            tau_Rdc_0, k1, sigma_cp, tau_Rdc_max, tau_rdc_min
        ),
        expected,
        rel_tol=1e-5,
    )


@pytest.mark.parametrize(
    'tau_Rdc_0, k1, sigma_cp, tau_Rdc_max, tau_rdc_min',
    [
        (0, 0.5, 0.1, 2, 0.3),  # tau_Rdc_0 zero
        (-1, 0.5, 0.1, 2, 0.3),  # tau_Rdc_0 negative
        (1, 0, 0.1, 2, 0.3),  # k1 zero
        (1, -0.5, 0.1, 2, 0.3),  # k1 negative
        (1, 0.5, -0.1, 2, 0.3),  # sigma_cp negative
        (1, 0.5, 0.1, 0, 0.3),  # tau_Rdc_max zero
        (1, 0.5, 0.1, -2, 0.3),  # tau_Rdc_max negative
        (1, 0.5, 0.1, 2, 0),  # tau_rdc_min zero
        (1, 0.5, 0.1, 2, -0.3),  # tau_rdc_min negative
    ],
)
def test_tau_Rdc_comp_value_errors(
    tau_Rdc_0, k1, sigma_cp, tau_Rdc_max, tau_rdc_min
):
    """Test tau_Rdc_comp raises ValueError."""
    with pytest.raises(ValueError):
        _section_8_2_shear.tau_Rdc_comp(
            tau_Rdc_0, k1, sigma_cp, tau_Rdc_max, tau_rdc_min
        )


@pytest.mark.parametrize(
    'a_cs_0, e_p, A_c, b_w, z, d, expected',
    [
        (1000, 50, 10000, 200, 500, 200, 0.000000429),
        (1200, 60, 12000, 250, 600, 200, 0.000000263),
        (1100, 55, 11000, 220, 550, 200, 0.0000003396),
    ],
)
def test_calculate_k1(a_cs_0, e_p, A_c, b_w, z, d, expected):
    """Test the calculation of the factor k1."""
    result = _section_8_2_shear.k1(a_cs_0, e_p, A_c, b_w, z, d)
    assert math.isclose(result, expected, rel_tol=1e-3)


@pytest.mark.parametrize(
    'a_cs_0, e_p, A_c, b_w, z, d',
    [
        (0, 50, 10000, 200, 500, 200),  # a_cs_0 zero
        (-1000, 50, 10000, 200, 500, 200),  # a_cs_0 negative
        (1000, -1, 10000, 200, 500, 200),  # e_p negative
        (1000, 50, 0, 200, 500, 200),  # A_c zero
        (1000, 50, -10000, 200, 500, 200),  # A_c negative
        (1000, 50, 10000, 0, 500, 200),  # b_w zero
        (1000, 50, 10000, -200, 500, 200),  # b_w negative
        (1000, 50, 10000, 200, 0, 200),  # z zero
        (1000, 50, 10000, 200, -500, 200),  # z negative
        (1000, 50, 10000, 200, 500, 0),  # d zero
        (1000, 50, 10000, 200, 500, -200),  # d negative
    ],
)
def test_k1_value_errors(a_cs_0, e_p, A_c, b_w, z, d):
    """Test k1 raises ValueError for non-positive or negative arguments."""
    with pytest.raises(ValueError):
        _section_8_2_shear.k1(a_cs_0, e_p, A_c, b_w, z, d)


@pytest.mark.parametrize(
    'tau_Rdc_0, a_cs_0, d, expected',
    [
        (1, 1000, 500, 2.4132),
        (0.9, 900, 450, 2.1719),
        (1.1, 1100, 550, 2.6546),
    ],
)
def test_calculate_tau_Rdc_max(tau_Rdc_0, a_cs_0, d, expected):
    """Test the calculation of the maximum shear stress resistance."""
    result = _section_8_2_shear.tau_Rdc_max(tau_Rdc_0, a_cs_0, d)
    assert math.isclose(
        result,
        expected,
        rel_tol=1e-3,
    )


@pytest.mark.parametrize(
    'tau_Rdc_0, a_cs_0, d',
    [
        (0, 1000, 500),  # tau_Rdc_0 zero
        (-1, 1000, 500),  # tau_Rdc_0 negative
        (1, 0, 500),  # a_cs_0 zero
        (1, -1000, 500),  # a_cs_0 negative
        (1, 1000, 0),  # d zero
        (1, 1000, -500),  # d negative
        (0, 0, 0),  # all zero
        (-1, -1000, -500),  # all negative
    ],
)
def test_tau_Rdc_max_value_errors(tau_Rdc_0, a_cs_0, d):
    """Test tau_Rdc_max raises ValueError for non-positive arguments."""
    with pytest.raises(ValueError):
        _section_8_2_shear.tau_Rdc_max(tau_Rdc_0, a_cs_0, d)


@pytest.mark.parametrize(
    'ds, As, dp, Ap, expected',
    [
        (500, 2000, 600, 1500, 545.45),
        (400, 1000, 450, 500, 418.0),
    ],
)
def test_d_eff_p(ds, As, dp, Ap, expected):
    """Test the calculation of effective depth."""
    result = _section_8_2_shear.d_eff_p(ds, As, dp, Ap)
    assert math.isclose(
        result,
        expected,
        rel_tol=1e-2,
    )


@pytest.mark.parametrize(
    'ds, As, dp, Ap',
    [
        (-1, 2000, 600, 1500),  # Negative ds
        (500, -2000, 600, 1500),  # Negative As
        (500, 2000, -600, 1500),  # Negative dp
        (500, 2000, 600, -1500),  # Negative Ap
        (-1, -2000, -600, -1500),  # All negative
        (0, 0, 0, 0),  # All zero (division by zero)
        (500, 0, 600, 0),  # As and Ap zero (division by zero)
        (0, 2000, 0, 1500),  # ds and dp zero (division by zero)
    ],
)
def test_d_eff_p_value_errors(ds, As, dp, Ap):
    """Test that d_eff_p raises ValueError for negative arguments or
    division by zero.
    """
    with pytest.raises(ValueError):
        _section_8_2_shear.d_eff_p(ds, As, dp, Ap)


@pytest.mark.parametrize(
    'ds, As, dp, Ap, bw, d, expected',
    [
        (500, 2000, 600, 1500, 300, 545.45, 0.02128),
        (400, 1000, 450, 500, 250, 425.00, 0.0138),
    ],
)
def test_calculate_reinforcement_ratio(ds, As, dp, Ap, bw, d, expected):
    """Test the calculation of reinforcement ratio."""
    result = _section_8_2_shear.rho_l_p(ds, As, dp, Ap, bw, d)
    assert math.isclose(
        result,
        expected,
        rel_tol=1e-2,
    )


@pytest.mark.parametrize(
    'ds, As, dp, Ap, bw, d',
    [
        (-1, 2000, 600, 1500, 300, 545.45),  # Negative ds
        (500, -2000, 600, 1500, 300, 545.45),  # Negative As
        (500, 2000, -600, 1500, 300, 545.45),  # Negative dp
        (500, 2000, 600, -1500, 300, 545.45),  # Negative Ap
        (500, 2000, 600, 1500, -300, 545.45),  # Negative bw
        (500, 2000, 600, 1500, 0, 545.45),  # Zero bw
        (500, 2000, 600, 1500, 300, -545.45),  # Negative d
    ],
)
def test_rho_l_p_value_errors(ds, As, dp, Ap, bw, d):
    """Test that rho_l_p raises ValueError for neg arguments or zero bw."""
    with pytest.raises(ValueError):
        _section_8_2_shear.rho_l_p(ds, As, dp, Ap, bw, d)


@pytest.mark.parametrize(
    'vEd_y, vEd_x, rho_l_x, rho_l_y, expected',
    [
        (20, 40, 0.005, 0.008, 0.005),  # vEd_y/vEd_x <= 0.5
        (40, 20, 0.005, 0.008, 0.008),  # vEd_y/vEd_x >= 2
        (30, 40, 0.005, 0.008, 0.00308),  # 0.5 < vEd_y/vEd_x < 2
        (20, 0, 0.005, 0.008, 0.008),  # vEd_x = 0 (returns rho_l_y)
        (25, 50, 0.005, 0.008, 0.005),  # Additional test to cover line 598
    ],
)
def test_rho_l_planar(vEd_y, vEd_x, rho_l_x, rho_l_y, expected):
    """Test the calculation of reinforcement ratio for planar members."""
    result = _section_8_2_shear.rho_l_planar(vEd_y, vEd_x, rho_l_x, rho_l_y)
    assert math.isclose(
        result,
        expected,
        rel_tol=1e-2,
    )


@pytest.mark.parametrize(
    'vEd_y, vEd_x, rho_l_x, rho_l_y',
    [
        (-1.0, 40, 0.005, 0.008),  # Negative vEd_y
        (20, -40, 0.005, 0.008),  # Negative vEd_x
        (20, 40, -0.005, 0.008),  # Negative rho_l_x
        (20, 40, 0.005, -0.008),  # Negative rho_l_y
    ],
)
def test_rho_l_planar_value_errors(vEd_y, vEd_x, rho_l_x, rho_l_y):
    """Test that rho_l_planar raises ValueError for negative inputs."""
    with pytest.raises(ValueError):
        _section_8_2_shear.rho_l_planar(vEd_y, vEd_x, rho_l_x, rho_l_y)


# Tests using pytest
def test_cot_theta_min():
    """Tests the function cot_theta_min with various scenarios."""
    Ac = 100000  # 100,000 mm2
    # No axial force (NEd = 0)
    assert _section_8_2_shear.cot_theta_min(0, 100, 10, 400, Ac) == 2.5

    # Tension case (NEd > 0)
    assert _section_8_2_shear.cot_theta_min(10, 100, 10, 400, Ac) == 2.49

    # Compression case with x < 0.25d (uses interpolation)
    # NEd = -100 kN, sigma_c = 1 MPa, interpolation: 2.5 + 1.0/6.0 = 2.666...
    assert _section_8_2_shear.cot_theta_min(
        -100, 100, 50, 400, Ac
    ) == pytest.approx(2.5 + 1.0 / 6.0, rel=1e-9)

    # Compression case with significant stress and x < 0.25d
    NEd = -300  # -300 kN (3 MPa stress: 300 * 1000 / 100000 = 3.0 MPa)
    assert _section_8_2_shear.cot_theta_min(NEd, 100, 50, 400, Ac) == 3.0

    # Compression case with very high stress (6 MPa) and x < 0.25d
    NEd_high = -600  # -600 kN (6 MPa stress: 600 * 1000 / 100000 = 6.0 MPa)
    assert _section_8_2_shear.cot_theta_min(NEd_high, 100, 50, 400, Ac) == 3.0

    # Compression case with intermediate stress (4.5 MPa) and x < 0.25d
    # sigma_c >= 3.0, so returns 3.0 (capped)
    NEd_intermediate = -450  # -450 kN (4.5 MPa stress)
    assert (
        _section_8_2_shear.cot_theta_min(NEd_intermediate, 100, 50, 400, Ac)
        == 3.0
    )

    # Compression case with significant stress but x >= 0.25d (tension formula)
    # sigma_c = 3.0 MPa, x = 120 >= 100, so: 2.5 - 0.1*(-300)/100 = 2.8
    assert _section_8_2_shear.cot_theta_min(NEd, 100, 120, 400, Ac) == 2.8

    # Compression case with low stress and x < 0.25d (uses interpolation)
    NEd_low = -100  # -100 kN (1 MPa stress: 100 * 1000 / 100000 = 1.0 MPa)
    assert _section_8_2_shear.cot_theta_min(
        NEd_low, 100, 50, 400, Ac
    ) == pytest.approx(2.5 + 1.0 / 6.0, rel=1e-9)

    # Test ductility class A (20% reduction)
    assert _section_8_2_shear.cot_theta_min(
        0, 100, 10, 400, Ac, apply_ductility_class_a_reduction=True
    ) == pytest.approx(2.0, rel=1e-9)
    assert _section_8_2_shear.cot_theta_min(
        10, 100, 10, 400, Ac, apply_ductility_class_a_reduction=True
    ) == pytest.approx(1.992, rel=1e-9)
    assert _section_8_2_shear.cot_theta_min(
        NEd_high, 100, 50, 400, Ac, apply_ductility_class_a_reduction=True
    ) == pytest.approx(2.4, rel=1e-9)


@pytest.mark.parametrize(
    'NEd, VEd, x, d, Ac, expected',
    [
        # No axial force cases
        (0, 100, 10, 400, 100000, 2.5),
        (0, 200, 20, 500, 100000, 2.5),
        # Tension cases (NEd > 0)
        (10, 100, 10, 400, 100000, 2.49),
        (50, 200, 20, 500, 100000, 2.475),
        (100, 100, 10, 400, 100000, 2.4),
        (200, 100, 10, 400, 100000, 2.3),
        (500, 100, 10, 400, 100000, 2.0),  # 2.5 - 0.1*500/100 = 2.0
        (1000, 100, 10, 400, 100000, 1.5),  # 2.5 - 0.1*1000/100 = 1.5
        # Compression cases with x < 0.25d (uses interpolation)
        (-10, 100, 10, 400, 100000, 2.5 + 0.1 / 6.0),  # 0.1 MPa
        (-50, 200, 20, 500, 100000, 2.5 + 0.5 / 6.0),  # 0.5 MPa
        (-100, 100, 10, 400, 100000, 2.5 + 1.0 / 6.0),  # 1.0 MPa
        (-200, 100, 10, 400, 100000, 2.5 + 2.0 / 6.0),  # 2.0 MPa
        # Compression cases with significant stress and x < 0.25d
        (-300, 100, 50, 400, 100000, 3.0),  # 3 MPa stress, x < 0.25d
        (-400, 200, 80, 500, 100000, 3.0),  # 4 MPa stress, x < 0.25d
        (-500, 100, 90, 400, 100000, 3.0),  # 5 MPa stress, x < 0.25d
        # Compression cases with significant stress but x >= 0.25d
        # (tension formula)
        (-300, 100, 120, 400, 100000, 2.8),  # 3 MPa, x >= 0.25d
        (-400, 200, 150, 500, 100000, 2.7),  # 4 MPa, x >= 0.25d
        # Compression cases with low stress and x < 0.25d (uses interpolation)
        (-100, 100, 50, 400, 100000, 2.5 + 1.0 / 6.0),  # 1 MPa stress
        (-200, 200, 80, 500, 100000, 2.5 + 2.0 / 6.0),  # 2 MPa stress
    ],
)
def test_cot_theta_min_comprehensive(NEd, VEd, x, d, Ac, expected):
    """Test cot_theta_min with comprehensive scenarios."""
    result = _section_8_2_shear.cot_theta_min(NEd, VEd, x, d, Ac)
    assert result == pytest.approx(expected, rel=1e-9)


def test_cot_theta_min_edge_cases():
    """Test cot_theta_min edge cases and boundary conditions."""
    # Test exactly at 3 MPa stress boundary
    Ac = 100000  # 100,000 mm2
    NEd_exact = -300  # Exactly 3 MPa: 300 * 1000 / 100000 = 3.0 MPa
    assert _section_8_2_shear.cot_theta_min(NEd_exact, 100, 50, 400, Ac) == 3.0

    # Test just below 3 MPa stress
    NEd_below = -299  # 2.99 MPa: 299 * 1000 / 100000 = 2.99 MPa
    # sigma_c = 2.99 MPa, x < 0.25d, so: 2.5 + 2.99/6.0 = 2.9983...
    assert _section_8_2_shear.cot_theta_min(
        NEd_below, 100, 50, 400, Ac
    ) == pytest.approx(2.5 + 2.99 / 6.0, rel=1e-9)

    # Test exactly at x = 0.25d boundary (x = 100, not < 100, so tension)
    x_exact = 100  # 0.25 * 400
    # sigma_c = 3.0 MPa, x = 100 (not < 100), so: 2.5 - 0.1*(-300)/100 = 2.8
    assert (
        _section_8_2_shear.cot_theta_min(NEd_exact, 100, x_exact, 400, Ac)
        == 2.8
    )

    # Test just above x = 0.25d boundary (tension formula)
    x_above = 101  # Just above 0.25 * 400
    assert (
        _section_8_2_shear.cot_theta_min(NEd_exact, 100, x_above, 400, Ac)
        == 2.8
    )

    # Test minimum tension value (1.0)
    NEd_high_tension = 1500  # High tension force
    result = _section_8_2_shear.cot_theta_min(
        NEd_high_tension, 100, 10, 400, Ac
    )
    assert result == 1.0

    # Test with zero VEd when NEd is zero (valid case)
    assert _section_8_2_shear.cot_theta_min(0, 0, 10, 400, Ac) == 2.5

    # Test with zero VEd when NEd is not zero (should raise ValueError)
    with pytest.raises(ValueError, match='VEd must not be zero'):
        _section_8_2_shear.cot_theta_min(10, 0, 10, 400, Ac)


def test_cot_theta_min_ductility_class_a_reduction():
    """Test cot_theta_min ductility class A reduction."""
    Ac = 100000  # 100,000 mm2
    # Test with and without reduction
    base_value = _section_8_2_shear.cot_theta_min(0, 100, 10, 400, Ac)
    reduced_value = _section_8_2_shear.cot_theta_min(
        0, 100, 10, 400, Ac, apply_ductility_class_a_reduction=True
    )

    assert reduced_value == base_value * 0.8
    assert reduced_value == 2.0


def test_cot_theta_min_interpolation():
    """Test cot_theta_min interpolation for intermediate stress values."""
    Ac = 100000  # 100,000 mm2

    # Test various stress levels from 0 to 3 MPa and above
    test_cases = [
        (0.0, 2.5),  # Exactly 0 MPa: 2.5 + 0/6.0 = 2.5
        (1.0, 2.5 + 1.0 / 6.0),  # 1 MPa: 2.5 + 1.0/6.0 = 2.666...
        (2.0, 2.5 + 2.0 / 6.0),  # 2 MPa: 2.5 + 2.0/6.0 = 2.833...
        (2.5, 2.5 + 2.5 / 6.0),  # 2.5 MPa: 2.5 + 2.5/6.0 = 2.916...
        (3.0, 3.0),  # Exactly 3 MPa (capped at 3.0)
        (4.0, 3.0),  # 4 MPa (capped at 3.0)
        (5.0, 3.0),  # 5 MPa (capped at 3.0)
        (6.0, 3.0),  # 6 MPa (capped at 3.0)
        (7.0, 3.0),  # 7 MPa (capped at 3.0)
    ]

    for stress_mpa, expected in test_cases:
        # Convert stress (MPa) to NEd (kN): sigma_c = abs(NEd) * 1000 / Ac
        # So: abs(NEd) = sigma_c * Ac / 1000
        NEd = -stress_mpa * Ac / 1000  # Convert to kN
        result = _section_8_2_shear.cot_theta_min(NEd, 100, 50, 400, Ac)
        assert result == pytest.approx(expected, rel=1e-3)


@pytest.mark.parametrize(
    'NEd, VEd, x, d, Ac',
    [
        (0, 100, -1, 400, 100000),  # x negative
        (0, 100, 10, 0, 100000),  # d zero
        (0, 100, -1, 0, 100000),  # x negative and d zero
        (0, 100, 10, -400, 100000),  # d negative
        (0, 100, -10, -400, 100000),  # both negative
        (-10, 100, 10, 0, 100000),  # d zero with NEd negative
        (10, 100, 10, 0, 100000),  # d zero with NEd positive
        (0, 100, 10, 400, 0),  # Ac zero
        (0, 100, 10, 400, -100000),  # Ac negative
        (10, 0, 10, 400, 100000),  # VEd zero with NEd positive
        (-10, 0, 10, 400, 100000),  # VEd zero with NEd negative
    ],
)
def test_cot_theta_min_value_error_all_cases(NEd, VEd, x, d, Ac):
    """Test cot_theta_min raises ValueError."""
    with pytest.raises(ValueError):
        _section_8_2_shear.cot_theta_min(NEd, VEd, x, d, Ac)


@pytest.mark.parametrize(
    'NEd, VEd, x, d, Ac',
    [
        (0, 100, 10, 0, 100000),  # d is zero
        (0, 100, -10, 400, 100000),  # x is negative
        (0, 100, -10, 0, 100000),  # x negative and d zero
        (0, 100, 10, -400, 100000),  # d negative
        (0, 100, -10, -400, 100000),  # both negative
        (0, 100, 10, 400, 0),  # Ac zero
        (0, 100, 10, 400, -100000),  # Ac negative
        (10, 0, 10, 400, 100000),  # VEd zero with NEd positive
        (-10, 0, 10, 400, 100000),  # VEd zero with NEd negative
    ],
)
def test_cot_theta_min_value_error(NEd, VEd, x, d, Ac):
    """Test cot_theta_min raises ValueError for invalid dimensions."""
    with pytest.raises(ValueError):
        _section_8_2_shear.cot_theta_min(NEd, VEd, x, d, Ac)


def test_tau_Rd_sy():
    """Tests the function tau_Rd_sy."""
    assert _section_8_2_shear.tau_Rd_sy(0.01, 500, 2, 2.5) == 10
    assert _section_8_2_shear.tau_Rd_sy(0.02, 400, 2.5, 3.0) == 20
    # Test boundary cases
    # cot_theta = 1 (minimum)
    assert _section_8_2_shear.tau_Rd_sy(0.01, 500, 1.0, 2.5) == 5
    # cot_theta = cot_theta_min (maximum)
    assert _section_8_2_shear.tau_Rd_sy(0.01, 500, 2.5, 2.5) == 12.5
    # cot_theta_min = 1 (minimum allowed)
    assert _section_8_2_shear.tau_Rd_sy(0.01, 500, 1.0, 1.0) == 5


@pytest.mark.parametrize(
    'rho_w, fywd, cot_theta, cot_theta_min',
    [
        (-0.01, 500, 2, 2.5),  # Negative rho_w
        (0.01, -500, 2, 2.5),  # Negative fywd
        (-0.01, -500, 2, 2.5),  # Both negative
        (0, 500, 2, 2.5),  # Zero rho_w
        (0.01, 0, 2, 2.5),  # Zero fywd
        (0, 0, 2, 2.5),  # Both zero
        (-0.01, 0, 2, 2.5),  # Negative rho_w, zero fywd
        (0, -500, 2, 2.5),  # Zero rho_w, negative fywd
        (0.01, 500, 2, 0),  # Zero cot_theta_min
        (0.01, 500, 2, -2.5),  # Negative cot_theta_min
        (0.01, 500, 2, 0.5),  # cot_theta_min < 1
        (0.01, 500, 2, 0.9),  # cot_theta_min < 1 (edge case)
        (0.01, 500, 0.5, 2.5),  # cot_theta < 1
        (0.01, 500, 3.0, 2.5),  # cot_theta > cot_theta_min
        (0.01, 500, 2.6, 2.5),  # cot_theta > cot_theta_min (edge case)
    ],
)
def test_tau_Rd_sy_value_errors(rho_w, fywd, cot_theta, cot_theta_min):
    """Test tau_Rd_sy raises ValueError for invalid arguments."""
    with pytest.raises(ValueError):
        _section_8_2_shear.tau_Rd_sy(rho_w, fywd, cot_theta, cot_theta_min)


def test_rho_w():
    """Tests the function rho_w."""
    assert _section_8_2_shear.rho_w(100, 200, 50) == 0.01
    assert _section_8_2_shear.rho_w(200, 200, 100) == 0.01


@pytest.mark.parametrize(
    'Asw, bw, s',
    [
        (-100, 200, 50),  # Negative Asw
        (100, -200, 50),  # Negative bw
        (100, 200, -50),  # Negative s
        (0, 200, 50),  # Zero Asw
        (100, 0, 50),  # Zero bw
        (100, 200, 0),  # Zero s
        (0, 0, 0),  # All zero
        (-100, -200, -50),  # All negative
    ],
)
def test_rho_w_value_errors(Asw, bw, s):
    """Test that rho_w raises ValueError for non-positive arguments."""
    with pytest.raises(ValueError):
        _section_8_2_shear.rho_w(Asw, bw, s)


def test_sigma_cd():
    """Tests the function sigma_cd."""
    # tan_theta = 1/cot_theta = 1/2 = 0.5
    # sigma_cd = 1 * (2 + 0.5) = 2.5, min(2.5, 0.5 * 20) = 2.5
    assert _section_8_2_shear.sigma_cd(1, 2, 2.5, 0.5, 20) == 2.5
    # Test with cot_theta = cot_theta_min
    assert _section_8_2_shear.sigma_cd(1, 2.5, 2.5, 0.5, 20) == 2.9
    # Test with cot_theta = 1
    assert _section_8_2_shear.sigma_cd(1, 1, 2.5, 0.5, 20) == 2.0
    # Test with cot_theta_min = 1 (boundary case)
    assert _section_8_2_shear.sigma_cd(1, 1, 1, 0.5, 20) == 2.0


@pytest.mark.parametrize(
    'tau_Ed, cot_theta, cot_theta_min, nu, f_cd',
    [
        (-1, 2, 2.5, 0.5, 20),  # Negative tau_Ed
        (1, 2, 0, 0.5, 20),  # Zero cot_theta_min
        (1, 2, -2.5, 0.5, 20),  # Negative cot_theta_min
        (1, 2, 0.5, 0.5, 20),  # cot_theta_min < 1
        (1, 2, 0.9, 0.5, 20),  # cot_theta_min < 1 (edge case)
        (1, 0.5, 2.5, 0.5, 20),  # cot_theta < 1
        (1, 3.0, 2.5, 0.5, 20),  # cot_theta > cot_theta_min
        (1, 2.6, 2.5, 0.5, 20),  # cot_theta > cot_theta_min (edge case)
        (1, 2, 2.5, -0.5, 20),  # Negative nu
        (1, 2, 2.5, 0.5, 0),  # Zero f_cd
        (1, 2, 2.5, 0.5, -20),  # Negative f_cd
        (-1, 0.5, 0.5, -0.5, -20),  # Multiple errors
    ],
)
def test_sigma_cd_value_errors(tau_Ed, cot_theta, cot_theta_min, nu, f_cd):
    """Test sigma_cd raises ValueError for invalid arguments."""
    with pytest.raises(ValueError):
        _section_8_2_shear.sigma_cd(tau_Ed, cot_theta, cot_theta_min, nu, f_cd)


def test_tau_Rd():
    """Tests the function tau_Rd."""
    assert _section_8_2_shear.tau_Rd(0.01, 500, 2, 2.5, 0.5, 50) == 10
    assert (
        _section_8_2_shear.tau_Rd(0.01, 500, 3, 3.0, 0.5, 20) == 5
    )  # Limited by the compression field
    # Test boundary cases
    assert _section_8_2_shear.tau_Rd(0.01, 500, 1, 2.5, 0.5, 50) == 5
    assert _section_8_2_shear.tau_Rd(0.01, 500, 2.5, 2.5, 0.5, 50) == 12.5


@pytest.mark.parametrize(
    'rho_w, fywd, cot_theta, cot_theta_min, nu, f_cd',
    [
        (-0.01, 500, 2, 2.5, 0.5, 50),  # Negative rho_w
        (0.01, -500, 2, 2.5, 0.5, 50),  # Negative fywd
        (0.01, 500, -2, 2.5, 0.5, 50),  # Negative cot_theta
        (0.01, 500, 2, 0.5, 0.5, 50),  # cot_theta_min < 1
        (0.01, 500, 2, 2.5, -0.5, 50),  # Negative nu
        (0.01, 500, 2, 2.5, 0.5, -50),  # Negative f_cd
        (-0.01, -500, -2, 0.5, -0.5, -50),  # All negative
        (0.01, 500, 0.5, 2.5, 0.5, 50),  # cot_theta < 1
        (0.01, 500, 3, 2.5, 0.5, 50),  # cot_theta > cot_theta_min
    ],
)
def test_tau_Rd_value_errors(rho_w, fywd, cot_theta, cot_theta_min, nu, f_cd):
    """Test tau_Rd raises ValueError for invalid arguments."""
    with pytest.raises(ValueError):
        _section_8_2_shear.tau_Rd(
            rho_w, fywd, cot_theta, cot_theta_min, nu, f_cd
        )


def test_cot_theta_simultaneous():
    """Tests the function cot_theta_simultaneous."""
    # Test basic calculation: sqrt((0.5 * 20) / (0.01 * 500) - 1) = sqrt(1) = 1
    assert (
        _section_8_2_shear.cot_theta_simultaneous(0.5, 20, 0.01, 500, 2) == 1
    )
    # Test: sqrt((0.5 * 40) / (0.01 * 500) - 1) = sqrt(3) ≈ 1.732
    assert _section_8_2_shear.cot_theta_simultaneous(
        0.5, 40, 0.01, 500, 2
    ) == pytest.approx(1.732, rel=1e-3)
    # Test clamping to minimum (1.0)
    # sqrt((0.5 * 10) / (0.01 * 500) - 1) = sqrt(0) = 0, clamped to 1
    assert (
        _section_8_2_shear.cot_theta_simultaneous(0.5, 10, 0.01, 500, 2.5)
        == 1.0
    )
    # Test clamping to cot_theta_min
    # sqrt((0.5 * 60) / (0.01 * 500) - 1) = sqrt(5) ≈ 2.236
    assert _section_8_2_shear.cot_theta_simultaneous(
        0.5, 60, 0.01, 500, 2.5
    ) == pytest.approx(2.236, rel=1e-3)

    # Test case where expression under square root is negative
    # (0.5 * 8) / (0.01 * 500) - 1 = 4 / 5 - 1 = -0.2 < 0
    with pytest.raises(ValueError, match='Expression under square root'):
        _section_8_2_shear.cot_theta_simultaneous(0.5, 8, 0.01, 500, 2.5)


@pytest.mark.parametrize(
    'nu, f_cd, rho_w, fywd, cot_theta_min',
    [
        (-0.1, 20, 0.01, 500, 2),  # Negative nu
        (0.5, -20, 0.01, 500, 2),  # Negative f_cd
        (0.5, 20, -0.01, 500, 2),  # Negative rho_w
        (0.5, 20, 0.01, -500, 2),  # Negative fywd
        (-0.5, -20, -0.01, -500, 2),  # All negative
        (0.5, 20, 0.01, 500, 0.5),  # cot_theta_min < 1
        (0.5, 20, 0, 500, 2),  # rho_w is zero
        (0.5, 20, 0.01, 0, 2),  # fywd is zero
        (0.5, 8, 0.01, 500, 2),  # Expression under square root is negative
    ],
)
def test_cot_theta_simultaneous_value_errors(
    nu, f_cd, rho_w, fywd, cot_theta_min
):
    """Test cot_theta_simultaneous raises ValueError for invalid arguments."""
    with pytest.raises(ValueError):
        _section_8_2_shear.cot_theta_simultaneous(
            nu, f_cd, rho_w, fywd, cot_theta_min
        )


@pytest.mark.parametrize(
    'Ftd, Es, Ast, expected', [(500, 210000, 1000, 0.002380952380952381)]
)
def test_epsilon_xt(Ftd, Es, Ast, expected):
    """Test εxt calculation."""
    assert _section_8_2_shear.epsilon_xt(Ftd, Es, Ast) == pytest.approx(
        expected
    )


@pytest.mark.parametrize(
    'Ftd, Es, Ast',
    [
        (500, -210000, 1000),  # Negative Es
        (500, 210000, -1000),  # Negative Ast
        (500, -210000, -1000),  # Both negative
        (500, 0, 1000),  # Es is zero
        (500, 210000, 0),  # Ast is zero
    ],
)
def test_epsilon_xt_value_errors(Ftd, Es, Ast):
    """Test epsilon_xt raises ValueError for negative or zero Es or Ast."""
    with pytest.raises(ValueError):
        _section_8_2_shear.epsilon_xt(Ftd, Es, Ast)


@pytest.mark.parametrize(
    'Fcd, Ec, Acc, expected', [(500, 30000, 1000, 0.016666666666666666)]
)
def test_epsilon_xc_compression(Fcd, Ec, Acc, expected):
    """Test εxc calculation for compression."""
    assert _section_8_2_shear.epsilon_xc_comp(Fcd, Ec, Acc) == pytest.approx(
        expected, rel=10e-3
    )


@pytest.mark.parametrize(
    'Fcd, Ec, Acc',
    [
        (500, -30000, 1000),  # Negative Ec
        (500, 30000, -1000),  # Negative Acc
        (500, -30000, -1000),  # Both negative
        (500, 0, 1000),  # Ec is zero
        (500, 30000, 0),  # Acc is zero
    ],
)
def test_epsilon_xc_comp_value_errors(Fcd, Ec, Acc):
    """Test epsilon_xc_comp raises ValueError for negative or zero Ec or
    Acc.
    """
    with pytest.raises(ValueError):
        _section_8_2_shear.epsilon_xc_comp(Fcd, Ec, Acc)


@pytest.mark.parametrize(
    'Fcd, Es, Asc, expected', [(500, 210000, 1000, 0.002380952380952381)]
)
def test_epsilon_xc_tension(Fcd, Es, Asc, expected):
    """Test εxc calculation for tension."""
    assert _section_8_2_shear.epsilon_xc_tens(Fcd, Es, Asc) == pytest.approx(
        expected, rel=10e-3
    )


@pytest.mark.parametrize(
    'Fcd, Es, Asc',
    [
        (500, -210000, 1000),  # Negative Es
        (500, 210000, -1000),  # Negative Asc
        (500, -210000, -1000),  # Both negative
        (500, 0, 1000),  # Es is zero
        (500, 210000, 0),  # Asc is zero
    ],
)
def test_epsilon_xc_tens_value_errors(Fcd, Es, Asc):
    """Test epsilon_xc_tens raises ValueError for negative or zero Es or
    Asc.
    """
    with pytest.raises(ValueError):
        _section_8_2_shear.epsilon_xc_tens(Fcd, Es, Asc)


@pytest.mark.parametrize(
    'epsilon_xt, epsilon_xc, expected', [(0.002, 0.001, 0.0015)]
)
def test_epsilon_x(epsilon_xt, epsilon_xc, expected):
    """Test εx calculation."""
    assert _section_8_2_shear.epsilon_x(
        epsilon_xt, epsilon_xc
    ) == pytest.approx(expected)


@pytest.mark.parametrize(
    'epsilon_x, cot_theta, cot_theta_min, expected',
    [(0.001, 2, 2.5, 0.5025)],
)
def test_nu(epsilon_x, cot_theta, cot_theta_min, expected):
    """Test ν calculation."""
    assert _section_8_2_shear.nu(
        epsilon_x, cot_theta, cot_theta_min
    ) == pytest.approx(expected, rel=10e-3)
    # Test boundary cases
    # For cot_theta = 1: nu = 1 / (1.0 + 110 * (0.001 + 0.002 * 1)) ≈ 0.7519
    assert _section_8_2_shear.nu(0.001, 1, 2.5) == pytest.approx(
        0.7519, rel=10e-3
    )
    # For cot_theta = 2.5: nu = 1 / (1.0 + 110 * (0.001 + 0.002 * 6.25))
    # ≈ 0.4024
    assert _section_8_2_shear.nu(0.001, 2.5, 2.5) == pytest.approx(
        0.4024, rel=10e-3
    )


@pytest.mark.parametrize(
    'epsilon_x, cot_theta, cot_theta_min',
    [
        (-0.001, 2, 2.5),  # Negative epsilon_x
        (-1.0, 2, 2.5),  # Negative epsilon_x
        (-0.5, 2, 2.5),  # Negative epsilon_x
        (0.001, 2, 0.5),  # cot_theta_min < 1
        (0.001, 0.5, 2.5),  # cot_theta < 1
        (0.001, 3, 2.5),  # cot_theta > cot_theta_min
    ],
)
def test_nu_value_errors(epsilon_x, cot_theta, cot_theta_min):
    """Test nu raises ValueError for invalid arguments."""
    with pytest.raises(ValueError):
        _section_8_2_shear.nu(epsilon_x, cot_theta, cot_theta_min)


@pytest.mark.parametrize(
    'VEd, cot_theta, cot_theta_min, expected',
    [
        (100, 2, 2.5, 200),
        (150, 1.5, 2.5, 225),
    ],
)
def test_calculate_nv(VEd, cot_theta, cot_theta_min, expected):
    """Test calculate_nv function with various inputs."""
    assert _section_8_2_shear.Nvd(
        VEd, cot_theta, cot_theta_min
    ) == pytest.approx(expected, rel=1e-6)
    # Test boundary cases
    assert _section_8_2_shear.Nvd(100, 1, 2.5) == pytest.approx(100, rel=1e-6)
    assert _section_8_2_shear.Nvd(100, 2.5, 2.5) == pytest.approx(
        250, rel=1e-6
    )


@pytest.mark.parametrize(
    'VEd, cot_theta, cot_theta_min',
    [
        (100, 2, 0.5),  # cot_theta_min < 1
        (100, 0.5, 2.5),  # cot_theta < 1
        (100, 3, 2.5),  # cot_theta > cot_theta_min
    ],
)
def test_Nvd_value_errors(VEd, cot_theta, cot_theta_min):
    """Test Nvd raises ValueError for invalid arguments."""
    with pytest.raises(ValueError):
        _section_8_2_shear.Nvd(VEd, cot_theta, cot_theta_min)


@pytest.mark.parametrize(
    'MEd, z, NVd, NE, expected',
    [
        (200, 400, 100, 50, 575.0),
        (300, 600, 150, 100, 625.0),
    ],
)
def test_calculate_ftd(MEd, z, NVd, NE, expected):
    """Test calculate_ftd function with various inputs."""
    assert _section_8_2_shear.Ftd(MEd, z, NVd, NE) == pytest.approx(
        expected, rel=1e-6
    )


@pytest.mark.parametrize(
    'MEd, z, NVd, NE, expected',
    [
        (200, 400, 100, 50, 425.0),
        (300, 600, 150, 100, 375.0),
    ],
)
def test_calculate_fcd(MEd, z, NVd, NE, expected):
    """Test calculate_fcd function with various inputs."""
    assert _section_8_2_shear.Fcd(MEd, z, NVd, NE) == pytest.approx(
        expected, rel=1e-6
    )


@pytest.mark.parametrize(
    'MEd, z, NVd, NE',
    [
        (200, 0, 100, 50),  # z is zero
        (200, -400, 100, 50),  # z is negative
    ],
)
def test_Ftd_value_errors(MEd, z, NVd, NE):
    """Test Ftd raises ValueError for invalid z."""
    with pytest.raises(ValueError):
        _section_8_2_shear.Ftd(MEd, z, NVd, NE)


@pytest.mark.parametrize(
    'MEd, z, NVd, NE',
    [
        (200, 0, 100, 50),  # z is zero
        (200, -400, 100, 50),  # z is negative
    ],
)
def test_Fcd_value_errors(MEd, z, NVd, NE):
    """Test Fcd raises ValueError for invalid z."""
    with pytest.raises(ValueError):
        _section_8_2_shear.Fcd(MEd, z, NVd, NE)


@pytest.mark.parametrize(
    'MEd_max, z, NEd, expected',
    [
        (200, 400, 50, 525.0),  # 200*1000/400 + 50/2 = 500 + 25 = 525
        (300, 600, 100, 550.0),  # 300*1000/600 + 100/2 = 500 + 50 = 550
    ],
)
def test_Ftd_max(MEd_max, z, NEd, expected):
    """Test Ftd_max calculation for direct intermediate support or concentrated
    loads.
    """
    assert _section_8_2_shear.Ftd_max(MEd_max, z, NEd) == pytest.approx(
        expected, rel=1e-6
    )


@pytest.mark.parametrize(
    'MEd_max, z, NEd',
    [
        (200, 0, 50),  # z is zero
        (200, -400, 50),  # z is negative
    ],
)
def test_Ftd_max_value_errors(MEd_max, z, NEd):
    """Test Ftd_max raises ValueError for invalid z."""
    with pytest.raises(ValueError):
        _section_8_2_shear.Ftd_max(MEd_max, z, NEd)


@pytest.mark.parametrize(
    'duct_material, is_grouted, wall_thickness, duct_diameter, expected',
    [
        ('steel', True, 3.0, 50.0, 0.5),  # Grouted steel ducts
        ('plastic', True, 1.0, 50.0, 0.8),  # Grouted plastic thin ducts
        ('plastic', True, 5.0, 50.0, 1.2),  # Grouted plastic thick ducts
        ('plastic', False, 1.0, 50.0, 1.2),  # Non-grouted plastic ducts
    ],
)
def test_calculate_k_duct(
    duct_material, is_grouted, wall_thickness, duct_diameter, expected
):
    """Test calculate_k_duct with various inputs."""
    assert (
        _section_8_2_shear.k_duct(
            duct_material, is_grouted, wall_thickness, duct_diameter
        )
        == expected
    )


@pytest.mark.parametrize(
    'duct_material, is_grouted, wall_thickness, duct_diameter',
    [
        ('steel', True, -1.0, 50.0),  # Negative wall_thickness
        ('steel', True, 2.0, -50.0),  # Negative duct_diameter
        ('plastic', False, -0.5, 40.0),  # Negative wall_thickness
        ('plastic', False, 1.0, -40.0),  # Negative duct_diameter
        ('steel', True, 0, 50.0),  # Zero wall_thickness
        ('steel', True, 2.0, 0),  # Zero duct_diameter
    ],
)
def test_k_duct_negative_dimensions(
    duct_material, is_grouted, wall_thickness, duct_diameter
):
    """Test k_duct raises ValueError for non-positive dimensions."""
    with pytest.raises(ValueError):
        _section_8_2_shear.k_duct(
            duct_material, is_grouted, wall_thickness, duct_diameter
        )


@pytest.mark.parametrize(
    'duct_material, is_grouted, wall_thickness, duct_diameter',
    [
        ('aluminum', True, 2.0, 50.0),  # Invalid material
        ('concrete', False, 1.0, 40.0),  # Invalid material
        ('', True, 1.0, 40.0),  # Empty string
        (None, True, 1.0, 40.0),  # None as material
        (123, True, 1.0, 40.0),  # Non-string type
    ],
)
def test_k_duct_invalid_material(
    duct_material, is_grouted, wall_thickness, duct_diameter
):
    """Test k_duct raises ValueError for invalid duct_material."""
    with pytest.raises(ValueError):
        _section_8_2_shear.k_duct(
            duct_material, is_grouted, wall_thickness, duct_diameter
        )


@pytest.mark.parametrize(
    'bw, duct_diameters, k_duct, expected',
    [
        (800, [50, 40], 0.5, 755.0),  # Grouted steel ducts
        (800, [50, 40], 0.8, 728.0),  # Grouted plastic thin ducts
        (800, [50, 40], 1.2, 692.0),  # Non-grouted or thick plastic ducts
    ],
)
def test_calculate_nominal_web_width(bw, duct_diameters, k_duct, expected):
    """Test calculate_nominal_web_width with various inputs."""
    assert _section_8_2_shear.bw_nom(
        bw, duct_diameters, k_duct
    ) == pytest.approx(expected, rel=1e-6)


@pytest.mark.parametrize(
    'bw, duct_diameters, k_duct',
    [
        (-800, [50, 40], 0.5),  # Negative bw
        (800, [-50, 40], 0.5),  # Negative duct diameter
        (800, [50, -40], 0.5),  # Negative duct diameter
        (800, [-50, -40], 0.5),  # All negative duct diameters
        (800, [800 / 8 + 1], 0.5),  # Sum of duct diameters exceeds bw/8
        (800, [100, 50, 50, 50, 50, 50, 50, 50, 50], 0.5),  # Sum exceeds bw/8
    ],
)
def test_bw_nom_value_errors(bw, duct_diameters, k_duct):
    """Test bw_nom raises ValueError for invalid arguments."""
    with pytest.raises(ValueError):
        _section_8_2_shear.bw_nom(bw, duct_diameters, k_duct)


@pytest.mark.parametrize(
    (
        'nu, f_cd, cot_theta, cot_theta_min, cot_beta_incl, rho_w, f_ywd, '
        'expected'
    ),
    [
        (0.6, 30, 1, 2.5, 1.5, 0.01, 500, 3),
        (0.5, 40, 1.5, 2.5, 1, 0.02, 400, 9.2307),
    ],
)
def test_tau_rd(
    nu, f_cd, cot_theta, cot_theta_min, cot_beta_incl, rho_w, f_ywd, expected
):
    """Test calculation of enhanced shear stress resistance τRd."""
    assert math.isclose(
        _section_8_2_shear.tau_rd(
            nu, f_cd, cot_theta, cot_theta_min, cot_beta_incl, rho_w, f_ywd
        ),
        expected,
        rel_tol=1e-2,
    )
    assert math.isclose(
        _section_8_2_shear.tau_rd(0.6, 30, 1, 2.5, 1.5, 0.01, 500),
        3,
        rel_tol=1e-2,
    )
    assert math.isclose(
        _section_8_2_shear.tau_rd(0.6, 30, 2.5, 2.5, 1.5, 0.01, 500),
        6.207,
        rel_tol=1e-2,
    )


@pytest.mark.parametrize(
    'nu, f_cd, cot_theta, cot_theta_min, cot_beta_incl, rho_w, f_ywd',
    [
        (-0.1, 30, 1, 2.5, 1.5, 0.01, 500),  # Negative nu
        (0.6, -30, 1, 2.5, 1.5, 0.01, 500),  # Negative f_cd
        (0.6, 0, 1, 2.5, 1.5, 0.01, 500),  # Zero f_cd
        (0.6, 30, 1, 2.5, 1.5, -0.01, 500),  # Negative rho_w
        (0.6, 30, 1, 2.5, 1.5, 0.01, -500),  # Negative f_ywd
        (0.6, 30, 1, 2.5, 1.5, 0.01, 0),  # Zero f_ywd
        (0.6, -30, 1, 2.5, 1.5, 0.01, -500),  # Both f_cd and f_ywd negative
        (0.6, 30, 1, 0.5, 1.5, 0.01, 500),  # cot_theta_min < 1
        (0.6, 30, 0.5, 2.5, 1.5, 0.01, 500),  # cot_theta < 1
        (0.6, 30, 3, 2.5, 1.5, 0.01, 500),  # cot_theta > cot_theta_min
    ],
)
def test_tau_rd_value_errors(
    nu, f_cd, cot_theta, cot_theta_min, cot_beta_incl, rho_w, f_ywd
):
    """Test tau_rd raises ValueError for invalid arguments."""
    with pytest.raises(ValueError):
        _section_8_2_shear.tau_rd(
            nu,
            f_cd,
            cot_theta,
            cot_theta_min,
            cot_beta_incl,
            rho_w,
            f_ywd,
        )


@pytest.mark.parametrize(
    'cot_beta_incl, cot_theta_min, expected',
    [
        (1.0, 2.5, 2.414213562373095),  # sqrt(2) ≈ 1.414
        (1.5, 2.5, 2.5),  # Capped at cot_theta_min
        (2.0, 3.0, 3.0),  # 2.0 + sqrt(5) ≈ 4.236, capped at 3.0
        (0.5, 2.0, 1.618033988749895),  # sqrt(1.25) ≈ 1.118
    ],
)
def test_cot_theta_max_shear_constant_nu(
    cot_beta_incl, cot_theta_min, expected
):
    """Test calculation of optimum cot_theta for constant nu."""
    result = _section_8_2_shear.cot_theta_max_shear_constant_nu(
        cot_beta_incl, cot_theta_min
    )
    assert math.isclose(result, expected, rel_tol=1e-2)


@pytest.mark.parametrize(
    'a, z, cot_theta_min, expected',
    [
        (1000, 500, 2.5, 2.5),  # 1.3 * 1000 / 500 = 2.6, capped at 2.5
        (1500, 500, 3.0, 3.0),  # 1.3 * 1500 / 500 = 3.9, capped at 3.0
        (1000, 1000, 3.0, 1.3),  # 1.3 * 1000 / 1000 = 1.3
        (500, 500, 2.5, 1.3),  # 1.3 * 500 / 500 = 1.3
    ],
)
def test_cot_theta_max_shear_variable_nu(a, z, cot_theta_min, expected):
    """Test calculation of optimum cot_theta for variable nu."""
    result = _section_8_2_shear.cot_theta_max_shear_variable_nu(
        a, z, cot_theta_min
    )
    assert math.isclose(result, expected, rel_tol=1e-2)


@pytest.mark.parametrize(
    'cot_beta_incl, cot_theta_min',
    [
        (1.0, 0.5),  # cot_theta_min < 1
        (-10.0, 1.0),  # Calculated cot_theta < 1
    ],
)
def test_cot_theta_max_shear_constant_nu_value_errors(
    cot_beta_incl, cot_theta_min
):
    """Test cot_theta_max_shear_constant_nu raises ValueError for invalid
    arguments.
    """
    with pytest.raises(ValueError):
        _section_8_2_shear.cot_theta_max_shear_constant_nu(
            cot_beta_incl, cot_theta_min
        )


@pytest.mark.parametrize(
    'a, z, cot_theta_min',
    [
        (0, 500, 2.5),  # a <= 0
        (-1000, 500, 2.5),  # a < 0
        (1000, 0, 2.5),  # z <= 0
        (1000, -500, 2.5),  # z < 0
        (1000, 500, 0.5),  # cot_theta_min < 1
        (
            100,
            10000,
            1.0,
        ),  # Calculated cot_theta < 1 (1.3 * 100 / 10000 = 0.013)
    ],
)
def test_cot_theta_max_shear_variable_nu_value_errors(a, z, cot_theta_min):
    """Test cot_theta_max_shear_variable_nu raises ValueError for invalid
    arguments.
    """
    with pytest.raises(ValueError):
        _section_8_2_shear.cot_theta_max_shear_variable_nu(a, z, cot_theta_min)


@pytest.mark.parametrize(
    'e_s, eps_x, f_ywd, cot_theta, cot_theta_min, expected',
    [
        (200000, 0.001, 500, 2, 2.5, 500),
        (210000, 0.002, 450, 1, 2.5, 420),
        (200000, 0.001, 500, 2.5, 2.5, 500),  # cot_theta = cot_theta_min
    ],
)
def test_sigma_swd(e_s, eps_x, f_ywd, cot_theta, cot_theta_min, expected):
    """Test calculation of stress σswd in shear reinforcement."""
    assert math.isclose(
        _section_8_2_shear.sigma_swd(
            e_s, eps_x, f_ywd, cot_theta, cot_theta_min
        ),
        expected,
        rel_tol=1e-2,
    )


@pytest.mark.parametrize(
    'Es, eps_x, f_ywd, cot_theta, cot_theta_min',
    [
        (-200000, 0.001, 500, 2, 2.5),  # Negative Es
        (0, 0.001, 500, 2, 2.5),  # Zero Es
        (200000, 0.001, -500, 2, 2.5),  # Negative f_ywd
        (200000, 0.001, 0, 2, 2.5),  # Zero f_ywd
        (-200000, 0.001, -500, 2, 2.5),  # Both negative
        (200000, 0.001, 500, 2, 0.5),  # cot_theta_min < 1
        (200000, 0.001, 500, 0.5, 2.5),  # cot_theta < 1
        (200000, 0.001, 500, 3, 2.5),  # cot_theta > cot_theta_min
    ],
)
def test_sigma_swd_value_errors(Es, eps_x, f_ywd, cot_theta, cot_theta_min):
    """Test sigma_swd raises ValueError for invalid arguments."""
    with pytest.raises(ValueError):
        _section_8_2_shear.sigma_swd(
            Es, eps_x, f_ywd, cot_theta, cot_theta_min
        )


@pytest.mark.parametrize(
    'tau_ed, rho_w, f_ywd, cot_theta, cot_theta_min, z, b_w, a, x, expected',
    [
        (0.5, 0.01, 500, 1, 2.5, 300, 200, 500, 250, 0),
        (0.4, 0.02, 400, 1.5, 2.5, 350, 250, 600, 200, -101.5),
        (
            0.5,
            0.01,
            500,
            2.5,
            2.5,
            300,
            200,
            500,
            250,
            0,
        ),  # cot_theta = cot_theta_min
    ],
)
def test_delta_m_ed(
    tau_ed, rho_w, f_ywd, cot_theta, cot_theta_min, z, b_w, a, x, expected
):
    """Test calculation of additional moment ΔMEd."""
    assert math.isclose(
        _section_8_2_shear.delta_MEd(
            tau_ed, rho_w, f_ywd, cot_theta, cot_theta_min, z, b_w, a, x
        ),
        expected,
        rel_tol=1e-2,
    )


@pytest.mark.parametrize(
    'tau_ed, rho_w, f_ywd, cot_theta, cot_theta_min, z, b_w, a, x',
    [
        (-0.5, 0.01, 500, 1, 2.5, 300, 200, 500, 250),  # Negative tau_ed
        (0.5, 0.01, 0, 1, 2.5, 300, 200, 500, 250),  # Zero f_ywd
        (0.5, 0.01, -500, 1, 2.5, 300, 200, 500, 250),  # Negative f_ywd
        (0.5, 0.01, 500, 1, 0.5, 300, 200, 500, 250),  # cot_theta_min < 1
        (0.5, 0.01, 500, 0.5, 2.5, 300, 200, 500, 250),  # cot_theta < 1
        (
            0.5,
            0.01,
            500,
            3,
            2.5,
            300,
            200,
            500,
            250,
        ),  # cot_theta > cot_theta_min
        (0.5, 0.01, 500, 1, 2.5, 0, 200, 500, 250),  # Zero z
        (0.5, 0.01, 500, 1, 2.5, -300, 200, 500, 250),  # Negative z
        (0.5, 0.01, 500, 1, 2.5, 300, 0, 500, 250),  # Zero b_w
        (0.5, 0.01, 500, 1, 2.5, 300, -200, 500, 250),  # Negative b_w
        (0.5, 0.01, 500, 1, 2.5, 300, 200, -500, 250),  # Negative a
        (0.5, 0.01, 500, 1, 2.5, 300, 200, 500, -250),  # Negative x
    ],
)
def test_delta_MEd_value_errors(
    tau_ed, rho_w, f_ywd, cot_theta, cot_theta_min, z, b_w, a, x
):
    """Test delta_MEd raises ValueError for invalid arguments."""
    with pytest.raises(ValueError):
        _section_8_2_shear.delta_MEd(
            tau_ed, rho_w, f_ywd, cot_theta, cot_theta_min, z, b_w, a, x
        )


@pytest.mark.parametrize(
    'cot_theta, alpha_w, cot_theta_min, expected',
    [
        (1, 45, 2.5, 1.0),  # tan(45/2) ≈ 0.414, 1.0 is within range
        (2, 60, 3.0, 2.0),  # tan(60/2) ≈ 0.577, so 2.0 is valid
        (0.5, 60, 3.0, 0.5773502691896257),  # Clamped to tan(30)
        (4, 60, 3.0, 3.0),  # Clamped to cot_theta_min
    ],
)
def test_cot_theta_inclined(cot_theta, alpha_w, cot_theta_min, expected):
    """Test calculation and validation of cot_theta for inclined
    reinforcement (Eq. 8.58).
    """
    result = _section_8_2_shear.cot_theta_inclined(
        cot_theta, alpha_w, cot_theta_min
    )
    assert math.isclose(result, expected, rel_tol=1e-2)


@pytest.mark.parametrize(
    'cot_theta, alpha_w, cot_theta_min',
    [
        (1, 44, 2.5),  # alpha_w < 45
        (1, 100, 2.5),  # alpha_w > 90
        (1, 60, 0.5),  # cot_theta_min < 1
    ],
)
def test_cot_theta_inclined_value_errors(cot_theta, alpha_w, cot_theta_min):
    """Test cot_theta_inclined raises ValueError for invalid arguments."""
    with pytest.raises(ValueError):
        _section_8_2_shear.cot_theta_inclined(
            cot_theta, alpha_w, cot_theta_min
        )


@pytest.mark.parametrize(
    'rho_w, f_ywd, cot_theta, alpha_w, cot_theta_min, expected',
    [
        (0.01, 500, 1, 45, 3.0, 7.071),
        (0.02, 450, 2, 60, 3.0, 20.088),
        (0.01, 500, 1, 90, 3.0, 5.0),  # alpha_w = 90°
    ],
)
def test_tau_Rd_sy_inclined(
    rho_w, f_ywd, cot_theta, alpha_w, cot_theta_min, expected
):
    """Test calculation of shear stress resistance τRd,sy for inclined
    reinforcement (Eq. 8.59).
    """
    assert math.isclose(
        _section_8_2_shear.tau_Rd_sy_inclined(
            rho_w, f_ywd, cot_theta, alpha_w, cot_theta_min
        ),
        expected,
        rel_tol=1e-2,
    )


@pytest.mark.parametrize(
    'rho_w, f_ywd, cot_theta, alpha_w, cot_theta_min',
    [
        (-0.01, 500, 2, 60, 3.0),  # Negative rho_w
        (0.01, 0, 2, 60, 3.0),  # Zero f_ywd
        (0.01, -500, 2, 60, 3.0),  # Negative f_ywd
        (0.01, 500, 2, 44, 3.0),  # alpha_w < 45
        (0.01, 500, 2, 100, 3.0),  # alpha_w > 90
        (0.01, 500, 2, 60, 0.5),  # cot_theta_min < 1
    ],
)
def test_tau_Rd_sy_inclined_value_errors(
    rho_w, f_ywd, cot_theta, alpha_w, cot_theta_min
):
    """Test tau_Rd_sy_inclined raises ValueError for invalid arguments."""
    with pytest.raises(ValueError):
        _section_8_2_shear.tau_Rd_sy_inclined(
            rho_w, f_ywd, cot_theta, alpha_w, cot_theta_min
        )


@pytest.mark.parametrize(
    'tau_ed, theta, alpha_w, nu, f_cd, cot_theta_min, expected',
    [
        (0.5, 1, 45, 0.6, 30, 2, 0.5),
        (0.7759, 2, 60, 0.5, 40, 2, 1.50522),
        (0.5, 1, 90, 0.6, 30, 2, 1.0),  # alpha_w = 90°
    ],
)
def test_sigma_cd_inclined(
    tau_ed, theta, alpha_w, nu, f_cd, cot_theta_min, expected
):
    """Test calculation of compression stress σcd for inclined reinforcement
    (Eq. 8.60).
    """
    assert math.isclose(
        _section_8_2_shear.sigma_cd_inclined(
            tau_ed, theta, alpha_w, nu, f_cd, cot_theta_min
        ),
        expected,
        rel_tol=1e-2,
    )


@pytest.mark.parametrize(
    'tau_ed, cot_theta, alpha_w, nu, f_cd, cot_theta_min',
    [
        (-0.5, 1, 60, 0.6, 30, 2),  # Negative tau_ed
        (0.5, 1, 60, -0.6, 30, 2),  # Negative nu
        (0.5, 1, 60, 0.6, 0, 2),  # Zero f_cd
        (0.5, 1, 60, 0.6, -30, 2),  # Negative f_cd
        (0.5, 1, 44, 0.6, 30, 2),  # alpha_w < 45
        (0.5, 1, 100, 0.6, 30, 2),  # alpha_w > 90
        (0.5, 1, 60, 0.6, 30, 0.5),  # cot_theta_min < 1
    ],
)
def test_sigma_cd_inclined_value_errors(
    tau_ed, cot_theta, alpha_w, nu, f_cd, cot_theta_min
):
    """Test sigma_cd_inclined raises ValueError for invalid arguments."""
    with pytest.raises(ValueError):
        _section_8_2_shear.sigma_cd_inclined(
            tau_ed, cot_theta, alpha_w, nu, f_cd, cot_theta_min
        )


@pytest.mark.parametrize(
    'v_ed, theta, alpha_w, cot_theta_min, expected',
    [
        (100, 0.3, 45, 1, -58.5786),
        (80, 0.1, 50, 3, -29.8233),
        (100, 1, 90, 2, 100.0),  # alpha_w = 90°
    ],
)
def test_n_vd(v_ed, theta, alpha_w, cot_theta_min, expected):
    """Test calculation of axial tensile force NVd."""
    assert math.isclose(
        _section_8_2_shear.NVds_inclined(v_ed, theta, alpha_w, cot_theta_min),
        expected,
        rel_tol=1e-2,
    )


@pytest.mark.parametrize(
    'VEd, cot_theta, alpha_w, cot_theta_min',
    [
        (100, 1, 44, 2),  # alpha_w < 45
        (100, 1, 100, 2),  # alpha_w > 90
        (100, 1, 60, 0.5),  # cot_theta_min < 1
    ],
)
def test_NVds_value_errors(VEd, cot_theta, alpha_w, cot_theta_min):
    """Test NVds_inclined raises ValueError for invalid arguments."""
    with pytest.raises(ValueError):
        _section_8_2_shear.NVds_inclined(
            VEd, cot_theta, alpha_w, cot_theta_min
        )


@pytest.mark.parametrize(
    'nu, f_cd, theta, beta_incl, rho_w, f_ywd, '
    + 'alpha_w, cot_theta_min, expected',
    [
        (0.6, 30, 1, 3, 0.01, 500, 45, 2, -3.8578),
        (0.5, 40, 0.7, 1, 0.02, 400, 60, 5, 6.9013),
        (0.6, 30, 1, 3, 0.01, 500, 90, 2, -3.0),  # alpha_w = 90°
    ],
)
def test_tau_rd_incl(
    nu, f_cd, theta, beta_incl, rho_w, f_ywd, alpha_w, cot_theta_min, expected
):
    """Test calculation of shear stress resistance τRd."""
    assert math.isclose(
        _section_8_2_shear.tau_Rd_inclined(
            nu, f_cd, theta, beta_incl, rho_w, f_ywd, alpha_w, cot_theta_min
        ),
        expected,
        rel_tol=1e-2,
    )


@pytest.mark.parametrize(
    'nu, f_cd, cot_theta, cot_beta_incl, rho_w, f_ywd, alpha_w, cot_theta_min',
    [
        (-0.1, 30, 1, 3, 0.01, 500, 45, 2),  # Negative nu
        (0.6, -30, 1, 3, 0.01, 500, 45, 2),  # Negative f_cd
        (0.6, 0, 1, 3, 0.01, 500, 45, 2),  # Zero f_cd
        (0.6, 30, 1, 3, -0.01, 500, 45, 2),  # Negative rho_w
        (0.6, 30, 1, 3, 0.01, -500, 45, 2),  # Negative f_ywd
        (0.6, 30, 1, 3, 0.01, 0, 45, 2),  # Zero f_ywd
        (0.6, 30, 1, 3, 0.01, 500, 44, 2),  # alpha_w < 45
        (0.6, 30, 1, 3, 0.01, 500, 100, 2),  # alpha_w > 90
        (0.6, 30, 1, 3, 0.01, 500, 45, 0.5),  # cot_theta_min < 1
    ],
)
def test_tau_rd_incl_value_errors(
    nu, f_cd, cot_theta, cot_beta_incl, rho_w, f_ywd, alpha_w, cot_theta_min
):
    """Test tau_rd_incl raises ValueError for invalid arguments."""
    with pytest.raises(ValueError):
        _section_8_2_shear.tau_Rd_inclined(
            nu,
            f_cd,
            cot_theta,
            cot_beta_incl,
            rho_w,
            f_ywd,
            alpha_w,
            cot_theta_min,
        )


@pytest.mark.parametrize(
    'Es, eps_x, cot_theta, alpha_w, f_ywd',
    [
        (-200000, 0.001, 2, 45, 500),  # Negative Es
        (0, 0.001, 2, 45, 500),  # Zero Es
        (200000, 0.001, 2, 45, -500),  # Negative f_ywd
        (200000, 0.001, 2, 45, 0),  # Zero f_ywd
        (200000, 0.001, 2, 44, 500),  # alpha_w < 45
        (200000, 0.001, 2, 100, 500),  # alpha_w > 90
        (-200000, 0.001, 2, 45, -500),  # Both negative
    ],
)
def test_sigma_swd_v2_value_errors(Es, eps_x, cot_theta, alpha_w, f_ywd):
    """Test sigma_swd_inclined raises ValueError for invalid arguments."""
    with pytest.raises(ValueError):
        _section_8_2_shear.sigma_swd_inclined(
            Es, eps_x, cot_theta, alpha_w, f_ywd
        )


# Test for sigma_swd_inclined with valid parameters (covers lines 1863-1874)
@pytest.mark.parametrize(
    'Es, eps_x, cot_theta, alpha_w, f_ywd, expected',
    [
        # Test with alpha_w = 90.0 to cover lines 1863-1864 (isclose check)
        (200000, 0.001, 1.5, 90.0, 500, 500),
        # Test with alpha_w = 45.0 to cover line 1866 (else branch)
        (200000, 0.001, 1.5, 45.0, 500, 500),
        # Test with alpha_w = 60.0 to cover line 1866 (else branch)
        (200000, 0.001, 1.5, 60.0, 500, 500),
        # Test with alpha_w = 89.0 to cover line 1866 (else branch)
        (200000, 0.001, 1.5, 89.0, 500, 500),
        # Test with smaller eps_x to get a value less than f_ywd
        (200000, 0.0001, 1.5, 90.0, 500, 295.0),
    ],
)
def test_sigma_swd_inclined_valid(
    Es, eps_x, cot_theta, alpha_w, f_ywd, expected
):
    """Test sigma_swd_inclined with valid parameters.

    Covers lines 1863-1874.
    """
    result = _section_8_2_shear.sigma_swd_inclined(
        Es, eps_x, cot_theta, alpha_w, f_ywd
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'tau_rd, m_ed, m_rd, expected',
    [
        (5.0, 2.0, 10.0, 4.0),
        (10.0, 3.0, 15.0, 8.0),
    ],
)
def test_shear_stress_resistance_reduced(tau_rd, m_ed, m_rd, expected):
    """Test calculation of reduced shear stress resistance."""
    assert _section_8_2_shear.tau_Rdm(tau_rd, m_ed, m_rd) == pytest.approx(
        expected, rel=10e-2
    )


@pytest.mark.parametrize(
    'tau_rd, m_ed, m_rd',
    [
        (-1.0, 2.0, 10.0),  # Negative tau_rd
        (5.0, -2.0, 10.0),  # Negative m_ed
        (5.0, 2.0, -10.0),  # Negative m_rd
        (5.0, 2.0, 0),  # Zero m_rd
        (-1.0, -2.0, -10.0),  # All negative
    ],
)
def test_tao_Rd_m_value_errors(tau_rd, m_ed, m_rd):
    """Test tau_Rdm raises ValueError for invalid arguments."""
    with pytest.raises(ValueError):
        _section_8_2_shear.tau_Rdm(tau_rd, m_ed, m_rd)


@pytest.mark.parametrize(
    'delta_fd, hf, delta_x, expected',
    [
        (100.0, 200.0, 1000.0, 0.5),
        (50.0, 150.0, 500.0, 0.666),
    ],
)
def test_longitudinal_shear_stress(delta_fd, hf, delta_x, expected):
    """Test calculation of longitudinal shear stress."""
    assert _section_8_2_shear.tau_Ed_flange(
        delta_fd, hf, delta_x
    ) == pytest.approx(expected, rel=10e-2)


@pytest.mark.parametrize(
    'delta_fd, hf, delta_x',
    [
        (-1.0, 200.0, 1000.0),  # Negative delta_fd
        (100.0, 0.0, 1000.0),  # Zero hf
        (100.0, -200.0, 1000.0),  # Negative hf
        (100.0, 200.0, 0.0),  # Zero delta_x
        (100.0, 200.0, -1000.0),  # Negative delta_x
        (-1.0, -200.0, -1000.0),  # All negative
    ],
)
def test_tao_Ed_flang_value_errors(delta_fd, hf, delta_x):
    """Test tao_Ed_flang raises ValueError for invalid arguments."""
    with pytest.raises(ValueError):
        _section_8_2_shear.tau_Ed_flange(delta_fd, hf, delta_x)


@pytest.mark.parametrize(
    'tau_ed, sf, hf, fyd, cot_theta_f, expected',
    [
        (3, 200.0, 200.0, 500.0, 1.0, 240),
        (2, 150.0, 150.0, 400.0, 1.5, 75.0),
    ],
)
def test_transverse_reinforcement_flange(
    tau_ed, sf, hf, fyd, cot_theta_f, expected
):
    """Test calculation of transverse reinforcement in the flange."""
    assert _section_8_2_shear.Asf_flange(
        tau_ed, sf, hf, fyd, cot_theta_f
    ) == pytest.approx(expected, rel=10e-2)


@pytest.mark.parametrize(
    'tau_ed, sf, hf, fyd, cot_theta_f',
    [
        (-1.0, 150.0, 200.0, 500.0, 1.0),  # Negative tau_ed
        (2.0, 0.0, 200.0, 500.0, 1.0),  # Zero sf
        (2.0, -150.0, 200.0, 500.0, 1.0),  # Negative sf
        (2.0, 150.0, 0.0, 500.0, 1.0),  # Zero hf
        (2.0, 150.0, -200.0, 500.0, 1.0),  # Negative hf
        (2.0, 150.0, 200.0, -500.0, 1.0),  # Negative fyd
        (2.0, 150.0, 200.0, 0, 1.0),  # Zero fyd
        (2.0, 150.0, 200.0, 500.0, 0.5),  # cot_theta_f < 1
    ],
)
def test_Asf_flang_value_errors(tau_ed, sf, hf, fyd, cot_theta_f):
    """Test Asf_flange raises ValueError for invalid arguments."""
    with pytest.raises(ValueError):
        _section_8_2_shear.Asf_flange(tau_ed, sf, hf, fyd, cot_theta_f)


@pytest.mark.parametrize(
    'tau_ed, cot_theta_f, fcd, nu, expected',
    [
        (5, 1.0, 30.0, 0.5, 10),  # cot(45°) = 1.0
        (3, 1.732, 25.0, 0.5, 6.9282),  # cot(30°) ≈ 1.732
    ],
)
def test_sigma_cd_flange(tau_ed, cot_theta_f, fcd, nu, expected):
    """Test check of compression field stress in the flange."""
    assert _section_8_2_shear.sigma_cd_flange(
        tau_ed, cot_theta_f, fcd, nu
    ) == pytest.approx(expected, rel=10e-2)


@pytest.mark.parametrize(
    'tau_ed, Ast_min, sf, hf, fyd, expected',
    [
        (
            1.0,
            100,
            200,
            200,
            500,
            True,
        ),  # tau_ed <= (100/(200*200))*500 = 1.25
        (
            2.0,
            100,
            200,
            200,
            500,
            False,
        ),  # tau_ed > (100/(200*200))*500 = 1.25
        (
            1.25,
            100,
            200,
            200,
            500,
            True,
        ),  # tau_ed == (100/(200*200))*500 = 1.25
    ],
)
def test_check_tau_Ed_flange_verification(
    tau_ed, Ast_min, sf, hf, fyd, expected
):
    """Test check_tau_Ed_flange_verification."""
    assert (
        _section_8_2_shear.check_tau_Ed_flange_verification(
            tau_ed, Ast_min, sf, hf, fyd
        )
        == expected
    )


@pytest.mark.parametrize(
    'tau_ed, Ast_min, sf, hf, fyd',
    [
        (-1.0, 100, 200, 200, 500),  # Negative tau_ed
        (1.0, -100, 200, 200, 500),  # Negative Ast_min
        (1.0, 100, 0, 200, 500),  # Zero sf
        (1.0, 100, -200, 200, 500),  # Negative sf
        (1.0, 100, 200, 0, 500),  # Zero hf
        (1.0, 100, 200, -200, 500),  # Negative hf
        (1.0, 100, 200, 200, 0),  # Zero fyd
        (1.0, 100, 200, 200, -500),  # Negative fyd
    ],
)
def test_check_tau_Ed_flange_verification_value_errors(
    tau_ed, Ast_min, sf, hf, fyd
):
    """Test check_tau_Ed_flange_verification raises ValueError."""
    with pytest.raises(ValueError):
        _section_8_2_shear.check_tau_Ed_flange_verification(
            tau_ed, Ast_min, sf, hf, fyd
        )


@pytest.mark.parametrize(
    'tau_ed, cot_theta_f, fcd, nu',
    [
        (-1.0, 1.0, 30.0, 0.5),  # Negative tau_ed
        (5.0, 0.5, 30.0, 0.5),  # cot_theta_f < 1
        (5.0, 1.0, -30.0, 0.5),  # Negative fcd
        (5.0, 1.0, 0, 0.5),  # Zero fcd
        (5.0, 1.0, 30.0, 0.0),  # nu zero
        (5.0, 1.0, 30.0, -0.1),  # nu negative
        (-1.0, 0.5, -30.0, -0.1),  # All invalid
    ],
)
def test_sigma_cd_flange_value_errors(tau_ed, cot_theta_f, fcd, nu):
    """Test sigma_cd_flange raises ValueError for invalid arguments."""
    with pytest.raises(ValueError):
        _section_8_2_shear.sigma_cd_flange(tau_ed, cot_theta_f, fcd, nu)


# Valid inputs
@pytest.mark.parametrize(
    'Ftd, Ast, Es, expected',
    [
        (100, 500, 200000, 0.001),
        (50, 250, 210000, 0.0009523809523809524),
        (200, 800, 205000, 0.0012195121951219512),
    ],
)
def test_eps_x_flang(Ftd, Ast, Es, expected):
    """Test eps_x_flang."""
    assert _section_8_2_shear.eps_x_flang(Ftd, Ast, Es) == pytest.approx(
        expected, rel=10e-4
    )


@pytest.mark.parametrize(
    'Ftd, Ast, Es',
    [
        (0, 500, 200000),  # Ftd zero
        (-100, 500, 200000),  # Ftd negative
        (100, 0, 200000),  # Ast zero
        (100, -500, 200000),  # Ast negative
        (100, 500, 0),  # Es zero
        (100, 500, -200000),  # Es negative
        (0, 0, 0),  # All zero
        (-1, -1, -1),  # All negative
    ],
)
def test_eps_x_flang_value_errors(Ftd, Ast, Es):
    """Test eps_x_flang raises ValueError for non-positive arguments."""
    with pytest.raises(ValueError):
        _section_8_2_shear.eps_x_flang(Ftd, Ast, Es)


@pytest.mark.parametrize(
    'VEdi, Ai, expected',
    [(1000, 200000, 5), (0, 100000, 0.0), (500, 250000, 2)],
)
def test_calculate_tau_edi(VEdi, Ai, expected):
    """Test the basic functionality of calculate_tau_edi."""
    assert _section_8_2_shear.tau_Edi(VEdi, Ai) == pytest.approx(
        expected, rel=10e-2
    )


@pytest.mark.parametrize(
    'VEdi, Ai',
    [
        (-10, 1000),  # Negative VEdi
        (100, -1000),  # Negative Ai
        (100, 0),  # Zero Ai
        (-5, -1000),  # Both negative
    ],
)
def test_tau_Edi_value_errors(VEdi, Ai):
    """Test tau_Edi raises ValueError for invalid arguments."""
    with pytest.raises(ValueError):
        _section_8_2_shear.tau_Edi(VEdi, Ai)


# Tests for calculate_tau_edi_composite
@pytest.mark.parametrize(
    'beta_new, VEd, z, bi, expected',
    [
        (0.5, 1000, 200, 1000, 2.5),
        (1.0, 1000, 200, 1000, 5.0),
    ],
)
def test_calculate_tau_edi_composite(beta_new, VEd, z, bi, expected):
    """Test the basic functionality of calculate_tau_edi_composite."""
    assert _section_8_2_shear.tau_Edi_composite(
        beta_new, VEd, z, bi
    ) == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'beta_new, VEd, z, bi',
    [
        (-0.1, 1000, 200, 1000),  # Negative beta_new
        (0.5, -1, 200, 1000),  # Negative VEd
        (0.5, 1000, 0, 1000),  # Zero z
        (0.5, 1000, 200, -1),  # Negative bi
        (0.5, 1000, 200, 0),  # Zero bi
    ],
)
def test_calculate_tau_edi_composite_errors(beta_new, VEd, z, bi):
    """Test that errors are raised for invalid inputs."""
    with pytest.raises(ValueError):
        _section_8_2_shear.tau_Edi_composite(beta_new, VEd, z, bi)


@pytest.mark.parametrize(
    'fck, sigma_n, Ai, Asi, fyd, alpha_deg, cv1, mu_v, gamma_c, expected',
    [
        (30, 5, 10000, 100, 500, 45, 0.1, 0.2, 1.5, 5.61),
        (35, 0, 20000, 200, 550, 90, 0.2, 0.25, 1.5, 2.1638),
        (40, 12, 15000, 150, 600, 135, 0.15, 0.3, 1.5, 1.2626),
        # alpha_deg = 35 (boundary, covers lines 2225-2228 and 2246)
        (30, 5, 10000, 100, 500, 35, 0.1, 0.2, 1.5, 6.0345),
        # alpha_deg = 135 (boundary, covers lines 2225-2228 and 2246)
        (30, 5, 10000, 100, 500, 135, 0.1, 0.2, 1.5, -1.4633),
    ],
)
def test_tau_rdi(
    fck, sigma_n, Ai, Asi, fyd, alpha_deg, cv1, mu_v, gamma_c, expected
):
    """Test the basic functionality and expected results of tau_Rdi."""
    # Calculate fcd from gamma_c (assuming eta_cc=1, k_tc=1 for test cases)
    fcd = fck / gamma_c
    assert math.isclose(
        _section_8_2_shear.tau_Rdi(
            fck, sigma_n, Ai, Asi, fyd, alpha_deg, cv1, mu_v, gamma_c, fcd
        ),
        expected,
        rel_tol=1e-3,
    )


@pytest.mark.parametrize(
    'fck, sigma_n, Ai, Asi, fyd, alpha_deg, cv1, mu_v, gamma_c, fcd',
    [
        (-30, 5, 10000, 100, 500, 45, 0.1, 0.2, 1.5, 20),  # Negative fck
        (0, 5, 10000, 100, 500, 45, 0.1, 0.2, 1.5, 0),  # Zero fck
        (30, 5, -10000, 100, 500, 45, 0.1, 0.2, 1.5, 20),  # Negative Ai
        (30, 5, 0, 100, 500, 45, 0.1, 0.2, 1.5, 20),  # Zero Ai
        (30, 5, 10000, -100, 500, 45, 0.1, 0.2, 1.5, 20),  # Negative Asi
        (30, 5, 10000, 100, -500, 45, 0.1, 0.2, 1.5, 20),  # Negative fyd
        (30, 5, 10000, 100, 0, 45, 0.1, 0.2, 1.5, 20),  # Zero fyd
        # alpha_deg validation check (line 2226)
        (30, 5, 10000, 100, 500, 20, 0.1, 0.2, 1.5, 20),  # alpha_deg < 35
        (
            30,
            5,
            10000,
            100,
            500,
            34,
            0.1,
            0.2,
            1.5,
            20,
        ),  # alpha_deg < 35 (edge)
        (30, 5, 10000, 100, 500, 140, 0.1, 0.2, 1.5, 20),  # alpha_deg > 135
        (
            30,
            5,
            10000,
            100,
            500,
            136,
            0.1,
            0.2,
            1.5,
            20,
        ),  # alpha_deg > 135 (edge)
        (30, 5, 10000, 100, 500, 45, 0.1, 0.2, -1.5, -20),  # Negative gamma_c
        (30, 5, 10000, 100, 500, 45, 0.1, 0.2, 0, 0),  # Zero gamma_c
        (30, 5, 10000, 100, 500, 45, 0.1, 0.2, 1.5, 0),  # Zero fcd
        (30, 5, 10000, 100, 500, 45, -0.1, 0.2, 1.5, 20),  # Negative cv1
        (30, 5, 10000, 100, 500, 45, 0.1, -0.2, 1.5, 20),  # Negative mu_v
        (-30, 5, -10000, -100, -500, 20, 0.1, 0.2, -1.5, 20),  # All negative
    ],
)
def test_tau_Rdi_value_errors(
    fck, sigma_n, Ai, Asi, fyd, alpha_deg, cv1, mu_v, gamma_c, fcd
):
    """Test tau_Rdi raises ValueError."""
    with pytest.raises(ValueError):
        _section_8_2_shear.tau_Rdi(
            fck, sigma_n, Ai, Asi, fyd, alpha_deg, cv1, mu_v, gamma_c, fcd
        )


# Test cv1 function
@pytest.mark.parametrize(
    'surface_roughness, tensile_stress, expected',
    [
        ('very smooth', False, 0.01),
        ('smooth', True, 0),  # Test with tensile stress
        ('rough', False, 0.15),
        ('very rough', False, 0.19),
        ('keyed', False, 0.37),
    ],
)
def test_cv1(surface_roughness, tensile_stress, expected):
    """Test cv1 coefficient."""
    assert (
        _section_8_2_shear.cv1(surface_roughness, tensile_stress) == expected
    )


def test_cv1_value_errors():
    """Test cv1 raises ValueError for unknown surface roughness."""
    with pytest.raises(ValueError):
        _section_8_2_shear.cv1('unknown', False)


# Test mu_v function
@pytest.mark.parametrize(
    'surface_roughness, expected',
    [
        ('very smooth', 0.5),
        ('smooth', 0.6),
        ('rough', 0.7),
        ('very rough', 0.9),
        ('keyed', 0.9),
    ],
)
def test_mu_v(surface_roughness, expected):
    """Test mu_v."""
    assert _section_8_2_shear.mu_v(surface_roughness) == expected


def test_mu_v_value_errors():
    """Test mu_v raises ValueError for unknown surface roughness."""
    with pytest.raises(ValueError):
        _section_8_2_shear.mu_v('unknown')


# Test cv2 function
@pytest.mark.parametrize(
    'surface_roughness, tensile_stress, expected',
    [
        ('very smooth', False, 0),
        ('smooth', True, 0),  # Test with tensile stress
        ('rough', False, 0.08),
        ('very rough', False, 0.15),
    ],
)
def test_cv2(surface_roughness, tensile_stress, expected):
    """Test cv2."""
    assert (
        _section_8_2_shear.cv2(surface_roughness, tensile_stress) == expected
    )


def test_cv2_value_errors():
    """Test cv2 raises ValueError for unknown surface roughness."""
    with pytest.raises(ValueError):
        _section_8_2_shear.cv2('unknown', False)


# Test kv function
@pytest.mark.parametrize(
    'surface_roughness, expected',
    [
        ('very smooth', 0),
        ('smooth', 0.5),
        ('rough', 0.5),
        ('very rough', 0.5),
    ],
)
def test_kv(surface_roughness, expected):
    """Test kv."""
    assert _section_8_2_shear.kv(surface_roughness) == expected


def test_kv_value_errors():
    """Test kv raises ValueError for unknown surface roughness."""
    with pytest.raises(ValueError):
        _section_8_2_shear.kv('unknown')


# Test kdowel function
@pytest.mark.parametrize(
    'surface_roughness, expected',
    [
        ('very smooth', 1.5),
        ('smooth', 1.1),
        ('rough', 0.9),
        ('very rough', 0.9),
    ],
)
def test_kdowel(surface_roughness, expected):
    """Test kdowel."""
    assert _section_8_2_shear.kdowel(surface_roughness) == expected


def test_kdowel_value_errors():
    """Test kdowel raises ValueError for unknown surface roughness."""
    with pytest.raises(ValueError):
        _section_8_2_shear.kdowel('unknown')


@pytest.mark.parametrize(
    'cv2, fck, gamma_c, mu_v, sigma_n, kv, rho_i, fyd, kdowel, expected',
    [
        (0.1, 30, 1.5, 0.2, 0.3, 0.1, 0.02, 500, 0.0, 0.6251),
    ],
)
def test_shear_stress_resistance(
    cv2, fck, gamma_c, mu_v, sigma_n, kv, rho_i, fyd, kdowel, expected
):
    """Test the shear stress resistance calculation."""
    # Calculate fcd from gamma_c (assuming eta_cc=1, k_tc=1 for test cases)
    fcd = fck / gamma_c
    result = _section_8_2_shear.tau_Rdi_no_yielding(
        cv2, fck, gamma_c, mu_v, sigma_n, kv, rho_i, fyd, kdowel, fcd
    )
    assert result == pytest.approx(expected, rel=10e-3)


@pytest.mark.parametrize(
    'cv2, fck, gamma_c, mu_v, sigma_n, kv, rho_i, fyd, kdowel, fcd',
    [
        (-0.1, 30, 1.5, 0.2, 0.3, 0.1, 0.02, 500, 0.0, 20),  # Negative cv2
        (0.1, -30, 1.5, 0.2, 0.3, 0.1, 0.02, 500, 0.0, 20),  # Negative fck
        (0.1, 0, 1.5, 0.2, 0.3, 0.1, 0.02, 500, 0.0, 0),  # Zero fck
        (
            0.1,
            30,
            -1.5,
            0.2,
            0.3,
            0.1,
            0.02,
            500,
            0.0,
            -20,
        ),  # Negative gamma_c
        (0.1, 30, 0, 0.2, 0.3, 0.1, 0.02, 500, 0.0, 0),  # Zero gamma_c
        (0.1, 30, 1.5, 0.2, 0.3, 0.1, 0.02, 500, 0.0, 0),  # Zero fcd
        # fyd validation check (line 2473)
        (
            0.1,
            30,
            1.5,
            0.2,
            0.3,
            0.1,
            0.02,
            -500,
            0.0,
            20,
        ),  # Negative fyd
        (0.1, 30, 1.5, 0.2, 0.3, 0.1, 0.02, 0, 0.0, 20),  # Zero fyd
        (0.1, 30, 1.5, -0.2, 0.3, 0.1, 0.02, 500, 0.0, 20),  # Negative mu_v
        # kv validation check (line 2479)
        (
            0.1,
            30,
            1.5,
            0.2,
            0.3,
            -0.1,
            0.02,
            500,
            0.0,
            20,
        ),  # Negative kv
        # rho_i validation check (line 2481)
        (
            0.1,
            30,
            1.5,
            0.2,
            0.3,
            0.1,
            -0.02,
            500,
            0.0,
            20,
        ),  # Negative rho_i
        # kdowel validation check (line 2483)
        (
            0.1,
            30,
            1.5,
            0.2,
            0.3,
            0.1,
            0.02,
            500,
            -0.1,
            20,
        ),  # Negative kdowel
        (
            0.1,
            -30,
            -1.5,
            0.2,
            0.3,
            -0.1,
            -0.02,
            -500,
            -0.1,
            20,
        ),  # All negative
    ],
)
def test_tau_Rdi_no_yielding_value_errors(
    cv2, fck, gamma_c, mu_v, sigma_n, kv, rho_i, fyd, kdowel, fcd
):
    """Test tau_Rdi_no_yielding raises ValueError."""
    with pytest.raises(ValueError):
        _section_8_2_shear.tau_Rdi_no_yielding(
            cv2, fck, gamma_c, mu_v, sigma_n, kv, rho_i, fyd, kdowel, fcd
        )


@pytest.mark.parametrize(
    'tmin, fctm, fyk, expected',
    [
        (200, 2.9, 500, 1.16),
    ],
)
def test_min_interface_reinforcement(tmin, fctm, fyk, expected):
    """Test the minimum interface reinforcement calculation."""
    result = _section_8_2_shear.as_min(tmin, fctm, fyk)
    assert result == pytest.approx(expected, rel=10e-3)


@pytest.mark.parametrize(
    'tmin, fctm, fyk',
    [
        (0, 2.9, 500),  # tmin zero
        (-1, 2.9, 500),  # tmin negative
        (200, 0, 500),  # fctm zero
        (200, -2.9, 500),  # fctm negative
        (200, 2.9, 0),  # fyk zero
        (200, 2.9, -500),  # fyk negative
        (0, 0, 0),  # all zero
        (-1, -2.9, -500),  # all negative
    ],
)
def test_as_min_value_errors(tmin, fctm, fyk):
    """Test as_min raises ValueError for non-positive arguments."""
    with pytest.raises(ValueError):
        _section_8_2_shear.as_min(tmin, fctm, fyk)
