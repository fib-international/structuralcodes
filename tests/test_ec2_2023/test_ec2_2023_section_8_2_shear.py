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
        (20, 0, 0.005, 0.008),  # vEd_x zero (division by zero)
    ],
)
def test_rho_l_planar_value_errors(vEd_y, vEd_x, rho_l_x, rho_l_y):
    """Test that rho_l_planar raises ValueError for invalid inputs."""
    with pytest.raises(ValueError):
        _section_8_2_shear.rho_l_planar(vEd_y, vEd_x, rho_l_x, rho_l_y)


# Tests using pytest
def test_cot_theta_min():
    """Tests the function cot_theta_min with various scenarios."""
    assert _section_8_2_shear.cot_theta_min(0, 100, 10, 400) == 2.5
    assert _section_8_2_shear.cot_theta_min(0, 100, 10, 400) == 2.5
    assert _section_8_2_shear.cot_theta_min(10, 100, 10, 400) == 2.49
    assert _section_8_2_shear.cot_theta_min(-100, 100, 50, 400) == 3.0


@pytest.mark.parametrize(
    'NEd, VEd, x, d',
    [
        (0, 100, -1, 400),  # x negative
        (0, 100, 10, 0),  # d zero
        (0, 100, -1, 0),  # x negative and d zero
        (0, 100, 10, -400),  # d negative
        (0, 100, -10, -400),  # both negative
        (-10, 100, 10, 0),  # d zero with NEd negative
        (10, 100, 10, 0),  # d zero with NEd positive
    ],
)
def test_cot_theta_min_value_error_all_cases(NEd, VEd, x, d):
    """Test cot_theta_min raises ValueError."""
    with pytest.raises(ValueError):
        _section_8_2_shear.cot_theta_min(NEd, VEd, x, d)


@pytest.mark.parametrize(
    'NEd, VEd, x, d',
    [
        (0, 100, 10, 0),  # d is zero
        (0, 100, -10, 400),  # x is negative
        (0, 100, -10, 0),  # x negative and d zero
        (0, 100, 10, -400),  # d negative
        (0, 100, -10, -400),  # both negative
    ],
)
def test_cot_theta_min_value_error(NEd, VEd, x, d):
    """Test cot_theta_min raises ValueError for invalid dimensions."""
    with pytest.raises(ValueError):
        _section_8_2_shear.cot_theta_min(NEd, VEd, x, d)


def test_tau_Rd_sy():
    """Tests the function tau_Rd_sy."""
    assert _section_8_2_shear.tau_Rd_sy(0.01, 500, 2) == 10
    assert _section_8_2_shear.tau_Rd_sy(0.02, 400, 2.5) == 20


@pytest.mark.parametrize(
    'rho_w, fywd, cot_theta',
    [
        (-0.01, 500, 2),  # Negative rho_w
        (0.01, -500, 2),  # Negative fywd
        (-0.01, -500, 2),  # Both negative
    ],
)
def test_tau_Rd_sy_value_errors(rho_w, fywd, cot_theta):
    """Test tau_Rd_sy raises ValueError for negative rho_w or fywd."""
    with pytest.raises(ValueError):
        _section_8_2_shear.tau_Rd_sy(rho_w, fywd, cot_theta)


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
    assert _section_8_2_shear.sigma_cd(1, 2, 0.5, 0.5, 20) == 2.5


@pytest.mark.parametrize(
    'tau_Ed, cot_theta, tan_theta, nu, f_cd',
    [
        (-1, 1, 1, 0.5, 20),  # Negative tau_Ed
        (1, -1, 1, 0.5, 20),  # Negative cot_theta
        (1, 1, -1, 0.5, 20),  # Negative tan_theta
        (1, 1, 1, -0.5, 20),  # Negative nu
        (1, 1, 1, 0.5, -20),  # Negative f_cd
        (-1, -1, -1, -0.5, -20),  # All negative
    ],
)
def test_sigma_cd_value_errors(tau_Ed, cot_theta, tan_theta, nu, f_cd):
    """Test sigma_cd raises ValueError for negative arguments."""
    with pytest.raises(ValueError):
        _section_8_2_shear.sigma_cd(tau_Ed, cot_theta, tan_theta, nu, f_cd)


def test_tau_Rd():
    """Tests the function tau_Rd."""
    assert _section_8_2_shear.tau_Rd(0.01, 500, 2, 0.5, 50) == 10
    assert (
        _section_8_2_shear.tau_Rd(0.01, 500, 3, 0.5, 20) == 5
    )  # Limited by the compression field


@pytest.mark.parametrize(
    'rho_w, fywd, cot_theta, nu, f_cd',
    [
        (-0.01, 500, 2, 0.5, 50),  # Negative rho_w
        (0.01, -500, 2, 0.5, 50),  # Negative fywd
        (0.01, 500, -2, 0.5, 50),  # Negative cot_theta
        (0.01, 500, 2, -0.5, 50),  # Negative nu
        (0.01, 500, 2, 0.5, -50),  # Negative f_cd
        (-0.01, -500, -2, -0.5, -50),  # All negative
    ],
)
def test_tau_Rd_value_errors(rho_w, fywd, cot_theta, nu, f_cd):
    """Test tau_Rd raises ValueError for negative arguments."""
    with pytest.raises(ValueError):
        _section_8_2_shear.tau_Rd(rho_w, fywd, cot_theta, nu, f_cd)


def test_cot_theta():
    """Tests the function cot_theta."""
    assert _section_8_2_shear.cot_theta(0.5, 20, 0.01, 500, 2) == 1
    assert _section_8_2_shear.cot_theta(0.5, 40, 0.01, 500, 2) == 2


@pytest.mark.parametrize(
    'nu, f_cd, rho_w, fywd, cot_theta_min',
    [
        (-0.1, 20, 0.01, 500, 2),  # Negative nu
        (0.5, -20, 0.01, 500, 2),  # Negative f_cd
        (0.5, 20, -0.01, 500, 2),  # Negative rho_w
        (0.5, 20, 0.01, -500, 2),  # Negative fywd
        (-0.5, -20, -0.01, -500, 2),  # All negative
    ],
)
def test_cot_theta_value_errors(nu, f_cd, rho_w, fywd, cot_theta_min):
    """Test cot_theta raises ValueError for negative arguments."""
    with pytest.raises(ValueError):
        _section_8_2_shear.cot_theta(nu, f_cd, rho_w, fywd, cot_theta_min)


@pytest.mark.parametrize(
    'Ftd, Est, Ast, expected', [(500, 210000, 1000, 0.002380952380952381)]
)
def test_epsilon_xt(Ftd, Est, Ast, expected):
    """Test εxt calculation."""
    assert _section_8_2_shear.epsilon_xt(Ftd, Est, Ast) == pytest.approx(
        expected
    )


@pytest.mark.parametrize(
    'Ftd, Est, Ast',
    [
        (500, -210000, 1000),  # Negative Est
        (500, 210000, -1000),  # Negative Ast
        (500, -210000, -1000),  # Both negative
    ],
)
def test_epsilon_xt_value_errors(Ftd, Est, Ast):
    """Test epsilon_xt raises ValueError for negative Est or Ast."""
    with pytest.raises(ValueError):
        _section_8_2_shear.epsilon_xt(Ftd, Est, Ast)


@pytest.mark.parametrize(
    'Fcd, Ecc, Acc, expected', [(500, 30000, 1000, 0.016666666666666666)]
)
def test_epsilon_xc_compression(Fcd, Ecc, Acc, expected):
    """Test εxc calculation for compression."""
    assert _section_8_2_shear.epsilon_xc_comp(Fcd, Ecc, Acc) == pytest.approx(
        expected, rel=10e-3
    )


@pytest.mark.parametrize(
    'Fcd, Ecc, Acc',
    [
        (500, -30000, 1000),  # Negative Ecc
        (500, 30000, -1000),  # Negative Acc
        (500, -30000, -1000),  # Both negative
    ],
)
def test_epsilon_xc_comp_value_errors(Fcd, Ecc, Acc):
    """Test epsilon_xc_comp raises ValueError for negative Ecc or Acc."""
    with pytest.raises(ValueError):
        _section_8_2_shear.epsilon_xc_comp(Fcd, Ecc, Acc)


@pytest.mark.parametrize(
    'Fcd, Esc, Asc, expected', [(500, 210000, 1000, 0.002380952380952381)]
)
def test_epsilon_xc_tension(Fcd, Esc, Asc, expected):
    """Test εxc calculation for tension."""
    assert _section_8_2_shear.epsilon_xc_tens(Fcd, Esc, Asc) == pytest.approx(
        expected, rel=10e-3
    )


@pytest.mark.parametrize(
    'Fcd, Esc, Asc',
    [
        (500, -210000, 1000),  # Negative Esc
        (500, 210000, -1000),  # Negative Asc
        (500, -210000, -1000),  # Both negative
    ],
)
def test_epsilon_xc_tens_value_errors(Fcd, Esc, Asc):
    """Test epsilon_xc_tens raises ValueError for negative Esc or Asc."""
    with pytest.raises(ValueError):
        _section_8_2_shear.epsilon_xc_tens(Fcd, Esc, Asc)


@pytest.mark.parametrize(
    'epsilon_xt, epsilon_xc, expected', [(0.002, 0.001, 0.0015)]
)
def test_epsilon_x(epsilon_xt, epsilon_xc, expected):
    """Test εx calculation."""
    assert _section_8_2_shear.epsilon_x(
        epsilon_xt, epsilon_xc
    ) == pytest.approx(expected)


@pytest.mark.parametrize('epsilon_x, theta, expected', [(0.001, 2, 0.5025)])
def test_nu(epsilon_x, theta, expected):
    """Test ν calculation."""
    assert _section_8_2_shear.nu(epsilon_x, theta) == pytest.approx(
        expected, rel=10e-3
    )


@pytest.mark.parametrize(
    'epsilon_x, cot_theta',
    [
        (-0.001, 2),  # Negative epsilon_x
        (-1.0, 0),  # Negative epsilon_x, zero cot_theta
        (-0.5, -2),  # Negative epsilon_x, negative cot_theta
    ],
)
def test_nu_value_errors(epsilon_x, cot_theta):
    """Test nu raises ValueError for negative epsilon_x."""
    with pytest.raises(ValueError):
        _section_8_2_shear.nu(epsilon_x, cot_theta)


@pytest.mark.parametrize(
    'VEd, theta, expected',
    [
        (100, 2, 200),
        (150, 1.5, 225),
    ],
)
def test_calculate_nv(VEd, theta, expected):
    """Test calculate_nv function with various inputs."""
    assert _section_8_2_shear.Nvd(VEd, theta) == pytest.approx(
        expected, rel=1e-6
    )


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
    ],
)
def test_k_duct_negative_dimensions(
    duct_material, is_grouted, wall_thickness, duct_diameter
):
    """Test k_duct raises ValueError."""
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
    'nu, f_cd, theta, beta_incl, rho_w, f_ywd, expected',
    [
        (0.6, 30, 1, 1.5, 0.01, 500, 3),
        (0.5, 40, 1.5, 1, 0.02, 400, 9.2307),
    ],
)
def test_tau_rd(nu, f_cd, theta, beta_incl, rho_w, f_ywd, expected):
    """Test calculation of enhanced shear stress resistance τRd."""
    assert math.isclose(
        _section_8_2_shear.tau_rd(nu, f_cd, theta, beta_incl, rho_w, f_ywd),
        expected,
        rel_tol=1e-2,
    )


@pytest.mark.parametrize(
    'nu, f_cd, cot_theta, cot_beta_incl, rho_w, f_ywd',
    [
        (0.6, -30, 1, 1.5, 0.01, 500),  # Negative f_cd
        (0.6, 30, 1, 1.5, 0.01, -500),  # Negative f_ywd
        (0.6, -30, 1, 1.5, 0.01, -500),  # Both negative
    ],
)
def test_tau_rd_value_errors(nu, f_cd, cot_theta, cot_beta_incl, rho_w, f_ywd):
    """Test tau_rd raises ValueError for negative f_cd or f_ywd."""
    with pytest.raises(ValueError):
        _section_8_2_shear.tau_rd(
            nu, f_cd, cot_theta, cot_beta_incl, rho_w, f_ywd
        )


@pytest.mark.parametrize(
    'e_s, eps_x, f_ywd, cot_theta, expected',
    [
        (200000, 0.001, 500, 2, 500),
        (210000, 0.002, 450, 1, 420),
    ],
)
def test_sigma_swd(e_s, eps_x, f_ywd, cot_theta, expected):
    """Test calculation of stress σswd in shear reinforcement."""
    assert math.isclose(
        _section_8_2_shear.sigma_swd(e_s, eps_x, f_ywd, cot_theta),
        expected,
        rel_tol=1e-2,
    )


@pytest.mark.parametrize(
    'Es, eps_x, f_ywd, cot_theta',
    [
        (-200000, 0.001, 500, 2),  # Negative Es
        (200000, 0.001, -500, 2),  # Negative f_ywd
        (-200000, 0.001, -500, 2),  # Both negative
    ],
)
def test_sigma_swd_value_errors(Es, eps_x, f_ywd, cot_theta):
    """Test sigma_swd raises ValueError for negative Es or f_ywd."""
    with pytest.raises(ValueError):
        _section_8_2_shear.sigma_swd(Es, eps_x, f_ywd, cot_theta)


@pytest.mark.parametrize(
    'tau_ed, rho_w, f_ywd, theta, z, b_w, a, x, expected',
    [
        (0.5, 0.01, 500, 1, 300, 200, 500, 250, 0),
        (0.4, 0.02, 400, 0.7, 350, 250, 600, 200, -45.5),
    ],
)
def test_delta_m_ed(tau_ed, rho_w, f_ywd, theta, z, b_w, a, x, expected):
    """Test calculation of additional moment ΔMEd."""
    assert math.isclose(
        _section_8_2_shear.delta_MEd(
            tau_ed, rho_w, f_ywd, theta, z, b_w, a, x
        ),
        expected,
        rel_tol=1e-2,
    )


@pytest.mark.parametrize(
    'tau_ed, rho_w, f_ywd, cot_theta, z, b_w, a, x',
    [
        (-0.5, 0.01, 500, 1, 300, 200, 500, 250),  # Negative tau_ed
        (0.5, 0.01, -500, 1, 300, 200, 500, 250),  # Negative f_ywd
        (0.5, 0.01, 500, 1, -300, 200, 500, 250),  # Negative z
        (0.5, 0.01, 500, 1, 300, -200, 500, 250),  # Negative b_w
        (0.5, 0.01, 500, 1, 300, 200, -500, 250),  # Negative a
        (0.5, 0.01, 500, 1, 300, 200, 500, -250),  # Negative x
        (-0.5, 0.01, -500, 1, -300, -200, -500, -250),  # All negative
    ],
)
def test_delta_MEd_value_errors(tau_ed, rho_w, f_ywd, cot_theta, z, b_w, a, x):
    """Test delta_MEd raises ValueError for negative arguments."""
    with pytest.raises(ValueError):
        _section_8_2_shear.delta_MEd(
            tau_ed, rho_w, f_ywd, cot_theta, z, b_w, a, x
        )


@pytest.mark.parametrize(
    'rho_w, f_ywd, theta, alpha_w, cot_theta_min, expected',
    [
        (0.01, 500, 1, 45, 1.5, 7.071),
        (0.02, 450, 2, 60, 3, 20.088),
    ],
)
def test_tau_rd_sy(rho_w, f_ywd, theta, alpha_w, cot_theta_min, expected):
    """Test calculation of shear stress resistance τRd,sy."""
    assert math.isclose(
        _section_8_2_shear.tau_rd_sy(
            rho_w, f_ywd, theta, alpha_w, cot_theta_min
        ),
        expected,
        rel_tol=1e-2,
    )


@pytest.mark.parametrize(
    'rho_w, f_ywd, cot_theta, alpha_w, cot_theta_min',
    [
        (0.01, -500, 2, 60, 3),  # Negative f_ywd
        (0.01, 500, 2, 44, 3),  # alpha_w < 45
        (0.01, 500, 2, 90, 3),  # alpha_w >= 90
        (0.01, 500, 2, 100, 3),  # alpha_w > 90
    ],
)
def test_tau_rd_sy_value_errors(
    rho_w, f_ywd, cot_theta, alpha_w, cot_theta_min
):
    """Test tau_rd_sy raises ValueError for invalid f_ywd or alpha_w."""
    with pytest.raises(ValueError):
        _section_8_2_shear.tau_rd_sy(
            rho_w, f_ywd, cot_theta, alpha_w, cot_theta_min
        )


@pytest.mark.parametrize(
    'tau_ed, theta, alpha_w, nu, f_cd, cot_theta_min, expected',
    [
        (0.5, 1, 45, 0.6, 30, 2, 0.5),
        (0.7759, 2, 60, 0.5, 40, 2, 1.50522),
    ],
)
def test_sigma_cd_s(tau_ed, theta, alpha_w, nu, f_cd, cot_theta_min, expected):
    """Test calculation of compression stress σcd."""
    assert math.isclose(
        _section_8_2_shear.sigma_cd_s(
            tau_ed, theta, alpha_w, nu, f_cd, cot_theta_min
        ),
        expected,
        rel_tol=1e-2,
    )


@pytest.mark.parametrize(
    'tau_ed, cot_theta, alpha_w, nu, f_cd, cot_theta_min',
    [
        (-0.5, 1, 60, 0.6, 30, 2),  # Negative tau_ed
        (0.5, 1, 60, 0.6, -30, 2),  # Negative f_cd
        (0.5, 1, 44, 0.6, 30, 2),  # alpha_w < 45
        (0.5, 1, 90, 0.6, 30, 2),  # alpha_w >= 90
        (0.5, 1, 100, 0.6, 30, 2),  # alpha_w > 90
        (-0.5, 1, 44, 0.6, -30, 2),  # All invalid
    ],
)
def test_sigma_cd_s_value_errors(
    tau_ed, cot_theta, alpha_w, nu, f_cd, cot_theta_min
):
    """Test sigma_cd_s raises ValueError for invalid arguments."""
    with pytest.raises(ValueError):
        _section_8_2_shear.sigma_cd_s(
            tau_ed, cot_theta, alpha_w, nu, f_cd, cot_theta_min
        )


@pytest.mark.parametrize(
    'v_ed, theta, alpha_w, cot_theta_min, expected',
    [
        (100, 0.3, 45, 1, -58.5786),
        (80, 0.1, 50, 3, -29.8233),
    ],
)
def test_n_vd(v_ed, theta, alpha_w, cot_theta_min, expected):
    """Test calculation of axial tensile force NVd."""
    assert math.isclose(
        _section_8_2_shear.NVds(v_ed, theta, alpha_w, cot_theta_min),
        expected,
        rel_tol=1e-2,
    )


@pytest.mark.parametrize(
    'VEd, cot_theta, alpha_w, cot_theta_min',
    [
        (-100, 1, 60, 2),  # Negative VEd
        (100, 1, 44, 2),  # alpha_w < 45
        (100, 1, 90, 2),  # alpha_w >= 90
        (100, 1, 100, 2),  # alpha_w > 90
    ],
)
def test_NVds_value_errors(VEd, cot_theta, alpha_w, cot_theta_min):
    """Test NVds raises ValueError for invalid arguments."""
    with pytest.raises(ValueError):
        _section_8_2_shear.NVds(VEd, cot_theta, alpha_w, cot_theta_min)


@pytest.mark.parametrize(
    'nu, f_cd, theta, beta_incl, rho_w, f_ywd, '
    + 'alpha_w, cot_theta_min, expected',
    [
        (0.6, 30, 1, 3, 0.01, 500, 45, 2, -3.8578),
        (0.5, 40, 0.7, 1, 0.02, 400, 60, 5, 6.9013),
    ],
)
def test_tau_rd_incl(
    nu, f_cd, theta, beta_incl, rho_w, f_ywd, alpha_w, cot_theta_min, expected
):
    """Test calculation of shear stress resistance τRd."""
    assert math.isclose(
        _section_8_2_shear.tau_rd_incl(
            nu, f_cd, theta, beta_incl, rho_w, f_ywd, alpha_w, cot_theta_min
        ),
        expected,
        rel_tol=1e-2,
    )


@pytest.mark.parametrize(
    'nu, f_cd, cot_theta, cot_beta_incl, rho_w, f_ywd, alpha_w, cot_theta_min',
    [
        (0.6, -30, 1, 3, 0.01, 500, 45, 2),  # Negative f_cd
        (0.6, 30, 1, 3, 0.01, -500, 45, 2),  # Negative f_ywd
        (0.6, 30, 1, 3, 0.01, 500, 44, 2),  # alpha_w < 45
        (0.6, 30, 1, 3, 0.01, 500, 90, 2),  # alpha_w >= 90
        (0.6, 30, 1, 3, 0.01, 500, 100, 2),  # alpha_w > 90
    ],
)
def test_tau_rd_incl_value_errors(
    nu, f_cd, cot_theta, cot_beta_incl, rho_w, f_ywd, alpha_w, cot_theta_min
):
    """Test tau_rd_incl raises ValueError for invalid arguments."""
    with pytest.raises(ValueError):
        _section_8_2_shear.tau_rd_incl(
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
    'e_s, eps_x, theta, alpha_w, f_ywd, expected',
    [
        (200000, 0.001, 1, 45, 500, 500),  # Example values
        (210000, 0.002, 2, 60, 450, 450),  # Example values
    ],
)
def test_sigma_swd_v2(e_s, eps_x, theta, alpha_w, f_ywd, expected):
    """Test calculation of stress σswd."""
    assert math.isclose(
        _section_8_2_shear.sigma_swd_v2(e_s, eps_x, theta, alpha_w, f_ywd),
        expected,
        rel_tol=1e-2,
    )


@pytest.mark.parametrize(
    'Es, eps_x, cot_theta, alpha_w, f_ywd',
    [
        (-200000, 0.001, 2, 45, 500),  # Negative Es
        (200000, 0.001, 2, 45, -500),  # Negative f_ywd
        (-200000, 0.001, 2, 45, -500),  # Both negative
    ],
)
def test_sigma_swd_v2_value_errors(Es, eps_x, cot_theta, alpha_w, f_ywd):
    """Test sigma_swd_v2 raises ValueError for negative Es or f_ywd."""
    with pytest.raises(ValueError):
        _section_8_2_shear.sigma_swd_v2(Es, eps_x, cot_theta, alpha_w, f_ywd)


@pytest.mark.parametrize(
    'tau_rd, m_ed, m_rd, expected',
    [
        (5.0, 2.0, 10.0, 4.0),
        (10.0, 3.0, 15.0, 8.0),
    ],
)
def test_shear_stress_resistance_reduced(tau_rd, m_ed, m_rd, expected):
    """Test calculation of reduced shear stress resistance."""
    assert _section_8_2_shear.tao_Rd_m(tau_rd, m_ed, m_rd) == pytest.approx(
        expected, rel=10e-2
    )


@pytest.mark.parametrize(
    'tau_rd, m_ed, m_rd',
    [
        (-1.0, 2.0, 10.0),  # Negative tau_rd
        (5.0, -2.0, 10.0),  # Negative m_ed
        (5.0, 2.0, -10.0),  # Negative m_rd
        (-1.0, -2.0, -10.0),  # All negative
    ],
)
def test_tao_Rd_m_value_errors(tau_rd, m_ed, m_rd):
    """Test tao_Rd_m raises ValueError for negative arguments."""
    with pytest.raises(ValueError):
        _section_8_2_shear.tao_Rd_m(tau_rd, m_ed, m_rd)


@pytest.mark.parametrize(
    'delta_fd, hf, delta_x, expected',
    [
        (100.0, 200.0, 1000.0, 0.5),
        (50.0, 150.0, 500.0, 0.666),
    ],
)
def test_longitudinal_shear_stress(delta_fd, hf, delta_x, expected):
    """Test calculation of longitudinal shear stress."""
    assert _section_8_2_shear.tao_Ed_flang(
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
        _section_8_2_shear.tao_Ed_flang(delta_fd, hf, delta_x)


@pytest.mark.parametrize(
    'tau_ed, sf, hf, fyd, expected',
    [
        (2.0, 150.0, 200.0, 500.0, 120.0),
        (1.5, 100.0, 100.0, 400.0, 37.5),
        (0.0, 200.0, 300.0, 600.0, 0.0),  # zero shear stress
    ],
)
def test_Ast_min_flang_valid(tau_ed, sf, hf, fyd, expected):
    """Test Ast_min_flang with valid inputs."""
    assert _section_8_2_shear.Ast_min_flang(
        tau_ed, sf, hf, fyd
    ) == pytest.approx(expected, rel=1e-6)


@pytest.mark.parametrize(
    'tau_ed, sf, hf, fyd',
    [
        (-1.0, 150.0, 200.0, 500.0),  # Negative tau_ed
        (2.0, 0.0, 200.0, 500.0),  # Zero sf
        (2.0, -150.0, 200.0, 500.0),  # Negative sf
        (2.0, 150.0, 0.0, 500.0),  # Zero hf
        (2.0, 150.0, -200.0, 500.0),  # Negative hf
        (2.0, 150.0, 200.0, -500.0),  # Negative fyd
    ],
)
def test_Ast_min_flang_value_errors(tau_ed, sf, hf, fyd):
    """Test Ast_min_flang raises ValueError for invalid arguments."""
    with pytest.raises(ValueError):
        _section_8_2_shear.Ast_min_flang(tau_ed, sf, hf, fyd)


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
    assert _section_8_2_shear.Asf_flang(
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
    ],
)
def test_Asf_flang_value_errors(tau_ed, sf, hf, fyd, cot_theta_f):
    """Test Asf_flang raises ValueError for invalid arguments."""
    with pytest.raises(ValueError):
        _section_8_2_shear.Asf_flang(tau_ed, sf, hf, fyd, cot_theta_f)


@pytest.mark.parametrize(
    'tau_ed, theta_f, fcd, nu, expected',
    [
        (5, 45.0, 30.0, 0.5, 10),
        (3, 30.0, 25.0, 0.5, 6.9282),
    ],
)
def test_sigma_cd_flang(tau_ed, theta_f, fcd, nu, expected):
    """Test check of compression field stress in the flange."""
    assert _section_8_2_shear.sigma_cd_flang(
        tau_ed, theta_f, fcd, nu
    ) == pytest.approx(expected, rel=10e-2)


@pytest.mark.parametrize(
    'tau_ed, theta_f, fcd, nu',
    [
        (-1.0, 45.0, 30.0, 0.5),  # Negative tau_ed
        (5.0, -45.0, 30.0, 0.5),  # Negative theta_f
        (5.0, 45.0, -30.0, 0.5),  # Negative fcd
        (5.0, 45.0, 30.0, 0.0),  # nu zero
        (5.0, 45.0, 30.0, -0.1),  # nu negative
        (-1.0, -45.0, -30.0, -0.1),  # All negative
    ],
)
def test_sigma_cd_flang_value_errors(tau_ed, theta_f, fcd, nu):
    """Test sigma_cd_flang raises ValueError for invalid arguments."""
    with pytest.raises(ValueError):
        _section_8_2_shear.sigma_cd_flang(tau_ed, theta_f, fcd, nu)


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
        (-5, -1000),  # Both negative
    ],
)
def test_tau_Edi_value_errors(VEdi, Ai):
    """Test tau_Edi raises ValueError for negative VEdi or Ai."""
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
        (-0.1, 1000, 200, 1000),
        (0.5, -1, 200, 1000),
        (0.5, 1000, 0, 1000),
        (0.5, 1000, 200, -1),
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
    ],
)
def test_tau_rdi(
    fck, sigma_n, Ai, Asi, fyd, alpha_deg, cv1, mu_v, gamma_c, expected
):
    """Test the basic functionality and expected results of tau_Rdi."""
    assert math.isclose(
        _section_8_2_shear.tau_Rdi(
            fck, sigma_n, Ai, Asi, fyd, alpha_deg, cv1, mu_v, gamma_c
        ),
        expected,
        rel_tol=1e-3,
    )


@pytest.mark.parametrize(
    'fck, sigma_n, Ai, Asi, fyd, alpha_deg, cv1, mu_v, gamma_c',
    [
        (-30, 5, 10000, 100, 500, 45, 0.1, 0.2, 1.5),  # Negative fck
        (30, 5, -10000, 100, 500, 45, 0.1, 0.2, 1.5),  # Negative Ai
        (30, 5, 10000, -100, 500, 45, 0.1, 0.2, 1.5),  # Negative Asi
        (30, 5, 10000, 100, -500, 45, 0.1, 0.2, 1.5),  # Negative fyd
        (30, 5, 10000, 100, 500, 20, 0.1, 0.2, 1.5),  # alpha_deg < 35
        (30, 5, 10000, 100, 500, 140, 0.1, 0.2, 1.5),  # alpha_deg > 135
        (30, 5, 10000, 100, 500, 45, 0.1, 0.2, -1.5),  # Negative gamma_c
        (-30, 5, -10000, -100, -500, 20, 0.1, 0.2, -1.5),  # All negative
    ],
)
def test_tau_Rdi_value_errors(
    fck, sigma_n, Ai, Asi, fyd, alpha_deg, cv1, mu_v, gamma_c
):
    """Test tau_Rdi raises ValueError."""
    with pytest.raises(ValueError):
        _section_8_2_shear.tau_Rdi(
            fck, sigma_n, Ai, Asi, fyd, alpha_deg, cv1, mu_v, gamma_c
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
    result = _section_8_2_shear.tau_Rdi_ny(
        cv2, fck, gamma_c, mu_v, sigma_n, kv, rho_i, fyd, kdowel
    )
    assert result == pytest.approx(expected, rel=10e-3)


@pytest.mark.parametrize(
    'cv2, fck, gamma_c, mu_v, sigma_n, kv, rho_i, fyd, kdowel',
    [
        (0.1, -30, 1.5, 0.2, 0.3, 0.1, 0.02, 500, 0.0),  # Negative fck
        (0.1, 30, -1.5, 0.2, 0.3, 0.1, 0.02, 500, 0.0),  # Negative gamma_c
        (0.1, 30, 1.5, 0.2, 0.3, -0.1, 0.02, 500, 0.0),  # Negative kv
        (0.1, 30, 1.5, 0.2, 0.3, 0.1, -0.02, 500, 0.0),  # Negative rho_i
        (0.1, 30, 1.5, 0.2, 0.3, 0.1, 0.02, -500, 0.0),  # Negative fyd
        (0.1, 30, 1.5, 0.2, 0.3, 0.1, 0.02, 500, -0.1),  # Negative kdowel
        (0.1, -30, -1.5, 0.2, 0.3, -0.1, -0.02, -500, -0.1),  # All negative
    ],
)
def test_tau_Rdi_ny_value_errors(
    cv2, fck, gamma_c, mu_v, sigma_n, kv, rho_i, fyd, kdowel
):
    """Test tau_Rdi_ny raises ValueError."""
    with pytest.raises(ValueError):
        _section_8_2_shear.tau_Rdi_ny(
            cv2, fck, gamma_c, mu_v, sigma_n, kv, rho_i, fyd, kdowel
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
