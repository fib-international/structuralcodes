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
    assert _section_8_2_shear.tao_Ed(VEd, bw, d) == pytest.approx(expected)


@pytest.mark.parametrize(
    'vEd, d, expected',
    [
        (50.0, 0.5, 50.0 / (0.9 * 0.5)),
        (100.0, 0.6, 100.0 / (0.9 * 0.6)),
    ],
)
def test_tao_Ed_planar(vEd, d, expected):
    """Test shear_stress_planar_members."""
    assert _section_8_2_shear.tao_Ed_planar(vEd, d) == pytest.approx(expected)


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
    result = _section_8_2_shear.tau_rdc_min(gamma_v, f_ck, f_yd, d, d_lower)
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
        _section_8_2_shear.d_eff_with_angle(dx, dy, vEd_x, vEd_y),
        expected,
        rel_tol=1e-9,
    )


@pytest.mark.parametrize(
    'gamma_v, rho_l, f_ck, d, d_g, tau_rdc_min, expected',
    [
        (1.5, 0.02, 30, 500, 16, 0.3, 0.5468),
        (1.4, 0.03, 40, 450, 20, 0.5, 0.82366),
        (1.6, 0.025, 35, 600, 18, 0.1, 0.5690),
    ],
)
def test_calculate_tau_Rdc(
    gamma_v, rho_l, f_ck, d, d_g, tau_rdc_min, expected
):
    """Test the calculation of the shear stress resistance."""
    result = _section_8_2_shear.tau_Rdc(
        gamma_v, rho_l, f_ck, d, d_g, tau_rdc_min
    )
    assert math.isclose(
        result,
        expected,
        rel_tol=1e-3,
    )


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
    'gamma_v, rho_l, f_ck, d, d_g, expected',
    [
        (1.5, 0.02, 30, 500, 16, 0.5468),
        (1.4, 0.03, 40, 450, 20, 0.8236),
        (1.6, 0.025, 35, 600, 18, 0.5690),
    ],
)
def test_calculate_tau_Rdc_0(gamma_v, rho_l, f_ck, d, d_g, expected):
    """Test the sh stress resistance wo/ axial force effects."""
    result = _section_8_2_shear.tau_Rdc_0(gamma_v, rho_l, f_ck, d, d_g)
    assert math.isclose(
        result,
        expected,
        rel_tol=1e-3,
    )


@pytest.mark.parametrize(
    'tau_Rdc_0, k1, sigma_cp, tau_Rdc_max, expected',
    [
        (1, 0.5, 0.1, 2, 0.95),
        (1, 0.6, 0.2, 2, 0.88),
        (1, 0.4, 0.3, 2, 0.88),
    ],
)
def test_calculate_tau_Rdc_comp(
    tau_Rdc_0, k1, sigma_cp, tau_Rdc_max, expected
):
    """Test the calculation of the shear considering comp normal forces."""
    assert math.isclose(
        _section_8_2_shear.tau_Rdc_comp(tau_Rdc_0, k1, sigma_cp, tau_Rdc_max),
        expected,
        rel_tol=1e-5,
    )


@pytest.mark.parametrize(
    'a_cs_0, e_p, A_c, b_w, z, d, expected',
    [
        (1000, 50, 10000, 200, 500, 200, 0.018),
        (1200, 60, 12000, 250, 600, 200, 0.0144),
        (1100, 55, 11000, 220, 550, 200, 0.01636),
    ],
)
def test_calculate_k1(a_cs_0, e_p, A_c, b_w, z, d, expected):
    """Test the calculation of the factor k1."""
    result = _section_8_2_shear.k1(a_cs_0, e_p, A_c, b_w, z, d)
    assert math.isclose(result, expected, rel_tol=1e-3)


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
