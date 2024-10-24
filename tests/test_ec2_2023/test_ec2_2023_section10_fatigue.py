"""Test for functions from Section 5 of EN 1992-1-1:2023."""

import pytest

from structuralcodes.codes.ec2_2023 import _section10_fatigue


@pytest.mark.parametrize(
    'Ap, phi, phi_p, xi, expected',
    [
        (500, 20, 15, 0.8, 516.397),
        (1000, 25, 20, 0.7, 935.414),
        (800, 22, 18, 0.9, 839.0470),
    ],
)
def test_Ae(Ap, phi, phi_p, xi, expected):
    """Test equivalent area of reinforcement
    calculation with different parameters.
    """
    result = _section10_fatigue.Ae(Ap, phi, phi_p, xi)
    assert pytest.approx(result, 0.01) == expected


@pytest.mark.parametrize(
    'Ap, type, phi_wire, expected',
    [
        (600, 'bundle', None, 39.191),
        (700, 'single_7_wire', 8, 14.0),
        (800, 'single_3_wire', 10, 12.0),
    ],
)
def test_phi_p_eq(Ap, type, phi_wire, expected):
    """Test prestressing steel diameter
    calculation with different parameters.
    """
    result = _section10_fatigue.phi_p_eq(Ap, type, phi_wire)
    assert pytest.approx(result, 0.01) == expected


@pytest.mark.parametrize(
    'theta_uls, expected',
    [
        (1, 1),
        (2, 1.414),
        (0.5, 1),
    ],
)
def test_cot_theta_fat(theta_uls, expected):
    """Test inclination of the compressive
    struts calculation with different ULS angles.
    """
    result = _section10_fatigue.cot_theta_fat(theta_uls)
    assert pytest.approx(result, 0.01) == expected


@pytest.mark.parametrize(
    'fck, state, prestressing_steel_type, expected',
    [
        (60, 'post', 'smooth', 0.225),
        (55, 'pre', 'strand', 0.6),
        (65, 'post', 'indented', 0.40),
        (75, 'pre', 'ribbed', 0.8),
        (50, 'post', 'strand', 0.60),
        (70, 'pre', 'strand', 0.6),
    ],
)
def test_xi_bond(fck, state, prestressing_steel_type, expected):
    """Test bond strength ratio calculation with different parameters."""
    result = _section10_fatigue.xi_bond(fck, state, prestressing_steel_type)
    assert pytest.approx(result, 0.01) == expected


@pytest.mark.parametrize(
    'fck, gamma_c, beta_cc, eta_cc, expected',
    [
        (30, 1.5, 1.0, 1.0, 16.0),
        (25, 1.4, 1.1, 0.9, 15.026),
        (40, 1.6, 0.9, 0.8, 15.3),
        (50, 1.5, 0.9, 1.0, 24),
        (60, 1.4, 1.2, 0.85, 37.157),
    ],
)
def test_fcd_fat(fck, gamma_c, beta_cc, eta_cc, expected):
    """Test calculation of design fatigue strength of concrete."""
    result = _section10_fatigue.fcd_fat(fck, gamma_c, beta_cc, eta_cc)
    assert pytest.approx(result, 0.01) == expected


@pytest.mark.parametrize(
    'sigma_cd_max, sigma_cd_min, fcd_fat, expected',
    [
        (15, 10, 17.0, False),
        (20, 5, 13.23, False),
        (5, 5, 14.4, True),
        (5, 15, 25.5, True),
        (25, 20, 30.6, False),
    ],
)
def test_check_concrete_fatigue_resistance(
    sigma_cd_max, sigma_cd_min, fcd_fat, expected
):
    """Test concrete fatigue resistance check for various stress conditions."""
    result = _section10_fatigue.fatigue_compression_check(
        sigma_cd_max, sigma_cd_min, fcd_fat
    )
    assert result == expected


@pytest.mark.parametrize(
    'tau_ed_max, tau_ed_min, tau_rd_c, expected',
    [
        (0.5, 0.1, 1.0, True),
        (0.9, 0.4, 1.0, False),
        (0.5, -0.2, 1.0, False),
        (0.5, -0.6, 1.0, False),
    ],
)
def test_shear_fatigue_check(tau_ed_max, tau_ed_min, tau_rd_c, expected):
    """Test the shear fatigue check function."""
    result = _section10_fatigue.no_reinf_shear_fatigue_check(
        tau_ed_max, tau_ed_min, tau_rd_c
    )
    assert result == pytest.approx(expected, rel=10e-2)


@pytest.mark.parametrize(
    'surface, expected',
    [
        ('rough', 0.075),
        ('very rough', 0.095),
        ('keyed', 0.185),
    ],
)
def test_cv_1_fat(surface, expected):
    """Test test_cv_1_fat."""
    result = _section10_fatigue.cv_1_fat(surface)
    assert result == expected


@pytest.mark.parametrize(
    'surface, expected',
    [
        ('rough', 0.7),
        ('very rough', 0.9),
        ('keyed', 0.9),
    ],
)
def test_mu_v_fat(surface, expected):
    """Test mu_v_fat."""
    result = _section10_fatigue.mu_v_fat(surface)
    assert result == expected


@pytest.mark.parametrize(
    'cv1_fat, mu_v_fat, sigma_n, f_ck, gamma_c, rho_i, f_yk, alpha, expected',
    [
        (0.075, 0.7, 5.0, 30.0, 1.5, 0.02, 500.0, 30.0, 3.77),
        (0.095, 0.9, 10.0, 40.0, 1.5, 0.03, 600.0, 45.0, 9.40),
    ],
)
def test_shear_strength_no_reinforcement(
    cv1_fat, mu_v_fat, sigma_n, f_ck, gamma_c, rho_i, f_yk, alpha, expected
):
    """Test tau_Rdi_fat_no_reinf."""
    result = _section10_fatigue.tau_Rdi_fat_no_reinf(
        cv1_fat, mu_v_fat, sigma_n, f_ck, gamma_c, rho_i, f_yk, alpha
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'mu_v_fat, sigma_n, rho, sigma_Rsk, gamma_s, alpha, expected',
    [
        (0.7, 5.0, 0.02, 50.0, 1.15, 30.0, 5.849),
        (0.9, 10.0, 0.03, 60.0, 1.15, 45.0, 13.673),
    ],
)
def test_shear_strength_with_reinforcement(
    mu_v_fat, sigma_n, rho, sigma_Rsk, gamma_s, alpha, expected
):
    """Test delta_tau_Rdi_fat_reinf."""
    result = _section10_fatigue.delta_tau_Rdi_fat_reinf(
        mu_v_fat, sigma_n, rho, sigma_Rsk, gamma_s, alpha
    )
    assert result == pytest.approx(expected, rel=1e-2)
