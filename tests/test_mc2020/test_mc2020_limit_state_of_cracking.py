"""Tests for fib MC2020 Chapter 30.5.2.4."""

import pytest

from structuralcodes.codes.mc2020 import _concrete_limit_state_of_cracking


@pytest.mark.parametrize(
    'exposure_class, expected',
    [
        ('X0', 0.4),
        ('XC', 0.3),
        ('XD', 0.3),
        ('XS', 0.3),
    ],
)
def test_wlim_rfc(exposure_class, expected):
    """Test wlim_rfc."""
    result = _concrete_limit_state_of_cracking.wlim_rfc(exposure_class)
    assert result == pytest.approx(expected)


@pytest.mark.parametrize(
    'exposure_class, stress_in_concrete, expected',
    [
        ('X0', 'PL2', 0.4),
        ('XC', 'PL3', 0.3),
        ('XC', 'PL1', 0.2),
        ('XS', 'PL1', 0.0),
    ],
)
def test_wlim_prestressed(exposure_class, stress_in_concrete, expected):
    """Test wlim_prestressed."""
    result = _concrete_limit_state_of_cracking.wlim_prestressed(
        exposure_class, stress_in_concrete
    )
    assert result == pytest.approx(expected)


def test_frp_reinforced_crack_width_limit():
    """Test wlim_frp."""
    result = _concrete_limit_state_of_cracking.wlim_frp()
    assert result == pytest.approx(0.7)


@pytest.mark.parametrize(
    'self_healing, expected',
    [
        (False, 0.15),
        (True, 0.2),
    ],
)
def test_wlim_fluid_tightness(self_healing, expected):
    """Test wlim_fluid_tightness."""
    result = _concrete_limit_state_of_cracking.wlim_fluid_tightness(
        self_healing
    )
    assert result == pytest.approx(expected)


@pytest.mark.parametrize(
    'Ac, fctm, alpha_e, rho_s, expected',
    [
        (30000, 2.9, 1.15, 0.015, 88.50),
        (25000, 3.0, 1.20, 0.02, 76.8),
        (40000, 2.5, 1.10, 0.012, 101.32),
        (35000, 3.2, 1.30, 0.018, 114.62),
        (45000, 2.8, 1.05, 0.014, 127.85),
    ],
)
def test_Nr_crack(Ac, fctm, alpha_e, rho_s, expected):
    """Test Nr_crack with various input scenarios."""
    result = _concrete_limit_state_of_cracking.Nr_crack(
        Ac, fctm, alpha_e, rho_s
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'Es, beta_TS, rho_s, f_ctm, expected',
    [
        (200000, 0.1, 0.015, 2.9, 0.00087),
        (210000, 0.12, 0.02, 3.0, 0.0006285),
        (190000, 0.15, 0.018, 2.7, 0.00067),
        (205000, 0.08, 0.017, 3.1, 0.000818),
        (200000, 0.05, 0.012, 2.5, 0.000989),
    ],
)
def test_eps_crack(Es, beta_TS, rho_s, f_ctm, expected):
    """Test eps_crack with various input scenarios."""
    result = _concrete_limit_state_of_cracking.eps_crack(
        Es, beta_TS, rho_s, f_ctm
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'h, d, x, expected',
    [
        (500, 450, 150, 1.16666),
        (600, 500, 100, 1.25),
    ],
)
def test_k1_r(h, d, x, expected):
    """Test k1_r."""
    result = _concrete_limit_state_of_cracking.k1_r(h, d, x)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'k1_over_r, sr_max, eps_sm_eps_cm, expected',
    [
        (1.0, 510, 0.0013, 0.663),
        (1.25, 600, 0.0019, 1.425),
    ],
)
def test_wcal(k1_over_r, sr_max, eps_sm_eps_cm, expected):
    """Test wcal."""
    result = _concrete_limit_state_of_cracking.wcal(
        k1_over_r, sr_max, eps_sm_eps_cm
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'h, xg, hc_ef, expected',
    [
        (500, 150, 50, 0.9285),
        (600, 200, 60, 0.925),
    ],
)
def test_calculate_kfl(h, xg, hc_ef, expected):
    """Test the calculation of kfl."""
    result = _concrete_limit_state_of_cracking.kfl(h, xg, hc_ef)
    assert result == pytest.approx(expected, rel=1e-4)


@pytest.mark.parametrize(
    'cracking_stage, c, kfl, tau_bms, fctm,'
    + ' rho_s_eff, casting_condition, kc, k_phi_rho, expected',
    [
        ('stabilized', 40, 0.85, 2.5, 2.9, 0.015, 'good', 1.50, 0.25, 127.143),
        ('formation', 50, 0.90, 2.2, 3.1, 0.020, 'poor', 1.60, 0.30, 205.654),
    ],
)
def test_sr_max(
    cracking_stage,
    c,
    kfl,
    tau_bms,
    fctm,
    rho_s_eff,
    casting_condition,
    kc,
    k_phi_rho,
    expected,
):
    """Test the calculation of maximum crack spacing sr,max."""
    result = _concrete_limit_state_of_cracking.sr_max(
        cracking_stage=cracking_stage,
        c=c,
        kfl=kfl,
        tau_bms=tau_bms,
        fctm=fctm,
        rho_s_eff=rho_s_eff,
        casting_condition=casting_condition,
        kc=kc,
        k_phi_rho=k_phi_rho,
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'rx, ry, phi, h, x, expected',
    [
        (30, 40, 16, 500, 100, 12600),
        (50, 60, 20, 600, 150, 24000),
    ],
)
def test_Ac_ef_bar(rx, ry, phi, h, x, expected):
    """Test Ac_ef_bar."""
    result = _concrete_limit_state_of_cracking.Ac_ef_bar(rx, ry, phi, h, x)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'rx, ry, phi, h, x, nl, sy, b, expected',
    [
        (30, 40, 16, 500, 100, 2, 50, 300, 51000),
        (50, 60, 20, 600, 150, 3, 60, 400, 112000),
    ],
)
def test_Ac_ef_group(rx, ry, phi, h, x, nl, sy, b, expected):
    """Test the calculation of Ac_ef_group."""
    result = _concrete_limit_state_of_cracking.Ac_ef_group(
        rx, ry, phi, h, x, nl, sy, b
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'As, Ac_ef, expected',
    [
        (1200, 16000, 0.075),
        (1800, 30000, 0.060),
    ],
)
def test_calculate_rho_s_ef(As, Ac_ef, expected):
    """Test the calculation of rho_s_ef."""
    result = _concrete_limit_state_of_cracking.rho_s_ef(As, Ac_ef)
    assert result == pytest.approx(expected, rel=1e-4)


@pytest.mark.parametrize(
    'sigma_s, sigma_sr_ef, Es, beta_TS, expected',
    [
        (500, 300, 200000, 0.6, 0.0016),
        (400, 200, 210000, 0.6, 0.00133),
    ],
)
def test_relative_mean_strain_direct_loads(
    sigma_s, sigma_sr_ef, Es, beta_TS, expected
):
    """Test relative mean strain for direct loads or imposed strains."""
    result = _concrete_limit_state_of_cracking.eps_sm_eps_cm(
        sigma_s, sigma_sr_ef, Es, beta_TS
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'cracking_stage, srm, expected',
    [
        ('stabilized', 0.0016, 0.00272),
        ('formation', 0.0016, 0.0032),
    ],
)
def test_sr_max_frc(cracking_stage, srm, expected):
    """Test test_sr_max_frc."""
    result = _concrete_limit_state_of_cracking.sr_max_frc(cracking_stage, srm)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'Rax, epsilon_free, beta_TS, fct_eff, Ec, expected',
    [
        (0.8, 0.002, 0.6, 2.5, 30000, 0.00155),
        (0.9, 0.003, 0.4, 2.8, 29000, 0.00266),
    ],
)
def test_relative_mean_strain_restrained(
    Rax, epsilon_free, beta_TS, fct_eff, Ec, expected
):
    """Test relative mean strain for restrained imposed strains."""
    result = _concrete_limit_state_of_cracking.eps_sm_eps_cm_restrained(
        Rax, epsilon_free, beta_TS, fct_eff, Ec
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'fct_eff, rho_s_ef, alpha_e, expected',
    [
        (2.5, 0.01, 8.0, 270.0),
        (3.0, 0.015, 7.0, 221.0),
    ],
)
def test_steel_stress_crack_formation(fct_eff, rho_s_ef, alpha_e, expected):
    """Test steel stress during crack formation stage."""
    result = _concrete_limit_state_of_cracking.sigma_sr_ef(
        fct_eff, rho_s_ef, alpha_e
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'epsilon_restr, epsilon_imp, expected',
    [
        (0.001, 0.002, 0.5),
        (0.002, 0.004, 0.5),
    ],
)
def test_Rax(epsilon_restr, epsilon_imp, expected):
    """Test Rax."""
    result = _concrete_limit_state_of_cracking.Rax(epsilon_restr, epsilon_imp)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'fctm_t, load_type, stage, expected',
    [
        (2.5, 'short-term', 'crack_formation', 4.5),
        (3.0, 'long-term', 'crack_formation', 4.05),
        (3.0, 'long-term', 'stabilized_cracking', 5.4),
    ],
)
def test_tau_bms(fctm_t, load_type, stage, expected):
    """Test tau_bms."""
    result = _concrete_limit_state_of_cracking.tau_bms(
        fctm_t, load_type, stage
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'load_type, stage, expected',
    [
        ('short-term', 'crack_formation', 0.6),
        ('long-term', 'crack_formation', 0.6),
        ('short-term', 'stabilized_cracking', 0.6),
        ('long-term', 'stabilized_cracking', 0.4),
    ],
)
def test_beta_TS(load_type, stage, expected):
    """Test beta_TS."""
    result = _concrete_limit_state_of_cracking.beta_TS(load_type, stage)
    assert result == expected


@pytest.mark.parametrize(
    'diameters, counts, expected',
    [
        ([12, 16], [4, 2], 13.6),
        ([10, 20, 25], [3, 2, 1], 18.157),
        ([8, 16, 20], [5, 1, 1], 12.842),
    ],
)
def test_equivalent_diameter(diameters, counts, expected):
    """Test calculation of equivalent bar diameter
    φeq for various bar diameters and counts.
    """
    result = _concrete_limit_state_of_cracking.phi_eq(diameters, counts)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'sr_max_x, sr_max_y, theta, expected',
    [
        (200, 300, 45, 169.705),
        (250, 250, 30, 183.012),
        (150, 350, 60, 172.18),
    ],
)
def test_crack_spacing_theta(sr_max_x, sr_max_y, theta, expected):
    """Test calculation of maximum crack spacing s_r,max,θ."""
    result = _concrete_limit_state_of_cracking.sr_max_theta(
        sr_max_x, sr_max_y, theta
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'sigma_x, sigma_y, tau_xy, rho_s_ef_x, rho_s_ef_y, expected',
    [
        (10, 5, 2, 0.01, 0.015, 34.299),
        (12, 7, 3, 0.012, 0.014, 38.657),
        (8, 4, 1.5, 0.015, 0.018, 36.576),
    ],
)
def test_theta_reinf(
    sigma_x, sigma_y, tau_xy, rho_s_ef_x, rho_s_ef_y, expected
):
    """Test theta_reinf."""
    result = _concrete_limit_state_of_cracking.theta_reinf(
        sigma_x, sigma_y, tau_xy, rho_s_ef_x, rho_s_ef_y
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'sigma_s_x, sigma_sr_ef_x, Es, beta_TS_x, expected',
    [
        (500, 300, 200000, 0.6, 0.0016),
        (400, 200, 210000, 0.6, 0.00133),
    ],
)
def test_eps_sm_x_eps_cm_x(sigma_s_x, sigma_sr_ef_x, Es, beta_TS_x, expected):
    """Test eps_sm_x_eps_cm_x."""
    result = _concrete_limit_state_of_cracking.eps_sm_x_eps_cm_x(
        sigma_s_x, sigma_sr_ef_x, Es, beta_TS_x
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'sigma_s_y, sigma_sr_ef_y, Es, beta_TS_y, expected',
    [
        (500, 300, 200000, 0.6, 0.0016),
        (400, 200, 210000, 0.6, 0.00133),
    ],
)
def test_eps_sm_y_eps_cm_y(sigma_s_y, sigma_sr_ef_y, Es, beta_TS_y, expected):
    """Test eps_sm_y_eps_cm_y."""
    result = _concrete_limit_state_of_cracking.eps_sm_y_eps_cm_y(
        sigma_s_y, sigma_sr_ef_y, Es, beta_TS_y
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'tau_xy, Ec, theta, expected',
    [
        (2, 30000, 45, -0.00013333),
        (1.5, 25000, 30, -0.0001385),
    ],
)
def test_eps_2(tau_xy, Ec, theta, expected):
    """Test principal compressive strain calculation."""
    result = _concrete_limit_state_of_cracking.eps_2(tau_xy, Ec, theta)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'sigma_x, tau_xy, rho_s_ef_x, theta, expected',
    [
        (10, 5, 0.01, 30, 1288.67),
        (8, 4, 0.012, 45, 1000.0),
    ],
)
def test_sigma_s_x(sigma_x, tau_xy, rho_s_ef_x, theta, expected):
    """Test sigma_s_x."""
    result = _concrete_limit_state_of_cracking.sigma_s_x(
        sigma_x, tau_xy, rho_s_ef_x, theta
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'sigma_y, tau_xy, rho_s_ef_y, theta, expected',
    [
        (10, 5, 0.01, 30, 1866.025),
        (8, 4, 0.012, 45, 1000.0),
    ],
)
def test_sigma_s_y(sigma_y, tau_xy, rho_s_ef_y, theta, expected):
    """Test sigma_s_y."""
    result = _concrete_limit_state_of_cracking.sigma_s_y(
        sigma_y, tau_xy, rho_s_ef_y, theta
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'beta_TS, sr_max_theta, theta, sr_max_x, expected',
    [
        (1.2, 200, 45, 150, 2.2627),
        (1.2, 180, 60, 130, 1.9185),
    ],
)
def test_beta_TS_x(beta_TS, sr_max_theta, theta, sr_max_x, expected):
    """Test beta_TS_x."""
    result = _concrete_limit_state_of_cracking.beta_TS_x(
        beta_TS, sr_max_theta, theta, sr_max_x
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'beta_TS, sr_max_theta, theta, sr_max_y, expected',
    [
        (1.2, 200, 45, 150, 2.2627),
        (1.2, 180, 60, 130, 3.3230),
    ],
)
def test_beta_TS_y(beta_TS, sr_max_theta, theta, sr_max_y, expected):
    """Test beta_TS_y."""
    result = _concrete_limit_state_of_cracking.beta_TS_y(
        beta_TS, sr_max_theta, theta, sr_max_y
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'eps_sm_x_eps_cm_x, eps_sm_y_eps_cm_y, eps_2, expected',
    [
        (0.0015, 0.0020, -0.001, 0.0005),
        (0.0020, 0.001, -0.0016, 0.0026),
    ],
)
def test_eps_sm_eps_cm_theta(
    eps_sm_x_eps_cm_x, eps_sm_y_eps_cm_y, eps_2, expected
):
    """Test eps_sm_eps_cm_theta."""
    result = _concrete_limit_state_of_cracking.eps_sm_eps_cm_theta(
        eps_sm_x_eps_cm_x, eps_sm_y_eps_cm_y, eps_2
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'rho_s_ef, rho_p_ef, xi_1, expected',
    [
        (0.01, 0.005, 0.8, 0.0132),
        (0.012, 0.004, 0.9, 0.01524),
    ],
)
def test_rho_s_p_ef(rho_s_ef, rho_p_ef, xi_1, expected):
    """Test calculation of equivalent reinforcement ratio ρs+p,ef."""
    result = _concrete_limit_state_of_cracking.rho_s_p_ef(
        rho_s_ef, rho_p_ef, xi_1
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'tau_bmp_tau_bms, phi, phi_p_eq, expected',
    [
        (1.25, 16, 12, 1.2909),
        (2, 20, 15, 1.6329),
    ],
)
def test_xi_1(tau_bmp_tau_bms, phi, phi_p_eq, expected):
    """Test xi_1."""
    result = _concrete_limit_state_of_cracking.xi_1(
        tau_bmp_tau_bms, phi, phi_p_eq
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'Ap, up, expected',
    [
        (150, 60, 10.0),
        (200, 80, 10.0),
    ],
)
def test_phi_p_eq(Ap, up, expected):
    """Test phi_p_eq."""
    result = _concrete_limit_state_of_cracking.phi_p_eq(Ap, up)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'Ap_i, type_, phi_wire, expected',
    [
        ('bundle', 200, None, 71.086),
        ('7-wire', None, 20, 109.956),
        ('3-wire', None, 12, 45.239),
    ],
)
def test_up_i(Ap_i, type_, phi_wire, expected):
    """Test up_i."""
    result = _concrete_limit_state_of_cracking.up_i(Ap_i, type_, phi_wire)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'surface_condition, tensioning_type, expected',
    [
        ('no bond', 'pretensioned', 0.0),
        ('smooth wire', 'pretensioned', 0.4),
        ('strand', 'pretensioned', 0.6),
        ('indented wire', 'pretensioned', 0.8),
        ('ribbed bar', 'post-tensioned', 1.0),
        ('strand', 'post-tensioned', 0.4),
        ('no bond', 'post-tensioned', 0.0),
        ('smooth wire', 'post-tensioned', 0.2),
        ('indented wire', 'post-tensioned', 0.6),
        ('smooth wire', 'no-bond', 0.0),
    ],
)
def test_tau_bmp_tau_bms(surface_condition, tensioning_type, expected):
    """Test tau_bmp_tau_bms."""
    result = _concrete_limit_state_of_cracking.tau_bmp_tau_bms(
        surface_condition, tensioning_type
    )
    assert result == expected


@pytest.mark.parametrize(
    'fctm, fFtsm, Act, sigma_s, kc, h, expected',
    [
        (2.5, 1.0, 5000, 400, 1.0, 300, 18.75),
        (3.0, 2.0, 6000, 500, 1.0, 500, 10.32),
        (2.0, 2.5, 5500, 450, 1.0, 800, 0.0),
    ],
)
def test_As_min_frc(fctm, fFtsm, Act, sigma_s, kc, h, expected):
    """Test As_min_frc."""
    result = _concrete_limit_state_of_cracking.As_min_frc(
        fctm, fFtsm, Act, sigma_s, h, kc
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'rho_s_ef, ry, d, w_lim_cal, beta_w, kf_fss, '
    + 'k_b, k1_r_simpl, sigma_s, Es, c, expected',
    [
        (0.01, 10, 200, 0.3, 1.0, 0.9, 1.0, 1.0, 400, 200000, 30, 56.777),
        (0.012, 12, 250, 0.25, 1.0, 0.85, 0.9, 0.95, 450, 210000, 40, 52.467),
    ],
)
def test_max_phi_crack(
    rho_s_ef,
    ry,
    d,
    w_lim_cal,
    beta_w,
    kf_fss,
    k_b,
    k1_r_simpl,
    sigma_s,
    Es,
    c,
    expected,
):
    """Test max_phi_crack."""
    result = _concrete_limit_state_of_cracking.max_phi_crack(
        rho_s_ef,
        ry,
        d,
        w_lim_cal,
        beta_w,
        kf_fss,
        k_b,
        k1_r_simpl,
        sigma_s,
        Es,
        c,
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'rho_s_ef, ry, d, w_lim_cal, beta_w, kf_fss, k_b, '
    + 'k1_r_simpl, sigma_s, Es, c, expected',
    [
        (0.01, 10, 200, 0.3, 1.0, 0.9, 1.0, 1.0, 400, 200000, 30, 1260.97),
        (0.012, 12, 250, 0.25, 1.0, 0.85, 0.9, 0.95, 450, 210000, 40, 717.85),
    ],
)
def test_max_bar_spacing(
    rho_s_ef,
    ry,
    d,
    w_lim_cal,
    beta_w,
    kf_fss,
    k_b,
    k1_r_simpl,
    sigma_s,
    Es,
    c,
    expected,
):
    """Test calculation of maximum allowable bar spacing sl."""
    result = _concrete_limit_state_of_cracking.max_st_crack(
        rho_s_ef,
        ry,
        d,
        w_lim_cal,
        beta_w,
        kf_fss,
        k_b,
        k1_r_simpl,
        sigma_s,
        Es,
        c,
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'r_y, h, tension, expected',
    [
        (10, 200, 'one-face', 0.825),
        (5, 300, 'two-faces', 1.0),
    ],
)
def test_kfl_simpl(r_y, h, tension, expected):
    """Test kfl_simpl."""
    result = _concrete_limit_state_of_cracking.kfl_simpl(r_y, h, tension)
    assert result == pytest.approx(expected, rel=1e-3)


@pytest.mark.parametrize(
    'rho_s_p_ef, h, d, state, expected',
    [
        (0.01, 250, 200, 'tension', 1.0),
        (0.015, 300, 250, 'bending', 1.0),
    ],
)
def test_k1_r_simpl(rho_s_p_ef, h, d, state, expected):
    """Test k1_r_simpl."""
    result = _concrete_limit_state_of_cracking.k1_r_simpl(
        rho_s_p_ef, h, d, state
    )
    assert result == pytest.approx(expected, rel=1e-3)
