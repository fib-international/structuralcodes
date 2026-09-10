"""Tests for fib MC2020 Chapter 14.6.1."""

import pytest

from structuralcodes.codes.mc2020 import _concrete_material_properties


@pytest.mark.parametrize(
    'fck, delta_f, expected',
    [
        (30, 8, 38),
        (40, 10, 50),
        (50, 5, 55),
    ],
)
def test_fcm(fck, delta_f, expected):
    """Test fcm."""
    result = _concrete_material_properties.fcm(fck, delta_f)
    assert result == expected


@pytest.mark.parametrize(
    'flck, delta_f, expected',
    [
        (25, 8, 33),
        (35, 10, 45),
        (45, 5, 50),
    ],
)
def test_flcm(flck, delta_f, expected):
    """Test flcm."""
    result = _concrete_material_properties.flcm(flck, delta_f)
    assert result == expected


@pytest.mark.parametrize(
    'grade, expected',
    [
        ('C20', 20),
        ('C30', 30),
        ('C50', 50),
        ('C100', 100),
    ],
)
def test_fck(grade, expected):
    """Test fck."""
    result = _concrete_material_properties.fck(grade)
    assert result == expected


@pytest.mark.parametrize(
    'grade, expected',
    [
        ('LC12', 13),
        ('LC25', 28),
        ('LC50', 55),
        ('LC80', 88),
    ],
)
def test_fck_cube(grade, expected):
    """Test fck_cube."""
    result = _concrete_material_properties.fck_cube(grade)
    assert result == expected


@pytest.mark.parametrize(
    'fck, expected',
    [
        (20, 2.292),
        (30, 3.022),
        (50, 3.942),
    ],
)
def test_calculate_fctm(fck, expected):
    """Test calculation of mean tensile strength fctm."""
    result = _concrete_material_properties.fctm(fck)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'fctm, expected',
    [
        (3.5, 2.45),
        (4.0, 2.8),
        (5.0, 3.5),
    ],
)
def test_fctk_min(fctm, expected):
    """Test fctk_min."""
    result = _concrete_material_properties.fctk_min(fctm)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'fctm, expected',
    [
        (3.5, 4.55),
        (4.0, 5.2),
        (5.0, 6.5),
    ],
)
def test_fctk_max(fctm, expected):
    """Test fctk_max."""
    result = _concrete_material_properties.fctk_max(fctm)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'fctm, density, expected',
    [
        (4.0, 1800, 3.563),
        (5.0, 1500, 4.045),
        (3.0, 2000, 2.836),
    ],
)
def test_flctm(fctm, density, expected):
    """Test flctm."""
    result = _concrete_material_properties.flctm(fctm, density)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'fctm, expected',
    [
        (3.5, 2.45),
        (4.0, 2.8),
        (5.0, 3.5),
    ],
)
def test_flctk_min(fctm, expected):
    """Test flctk_min."""
    result = _concrete_material_properties.flctk_min(fctm)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'fctm, expected',
    [
        (3.5, 4.55),
        (4.0, 5.2),
        (5.0, 6.5),
    ],
)
def test_flctk_max(fctm, expected):
    """Test flctk_max."""
    result = _concrete_material_properties.flctk_max(fctm)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'fctm_sp, alpha_sp, expected_fctm',
    [
        (2.5, 1.0, 2.5),
        (3.0, 1.0, 3.0),
        (4.0, 1.0, 4.0),
    ],
)
def test_fctm_from_splitting(fctm_sp, alpha_sp, expected_fctm):
    """Test fctm_from_splitting."""
    result = _concrete_material_properties.fctm_from_splitting(
        fctm_sp, alpha_sp
    )
    assert result == expected_fctm


@pytest.mark.parametrize(
    'fctm_fl, hb, expected_fctm',
    [
        (4.0, 500, 3.292),
        (3.5, 300, 2.676),
        (5.0, 800, 4.329),
    ],
)
def test_fctm_from_flexural(fctm_fl, hb, expected_fctm):
    """Test fctm_from_flexural."""
    result = _concrete_material_properties.fctm_from_flexural(fctm_fl, hb)
    assert result == pytest.approx(expected_fctm, rel=1e-2)


@pytest.mark.parametrize(
    'fck, expected_GF',
    [
        (20, 133.220),
        (30, 141.575),
        (50, 152.849),
    ],
)
def test_GF(fck, expected_GF):
    """Test GF."""
    result = _concrete_material_properties.GF(fck)
    assert result == pytest.approx(expected_GF, rel=1e-2)


@pytest.mark.parametrize(
    'flctm, use_normal_weight_sand, expected_GF_l',
    [
        (3.5, True, 80),
        (4.0, False, 64),
        (5.0, True, 104),
    ],
)
def test_GF_l(flctm, use_normal_weight_sand, expected_GF_l):
    """Test GF_l."""
    result = _concrete_material_properties.GF_l(flctm, use_normal_weight_sand)
    assert result == expected_GF_l


@pytest.mark.parametrize(
    'sigma_1, sigma_2, sigma_3, fcm, fctm, fc2cm, sigma_com, '
    + 'tau_com, is_lightweight, expected',
    [
        # Test case 1: Normal weight concrete with given sigma_com and tau_com
        (30, 20, 10, 40, 4, 35, None, None, False, -16.71),
        (25, 15, 5, 30, 3.5, 28, -240, None, False, 115.88),
        (20, 15, 10, 25, 2.5, 22, None, None, True, 5.51),
        (15, 10, 5, 20, 2, 18, -60, 75, True, 4.99),
    ],
)
def test_multiaxial_stress_equation(
    sigma_1,
    sigma_2,
    sigma_3,
    fcm,
    fctm,
    fc2cm,
    sigma_com,
    tau_com,
    is_lightweight,
    expected,
):
    """Test multiaxial_stress_equation."""
    result = _concrete_material_properties.multiaxial_stress_equation(
        sigma_1=sigma_1,
        sigma_2=sigma_2,
        sigma_3=sigma_3,
        fcm=fcm,
        fctm=fctm,
        fc2cm=fc2cm,
        sigma_com=sigma_com,
        tau_com=tau_com,
        is_lightweight=is_lightweight,
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'fck, alpha_E, expected',
    [
        (30, 1.0, 31008.365),
        (40, 1.2, 40954.947),
        (50, 0.8, 29411.586),
    ],
)
def test_E_ci(fck, alpha_E, expected):
    """Test Eci."""
    result = _concrete_material_properties.E_ci(fck, alpha_E)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'fcm, density, alpha_E, expected',
    [
        (30, 1800, 1.0, 20757.67),
        (40, 1600, 1.2, 21662.12),
        (50, 2000, 0.8, 24307.10),
    ],
)
def test_El_ci(fcm, density, alpha_E, expected):
    """Test El_ci."""
    result = _concrete_material_properties.El_ci(fcm, density, alpha_E)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'Eci, fcm, expected',
    [
        (30000, 30, 26045.45),
        (35000, 40, 31181.82),
        (40000, 50, 36545.45),
    ],
)
def test_E_ci_red_el(Eci, fcm, expected):
    """Test E_ci_red_el."""
    result = _concrete_material_properties.E_ci_red_el(Eci, fcm)
    assert result == pytest.approx(expected, rel=1e-2)


def test_nu_c():
    """Test test_nu_c."""
    result = _concrete_material_properties.nu_c()
    assert result == 0.20


@pytest.mark.parametrize(
    'concrete_grade, expected',
    [
        ('C20', 13300),
        ('C30', 16500),
        ('C50', 23200),
        ('C100', 36000),
    ],
)
def test_Ec1(concrete_grade, expected):
    """Test Ec1."""
    result = _concrete_material_properties.Ec1(concrete_grade)
    assert result == expected


@pytest.mark.parametrize(
    'concrete_grade, expected',
    [
        ('C20', 0.0021),
        ('C30', 0.0023),
        ('C50', 0.0026),
        ('C100', 0.003),
    ],
)
def test_eps_c1(concrete_grade, expected):
    """Test eps_c1."""
    result = _concrete_material_properties.eps_c1(concrete_grade)
    assert result == pytest.approx(expected, rel=1e-1)


@pytest.mark.parametrize(
    'concrete_grade, expected',
    [
        ('C20', 0.0035),
        ('C30', 0.0035),
        ('C50', 0.0034),
        ('C100', 0.003),
    ],
)
def test_eps_c_lim(concrete_grade, expected):
    """Test eps_c_lim."""
    result = _concrete_material_properties.eps_c_lim(concrete_grade)
    assert result == pytest.approx(expected, rel=1e-1)


@pytest.mark.parametrize(
    'f_lck, E_lc, kappa_lc, expected',
    [
        (25, 18000, 1.1, 0.0024),
        (30, 20000, 1.3, 0.00247),
    ],
)
def test_eps_lc1(f_lck, E_lc, kappa_lc, expected):
    """Test calculation of strain εlc1 for lightweight aggregate concrete."""
    result = _concrete_material_properties.eps_lc1(f_lck, E_lc, kappa_lc)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'epsilon_ct, E_ci, f_ctm, expected_sigma_ct',
    [
        (0.0001, 30000, 3.0, 2.75),
        (0.00005, 30000, 3.0, 1.5),
        (0.0001, 30000, 3.5, 3.0),
    ],
)
def test_sigma_ct_uncracked(epsilon_ct, E_ci, f_ctm, expected_sigma_ct):
    """Test csigma_ct_uncracked."""
    result = _concrete_material_properties.sigma_ct_uncracked(
        epsilon_ct, E_ci, f_ctm
    )
    assert result == pytest.approx(expected_sigma_ct, rel=1e-2)


@pytest.mark.parametrize(
    'w, GF, f_ctm, expected_sigma_ct',
    [
        (0.001, 0.15, 3.0, 2.95),
        (0.01, 0.15, 3.0, 2.5),
        (0.05, 0.15, 3.0, 0.6),
    ],
)
def test_sigma_ct_cracked(w, GF, f_ctm, expected_sigma_ct):
    """Test sigma_ct_cracked."""
    result = _concrete_material_properties.sigma_ct_cracked(w, GF, f_ctm)
    assert result == pytest.approx(expected_sigma_ct, rel=1e-2)


@pytest.mark.parametrize(
    'delta, w, f_cm, C, expected',
    [
        (0.5, 0.2, 30, 1.0, 11.825),
        (0.3, 0.1, 40, 1.0, 16.368),
    ],
)
def test_calculate_mean_shear_stress(delta, w, f_cm, C, expected):
    """Test calculation of mean shear stress τ in an open crack."""
    result = _concrete_material_properties.tau_crack_friction(
        delta, w, f_cm, C
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'delta, w, f_cm, C, expected',
    [
        (0.5, 0.2, 30, 1.0, 6.007),
        (0.3, 0.1, 40, 1.0, 7.351),
    ],
)
def test_calculate_mean_normal_stress(delta, w, f_cm, C, expected):
    """Test calculation of mean normal stress σ in an open crack."""
    result = _concrete_material_properties.sigma_crack_friction(
        delta, w, f_cm, C
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'f_ck, strength_class, expected',
    [
        (30, 'CR', 0.3),
        (50, 'CN', 0.4),
        (70, 'CS', 0.4),
    ],
)
def test_get_sc_coefficient(f_ck, strength_class, expected):
    """Test sC."""
    result = _concrete_material_properties.sC(f_ck, strength_class)
    assert result == expected


@pytest.mark.parametrize(
    't, t_ref, s_C, expected_beta_cc',
    [
        (14, 28, 0.5, 0.8129),
        (7, 28, 0.3, 0.7408),
    ],
)
def test_beta_cc(t, t_ref, s_C, expected_beta_cc):
    """Test beta_cc."""
    result = _concrete_material_properties.beta_cc(t, t_ref, s_C)
    assert result == pytest.approx(expected_beta_cc, rel=1e-2)


@pytest.mark.parametrize(
    'f_cm, beta_cc, expected',
    [
        (38, 0.7071, 26.67),
        (48, 0.5443, 26.13),
    ],
)
def test_fcm_t(f_cm, beta_cc, expected):
    """Test fcm_t."""
    result = _concrete_material_properties.fcm_t(f_cm, beta_cc)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'aggregate_strength, expected_slc',
    [
        ('high', 0.05),
        ('low', 0.25),
    ],
)
def test_slC(aggregate_strength, expected_slc):
    """Test slC."""
    result = _concrete_material_properties.slC(aggregate_strength)
    assert result == expected_slc


@pytest.mark.parametrize(
    't, t_ref, s_lc, expected',
    [
        (14, 28, 0.05, 0.979),
        (7, 28, 0.25, 0.779),
    ],
)
def test_beta_lcc(t, t_ref, s_lc, expected):
    """Test beta_lcc."""
    result = _concrete_material_properties.beta_lcc(t, t_ref, s_lc)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'f_lcm, beta_lcc, expected',
    [
        (30, 0.747, 22.41),
        (35, 0.657, 22.995),
    ],
)
def test_flcm_t(f_lcm, beta_lcc, expected):
    """Test flcm_t."""
    result = _concrete_material_properties.flcm_t(f_lcm, beta_lcc)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'beta_cc, expected',
    [
        (0.75, 0.909),
        (0.8, 0.929),
    ],
)
def test_beta_E(beta_cc, expected):
    """Test beta_E."""
    result = _concrete_material_properties.beta_E(beta_cc)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'E_ci, beta_E_t, expected',
    [
        (33000, 0.9005, 29716.5),
        (35000, 0.8178, 28623),
    ],
)
def test_E_ci_t(E_ci, beta_E_t, expected):
    """Test E_ci_t."""
    result = _concrete_material_properties.E_ci_t(E_ci, beta_E_t)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    't0, expected',
    [
        (7, 0.6621),
        (28, 0.6743),
    ],
)
def test_beta_t0(t0, expected):
    """Test beta_t0."""
    result = _concrete_material_properties.beta_t0(t0)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    't, t0, beta_T0, expected',
    [
        (180, 28, 0.9, 0.933),
        (365, 28, 0.85, 0.896),
    ],
)
def test_beta_c_sus(t, t0, beta_T0, expected):
    """Test beta_c_sus."""
    result = _concrete_material_properties.beta_c_sus(t, t0, beta_T0)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'f_cm, beta_cc_t, beta_c_sus_t_t0, expected',
    [
        (30, 0.8, 0.936, 22.464),
        (40, 0.9, 0.922, 33.192),
    ],
)
def test_fcm_sus(f_cm, beta_cc_t, beta_c_sus_t_t0, expected):
    """Test fcm_sus."""
    result = _concrete_material_properties.fcm_sus(
        f_cm, beta_cc_t, beta_c_sus_t_t0
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'f_ctk, alpha, expected',
    [
        (3, 'normal', 1.8),
        (4, 'high-strength', 3.0),
    ],
)
def test_fctk_sus(f_ctk, alpha, expected):
    """Test fctk_sus."""
    result = _concrete_material_properties.fctk_sus(f_ctk, alpha)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'sigma_c_t0, phi_t_t0, E_ci, expected',
    [
        (10, 2.0, 30000, 0.0006667),
        (15, 1.5, 28000, 0.0008036),
    ],
)
def test_eps_cc(sigma_c_t0, phi_t_t0, E_ci, expected):
    """Test eps_cc."""
    result = _concrete_material_properties.eps_cc(sigma_c_t0, phi_t_t0, E_ci)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'sigma_c_t0, E_ci_t0, E_ci, phi_t_t0, expected',
    [
        (10, 25000, 28000, 2.0, 0.00111),
        (15, 24500, 29000, 1.5, 0.00139),
    ],
)
def test_eps_c_sigma(sigma_c_t0, E_ci_t0, E_ci, phi_t_t0, expected):
    """Test eps_c_sigma."""
    result = _concrete_material_properties.eps_c_sigma(
        sigma_c_t0, E_ci_t0, E_ci, phi_t_t0
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'phi_bc_t_t0, phi_dc_t_t0, expected_phi',
    [
        (0.5, 0.3, 0.8),
        (0.0, 0.0, 0.0),
        (1.2, 0.8, 2.0),
    ],
)
def test_phi_t_t0(phi_bc_t_t0, phi_dc_t_t0, expected_phi):
    """Test calculation of total creep coefficient φ(t, t0)."""
    result = _concrete_material_properties.phi_t_t0(phi_bc_t_t0, phi_dc_t_t0)
    assert result == pytest.approx(expected_phi, rel=1e-2)


@pytest.mark.parametrize(
    'fcm, expected',
    [
        (30, 0.166),
        (40, 0.136),
        (50, 0.116),
    ],
)
def test_beta_bc_fcm(fcm, expected):
    """Test beta_bc_fcm."""
    result = _concrete_material_properties.beta_bc_fcm(fcm)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    't, t0, t0_adj, expected',
    [
        (180, 28, 14, 6.55),
        (365, 28, 28, 5.97),
    ],
)
def test_beta_bc_t_t0(t, t0, t0_adj, expected):
    """Test beta_bc_t_t0."""
    result = _concrete_material_properties.beta_bc_t_t0(t, t0, t0_adj)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'beta_bc_t_t0, beta_bc_fcm, expected',
    [
        (2.0, 0.22, 0.44),
        (1.5, 0.30, 0.45),
    ],
)
def test_phi_bc_t_t0(beta_bc_t_t0, beta_bc_fcm, expected):
    """Test phi_bc_t_t0."""
    result = _concrete_material_properties.phi_bc_t_t0(
        beta_bc_t_t0, beta_bc_fcm
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'fcm, expected',
    [
        (30, 3.523),
        (40, 2.355),
    ],
)
def test_beta_dc_fcm(fcm, expected):
    """Test beta_dc_fcm."""
    result = _concrete_material_properties.beta_dc_fcm(fcm)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'RH, h, expected',
    [
        (50, 200, 0.854),
        (75, 100, 0.538),
    ],
)
def test_beta_RH(RH, h, expected):
    """Test beta_RH."""
    result = _concrete_material_properties.beta_RH_c(RH, h)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    't0_adj, expected',
    [
        (28, 0.488),
        (14, 0.557),
    ],
)
def test_beta_dc_t0(t0_adj, expected):
    """Test beta_dc_t0."""
    result = _concrete_material_properties.beta_dc_t0(t0_adj)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    't0_adj, expected',
    [
        (28, 0.338),
        (14, 0.309),
    ],
)
def test_gamma_t0(t0_adj, expected):
    """Test gamma_t0."""
    result = _concrete_material_properties.gamma_t0(t0_adj)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'h, fcm, expected',
    [
        (200, 30, 570.030),
        (100, 40, 383.85),
    ],
)
def test_beta_h(h, fcm, expected):
    """Test beta_h."""
    result = _concrete_material_properties.beta_h(h, fcm)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    't, t0, gamma_t0, beta_h, expected',
    [
        (180, 28, 0.191, 875, 0.694),
        (365, 28, 0.228, 488, 0.815),
    ],
)
def test_beta_dc_t_t0(t, t0, gamma_t0, beta_h, expected):
    """Test beta_dc_t_t0."""
    result = _concrete_material_properties.beta_dc_t_t0(
        t, t0, gamma_t0, beta_h
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'beta_dc_fm, beta_RH, beta_dc_t0, beta_dc_t_t0, expected',
    [
        (3.67, 0.0886, 0.754, 0.143, 0.035),
        (2.58, 0.114, 0.852, 0.051, 0.0127),
    ],
)
def test_phi_dc_t_t0(beta_dc_fm, beta_RH, beta_dc_t0, beta_dc_t_t0, expected):
    """Test phi_dc_t_t0."""
    result = _concrete_material_properties.phi_dc_t_t0(
        beta_dc_fm, beta_RH, beta_dc_t0, beta_dc_t_t0
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'phi_t_t0, rho, concrete_grade, expected',
    [
        (2.0, 1800, 'LC12', 1.740),
        (2.0, 1800, 'LC16', 1.740),
        (2.0, 1800, 'LC20', 1.339),
        (1.5, 2000, 'LC16', 1.612),
    ],
)
def test_phi_l_t_t0(phi_t_t0, rho, concrete_grade, expected):
    """Test phi_l_t_t0."""
    result = _concrete_material_properties.phi_l_t_t0(
        phi_t_t0, rho, concrete_grade
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    't0_T, strength_class, expected',
    [
        (28, 'CS', 24.15),
        (28, 'CN', 28.0),
        (28, 'CR', 32.46),
        (14, 'CN', 14.0),
    ],
)
def test_t0_adj(t0_T, strength_class, expected):
    """Test t0_adj."""
    result = _concrete_material_properties.t0_adj(t0_T, strength_class)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'phi_t_t0, sigma_c, f_cm_t0, expected_phi_sigma',
    [
        (1.2, 18, 30, 1.620),
        (1.5, 15, 25, 2.025),
    ],
)
def test_phi_sigma_t_t0(phi_t_t0, sigma_c, f_cm_t0, expected_phi_sigma):
    """Test phi_sigma_t_t0."""
    result = _concrete_material_properties.phi_sigma_t_t0(
        phi_t_t0, sigma_c, f_cm_t0
    )
    assert result == pytest.approx(expected_phi_sigma, rel=1e-2)


@pytest.mark.parametrize(
    'eps_cbs_t, eps_cds_t_ts, expected',
    [
        (0.00024, 0.000168, 0.000408),
        (0.0002, 0.0002, 0.0004),
    ],
)
def test_eps_cs_t_ts(eps_cbs_t, eps_cds_t_ts, expected):
    """Test eps_cs_t_ts."""
    result = _concrete_material_properties.eps_cs_t_ts(eps_cbs_t, eps_cds_t_ts)
    assert result == expected


@pytest.mark.parametrize(
    'eps_cbs0, beta_bs_t, expected',
    [
        (0.0003, 0.8, 0.00024),
        (0.0002, 1.0, 0.0002),
    ],
)
def test_eps_cbs_t(eps_cbs0, beta_bs_t, expected):
    """Test eps_cbs_t."""
    result = _concrete_material_properties.eps_cbs_t(eps_cbs0, beta_bs_t)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'eps_cds0, beta_RH, beta_ds_t_ts, expected',
    [
        (0.0004, 0.7, 0.6, 0.000168),
        (0.0005, 0.8, 0.5, 0.0002),
    ],
)
def test_calculate_drying_shrinkage_strain(
    eps_cds0, beta_RH, beta_ds_t_ts, expected
):
    """Test eps_cds_t_ts."""
    result = _concrete_material_properties.eps_cds_t_ts(
        eps_cds0, beta_RH, beta_ds_t_ts
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'f_cm, alpha_bs, expected',
    [
        (30, 700, -0.000449),
        (25, 800, -0.000375),
    ],
)
def test_eps_cbs_0(f_cm, alpha_bs, expected):
    """Test eps_cbs_0."""
    result = _concrete_material_properties.eps_cbs_0(f_cm, alpha_bs)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    't, expected',
    [
        (28, 0.653),
        (180, 0.932),
    ],
)
def test_beta_bs_t(t, expected):
    """Test beta_bs_t."""
    result = _concrete_material_properties.beta_bs_t(t)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'strength_class, expected',
    [
        ('CS', 800),
        ('CN', 700),
        ('CR', 600),
        ('cs', 800),
        ('Cr', 600),
    ],
)
def test_alpha_bs(strength_class, expected):
    """Test alpha_bs."""
    result = _concrete_material_properties.alpha_bs(strength_class)
    assert result == expected


@pytest.mark.parametrize(
    'strength_class, expected',
    [
        ('CS', 3),
        ('CN', 4),
        ('CR', 6),
        ('cs', 3),
        ('Cr', 6),
    ],
)
def test_alpha_ds_1(strength_class, expected):
    """Test alpha_ds_1."""
    result = _concrete_material_properties.alpha_ds_1(strength_class)
    assert result == expected


@pytest.mark.parametrize(
    'strength_class, expected',
    [
        ('CS', 0.013),
        ('CN', 0.012),
        ('CR', 0.012),
        ('cs', 0.013),
        ('Cr', 0.012),
    ],
)
def test_alpha_ds_2(strength_class, expected):
    """Test alpha_ds_2."""
    result = _concrete_material_properties.alpha_ds_2(strength_class)
    assert result == expected


@pytest.mark.parametrize(
    'f_cm, alpha_ds1, alpha_ds2, expected',
    [
        (30, 4, 0.012, 0.00460),
        (25, 3, 0.013, 0.00397),
    ],
)
def test_eps_cds_0_fcm(f_cm, alpha_ds1, alpha_ds2, expected):
    """Test eps_cds_0_fcm."""
    result = _concrete_material_properties.eps_cds_0_fcm(
        f_cm, alpha_ds1, alpha_ds2
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'RH, RH_eq, expected',
    [
        (80, 99, -0.7320),
        (98, 99, -0.0464),
        (100, 99, 0.2815),
    ],
)
def test_beta_RH_s(RH, RH_eq, expected):
    """Test beta_RH_s."""
    result = _concrete_material_properties.beta_RH_s(RH, RH_eq)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    't, ts, h, expected',
    [
        (180, 28, 150, 0.8116),
        (365, 0, 200, 0.8502),
    ],
)
def test_beta_ds_t_ts(t, ts, h, expected):
    """Test beta_ds_t_ts."""
    result = _concrete_material_properties.beta_ds_t_ts(t, ts, h)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'fcm, expected',
    [
        (25, 99),
        (30, 99),
        (50, 95.5311),
    ],
)
def test_RH_eq(fcm, expected):
    """Test RH_eq."""
    result = _concrete_material_properties.RH_eq(fcm)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'eps_cs, concrete_grade, expected',
    [
        (0.0005, 'LC8', 0.00075),
        (0.0005, 'LC12', 0.00075),
        (0.0005, 'LC16', 0.00075),
        (0.0005, 'LC20', 0.0006),
        (0.0005, 'LC25', 0.0006),
    ],
)
def test_eps_lcs_t_ts(eps_cs, concrete_grade, expected):
    """Test eps_lcs_t_ts."""
    result = _concrete_material_properties.eps_lcs_t_ts(eps_cs, concrete_grade)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'temperature_intervals, expected',
    [
        (
            [(5, 20), (10, 15), (7, 25)],
            21.651,
        ),
        (
            [(1, 0), (2, 5), (1, -5)],
            1.602,
        ),
        (
            [(10, 30)],
            15.662,
        ),
    ],
)
def test_t_T_maturity(temperature_intervals, expected):
    """Test t_T_maturity."""
    result = _concrete_material_properties.t_T_maturity(temperature_intervals)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'delta_T, concrete_type, expected',
    [
        (20, 10e-6, 0.0002),
        (15, 8e-6, 0.00012),
        (10, 10e-6, 0.0001),
    ],
)
def test_eps_c_T(delta_T, concrete_type, expected):
    """Test eps_c_T."""
    result = _concrete_material_properties.eps_c_T(delta_T, concrete_type)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'concrete_type, expected',
    [
        ('normal', 10e-6),
        ('lightweight', 8e-6),
    ],
)
def test_alpha_T(concrete_type, expected):
    """Test alpha_T."""
    result = _concrete_material_properties.alpha_T(concrete_type)
    assert result == expected


@pytest.mark.parametrize(
    'f_cm, T, expected',
    [
        (30, 40, 28.2),
        (35, 60, 30.8),
    ],
)
def test_fcm_T(f_cm, T, expected):
    """Test fcm_T."""
    result = _concrete_material_properties.fcm_T(f_cm, T)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'f_lcm, T, expected',
    [
        (25, 40, 24.0),
        (20, 60, 18.4),
    ],
)
def test_flcm_T(f_lcm, T, expected):
    """Test flcm_T."""
    result = _concrete_material_properties.flcm_T(f_lcm, T)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'f_cm_ref, D_fc, S_fc, expected',
    [
        (30, 0.01, 0.04, 31.5),
    ],
)
def test_fcm_hT_T(f_cm_ref, D_fc, S_fc, expected):
    """Test fcm_hT_T."""
    result = _concrete_material_properties.fcm_hT_T(f_cm_ref, D_fc, S_fc)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'h_T, expected',
    [
        (0.8, 0.319),
    ],
)
def test_calculate_D_fc(h_T, expected):
    """Test D_fc."""
    result = _concrete_material_properties.D_fc(h_T)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'h_ref, K_hT, T, T_ref, expected',
    [
        (1.0, 0.01, 25, 20, 0.95),
    ],
)
def test_hT_Tref(h_ref, K_hT, T, T_ref, expected):
    """Test hT_Tref."""
    result = _concrete_material_properties.hT_Tref(h_ref, K_hT, T, T_ref)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'h_ref, expected',
    [
        (0.05, 0.00053),
    ],
)
def test_K_hT(h_ref, expected):
    """Test K_hT."""
    result = _concrete_material_properties.K_hT(h_ref)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'T, alpha_f, expected',
    [
        (50, 0.2, 0.0915),
    ],
)
def test_S_fc(T, alpha_f, expected):
    """Test S_fc."""
    result = _concrete_material_properties.S_fc(T, alpha_f)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'h_ref, expected',
    [
        (1.0, -0.589),
    ],
)
def test_calculate_alpha_f(h_ref, expected):
    """Test alpha_fc."""
    result = _concrete_material_properties.alpha_fc(h_ref)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'f_ctm, T, expected',
    [
        (3.0, 40, 2.52),
        (4.0, 60, 2.72),
    ],
)
def test_fctm_T(f_ctm, T, expected):
    """Test fctm_T."""
    result = _concrete_material_properties.fctm_T(f_ctm, T)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'G_F, T, concrete_type, expected',
    [
        (150, 50, 'dry', 136.5),
        (180, 60, 'mass', 136.8),
    ],
)
def test_GF_T(G_F, T, concrete_type, expected):
    """Test GF_T."""
    result = _concrete_material_properties.GF_T(G_F, T, concrete_type)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'f_ctm_ref, D_fct, S_fct, expected',
    [
        (3.0, 0.02, 0.05, 3.21),
    ],
)
def test_fctm_hT_T(f_ctm_ref, D_fct, S_fct, expected):
    """Test fctm_hT_T."""
    result = _concrete_material_properties.fctm_hT_T(f_ctm_ref, D_fct, S_fct)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'T, h_T, expected',
    [
        (50, 0.8, 0.2598),
    ],
)
def test_D_fct(T, h_T, expected):
    """Test D_fct."""
    result = _concrete_material_properties.D_fct(T, h_T)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'T, alpha_f, expected',
    [
        (50, 0.2, 0.1005),
    ],
)
def test_S_fct(T, alpha_f, expected):
    """Test S_fct."""
    result = _concrete_material_properties.S_fct(T, alpha_f)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'h_ref, expected',
    [
        (1.0, -0.7900),
    ],
)
def test_alpha_fct(h_ref, expected):
    """Test alpha_fct."""
    result = _concrete_material_properties.alpha_fct(h_ref)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'E_cm_ref, D_Ec, S_Ec, expected',
    [
        (35000, -0.01, -0.02, 33950.0),
    ],
)
def test_E_cm_hT_T(E_cm_ref, D_Ec, S_Ec, expected):
    """Test E_cm_hT_T."""
    result = _concrete_material_properties.E_cm_hT_T(E_cm_ref, D_Ec, S_Ec)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'T, h_T, expected',
    [
        (50, 0.8, -0.1417),
    ],
)
def test_calculate_D_Ec(T, h_T, expected):
    """Test D_Ec."""
    result = _concrete_material_properties.D_Ec(T, h_T)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'h_ref, expected',
    [
        (1.0, -0.2163),
    ],
)
def test_alpha_Ec(h_ref, expected):
    """Test alpha_Ec."""
    result = _concrete_material_properties.alpha_Ec(h_ref)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'T, alpha_Ec, expected',
    [
        (50, 0.2, 0.06),
    ],
)
def test_calculate_S_Ec(T, alpha_Ec, expected):
    """Test S_Ec."""
    result = _concrete_material_properties.S_Ec(T, alpha_Ec)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'beta_h, T, expected',
    [
        (1.0, 50, 112.029),
    ],
)
def test_beta_h_T(beta_h, T, expected):
    """Test beta_h_T."""
    result = _concrete_material_properties.beta_h_T(beta_h, T)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'phi_bc, T, expected',
    [
        (1.5, 50, 2.3524),
    ],
)
def test_phi_bc_T(phi_bc, T, expected):
    """Test phi_bc_T."""
    result = _concrete_material_properties.phi_bc_T(phi_bc, T)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'phi_dc, T, expected',
    [
        (0.5, 50, 0.8580),
    ],
)
def test_phi_dc_T(phi_dc, T, expected):
    """Test phi_dc_T."""
    result = _concrete_material_properties.phi_dc_T(phi_dc, T)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'phi_t_t0, T, expected',
    [
        (0.004, 50, 0.3640),
    ],
)
def test_phi_t_t0_T(phi_t_t0, T, expected):
    """Test phi_t_t0_T."""
    result = _concrete_material_properties.phi_t_t0_T(phi_t_t0, T)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'T, h, expected',
    [
        (50, 1.0, 0.00578),
    ],
)
def test_alpha_sT(T, h, expected):
    """Test alpha_sT."""
    result = _concrete_material_properties.alpha_sT(T, h)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'beta_RH, beta_sT, expected',
    [
        (1.2, 1.1, 1.32),
    ],
)
def test__(beta_RH, beta_sT, expected):
    """Test beta_RH_T."""
    result = _concrete_material_properties.beta_RH_T(beta_RH, beta_sT)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'RH, T, expected',
    [
        (50, 50, 1.0566),
    ],
)
def test_beta_sT(RH, T, expected):
    """Test beta_sT."""
    result = _concrete_material_properties.beta_sT(RH, T)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'RH, RH_T, expected',
    [
        (70, 98, -1.018),
        (98, 95, 0.25),
    ],
)
def test_beta_RH_(RH, RH_T, expected):
    """Test beta_RH."""
    result = _concrete_material_properties.beta_RH(RH, RH_T)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'fcm, T, expected',
    [
        (25, 60, 94.904),
        (40, 70, 89.687),
    ],
)
def test_RH_T(fcm, T, expected):
    """Test RH_T."""
    result = _concrete_material_properties.RH_T(fcm, T)
    assert result == pytest.approx(expected, rel=1e-2)
