"""Tests for _section5_materials module"""

import math

import pytest

from structuralcodes.codes.ec2_2023 import _section5_materials


@pytest.mark.parametrize(
    "fck, delta, expected",
    [(12, 8, 20), (-16, 8, 24), (25, -10, 35), (-55, -15, 70)],
)
def test_fck(fck, delta, expected):
    """Test the fck function"""
    assert math.isclose(
        _section5_materials.fcm(fck, delta_f=delta), expected, rel_tol=0.005
    )


@pytest.mark.parametrize(
    "test_input1, expected",
    [
        (12, 1.572),
        (16, 1.905),
        (20, 2.21),
        (25, 2.565),
        (30, 2.896),
        (35, 3.21),
        (40, 3.509),
        (45, 3.795),
        (50, 4.072),
        (55, 4.183),
        (60, 4.306),
        (70, 4.533),
        (80, 4.74),
        (90, 4.93),
        (100, 5.106),
    ],
)
def test_fctm(test_input1, expected):
    """Test the fctm function."""
    assert math.isclose(
        _section5_materials.fctm(test_input1), expected, rel_tol=0.05
    )


@pytest.mark.parametrize(
    "fctm, expected",
    [
        (1.6, 1.12),
        (-2.2, 1.54),
        (-4.1, 2.87),
    ],
)
def test_fctk_5(fctm, expected):
    """Test the fctk_5 function"""
    assert math.isclose(
        _section5_materials.fctk_5(fctm), expected, rel_tol=0.005
    )


@pytest.mark.parametrize(
    "fctm, expected",
    [
        (1.6, 2.08),
        (-2.2, 2.86),
        (-4.1, 5.33),
    ],
)
def test_fctm_95(fctm, expected):
    """Test the fctk_95 function"""
    assert math.isclose(
        _section5_materials.fctk_95(fctm), expected, rel_tol=0.005
    )


@pytest.mark.parametrize(
    "fcm, kE, expected",
    [
        (20, 9500, 25786.9673576516),
        (24, 5000, 14422.4957030741),
        (43, 13000, 45544.1747850274),
        (53, 10000, 37562.8575422107),
        (88, 6000, 26687.7610868318),
    ],
)
def test_Ecm(fcm, kE, expected):
    """Test the Ecm function"""
    assert math.isclose(
        _section5_materials.Ecm(fcm, kE=kE), expected, rel_tol=10e-5
    )


@pytest.mark.parametrize(
    "fcm, kE",
    [
        (-20, 9500),
        (24, -5000),
        (43, 50000),
        (53, 3000),
    ],
)
def test_Ecm_raises_errors(fcm, kE):
    """Test the Ecm function raises errors"""
    with pytest.raises(ValueError):
        _section5_materials.Ecm(fcm, kE)


@pytest.mark.parametrize(
    "Ac, u, expected",
    [
        (400, 20, 40),
        (500, 30, 33.3333333333333),
        (200, 15, 26.6666666666667),
        (100, 40, 5),
    ],
)
def test_hn(Ac, u, expected):
    """Test the hn function"""
    assert math.isclose(_section5_materials.hn(Ac, u), expected, rel_tol=10e-5)


@pytest.mark.parametrize("Ac, u", [(-40, 20), (40, -20), (-40, -20)])
def test_hn_raises_errors(Ac, u):
    """Test the hn function raises errors"""
    with pytest.raises(ValueError):
        _section5_materials.hn(Ac, u)


@pytest.mark.parametrize(
    "hn, atm_conditions, expected",
    [
        (100, "dry", 0.82),
        (150, "dry", 0.805),
        (800, "dry", 0.732),
        (200, "humid", 0.68),
        (400, "humid", 0.66666667),
        (800, "humid", 0.648),
    ],
)
def test_A_phi_correction_exp(hn, atm_conditions, expected):
    """Test the A_phi_correction_exp function"""
    assert math.isclose(
        _section5_materials.A_phi_correction_exp(hn, atm_conditions),
        expected,
        rel_tol=10e-5,
    )


@pytest.mark.parametrize(
    "hn, atm_conditions", [(100, "1231"), (80, "dry"), (1100, "humid")]
)
def test_A_phi_correction_exp_raises_errors(hn, atm_conditions):
    """Test A_phi_correction_exp raises errors"""
    with pytest.raises(ValueError):
        _section5_materials.A_phi_correction_exp(hn, atm_conditions)


@pytest.mark.parametrize(
    "t0, atm_conditions, _hn, concrete_class, expected",
    [
        (10, "dry", 500, "CS", 2.5),
        (28, "humid", 200, "CN", 1.6),
        (91, "dry", 750, "CR", 1.45),
        (60, "humid", 600, "CS", 1.41016),
    ],
)
def test_phi_50y_t0(t0, atm_conditions, _hn, concrete_class, expected):
    """Test phi_50y_t0 function"""
    assert math.isclose(
        _section5_materials.phi_50y_t0(
            t0, atm_conditions, _hn, concrete_class
        ),
        expected,
        rel_tol=10e-5,
    )


@pytest.mark.parametrize(
    "t0, atm_conditions, _hn, concrete_class",
    [
        (-1, "dry", 500, "CS"),
        (50, "ASDF", 500, "CS"),
        (50, "dry", 50, "CS"),
        (50, "dry", 1500, "CS"),
        (50, "dry", 500, "ASD"),
        (1, "dry", 100, "CS"),
    ],
)
def test_phi_50y_t0_raises_errors(t0, atm_conditions, _hn, concrete_class):
    """Test phi_50y_t0 raises errors"""
    with pytest.raises(ValueError):
        _section5_materials.phi_50y_t0(t0, atm_conditions, _hn, concrete_class)


@pytest.mark.parametrize(
    "fck_28, atm_conditions, _hn, concrete_class, expected",
    [
        (35, "dry", 500, "CS", 0.45),
        (50, "humid", 1000, "CN", 0.23),
        (80, "dry", 200, "CR", 0.54),
        (40, "humid", 300, "CS", 0.29),
        (25, "dry", 800, "CN", 0.46333),
    ],
)
def test_eps_cs_50y(fck_28, atm_conditions, _hn, concrete_class, expected):
    """Test phi_50y_t0 function"""
    assert math.isclose(
        _section5_materials.eps_cs_50y(
            fck_28, atm_conditions, _hn, concrete_class
        ),
        expected,
        rel_tol=10e-5,
    )


@pytest.mark.parametrize(
    "fck_28, atm_conditions, _hn, concrete_class",
    [
        (15, "dry", 500, "CS"),
        (90, "dry", 500, "CS"),
        (50, "ASDF", 500, "CS"),
        (50, "dry", 50, "CS"),
        (50, "dry", 1500, "CS"),
        (50, "dry", 500, "ASD"),
        (1, "dry", 100, "CS"),
    ],
)
def test_eps_cs_50y_raises_errors(fck_28, atm_conditions, _hn, concrete_class):
    """Test eps_cs_50y raises errors"""
    with pytest.raises(ValueError):
        _section5_materials.eps_cs_50y(
            fck_28, atm_conditions, _hn, concrete_class
        )


@pytest.mark.parametrize(
    "fck, fck_ref, expected",
    [(60, 40, 0.873580464736299), (40, 45, 1), (60, 50, 0.941036028881029)],
)
def test_eta_cc(fck, fck_ref, expected):
    """Test eta_cc function"""
    assert math.isclose(
        _section5_materials.eta_cc(fck, fck_ref),
        expected,
        rel_tol=10e-5,
    )


@pytest.mark.parametrize("fck, fck_ref", [(-10, 40), (0, 20), (30, -10)])
def test_eta_cc_raises_errors(fck, fck_ref):
    """Test eta_cc raises errors"""
    with pytest.raises(ValueError):
        _section5_materials.eta_cc(fck, fck_ref)


@pytest.mark.parametrize(
    "t_ref, t0, concrete_class, expected",
    [
        (20, 40, "CR", 0.85),
        (20, 40, "cr", 0.85),
        (20, 40, "rapid", 0.85),
        (27, 180, "CR", 1),
        (57, 180, "CS", 0.85),
        (55, 180, "CS", 1),
        (55, 180, "slow", 1),
    ],
)
def test_k_tc(t_ref, t0, concrete_class, expected):
    """Test k_tc function"""
    assert math.isclose(
        _section5_materials.k_tc(t_ref, t0, concrete_class),
        expected,
        rel_tol=10e-5,
    )


@pytest.mark.parametrize(
    "t_ref, t0, concrete_class",
    [(-3, 20, "CS"), (10, -5, "CS"), (10, 10, "aadsf")],
)
def test_k_tc_raises_errors(t_ref, t0, concrete_class):
    """Test k_tc taises errors"""
    with pytest.raises(ValueError):
        _section5_materials.k_tc(t_ref, t0, concrete_class)


@pytest.mark.parametrize(
    "fck, eta_cc, k_tc, gamma_C, expected",
    [
        (40, 0.9, 0.85, 1.35, 22.6666666666667),
        (35, 1, 1, 1.5, 23.3333333333333),
        (50, 0.85, 0.85, 1.4, 25.8035714285714),
    ],
)
def test_fcd(fck, eta_cc, k_tc, gamma_C, expected):
    """Test fcd function"""
    assert math.isclose(
        _section5_materials.fcd(fck, eta_cc, k_tc, gamma_C),
        expected,
        rel_tol=10e-5,
    )


@pytest.mark.parametrize(
    "fck, eta_cc, k_tc, gamma_C",
    [
        (-10, 0.5, 0.85, 1.5),
        (40, -1, 0.85, 1.5),
        (40, 2, 1, 1.5),
        (40, 0.9, 1.0, -2),
    ],
)
def test_fcd_raises_errors(fck, eta_cc, k_tc, gamma_C):
    """Test fcd function raises errors"""
    with pytest.raises(ValueError):
        _section5_materials.fcd(fck, eta_cc, k_tc, gamma_C)


@pytest.mark.parametrize(
    "t_ref, concrete_class, expected",
    [(25, "CR", 0.8), (30, "CR", 0.7), (48, "cs", 0.8), (70, "CS", 0.7)],
)
def test_k_tt(t_ref, concrete_class, expected):
    """Test k_tt function"""
    assert math.isclose(
        _section5_materials.k_tt(t_ref, concrete_class),
        expected,
        rel_tol=10e-5,
    )


@pytest.mark.parametrize("t_ref, concrete_class", [(-10, "CR"), (20, "ADSF")])
def test_k_tt_raises_errors(t_ref, concrete_class):
    """Test k_tt raises errors"""
    with pytest.raises(ValueError):
        _section5_materials.k_tt(t_ref, concrete_class)


@pytest.mark.parametrize(
    "fctk_5, k_tt, gamma_C", [(-20, 0.8, 1.5), (40, 0.7, -1)]
)
def test_fctd_raises_errors(fctk_5, k_tt, gamma_C):
    """Test fctd raises errors"""
    with pytest.raises(ValueError):
        _section5_materials.fctd(fctk_5, k_tt, gamma_C)


@pytest.mark.parametrize(
    "fctk_5, k_tt, gamma_C, expected",
    [
        (2, 0.7, 1.5, 0.933333333333333),
        (3, 0.8, 1.25, 1.92),
        (1.8, 0.7, 1.3, 0.969230769230769),
    ],
)
def test_fctd(fctk_5, k_tt, gamma_C, expected):
    """Test fctd"""
    assert math.isclose(
        _section5_materials.fctd(fctk_5, k_tt, gamma_C),
        expected,
        rel_tol=10e-5,
    )


@pytest.mark.parametrize(
    "fcm, expected",
    [(30, 0.0021750627541677), (50, 0.00257882204904827), (80, 0.0028)],
)
def test_eps_c1(fcm, expected):
    """Test eps_c1"""
    assert math.isclose(
        _section5_materials.eps_c1(fcm),
        expected,
        rel_tol=10e-5,
    )


@pytest.mark.parametrize("fcm", [(-20)])
def test_eps_c1_raises_errors(fcm):
    """Test eps_c1 raises errors"""
    with pytest.raises(ValueError):
        _section5_materials.eps_c1(fcm)


@pytest.mark.parametrize(
    "fcm, expected",
    [(70, 0.00301456920899968), (50, 0.0035), (80, 0.00286325067128806)],
)
def test_eps_cu1(fcm, expected):
    """Test eps_cu1"""
    assert math.isclose(
        _section5_materials.eps_cu1(fcm),
        expected,
        rel_tol=10e-5,
    )


@pytest.mark.parametrize("fcm", [(-20)])
def test_eps_cu1_raises_errors(fcm):
    """Test eps_cu1 raises errors"""
    with pytest.raises(ValueError):
        _section5_materials.eps_cu1(fcm)


@pytest.mark.parametrize(
    "Ecm, fcm, eps_c, eps_c1, eps_cu1, expected",
    [
        (37562, 50, 0.0025, 0.00257882204904827, 0.0035, 49.9547867728839),
        (25000, 80, 0.002, 0.0028, 0.00286325067128806, 51.3165266106443),
    ],
)
def test_sigma_c(Ecm, fcm, eps_c, eps_c1, eps_cu1, expected):
    """Test sigma_c function"""
    assert math.isclose(
        _section5_materials.sigma_c(Ecm, fcm, eps_c, eps_c1, eps_cu1),
        expected,
        rel_tol=10e-5,
    )


@pytest.mark.parametrize(
    "Ecm, fcm, eps_c, eps_c1, eps_cu1",
    [
        (-37562, 50, 0.0025, 0.00257882204904827, 0.0035),
        (37562, -50, 0.0025, 0.00257882204904827, 0.0035),
        (37562, 50, -0.0025, 0.00257882204904827, 0.0035),
        (37562, 50, 0.0025, 0, 0.0035),
        (37562, 50, 0.0025, -0.00257882204904827, 0.0035),
        (37562, 50, 0.0025, 0.00257882204904827, -0.0035),
        (37562, 50, 0.02, 0.00257882204904827, 0.0035),
    ],
)
def test_sigma_c_raises_errors(Ecm, fcm, eps_c, eps_c1, eps_cu1):
    """Test sigma_c_raises_errors"""
    with pytest.raises(ValueError):
        _section5_materials.sigma_c(Ecm, fcm, eps_c, eps_c1, eps_cu1)


@pytest.mark.parametrize("concrete_type, expected", [("nc", 25), ("npc", 24)])
def test_concrete_mean_unit_weight(concrete_type, expected):
    """Test concrete_mean_weight function"""
    assert math.isclose(
        _section5_materials.weight_c(concrete_type),
        expected,
        rel_tol=10e-5,
    )


@pytest.mark.parametrize(
    "concrete_type",
    [
        ("ADSF"),
    ],
)
def test_weight_c_raises_errors(concrete_type):
    """Test weight_c raises errors"""
    with pytest.raises(ValueError):
        _section5_materials.weight_c(concrete_type)


def test_alpha_c_th():
    """Test alpha_c_th function"""
    assert math.isclose(
        _section5_materials.alpha_c_th(), 10 * 10e-6, rel_tol=10e-5
    )


@pytest.mark.parametrize(
    "ductility_class, expected_k, expected_eps_uk",
    [("A", 1.05, 0.025), ("B", 1.08, 0.05), ("C", 1.25, 0.075)],
)
def test_r_steel_stress_strain_params(
    ductility_class, expected_k, expected_eps_uk
):
    """Test r_steel_stress_strain retuns proper values"""
    params = _section5_materials.r_steel_stress_strain_params(ductility_class)
    assert math.isclose(params[0], expected_k, rel_tol=10e-5)
    assert math.isclose(params[1], expected_eps_uk, rel_tol=10e-5)


def test_r_steel_stress_strain_params_raises_errors():
    """Test r_steel_stress_strain_params raises errors"""
    with pytest.raises(ValueError):
        _section5_materials.r_steel_stress_strain_params("asdf")


def test_Es():
    """Test Es function"""
    assert math.isclose(_section5_materials.Es(), 200000, rel_tol=10e-5)


def test_alpha_s_th():
    """Test alpha_s_th function"""
    assert math.isclose(
        _section5_materials.alpha_s_th(), 10 * 10e-6, rel_tol=10e-5
    )


def test_weight_s():
    """Test weight_s function"""
    assert math.isclose(_section5_materials.weight_s(), 78.5, rel_tol=10e-5)


@pytest.mark.parametrize(
    "fyk, gamma_S, expected",
    [(500, 1.15, 434.782609), (600, 1.15, 521.739130), (400, 1, 400)],
)
def test_fyd(fyk, gamma_S, expected):
    """Test fyd function"""
    assert math.isclose(
        _section5_materials.fyd(fyk, gamma_S), expected, rel_tol=10e-5
    )


@pytest.mark.parametrize("fyk, gamma_S", [(-300, 1.15), (600, -3)])
def test_fyd_raises_errors(fyk, gamma_S):
    """Test fyd function raises errors"""
    with pytest.raises(ValueError):
        _section5_materials.fyd(fyk, gamma_S)


@pytest.mark.parametrize(
    "eps_uk, gamma_S, expected",
    [
        (0.025, 1.15, 0.02174),
        (0.05, 1, 0.05),
    ],
)
def test_eps_ud(eps_uk, gamma_S, expected):
    """Test eps_ud function"""
    assert math.isclose(
        _section5_materials.eps_ud(eps_uk, gamma_S), expected, rel_tol=10e-5
    )


@pytest.mark.parametrize("eps_uk, gamma_S", [(-0.025, 1.15), (0.025, -1)])
def test_eps_ud_raises_errors(eps_uk, gamma_S):
    """Test eps_ud raises errors"""
    with pytest.raises(ValueError):
        _section5_materials.eps_ud(eps_uk, gamma_S)


@pytest.mark.parametrize(
    "eps, fy, k, eps_u, expected",
    [
        (0.015, 500, 1.05, 0.025, 513.888888888889),
        (0.043, 521, 1, 0.0435, 521),
        (0.001, 400, 1.25, 0.075, 200),
        (0.043, 521, 1.08, 0.0435, 562.170402249664),
    ],
)
def test_sigma_s(eps, fy, k, eps_u, expected):
    """Test sigma_s function"""
    assert math.isclose(
        _section5_materials.sigma_s(eps, fy, k, eps_u), expected, rel_tol=10e-5
    )


@pytest.mark.parametrize(
    "fp01k, gamma_S, expected",
    [(1560, 1.15, 1356.521739), (1740, 1, 1740), (950, 1.15, 826.086956)],
)
def test_fpd(fp01k, gamma_S, expected):
    """Test fpd function"""
    assert math.isclose(
        _section5_materials.fpd(fp01k, gamma_S), expected, rel_tol=10e-5
    )


@pytest.mark.parametrize(
    "prestress_class, element, expected",
    [
        ("Y1770", "W", (1550, 1770)),
        ("Y1960", "S", (1740, 1960)),
        ("Y1030", "b", (835, 1030)),
    ],
)
def test_p_steel_strain_params(prestress_class, element, expected):
    """Test p_steel_strain_params"""
    assert (
        _section5_materials.p_steel_stress_params(prestress_class, element)
        == expected
    )


@pytest.mark.parametrize(
    "prestress_class, element", [("Y1770", "o"), ("ASDF", "S")]
)
def test_p_steel_strain_params_raises_errors(prestress_class, element):
    """Test p_steel_strain_params raises errors"""
    with pytest.raises(ValueError):
        _section5_materials.p_steel_stress_params(prestress_class, element)


@pytest.mark.parametrize("fp01k, gamma_S", [(-1560, 1), (1560, 0), (1560, -1)])
def test_fpd_raises_errors(fp01k, gamma_S):
    """Test fpd raises errors"""
    with pytest.raises(ValueError):
        _section5_materials.fpd(fp01k, gamma_S)


@pytest.mark.parametrize(
    "eps, fpy, fpu, eps_u, Ep, expected",
    [
        (0.03, 1380, 1570, 0.035, 200000, 1536.19217081851),
        (0.004, 1740, 1960, 0.0304347826086957, 200000, 800),
        (0.02, 1650, 1860, 0.035, 200000, 1742.24299065421),
        (0.032, 900, 900, 0.0304347826086957, 200000, 900),
    ],
)
def test_steel_p(eps, fpy, fpu, eps_u, Ep, expected):
    """Test steel_p function"""
    assert math.isclose(
        _section5_materials.sigma_p(eps, fpy, fpu, eps_u, Ep),
        expected,
        rel_tol=10e-5,
    )


@pytest.mark.parametrize(
    "eps, fpy, fpu, eps_u, Ep",
    [
        (-0.03, 1380, 1570, 0.035, 200000),
        (0.035, 1380, 1570, 0.03, 200000),
        (0.02, -1380, 1570, 0.035, 200000),
        (0.02, 1380, 1200, 0.025, 200000),
        (0.02, 1380, -1570, 0.025, 200000),
        (0.02, 1380, 1570, -0.025, 200000),
        (0.02, 1380, 1570, 0.025, -200000),
    ],
)
def test_steel_p_raises_errors(eps, fpy, fpu, eps_u, Ep):
    """Test steel_p raises errors"""
    with pytest.raises(ValueError):
        _section5_materials.sigma_p(eps, fpy, fpu, eps_u, Ep)
