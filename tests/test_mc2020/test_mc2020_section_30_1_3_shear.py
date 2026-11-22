"""Tests for the Section 30.1.3 of fib Model Code 2020."""

import math

import pytest

from structuralcodes.codes.mc2020 import _section_30_1_3_shear


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                           Section 30.1.3.1                          #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
@pytest.mark.parametrize(
    'V_Rd, V_Ed, expected',
    [
        (500, 400, 0.8),
        (500, 600, 1.2),
        (500, 500, 1.0),
        (500, -250, 0.5),
    ],
)
def test_UC_V(V_Rd, V_Ed, expected):
    """Test function unity_check_V."""
    assert math.isclose(
        _section_30_1_3_shear.unity_check_V(V_Rd, V_Ed), expected, rel_tol=1e-9
    )


@pytest.mark.parametrize(
    'a_v, d, expected',
    [
        (2500, 1000, 1.0),
        (1500, 1000, 0.75),
        (1200, 1000, 0.60),
        (200, 1000, 0.50),
    ],
)
def test_k_dir(a_v, d, expected):
    """Test function k_dir."""
    assert math.isclose(
        _section_30_1_3_shear.k_dir(a_v, d), expected, rel_tol=1e-9
    )


@pytest.mark.parametrize(
    'z_s, z_p, A_s, A_p, h, expected',
    [
        (800, 750, 500, 250, 1000, 784.042553),
        (800, 600, 500, 250, 1000, 745.4545455),
        (800, 750, 600, 250, 1000, 785.9550562),
        (800, 750, 500, 400, 1000, 778.5714286),
        (800, 750, 500, 250, 800, 784.042553),
        (600, 200, 500, 250, 1000, 720),  # check max. condition
    ],
)
def test_z_v(z_s, z_p, A_s, A_p, h, expected):
    """Test function z_v."""
    assert math.isclose(
        _section_30_1_3_shear.z_v(z_s, z_p, A_s, A_p, h),
        expected,
        rel_tol=1e-9,
    )


@pytest.mark.parametrize(
    'M_Ed, V_Ed, N_Ed, E_s, A_s, z_v, cot_theta, delta_e, E_c, A_c_ten, '
    'expected',
    [
        (
            2.0e03,
            1.0e03,
            500,
            2.0e09,
            1.0e-03,
            0.85,
            2,
            0.2,
            3.0e07,
            0.04,
            8.71323529e-04,
        ),
        (
            2.0e03,
            1.2e03,
            500,
            2.0e09,
            1.0e-03,
            0.85,
            2,
            0.2,
            3.0e07,
            0.04,
            9.21323529e-04,
        ),
        (
            2.0e03,
            1.0e03,
            -3.0e004,
            2.0e09,
            1.0e-03,
            0.85,
            2,
            0.2,
            3.0e07,
            0.04,
            0,
        ),
        (
            2.0e03,
            1.0e03,
            500,
            2.5e09,
            1.0e-03,
            0.85,
            2,
            0.2,
            3.0e07,
            0.04,
            6.97058824e-04,
        ),
        (
            2.0e03,
            1.0e03,
            500,
            2.0e09,
            2.0e-03,
            0.85,
            2,
            0.2,
            3.0e07,
            0.04,
            4.35661765e-04,
        ),
        (
            2.0e03,
            1.0e03,
            500,
            2.0e09,
            1.0e-03,
            0.7,
            2,
            0.2,
            3.0e07,
            0.04,
            9.91071429e-04,
        ),
        (
            2.0e03,
            1.0e03,
            500,
            2.0e09,
            1.0e-03,
            0.85,
            1.5,
            0.2,
            3.0e07,
            0.04,
            8.08823529e-04,
        ),
        (
            2.0e03,
            1.0e03,
            500,
            2.0e09,
            1.0e-03,
            0.85,
            2,
            0.15,
            3.0e07,
            0.04,
            8.78676471e-04,
        ),
        (
            2.0e03,
            1.0e03,
            500,
            2.0e09,
            1.0e-03,
            0.85,
            2,
            0.2,
            1.0e07,
            0.04,
            8.71323529e-04,
        ),
        (
            2.0e03,
            1.0e03,
            500,
            2.0e09,
            1.0e-03,
            0.85,
            2,
            0.2,
            3.0e07,
            0.08,
            8.71323529e-04,
        ),
    ],
)
def test_eps_x(
    M_Ed,
    V_Ed,
    N_Ed,
    E_s,
    A_s,
    z_v,
    cot_theta,
    delta_e,
    E_c,
    A_c_ten,
    expected,
):
    """Test function eps_x."""
    assert math.isclose(
        _section_30_1_3_shear.eps_x(
            M_Ed,
            V_Ed,
            N_Ed,
            E_s,
            A_s,
            z_v,
            cot_theta,
            delta_e,
            E_c,
            A_c_ten,
        ),
        expected,
        rel_tol=1e-9,
    )


@pytest.mark.parametrize(
    'M_Ed, V_Ed, N_Ed, E_s, E_p, A_s, A_p, z_v, z_s, z_p, cot_theta, e_p, '
    'expected',
    [
        (
            2.0e03,
            1.0e03,
            500,
            2.0e09,
            1.95e09,
            1.0e-03,
            5.0e-04,
            0.85,
            0.9,
            0.9,
            2,
            0.45,
            5.74229692e-04,
        ),
        (
            2.0e03,
            1.2e03,
            500,
            2.0e09,
            1.95e09,
            1.0e-03,
            5.0e-04,
            0.85,
            0.9,
            0.9,
            2,
            0.45,
            6.05975724e-04,
        ),
        (
            2.0e03,
            1.0e03,
            -3.0e4,
            2.0e09,
            1.95e09,
            1.0e-03,
            5.0e-04,
            0.85,
            0.9,
            0.9,
            2,
            0.45,
            0,
        ),
        (
            2.0e03,
            1.0e03,
            500,
            1.50e09,
            1.95e09,
            1.0e-03,
            5.0e-04,
            0.85,
            0.9,
            0.9,
            2,
            0.45,
            6.90235690e-04,
        ),
        (
            2.0e03,
            1.0e03,
            500,
            2.0e09,
            3.0e09,
            1.0e-03,
            5.0e-04,
            0.85,
            0.9,
            0.9,
            2,
            0.45,
            4.88095238e-04,
        ),
        (
            2.0e03,
            1.0e03,
            500,
            2.0e09,
            1.95e09,
            3.0e-03,
            5.0e-04,
            0.85,
            0.9,
            0.9,
            2,
            0.45,
            2.44922342e-04,
        ),
        (
            2.0e03,
            1.0e03,
            500,
            2.0e09,
            1.95e09,
            1.0e-03,
            5.0e-03,
            0.85,
            0.9,
            0.9,
            2,
            0.45,
            1.45390071e-04,
        ),
        (
            2.0e03,
            1.0e03,
            500,
            2.0e09,
            1.95e09,
            1.0e-03,
            5.0e-04,
            0.7,
            0.9,
            0.9,
            2,
            0.45,
            5.46218487e-04,
        ),
        (
            2.0e03,
            1.0e03,
            500,
            2.0e09,
            1.95e09,
            1.0e-03,
            5.0e-04,
            0.85,
            0.6,
            0.9,
            2,
            0.45,
            7.40072202e-04,
        ),
        (
            2.0e03,
            1.0e03,
            500,
            2.0e09,
            1.95e09,
            1.0e-03,
            5.0e-04,
            0.85,
            0.9,
            1.3,
            2,
            0.45,
            5.33822331e-04,
        ),
        (
            2.0e03,
            1.0e03,
            500,
            2.0e09,
            1.95e09,
            1.0e-03,
            5.0e-04,
            0.85,
            0.9,
            0.9,
            2.5,
            0.45,
            6.13912232e-04,
        ),
        (
            2.0e03,
            1.0e03,
            500,
            2.0e09,
            1.95e09,
            1.0e-03,
            5.0e-04,
            0.85,
            0.9,
            0.9,
            2,
            0.3,
            5.88235294e-04,
        ),
    ],
)
def test_eps_x(
    M_Ed,
    V_Ed,
    N_Ed,
    E_s,
    E_p,
    A_s,
    A_p,
    z_v,
    z_s,
    z_p,
    cot_theta,
    e_p,
    expected,
):
    """Test function eps_x_bond."""
    assert math.isclose(
        _section_30_1_3_shear.eps_x_bond(
            M_Ed,
            V_Ed,
            N_Ed,
            E_s,
            E_p,
            A_s,
            A_p,
            z_v,
            z_s,
            z_p,
            cot_theta,
            e_p,
        ),
        expected,
        rel_tol=1e-5,
    )


@pytest.mark.parametrize(
    'M_Ed0, M_Pd, expected',
    [
        (5000, -300, 4700),
        (5000, 300, 5300),
        (2000, -300, 1700),
    ],
)
def test_M_Ed(
    M_Ed0,
    M_Pd,
    expected,
):
    """Test function M_Ed."""
    assert math.isclose(
        _section_30_1_3_shear.M_Ed(
            M_Ed0,
            M_Pd,
        ),
        expected,
        rel_tol=1e-5,
    )


@pytest.mark.parametrize(
    'N_Ed0, F_p, delta_p, expected',
    [
        (-5_000, 1000, 0.523598776, -5866.025404),
        (5_000, 1000, 0.523598776, 4133.974596),  # Tensile N
        (-5_000, 5000, 0.523598776, -9330.127019),
        (-5_000, 1000, 0.785398163, -5707.106781),
    ],
)
def test_N_Ed(
    N_Ed0,
    F_p,
    delta_p,
    expected,
):
    """Test function N_Ed."""
    assert math.isclose(
        _section_30_1_3_shear.N_Ed(
            N_Ed0,
            F_p,
            delta_p,
        ),
        expected,
        rel_tol=1e-9,
    )


@pytest.mark.parametrize(
    'V_Ed0, F_p, delta_p, expected',
    [
        (5000, 1000, 0.523598776, 4500),
        (5000, 5000, 0.523598776, 2500),
        (5000, 1000, 0.785398163, 4292.893219),
    ],
)
def test_V_Ed(
    V_Ed0,
    F_p,
    delta_p,
    expected,
):
    """Test function V_Ed."""
    assert math.isclose(
        _section_30_1_3_shear.V_Ed(
            V_Ed0,
            F_p,
            delta_p,
        ),
        expected,
        rel_tol=1e-5,
    )


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                           Section 30.1.3.4                          #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
@pytest.mark.parametrize(
    'k_s, f_FTU, gamma_F, cot_theta, z_v, b_w, expected',
    [
        (0.64, 8, 1.6, 1.3458876, 700, 1000, 3014788),
        (0.64, 8, 1.3, 1.3458876, 700, 1000, 3710509),
        (0.64, 8, 1.6, 1.2348972, 700, 1000, 2766170),
        (0.64, 8, 1.6, 1.3458876, 600, 1000, 2584104),
        (0.64, 8, 1.6, 1.3458876, 700, 420, 1266211),
    ],
)
def test_V_Rd_FRC(
    k_s,
    f_FTU,
    gamma_F,
    cot_theta,
    z_v,
    b_w,
    expected,
):
    """Test function V_Rd_FRC."""
    assert math.isclose(
        _section_30_1_3_shear.V_Rd_FRC(
            k_s,
            f_FTU,
            gamma_F,
            cot_theta,
            z_v,
            b_w,
        ),
        expected,
        rel_tol=1e-5,
    )


@pytest.mark.parametrize(
    'f_yd, E_s, zeta, z_v, expected',
    [
        (435, 200000, 1.2, 500, 0.083523654),
        (435, 100000, 1.2, 500, 0.058581),
        (435, 200000, 0, 500, 0.121615),
        (435, 200000, 1.2, 800, 0.067863),
    ],
)
def test_k_v_FRC_LOA_1(
    f_yd,
    E_s,
    zeta,
    z_v,
    expected,
):
    """Test function k_v_FRC_LOA_1."""
    assert math.isclose(
        _section_30_1_3_shear.k_v_FRC_LOA_1(
            f_yd,
            E_s,
            zeta,
            z_v,
        ),
        expected,
        rel_tol=1e-5,
    )


@pytest.mark.parametrize(
    'f_yd, E_s, f_FTU, f_ck, expected',
    [
        (435, 200000, 8, 30, 1.27333333),
        (435, 100000, 8, 30, 0),
        (435, 200000, 7, 30, 1.11416667),
        (435, 200000, 8, 20, 1.91000000),
    ],
)
def test_zeta_LOA_1(
    f_yd,
    E_s,
    f_FTU,
    f_ck,
    expected,
):
    """Test function zeta_LOA_1."""
    assert math.isclose(
        _section_30_1_3_shear.zeta_LOA_1(
            f_yd,
            E_s,
            f_FTU,
            f_ck,
        ),
        expected,
        rel_tol=1e-6,
    )


@pytest.mark.parametrize(
    'f_yd, E_s, expected',
    [
        (435, 200000, 36.6125),
        (280, 200000, 33.9),
        (435, 100000, 44.225),
    ],
)
def test_theta_F_LOA_1(
    f_yd,
    E_s,
    expected,
):
    """Test function theta_F_LOA_1."""
    assert math.isclose(
        _section_30_1_3_shear.theta_F_LOA_1(
            f_yd,
            E_s,
        ),
        expected,
        rel_tol=1e-9,
    )


@pytest.mark.parametrize(
    'eta, gamma_c, rho_l, f_ck, d_dg, d, f_Ftud, tau_Rdc_min, expected',
    [
        (0.9, 1.5, 0.02, 30, 40, 400, 3, 0.5, 3.65416341),
        (0.6, 1.5, 0.02, 30, 40, 400, 3, 0.5, 3.43610894),
        (0.9, 1.1, 0.02, 30, 40, 400, 3, 0.5, 3.89204102),
        (0.9, 1.5, 0.01, 30, 40, 400, 3, 0.5, 3.51920985),
        (0.9, 1.5, 0.02, 55, 40, 400, 3, 0.5, 3.80063283),
        (0.9, 1.5, 0.02, 30, 20, 400, 3, 0.5, 3.51920985),
        (0.9, 1.5, 0.02, 30, 40, 1200, 3, 0.5, 3.45357158),
        (0.9, 1.5, 0.02, 30, 40, 400, 10, 0.5, 10.65416341),
        (0.9, 1.5, 0.02, 30, 40, 400, 3, 10, 12),
    ],
)
def test_tau_RD_cF_LOA_1(
    eta,
    gamma_c,
    rho_l,
    f_ck,
    d_dg,
    d,
    f_Ftud,
    tau_Rdc_min,
    expected,
):
    """Test function tau_RD_cF_LOA_1."""
    assert math.isclose(
        _section_30_1_3_shear.tau_RD_cF_LOA_1(
            eta,
            gamma_c,
            rho_l,
            f_ck,
            d_dg,
            d,
            f_Ftud,
            tau_Rdc_min,
        ),
        expected,
        rel_tol=1e-7,
    )


@pytest.mark.parametrize(
    'f_Ftuk, expected',
    [
        (1, 0.7),
        (2, 0.4),
        (0.2, 1),
    ],
)
def test_eta_LOA_1(
    f_Ftuk,
    expected,
):
    """Test function eta_LOA_1."""
    assert math.isclose(
        _section_30_1_3_shear.eta_LOA_1(
            f_Ftuk,
        ),
        expected,
        rel_tol=1e-9,
    )


@pytest.mark.parametrize(
    'f_R3k, f_R1k, expected',
    [
        (1, 1, 0.33),
        (2, 0.8, 0.816),
        (1, 50, 0),
    ],
)
def test_f_FTU_alt(
    f_R3k,
    f_R1k,
    expected,
):
    """Test function f_FTU_alt."""
    assert math.isclose(
        _section_30_1_3_shear.f_FTU_alt(
            f_R3k,
            f_R1k,
        ),
        expected,
        rel_tol=1e-9,
    )


@pytest.mark.parametrize(
    'eps_x, zeta, k_dg, z_v, expected',
    [
        (0.002, 1.3, 1, 600, 0.061320755),
        (0.002, 0.8, 1, 600, 0.067708333),
        (0.002, 1.3, 1.5, 600, 0.05163853),
        (0.002, 1.3, 1, 2000, 0.032704403),
    ],
)
def test_k_v_FRC_LOA_2(
    eps_x,
    zeta,
    k_dg,
    z_v,
    expected,
):
    """Test function k_v_FRC_LOA_2."""
    assert math.isclose(
        _section_30_1_3_shear.k_v_FRC_LOA_2(
            eps_x,
            zeta,
            k_dg,
            z_v,
        ),
        expected,
        rel_tol=1e-7,
    )


@pytest.mark.parametrize(
    'eps_x, f_FTU, f_ck, expected',
    [
        (0.0002, 1.8, 30, 1.032),
        (0.0006, 1.8, 30, 0.696),
        (0.0002, 1.2, 30, 0.688),
        (0.0002, 1.8, 44, 0.703636),
        (0.01, 1.8, 30, 0),
    ],
)
def test_zeta_LOA_2(
    eps_x,
    f_FTU,
    f_ck,
    expected,
):
    """Test function zeta_LOA_2."""
    assert math.isclose(
        _section_30_1_3_shear.zeta_LOA_2(
            eps_x,
            f_FTU,
            f_ck,
        ),
        expected,
        rel_tol=1e-6,
    )


@pytest.mark.parametrize(
    'eps_x, expected',
    [
        (0.003, 50),
        (0.0003, 31.1),
    ],
)
def test_theta_F_LOA_2_3(
    eps_x,
    expected,
):
    """Test function theta_F_LOA_2."""
    assert math.isclose(
        _section_30_1_3_shear.theta_F_LOA_2_3(
            eps_x,
        ),
        expected,
        rel_tol=1e-6,
    )


@pytest.mark.parametrize(
    'eps_x, k_dg, z_v, expected',
    [
        (0.003, 0.8, 600, 0.063882064),
        (0.0006, 0.8, 600, 0.184921764),
        (0.003, 0.9, 600, 0.061393152),
        (0.003, 0.8, 2000, 0.036363636),
    ],
)
def test_k_v_FRC_LOA_3_1(
    eps_x,
    k_dg,
    z_v,
    expected,
):
    """Test function theta_F_LOA_3_1."""
    assert math.isclose(
        _section_30_1_3_shear.k_v_FRC_LOA_3_1(
            eps_x,
            k_dg,
            z_v,
        ),
        expected,
        rel_tol=1e-6,
    )


@pytest.mark.parametrize(
    'eps_x, k_dg, d_v, theta_F, expected',
    [
        (0.003, 0.8, 600, 30, 4.206662884),
        (0.0003, 0.8, 600, 30, 0.657291076),
        (0.003, 0.9, 600, 30, 4.377203272),
        (0.003, 0.8, 2000, 30, 7.390083446),
        (0.003, 0.8, 600, 70, 10.65164434),
        (0.0001, 0.8, 600, 30, 0.5),  # Check max. cond.
    ],
)
def test_w_u_LOA_3_1(
    eps_x,
    k_dg,
    d_v,
    theta_F,
    expected,
):
    """Test function w_u_LOA_3_1."""
    assert math.isclose(
        _section_30_1_3_shear.w_u_LOA_3_1(
            eps_x,
            k_dg,
            d_v,
            theta_F,
        ),
        expected,
        rel_tol=1e-9,
    )
