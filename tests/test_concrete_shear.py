import math

import pytest

from structuralcodes.codes.mc2010 import _concrete_shear


@pytest.mark.parametrize(
    'E_s, As, Med, Ved, Ned, z, deltaE, expected',
    [
        (200000, 1000, 50000000, 10000, 2000, 160, 50, 8.1e-4),
        (210000, 1000, 50000000, 10000, 2000, 160, 50, 7.7e-4),
        (210000, 5000, 50000000, 10000, 2000, 160, 50, 1.5e-4),
        (210000, 2000, 50000000, 10000, 2000, 160, 50, 3.9e-4),
        (210000, 2000, 40000000, 20000, 2000, 160, 50, 3.2e-4),
        (210000, 2000, 40000000, 20000, 1000, 160, 50, 3.2e-4),
        (210000, 2000, 40000000, 20000, 1000, 140, 50, 3.64965e-4),
        (210000, 2000, 40000000, 20000, 1000, 180, 50, 2.9e-4),
    ],
)
def test_epsilon_x(E_s, As, Med, Ved, Ned, z, deltaE, expected):
    """Test the epsilon_x function."""
    assert math.isclose(_concrete_shear.epsilon_x(
        E_s, As, Med, Ved, Ned, z, deltaE), expected, rel_tol=0.05
    )


@pytest.mark.parametrize(
    'approx_lvl_s, fck, bw, theta, z, E_s, As, Med, Ved, Ned, delta_e, alfa, gamma_c, expected',
    [
        (1, 30, 50, 20, 200, 210000, 1000, 200e6, 50e3, 10e3, 50, 20, 1.5, 70707),
        (2, 30, 50, 20, 200, 210000, 1000, 200e6, 50e3, 10e3, 50, 20, 1.5, 39997),
        (2, 30, 50, 20, 200, 210000, 1000, 50e6, 10e3, 10e3, 50, 20, 1.5, 55179.55),
        (2, 30, 50, 45, 200, 210000, 1000, 0, 0, 0, 50, 20, 1.5, 243586),
        (2, 30, 50, 45, 200, 210000, 1000, 0, 0, 0, 50, 45, 1.5, 130000),
        (3, 30, 50, 20, 200, 210000, 1000, 50e6, 10e3, 10e3, 50, 20, 1.5, 102995),

    ],
)
def test_vrd_max(
    approx_lvl_s, fck, bw, theta, z, E_s, As, Med, Ved, Ned,
    delta_e, alfa, gamma_c, expected
):
    """Test the v_rd_max function."""
    assert math.isclose(_concrete_shear.v_rd_max(
        approx_lvl_s, fck, bw, theta, z, E_s, As, Med, Ved, Ned, delta_e, alfa, gamma_c), expected, rel_tol=0.5
    )


@pytest.mark.parametrize(
    '''approx_lvl_c, approx_lvl_s, fck, z, bw, dg, E_s, As, Med,
     Ved, Ned, delta_e, alfa, gamma_c, expected''',
    [
        (1, 0, 35, 180, 300, 0, 0, 0, 0, 0, 0, 0, 0, 1.5, 31294),
        (1, 0, 35, 200, 300, 0, 0, 0, 0, 0, 0, 0, 0, 1.5, 34077),
        (1, 1, 35, 200, 300, 0, 0, 0, 0, 0, 0, 0, 0, 1.5, 34077),
        (2, 1, 35, 140, 300, 16, 21e4, 2000, 40e6, 2e4, 1000, 50, 0, 1.5, 48828),
        (2, 1, 35, 140, 300, 32, 21e4, 2000, 40e6, 2e4, 1000, 50, 0, 1.5, 50375),
        (0, 3, 35, 200, 300, 32, 21e4, 2000, 40e6, 2e4, 1000, 50, 1.5, 1.5, 67566),
        (0, 3, 35, 200, 300, 32, 21e4, 2000, 40e6, 2000e4, 1000, 50, 1.5, 1.5, 0),
    ],
)
def test_v_rdc(
    approx_lvl_c, approx_lvl_s, fck, z, bw, dg, E_s, As, Med,
    Ved, Ned, delta_e, alfa, gamma_c, expected
):

    """Test the v_rdc function."""
    assert math.isclose(_concrete_shear.v_rdc(
                approx_lvl_c, approx_lvl_s, fck, z, bw, dg, E_s, As, Med, Ved,
                Ned, delta_e, alfa, gamma_c
                ),
        expected, rel_tol=0.001)

@pytest.mark.parametrize(
    '''asw, sw, z, fywd, theta, alpha, expected''',
    [
        (1600, 50, 200, 355, 25, 30, 4403769),
        (2000, 50, 200, 355, 25, 30, 5504711),
        (1600, 50, 200, 355, 25, 30, 4403769),
        (1600, 100, 200, 355, 25, 30, 2201884),
        (1600, 50, 200, 275, 25, 30, 3411370),
        (1600, 50, 200, 355, 22, 30, 4779308),
        (1600, 50, 200, 355, 25, 25, 4118262),
    ],
)
def test_v_rds(asw, sw, z, fywd, theta, alpha, expected):

    """Test the v_rds function."""
    assert math.isclose(_concrete_shear.v_rds(
                asw, sw, z, fywd, theta, alpha
                ),
        expected, rel_tol=0.001)


@pytest.mark.parametrize(
    '''approx_lvl_h, f_ctd, i_c, s_c, b_w, sigma_cp, l_x, l_bd0, S_cy,
    b_wy, y, y_c, A_c, A_cy, y_pt, f_p_lx, f_p_lx_dx, expected''',
    [
        (1, 2.6, 6e8, 6e5, 50, 150, 40, 30, 3e6, 3e5, 100, 200, 2000, 1000, 80, 1000e3, 200e3, 839136),
        (1, 3.5, 6e8, 6e5, 50, 150, 40, 30, 3e6, 3e5, 100, 200, 2000, 1000, 80, 1000e3, 200e3, 976183),
        (1, 2.6, 5e8, 6e5, 50, 150, 40, 30, 3e6, 3e5, 100, 200, 2000, 1000, 80, 1000e3, 200e3, 699280),
        (1, 2.6, 6e8, 5e5, 50, 150, 40, 30, 3e6, 3e5, 100, 200, 2000, 1000, 80, 1000e3, 200e3, 1006963),
        (1, 2.6, 6e8, 6e5, 40, 150, 40, 30, 3e6, 3e5, 100, 200, 2000, 1000, 80, 1000e3, 200e3, 671309),
        (1, 2.6, 6e8, 6e5, 50, 180, 40, 30, 3e6, 3e5, 100, 200, 2000, 1000, 80, 1000e3, 200e3, 918050),
        (1, 2.6, 6e8, 6e5, 50, 150, 35, 30, 3e6, 3e5, 100, 200, 2000, 1000, 80, 1000e3, 200e3, 785800),
        (1, 2.6, 6e8, 6e5, 50, 150, 40, 25, 3e6, 3e5, 100, 200, 2000, 1000, 80, 1000e3, 200e3, 918050),
        (1, 2.6, 6e8, 6e5, 50, 150, 40, 30, 2e6, 3e5, 100, 200, 2000, 1000, 80, 1000e3, 200e3, 839136),
        (1, 2.6, 6e8, 6e5, 50, 150, 40, 30, 3e6, 2e5, 100, 200, 2000, 1000, 80, 1000e3, 200e3, 839136),
        (1, 2.6, 6e8, 6e5, 50, 150, 40, 30, 3e6, 3e5, 80, 200, 2000, 1000, 80, 1000e3, 200e3, 839136),
        (1, 2.6, 6e8, 6e5, 50, 150, 40, 30, 3e6, 3e5, 100, 180, 2000, 1000, 80, 1000e3, 200e3, 839136),
        (1, 2.6, 6e8, 6e5, 50, 150, 40, 30, 3e6, 3e5, 100, 200, 1800, 1000, 80, 1000e3, 200e3, 839136),
        (1, 2.6, 6e8, 6e5, 50, 150, 40, 30, 3e6, 3e5, 100, 200, 2000, 1200, 80, 1000e3, 200e3, 839136),
        (1, 2.6, 6e8, 6e5, 50, 150, 40, 30, 3e6, 3e5, 100, 200, 2000, 1000, 60, 1000e3, 200e3, 839136),
        (1, 2.6, 6e8, 6e5, 50, 150, 40, 30, 3e6, 3e5, 100, 200, 2000, 1000, 80, 800e3, 200e3, 839136),
        (1, 2.6, 6e8, 6e5, 50, 150, 40, 30, 3e6, 3e5, 100, 200, 2000, 1000, 80, 1000e3, 250e3, 839136),

        (2, 2.6, 6e8, 6e5, 50, 150, 40, 30, 3e6, 3e5, 100, 200, 2000, 1000, 80, 1000e3, 200e3, 160002777),
        (2, 3.5, 6e8, 6e5, 50, 150, 40, 30, 3e6, 3e5, 100, 200, 2000, 1000, 80, 1000e3, 200e3, 214002777),
        (2, 2.6, 5e8, 6e5, 50, 150, 40, 30, 3e6, 3e5, 100, 200, 2000, 1000, 80, 1000e3, 200e3, 137336111),
        (2, 2.6, 6e8, 6e5, 50, 150, 40, 30, 3e6, 3e5, 100, 200, 2000, 1000, 80, 1000e3, 200e3, 160002777),
        (2, 2.6, 6e8, 6e5, 40, 150, 40, 30, 3e6, 3e5, 100, 200, 2000, 1000, 80, 1000e3, 200e3, 160002777),
        (2, 2.6, 6e8, 6e5, 50, 180, 40, 30, 3e6, 3e5, 100, 200, 2000, 1000, 80, 1000e3, 200e3, 160002777),
        (2, 2.6, 6e8, 6e5, 50, 150, 35, 30, 3e6, 3e5, 100, 200, 2000, 1000, 80, 1000e3, 200e3, 160002430),
        (2, 2.6, 6e8, 6e5, 50, 150, 40, 25, 3e6, 3e5, 100, 200, 2000, 1000, 80, 1000e3, 200e3, 160003333),
        (2, 2.6, 6e8, 6e5, 50, 150, 40, 30, 2e6, 3e5, 100, 200, 2000, 1000, 80, 1000e3, 200e3, 228004166),
        (2, 2.6, 6e8, 6e5, 50, 150, 40, 30, 3e6, 2e5, 100, 200, 2000, 1000, 80, 1000e3, 200e3, 108001851),
        (2, 2.6, 6e8, 6e5, 50, 150, 40, 30, 3e6, 3e5, 80, 200, 2000, 1000, 80, 1000e3, 200e3, 160003333),
        (2, 2.6, 6e8, 6e5, 50, 150, 40, 30, 3e6, 3e5, 100, 180, 2000, 1000, 80, 1000e3, 200e3, 156002222),
        (2, 2.6, 6e8, 6e5, 50, 150, 40, 30, 3e6, 3e5, 100, 200, 1800, 1000, 80, 1000e3, 200e3, 157780864),
        (2, 2.6, 6e8, 6e5, 50, 150, 40, 30, 3e6, 3e5, 100, 200, 2000, 1200, 80, 1000e3, 200e3, 156002777),
        (2, 2.6, 6e8, 6e5, 50, 150, 40, 30, 3e6, 3e5, 100, 200, 2000, 1000, 60, 1000e3, 200e3, 164002777),
        (2, 2.6, 6e8, 6e5, 50, 150, 40, 30, 3e6, 3e5, 100, 200, 2000, 1000, 80, 800e3, 200e3, 160002222),
        (2, 2.6, 6e8, 6e5, 50, 150, 40, 30, 3e6, 3e5, 100, 200, 2000, 1000, 80, 1000e3, 250e3, 161002777),
        
    ],
)
def test_v_rd_ct(
    approx_lvl_h, f_ctd, i_c, s_c, b_w, sigma_cp, l_x, l_bd0, S_cy,
    b_wy, y, y_c, A_c, A_cy, y_pt, f_p_lx, f_p_lx_dx, expected
):

    """Test the v_rd_ct function."""
    assert math.isclose(_concrete_shear.v_rd_ct(
        approx_lvl_h, f_ctd, i_c, s_c, b_w, sigma_cp, l_x, l_bd0, S_cy,
        b_wy, y, y_c, A_c, A_cy, y_pt, f_p_lx, f_p_lx_dx
    ),
        expected, rel_tol=0.001)


@pytest.mark.parametrize(
    '''beta, v_ed, z, b_i, expected''',
    [
        (0.7, 50e3, 200, 50, 3.5),
        (0.75, 50e3, 200, 50, 3.75),
        (0.7, 40e3, 200, 50, 2.8),
        (0.7, 50e3, 180, 50, 3.888),
        (0.7, 50e3, 200, 60, 2.916),
    ],
)
def test_tau_edi(beta, v_ed, z, b_i, expected):

    """Test the tau_edi function."""
    assert math.isclose(_concrete_shear.tau_edi(
        beta, v_ed, z, b_i
    ),
        expected, rel_tol=0.001)


@pytest.mark.parametrize(
    '''c_a, f_ctd, mu, sigma_n, f_ck, f_cd, expected''',
    [
        (0.2, 2.6, 0.6, 100, 30, 17, 4.675),
        (0.2, 3.5, 0.6, 100, 30, 17, 4.675),
        (0.2, 2.6, 0.7, 100, 30, 17, 4.675),
        (0.2, 2.6, 0.6, 80, 30, 17, 4.675),
        (0.2, 2.6, 0.6, 100, 20, 11.3, 4.675),
        (0.2, 2.6, 0.6, 100, 30, 17, 4.675),
    ],
)
def test_tau_rdi_without_reinforceent(
    c_a, f_ctd, mu, sigma_n, f_ck, f_cd, expected
):

    """Test the tau_edi function."""
    assert math.isclose(_concrete_shear.tau_rdi_without_reinforceent(
        c_a, f_ctd, mu, sigma_n, f_ck, f_cd
    ),
        expected, rel_tol=0.001)


@pytest.mark.parametrize(
    '''approx_lvl_c, approx_lvl_s, reinforcment, fck, z,
    bw, dg, E_s, As, Med, Ved, Ned, delta_e, alfa, gamma_c,
    asw, sw, f_ywd, theta, expected''',
    [
        (1, 0, False, 35, 180, 200, 16, 200000, 2000, 0, 2000, 0, 20, 90, 1.5, 0, 0, 434, 40, 20863),
        (2, 0, False, 35, 180, 200, 16, 200000, 2000, 0, 2000, 0, 20, 90, 1.5, 0, 0, 434, 40, 62336),
        (2, 3, True, 35, 180, 200, 16, 200000, 2000, 0, 2000, 0, 20, 90, 1.5, 500, 200, 434, 40, 62336),
    ],
)
def test_v_rd(
    approx_lvl_c, approx_lvl_s, reinforcment, fck, z,
    bw, dg, E_s, As, Med, Ved, Ned, delta_e, alfa, gamma_c,
    asw, sw, f_ywd, theta, expected
):

    """Test the tau_edi function."""
    assert math.isclose(_concrete_shear.v_rd(
        approx_lvl_c, approx_lvl_s, reinforcment, fck, z,
        bw, dg, E_s, As, Med, Ved, Ned, delta_e, alfa, gamma_c,
        asw, sw, f_ywd, theta
    ),
        expected, rel_tol=0.001)
    