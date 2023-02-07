import math

import pytest

from structuralcodes.code.mc2010 import _concrete_shear


# @pytest.mark.parametrize(
#     'test_input1, testinput2, test, expected',
#     [
#         (12, 150, 50, 2625.424),
#         (14, 150, 50, 2835.782),
#         (16, 150, 50, 3031.579),
#         (20, 150, 50, 3389.408),
#         (25, 150, 50, 3789.474),
#         (30, 150, 50, 4151.160)
#     ]
# )
# def test_fcm(test_input1, test_input2, test_input3, expected):
#     """Test the fcm function."""
#     assert math.isclose(
#         _concrete_shear.vrdc(test_input1, test_input2, test_input3, expected)
#     )


@pytest.mark.parametrize(
    'E, As, Med, Ved, Ned, z, deltaE, expected',
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
def test_epsilon_x(E, As, Med, Ved, Ned, z, deltaE, expected):
    """Test the epsilon_x function."""
    assert math.isclose(_concrete_shear.epsilon_x(
        E, As, Med, Ved, Ned, z, deltaE), expected, rel_tol=0.05
    )


@pytest.mark.parametrize(
    'approx_lvl_s, fck, bw, theta, z, E, As, Med, Ved, Ned, delta_e, alfa, gamma_c, expected',
    [
        (1, 30, 50, 20, 200, 210000, 1000, 200e6, 50e3, 10e3, 50, 20, 1.5, 70707),
        (2, 30, 50, 20, 200, 210000, 1000, 200e6, 50e3, 10e3, 50, 20, 1.5, 39997),
        (2, 30, 50, 20, 200, 210000, 1000, 50e6, 10e3, 10e3, 50, 20, 1.5, 39997),
    ],
)
def test_vrd_max(approx_lvl_s, fck, bw, theta, z, E, As, Med, Ved, Ned,
    delta_e, alfa, gamma_c, expected):
    """Test the v_rd_max function."""
    assert math.isclose(_concrete_shear.v_rd_max(
        approx_lvl_s, fck, bw, theta, z, E, As, Med, Ved, Ned, delta_e, alfa, gamma_c), expected, rel_tol=0.05
    )


@pytest.mark.parametrize(
    '''approx_lvl_c, approx_lvl_s, fck, z, bw, dg, E, As, Med,
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
def test_v_rdc(approx_lvl_c, approx_lvl_s, fck, z, bw, dg, E, As, Med,
    Ved, Ned, delta_e, alfa, gamma_c, expected):

    """Test the v_rdc function."""
    assert math.isclose(_concrete_shear.v_rdc(
                approx_lvl_c, approx_lvl_s, fck, z, bw, dg, E, As, Med, Ved,
                Ned, delta_e, alfa, gamma_c
                ),
        expected, rel_tol=0.001)
        