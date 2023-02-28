"""Test for the function of _concrete_torsion"""
import math

import pytest

from structuralcodes.codes.mc2010 import _concrete_torsion


@pytest.mark.parametrize(
    '''t_ed, a_k, z_i, expected''',
    [
        (2000000, 2000, 300, 150000),
    ],
)
def test_ved_ti(t_ed, a_k, z_i, expected):

    """Test the ved_ti function."""
    assert math.isclose(_concrete_torsion.v_ed_ti(t_ed, a_k, z_i),
                        expected, rel_tol=0.001)


@pytest.mark.parametrize(
    '''f_ck, gamma_c, d_k, a_k, theta, approx_lvl_s, E_s, As,
        Med, Ved, Ned, z, delta_e, expected''',
    [
        (35, 1.5, 150, 50000, 40, 1, 200000, 2000, 0, 2000, 0, 180, 20,
         11256044),
        (35, 1.5, 150, 50000, 40, 2, 200000, 2000, 0, 2000, 0, 180, 20,
         13301400),
        (35, 1.5, 150, 50000, 40, 3, 200000, 2000, 0, 2000, 0, 180, 20,
         10084000),

    ],
)
def test_t_rd_max(
        f_ck, gamma_c, d_k, a_k, theta, approx_lvl_s, E_s, As,
        Med, Ved, Ned, z, delta_e, expected):

    """Test the t_rd_max function."""
    assert math.isclose(_concrete_torsion.t_rd_max(
            f_ck, gamma_c, d_k, a_k, theta, approx_lvl_s, E_s, As,
            Med, Ved, Ned, z, delta_e), expected, rel_tol=0.001)


@pytest.mark.parametrize(
    '''t_ed, approx_lvl, fck, bw, theta, z, E_s, As,
        Med, Ved, Ned, delta_e, alfa, d_k, a_k, gamma_c, expected''',
    [
        (100e3, 3, 35, 200, 40, 180, 200000, 2000, 0, 10e3, 10e3, 20,
         90, 150, 50000, 1.5, True),
        (100e3, 3, 35, 200, 40, 180, 200000, 2000, 0, 1000e3, 10e3, 20,
         90, 150, 50000, 1.5, False),
        (100e3, 3, 35, 200, 40, 180, 200000, 2000, 0, 10e3, 10e3, 20,
         90, 150, 50000, 1.5, True),
        (10000e3, 3, 35, 200, 40, 180, 200000, 2000, 0, 10e3, 10e3, 20,
         90, 150, 50000, 1.5, False),
    ],
)
def test_t_rd(
        t_ed, approx_lvl, fck, bw, theta, z, E_s, As,
        Med, Ved, Ned, delta_e, alfa, d_k, a_k, gamma_c, expected
        ):

    """Test the t_rd function."""
    assert math.isclose(_concrete_torsion.t_rd(
        t_ed, approx_lvl, fck, bw, theta, z, E_s, As,
        Med, Ved, Ned, delta_e, alfa, d_k, a_k, gamma_c), expected)

