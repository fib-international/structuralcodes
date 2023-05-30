"""Test for the function of _concrete_torsion"""
import math

import pytest

from structuralcodes.codes.mc2010 import _concrete_torsion


def create_load_dict(
    Med: float, Ved: float, Ned: float, delta_e: float
) -> dict:
    """returns dictionary assosiated with loads"""
    dictionary = {'Med': Med, 'Ved': Ved, 'Ned': Ned, 'delta_e': delta_e}
    return dictionary


@pytest.mark.parametrize(
    '''t_ed, a_k, z_i, expected''',
    [
        (2000000, 2000, 300, 150000),
    ],
)
def test_ved_ti(t_ed, a_k, z_i, expected):
    """Test the ved_ti function."""
    assert math.isclose(
        _concrete_torsion.v_ed_ti(t_ed, a_k, z_i), expected, rel_tol=0.001
    )


@pytest.mark.parametrize(
    '''f_ck, d_k, a_k, theta, approx_lvl, z, E_s, As,
        loads, gamma_c, expected''',
    [
        (
            35,
            150,
            50000,
            40,
            1,
            180,
            200000,
            2000,
            create_load_dict(0, 2000, 0, 20),
            1.5,
            10084252,
        ),
        (
            35,
            150,
            50000,
            40,
            2,
            180,
            200000,
            2000,
            create_load_dict(0, 2000, 0, 20),
            1.5,
            13301397,
        ),
        (
            35,
            150,
            50000,
            40,
            3,
            180,
            200000,
            2000,
            create_load_dict(0, 2000, 0, 20),
            1.5,
            10084000,
        ),
    ],
)
def test_t_rd_max(
    f_ck, d_k, a_k, theta, approx_lvl, z, E_s, As, loads, gamma_c, expected
):
    """Test the t_rd_max function."""
    assert math.isclose(
        _concrete_torsion.t_rd_max(
            f_ck, d_k, a_k, theta, approx_lvl, z, E_s, As, loads, gamma_c
        ),
        expected,
        rel_tol=0.001,
    )


@pytest.mark.parametrize(
    '''t_ed, approx_lvl, fck, bw, theta, z, E_s, As,
        loads, d_k, a_k, alpha, gamma_c, expected''',
    [
        (
            100e3,
            3,
            35,
            200,
            40,
            180,
            200000,
            2000,
            create_load_dict(0, 10e3, 10e3, 20),
            150,
            50000,
            90,
            1.5,
            True,
        ),
        (
            100e3,
            3,
            35,
            200,
            40,
            180,
            200000,
            2000,
            create_load_dict(0, 1000e3, 10e3, 20),
            150,
            50000,
            90,
            1.5,
            False,
        ),
        (
            100e3,
            3,
            35,
            200,
            40,
            180,
            200000,
            2000,
            create_load_dict(0, 10e3, 10e3, 20),
            150,
            50000,
            90,
            1.5,
            True,
        ),
        (
            15000e3,
            3,
            35,
            200,
            40,
            180,
            200000,
            2000,
            create_load_dict(0, 10e3, 10e3, 20),
            150,
            50000,
            90,
            1.5,
            False,
        ),
    ],
)
def test_t_rd(
    t_ed,
    approx_lvl,
    fck,
    bw,
    theta,
    z,
    E_s,
    As,
    loads,
    d_k,
    a_k,
    alpha,
    gamma_c,
    expected,
):
    """Test the t_rd function."""
    assert math.isclose(
        _concrete_torsion.t_rd(
            t_ed,
            approx_lvl,
            fck,
            bw,
            theta,
            z,
            E_s,
            As,
            loads,
            d_k,
            a_k,
            alpha,
            gamma_c,
        ),
        expected,
    )
