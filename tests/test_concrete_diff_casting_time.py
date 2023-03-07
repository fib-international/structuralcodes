"""Test for the function of _interface_different_casting_time"""
import math

import pytest

from structuralcodes.codes.mc2010 import _concrete_interface_different_casting_times

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
    assert math.isclose(_concrete_interface_different_casting_times.tau_edi(
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
        (0.2, 2.6, 0.6, 100, 20, 11.3, 3.1075),
        (0.2, 2.6, 0.6, 100, 30, 17, 4.675),
    ],
)
def test_tau_rdi_without_reinforcement(
    c_a, f_ctd, mu, sigma_n, f_ck, f_cd, expected
):

    """Test the tau_rdi_without_reinforcement function."""
    assert math.isclose(_concrete_interface_different_casting_times.tau_rdi_without_reinforcement(
        c_a, f_ctd, mu, sigma_n, f_ck, f_cd
    ),
        expected, rel_tol=0.001)


@pytest.mark.parametrize(
    '''c_r, k1, k2, mu, ro, sigma_n, alpha,
    beta_c, f_ck, f_yd, f_cd, expected''',
    [
        (0.1, 0.5, 0.9, 0.7, 0.05, 100, 15, 0.5, 30, 434, 17, 4.675),
    ],
)
def test_tau_rdi_with_reinforcement(
    c_r, k1, k2, mu, ro, sigma_n, alpha, beta_c, f_ck, f_yd, f_cd, expected
):

    """Test the tau_rdi_with_reinforcement function."""
    assert math.isclose(_concrete_interface_different_casting_times.tau_rdi_with_reinforcement(
        c_r, k1, k2, mu, ro, sigma_n, alpha, beta_c, f_ck, f_yd, f_cd
    ),
        expected, rel_tol=0.001)
