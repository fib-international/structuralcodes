"""Test module for AASHTO LRFD 2024 Crack Control."""

import math

import pytest

from structuralcodes.codes.aashto_2020 import _crack_control


@pytest.mark.parametrize(
    'h, dc, expected',
    [
        (13.7795, 1.37795, 1.15873),
        (13.7795, 0.98425, 1.10989),
        (13.7795, 2.36220, 1.29557),
    ],
)
def test_beta_s(h, dc, expected):
    """Test the besta_s function."""
    assert math.isclose(_crack_control.beta_s(h, dc), expected, rel_tol=0.005)


@pytest.mark.parametrize('h, dc', [(-12, 1.45), (15, -1.23)])
def test_beta_s_errors(h, dc):
    """Test beta_s errors."""
    with pytest.raises(ValueError):
        _crack_control.beta_s(h, dc)


@pytest.mark.parametrize(
    'fyk, beta_s, gamma_e, dc, expected',
    [
        (72.5, 1.15873, 0.4632, 1.37795, 3.676),
        (72.5, 1.15873, 0.6948, 1.37795, 6.893),
        (72.5, 1.29557, 0.6948, 2.36220, 3.905),
    ],
)
def test_s(fyk, beta_s, gamma_e, dc, expected):
    """Test s function."""
    assert math.isclose(
        _crack_control.s(fyk, beta_s, gamma_e, dc), expected, rel_tol=0.005
    )
