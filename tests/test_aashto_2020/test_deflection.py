"""Test module for AASHTO LRFD 2024 Deflections."""

import math

import pytest

from structuralcodes.codes.aashto_2020 import _deflections


@pytest.mark.parametrize(
    'fc_prime, fy, rho, b, h, d, L, phi, expected',
    [
        (
            3.625,
            65.25,
            0.00115,
            39.37,
            27.559,
            25.591,
            307.086,
            0.9,
            15609.715,
        ),
        (
            3.625,
            65.25,
            0.00115,
            39.37,
            27.559,
            25.591,
            358.267,
            0.9,
            15609.715,
        ),
        (
            3.625,
            65.25,
            0.00115,
            39.37,
            27.559,
            25.591,
            383.858,
            0.9,
            15609.715,
        ),
    ],
)
def test_Mu(fc_prime, fy, rho, b, h, d, L, phi, expected):
    """Test the Mu function."""
    assert math.isclose(
        _deflections.Mu(fc_prime, fy, rho, b, h, d, L, phi),
        expected,
        rel_tol=0.005,
    )


@pytest.mark.parametrize(
    'fc_prime, b, h, d, expected',
    [
        (3.625, 39.37, 27.559, 25.5905, 2656.765),
        (7.250, 39.37, 27.559, 25.5905, 3757.233),
        (14.50, 39.37, 27.559, 25.5905, 5313.530),
    ],
)
def test_Mcr(fc_prime, b, h, d, expected):
    """Test the Mcr function."""
    assert math.isclose(
        _deflections.Mcr(fc_prime, b, h, d),
        expected,
        rel_tol=0.005,
    )


@pytest.mark.parametrize(
    'fc_prime, fy, rho, b, h, d, L, phi, Es, Ec, expected',
    [
        (
            3.625,
            65.25,
            0.00115,
            39.37,
            27.559,
            25.591,
            307.086,
            0.9,
            29000,
            3823.945,
            68671.81,
        ),
        (
            7.25,
            65.25,
            0.00115,
            39.37,
            27.559,
            25.591,
            307.086,
            0.9,
            29000,
            4806.75,
            68671.81,
        ),
        (
            14.50,
            65.25,
            0.00115,
            39.37,
            27.559,
            25.591,
            307.086,
            0.9,
            29000,
            6042.149,
            83788.31,
        ),
    ],
)
def test_Ie(fc_prime, fy, rho, b, h, d, L, phi, Es, Ec, expected):
    """Test the Ie function."""
    assert math.isclose(
        _deflections.Ie(fc_prime, fy, rho, b, h, d, L, phi, Es, Ec),
        expected,
        rel_tol=0.005,
    )


@pytest.mark.parametrize(
    'qqp, L, Ec, Ieff, expected',
    [
        (1.337, 307.086, 3823.945, 35648.81, 1.136),
        (2.673, 307.086, 4806.750, 54188.30, 1.188),
        (5.346, 307.086, 6042.149, 83788.31, 1.223),
    ],
)
def test_delta(qqp, L, Ec, Ieff, expected):
    """Test delta function."""
    assert math.isclose(
        _deflections.delta(qqp, L, Ec, Ieff), expected, rel_tol=0.005
    )


@pytest.mark.parametrize(
    'deflection, factor, rho_prime, expected',
    [
        (1.136, 2, 0.003776, 1.910),
        (1.188, 2, 0.003776, 1.999),
        (1.223, 2, 0.003776, 2.057),
    ],
)
def test_time_delta(deflection, factor, rho_prime, expected):
    """Test time_delta function."""
    assert math.isclose(
        _deflections.time_delta(deflection, factor, rho_prime),
        expected,
        rel_tol=0.005,
    )
