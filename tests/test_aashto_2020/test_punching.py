"""Test module for testing the AASHTO 2020 punching design."""

import math

import pytest

from structuralcodes.codes.aashto_2020 import _punching


@pytest.mark.parametrize(
    'W, L, df, expected',
    [
        (23.622, 19.685, 11.811, 86.614),
    ],
)
def test_b0_edge(W, L, df, expected):
    """Test the b0_edge function."""
    assert math.isclose(_punching.b0_edge(W, L, df), expected, rel_tol=0.005)


@pytest.mark.parametrize(
    'W, L, df, c, expected',
    [
        (23.622, 19.685, 11.811, 25.685, 68.992),
    ],
)
def test_b0_corner(W, L, df, c, expected):
    """Test the b0_corner function."""
    assert math.isclose(
        _punching.b0_corner(W, L, df, c), expected, rel_tol=0.005
    )


@pytest.mark.parametrize(
    'W, L, df, expected',
    [
        (23.622, 19.685, 11.811, 133.858),
    ],
)
def test_b0_interior(W, L, df, expected):
    """Test the b0_interior function."""
    assert math.isclose(
        _punching.b0_interior(W, L, df), expected, rel_tol=0.005
    )
