"""Tests for reinforcement steel material properties."""

import math

import pytest

from structuralcodes.codes.ec2_2004 import _reinforcement_material_properties


@pytest.mark.parametrize(
    'fyk, gamma_S, expected',
    [(500, 1.15, 434.782609), (600, 1.15, 521.739130), (400, 1, 400)],
)
def test_fyd(fyk, gamma_S, expected):
    """Test fyd function."""
    assert math.isclose(
        _reinforcement_material_properties.fyd(fyk, gamma_S),
        expected,
        rel_tol=10e-5,
    )


@pytest.mark.parametrize('fyk, gamma_S', [(-300, 1.15), (600, -3)])
def test_fyd_raises_errors(fyk, gamma_S):
    """Test fyd function raises errors."""
    with pytest.raises(ValueError):
        _reinforcement_material_properties.fyd(fyk, gamma_S)
