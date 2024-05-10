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


def test_duct_class_not_existing():
    """Test getting ductility properties for a ductility class that is not
    available.
    """
    # Arrange
    fyk = 500
    ductility_class = 'not a class'
    # Assert
    with pytest.raises(ValueError):
        _reinforcement_material_properties.reinforcement_duct_props(
            fyk=fyk, ductility_class=ductility_class
        )


@pytest.mark.parametrize(
    'ductility_class, exp_ratio, exp_strain',
    [
        ('a', 1.05, 2.5e-2),
        ('b', 1.08, 5e-2),
        ('c', 1.15, 7.5e-2),
    ],
)
def test_duct_class_props(ductility_class, exp_ratio, exp_strain):
    """Test getting ductility class properties."""
    # Arrange
    fyk = 500

    # Act
    props = _reinforcement_material_properties.reinforcement_duct_props(
        fyk=fyk, ductility_class=ductility_class
    )

    # Assert
    assert math.isclose(props['ftk'] / fyk, exp_ratio)
    assert math.isclose(props['epsuk'], exp_strain)
