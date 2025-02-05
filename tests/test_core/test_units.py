"""Tests for the classes that handle units."""

import math

import pytest

from structuralcodes.core import UnitConverter, UnitSet
from structuralcodes.core._units import _UNITS


@pytest.mark.parametrize(
    'force_unit',
    list(_UNITS['force'].keys()),
)
@pytest.mark.parametrize(
    'length_unit',
    list(_UNITS['length'].keys()),
)
def test_existing_units(length_unit: str, force_unit: str):
    """Test if an existing unit is applied correctly."""
    # Arrange
    units = UnitSet(length=length_unit, force=force_unit)
    expected_length = _UNITS['length'].get(length_unit)
    expected_force = _UNITS['force'].get(force_unit)
    expected_stress = expected_force / expected_length**2

    # Act
    current_length = units.length_unit
    current_force = units.force_unit
    current_stress = units.stress_unit

    # Assert
    assert current_length == expected_length
    assert current_force == expected_force
    assert current_stress == expected_stress


@pytest.mark.parametrize(
    'force, expected_force',
    [
        ('n', 'N'),
        ('KN', 'kN'),
        ('mn', 'MN'),
    ],
)
@pytest.mark.parametrize(
    'length, expected_length',
    [
        ('M', 'm'),
        ('MM', 'mm'),
        ('Inch', 'inch'),
        ('FoOt', 'foot'),
    ],
)
def test_existing_unit_not_matching_case(
    length: str, expected_length: str, force: str, expected_force: str
):
    """Test if the unitset is not case sensitive."""
    # Act
    units = UnitSet(length=length, force=force)

    # Assert
    assert units.length == expected_length
    assert units.force == expected_force


@pytest.mark.parametrize(
    'force',
    [
        'not valid force unit',
    ],
)
@pytest.mark.parametrize(
    'length',
    [
        'not valid length',
    ],
)
def test_not_existing_unit(length: str, force: str):
    """Test if a ValueError is raised if an invalid unit is input."""
    # Assert
    with pytest.raises(ValueError):
        UnitSet(length=length, force=force)


@pytest.mark.parametrize(
    'force_from, force_to, scale',
    [
        ('N', 'N', 1),
        ('N', 'kN', 1e-3),
        ('N', 'MN', 1e-6),
    ],
)
def test_converting_force(force_from: str, force_to: str, scale: float):
    """Test converting force from unit to unit."""
    # Arrange
    from_units = UnitSet(force=force_from)
    to_units = UnitSet(force=force_to)
    converter = UnitConverter(from_units=from_units, to_units=to_units)

    # Act
    converted_force_forwards = converter.convert_force_forwards(1)
    converted_force_backwards = converter.convert_force_backwards(1)

    # Assert
    assert math.isclose(converted_force_forwards, scale)
    assert math.isclose(converted_force_backwards, 1 / scale)


@pytest.mark.parametrize(
    'length_from, length_to, scale',
    [
        ('mm', 'mm', 1),
        ('mm', 'm', 1e-3),
    ],
)
def test_converting_length(length_from: str, length_to: str, scale: float):
    """Test converting length from unit to unit."""
    # Arrange
    from_units = UnitSet(length=length_from)
    to_units = UnitSet(length=length_to)
    converter = UnitConverter(from_units=from_units, to_units=to_units)

    # Act
    converted_length_forwards = converter.convert_length_forwards(1)
    converted_length_backwards = converter.convert_length_backwards(1)

    # Assert
    assert math.isclose(converted_length_forwards, scale)
    assert math.isclose(converted_length_backwards, 1 / scale)


@pytest.mark.parametrize(
    'force_from, force_to, length_from, length_to, scale',
    [
        ('N', 'kN', 'mm', 'm', 1e3),
        ('N', 'MN', 'mm', 'm', 1),
        ('N', 'MN', 'm', 'm', 1e-6),
        ('N', 'N', 'm', 'mm', 1e-6),
    ],
)
def test_converting_stress(
    force_from: str,
    force_to: str,
    length_from: str,
    length_to: str,
    scale: float,
):
    """Test converting stress from unit to unit."""
    # Arrange
    from_units = UnitSet(force=force_from, length=length_from)
    to_units = UnitSet(force=force_to, length=length_to)
    converter = UnitConverter(from_units=from_units, to_units=to_units)

    # Act
    converted_stress_forwards = converter.convert_stress_forwards(1)
    converted_stress_backwards = converter.convert_stress_backwards(1)

    # Assert
    assert math.isclose(converted_stress_forwards, scale)
    assert math.isclose(converted_stress_backwards, 1 / scale)
