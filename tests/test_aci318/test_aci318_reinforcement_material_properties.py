"""Tests for reinforcement material properties of ACI 318-19."""

import math

import pytest

from structuralcodes.codes.aci318 import _reinforcement_material_properties


def test_Es():
    """Test modulus of elasticity of reinforcement."""
    assert _reinforcement_material_properties.Es() == 200000.0


@pytest.mark.parametrize(
    'fy, expected',
    [
        (280, 280),
        (420, 420),
        (550, 550),
        (690, 690),
    ],
)
def test_fy_design_default(fy, expected):
    """Test design yield strength with default phi=1.0."""
    assert math.isclose(
        _reinforcement_material_properties.fy_design(fy), expected
    )


def test_fy_design_with_phi():
    """Test design yield strength with explicit phi."""
    assert math.isclose(
        _reinforcement_material_properties.fy_design(420, phi=0.9), 378
    )


def test_fy_design_invalid_fy():
    """Test fy_design raises for non-positive fy."""
    with pytest.raises(ValueError):
        _reinforcement_material_properties.fy_design(-1)


def test_fy_design_invalid_phi():
    """Test fy_design raises for phi outside (0, 1]."""
    with pytest.raises(ValueError):
        _reinforcement_material_properties.fy_design(420, phi=1.5)


@pytest.mark.parametrize(
    'fy, expected',
    [
        (420, 420 / 200000),
        (280, 280 / 200000),
        (550, 550 / 200000),
    ],
)
def test_epsyd(fy, expected):
    """Test yield strain."""
    assert math.isclose(_reinforcement_material_properties.epsyd(fy), expected)


def test_epsyd_invalid_fy():
    """Test epsyd raises for non-positive fy."""
    with pytest.raises(ValueError):
        _reinforcement_material_properties.epsyd(-1)


@pytest.mark.parametrize(
    'grade, exp_fy, exp_fu',
    [
        ('40', 280, 420),
        ('60', 420, 550),
        ('80', 550, 690),
        ('100', 690, 860),
    ],
)
def test_reinforcement_grade_props(grade, exp_fy, exp_fu):
    """Test reinforcement grade property lookup."""
    props = _reinforcement_material_properties.reinforcement_grade_props(grade)
    assert math.isclose(props['fy'], exp_fy)
    assert math.isclose(props['fu'], exp_fu)


def test_reinforcement_grade_props_invalid():
    """Test grade lookup raises for unknown grade."""
    with pytest.raises(ValueError):
        _reinforcement_material_properties.reinforcement_grade_props('999')
