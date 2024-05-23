"""Tests for the concrete EC2 2004."""

import math

import pytest

from structuralcodes.codes import ec2_2004
from structuralcodes.materials.concrete import (
    ConcreteEC2_2004,
    create_concrete,
)


@pytest.mark.parametrize(
    'design_code_to_set',
    ['ec2_2004', 'EC2_2004', 'eC2_2004'],
)
@pytest.mark.parametrize(
    'fck, expected_name',
    [(20, 'C20'), (25, 'C25'), (30, 'C30'), (35, 'C35'), (40, 'C40')],
)
def test_create_concrete(design_code_to_set, fck, expected_name):
    """Test creating a concrete with EC2 2004."""
    # Arrange
    expected_density = 2400

    # Act
    c = create_concrete(fck=fck, design_code=design_code_to_set)

    # Assert
    assert isinstance(c, ConcreteEC2_2004)
    assert c.name == expected_name
    assert c.density == expected_density


@pytest.mark.parametrize(
    'fck, fcm',
    [(12, 20), (35, 43), (55, 63), (90, 98), (120, 128)],
)
def test_update_attributes(fck, fcm):
    """Test update_attributes function."""
    c = create_concrete(fck=fck, design_code='ec2_2004')
    c.update_attributes({'fcm': fcm})
    # Test a warning is raised when a not valid key is inputted
    with pytest.warns(UserWarning):
        c.update_attributes({'not_valid_key': fcm})

    assert c.fcm is not None
    assert c.fcm == fcm


fck_parametrized = pytest.mark.parametrize('fck', [20, 25, 30, 35, 40])


@fck_parametrized
def test_fck_getter(fck):
    """Test fck getter."""
    c = ConcreteEC2_2004(fck)

    assert c.fck == fck


@fck_parametrized
def test_fck_setter(fck):
    """Test fck setter."""
    c = ConcreteEC2_2004(fck)
    c.fck = fck + 5

    assert c.fck == fck + 5


def test_properties_initialized_to_none():
    """Test if a ConcreteEC2_2004 has the attributes set to None."""
    c = ConcreteEC2_2004(fck=25)

    assert c._fcm is None
    assert c._fctm is None
    assert c._fctk_5 is None
    assert c._fctk_95 is None
    assert c._Ecm is None


def test_reset_properties():
    """Test _reset_attributes function."""
    c = ConcreteEC2_2004(fck=25)

    _ = c.fcm
    _ = c.fctm
    _ = c.fctk_5
    _ = c.fctk_95
    _ = c.Ecm

    c._reset_attributes()

    assert c._fcm is None
    assert c._fctm is None
    assert c._fctk_5 is None
    assert c._fctk_95 is None
    assert c._Ecm is None


@fck_parametrized
def test_fcm_getter(fck):
    """Test the fcm getter."""
    c = ConcreteEC2_2004(fck=fck)
    expected = ec2_2004.fcm(fck)
    assert math.isclose(c.fcm, expected)


fcm_parametrized = pytest.mark.parametrize(
    'test_input, expected',
    [(12, 20), (35, 43), (55, 63), (90, 98)],
)


@fcm_parametrized
def test_fcm_setter(test_input, expected):
    """Test the fcm setter."""
    c = ConcreteEC2_2004(fck=test_input)
    c.fcm = expected

    assert math.isclose(c.fcm, expected)


@pytest.mark.parametrize(
    'test_input',
    [12, 35, 55, 90, 120],
)
def test_fcm_setter_exception(test_input):
    """Test the fcm setter with a wrong value."""
    c = ConcreteEC2_2004(fck=test_input)
    with pytest.raises(ValueError):
        c.fcm = test_input - 1


@fck_parametrized
def test_fctm_getter(fck):
    """Test the fctm getter function."""
    c = ConcreteEC2_2004(fck=fck)
    expected = ec2_2004.fctm(fck)
    assert math.isclose(c.fctm, expected, rel_tol=1e-6)


fctm_parametrized = pytest.mark.parametrize(
    'test_input, expected',
    [(12, 1.6), (35, 3.2), (55, 4.2), (90, 5.0)],
)


@fctm_parametrized
def test_fctm_setter(test_input, expected):
    """Test the fctm setter."""
    c = ConcreteEC2_2004(fck=test_input)
    c.fctm = expected

    assert math.isclose(c.fctm, expected)


@pytest.mark.parametrize(
    'test_input',
    [12, 35, 55, 90],
)
def test_fctm_setter_warning(test_input):
    """Test the fctm setter with a wrong value."""
    c = ConcreteEC2_2004(fck=test_input)
    with pytest.warns(UserWarning):
        c.fctm = test_input * 0.5 + 1


@fck_parametrized
def test_fctk5_getter(fck):
    """Test the fctk_5 getter function."""
    c = ConcreteEC2_2004(fck=fck)
    expected = ec2_2004.fctk_5(ec2_2004.fctm(fck))
    assert math.isclose(c.fctk_5, expected, rel_tol=1e-6)


fctk5_parametrized = pytest.mark.parametrize(
    'test_input, expected',
    [(12, 1.1), (35, 2.2), (55, 3.0), (90, 3.5)],
)


@fctk5_parametrized
def test_fctk5_setter(test_input, expected):
    """Test the fctk_5 setter."""
    c = ConcreteEC2_2004(fck=test_input)
    c.fctk_5 = expected

    assert math.isclose(c.fctk_5, expected)


@fck_parametrized
def test_fctk95_getter(fck):
    """Test the fctk_95 getter function."""
    c = ConcreteEC2_2004(fck=fck)
    expected = ec2_2004.fctk_95(ec2_2004.fctm(fck))
    assert math.isclose(c.fctk_95, expected, rel_tol=1e-6)


fctk95_parametrized = pytest.mark.parametrize(
    'test_input, expected',
    [(12, 2.0), (35, 4.2), (55, 5.5), (90, 6.6)],
)


@fctk95_parametrized
def test_fctk95_setter(test_input, expected):
    """Test the fctk_95 setter."""
    c = ConcreteEC2_2004(fck=test_input)
    c.fctk_95 = expected

    assert math.isclose(c.fctk_95, expected)


@fck_parametrized
def test_Ecm_getter(fck):
    """Test the Ecm getter function."""
    c = ConcreteEC2_2004(fck=fck)
    expected = ec2_2004.Ecm(ec2_2004.fcm(fck))
    assert math.isclose(c.Ecm, expected, rel_tol=1e-6)


Ecm_parametrized = pytest.mark.parametrize(
    'test_input, expected',
    [(12, 27000), (35, 34000), (55, 38000), (90, 44000)],
)


@Ecm_parametrized
def test_Ecm_setter(test_input, expected):
    """Test the Ecm setter."""
    c = ConcreteEC2_2004(fck=test_input)
    c.Ecm = expected

    assert math.isclose(c.Ecm, expected)


@pytest.mark.parametrize(
    'fck, alpha_cc, fcd',
    [
        (35, 0.85, 19.833),
        (45, 0.85, 25.5),
        (70, 0.85, 39.667),
        (90, 0.85, 51.0),
    ],
)
def test_fcd(fck, alpha_cc, fcd):
    """Test the fcd method on the concrete class."""
    # Arrange
    concrete = ConcreteEC2_2004(fck=fck)

    # Act and assert
    assert math.isclose(
        concrete.fcd(alpha_cc=alpha_cc),
        fcd,
        rel_tol=1e-4,
    )
