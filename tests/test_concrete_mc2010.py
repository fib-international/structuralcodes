"""Tests for the concrete mc2010"""
import math

import pytest

from structuralcodes.materials.concrete import create_concrete, ConcreteMC2010

# Series of tests using the factory function


@pytest.mark.parametrize(
    'design_code_to_set',
    ['mc2010', 'MC2010', 'mC2010'],
)
@pytest.mark.parametrize(
    'fck, expected_name',
    [(20, 'C20'), (25, 'C25'), (30, 'C30'), (35, 'C35'), (40, 'C40')],
)
def test_create_concrete(design_code_to_set, fck, expected_name):
    """Test creating a concrete with MC2010."""
    # Arrange
    expected_density = 2400

    # Act
    c = create_concrete(fck=fck, design_code=design_code_to_set)

    # Assert
    assert isinstance(c, ConcreteMC2010)
    assert c.name == expected_name
    assert c.density == expected_density


def test_create_concrete_wrong_code():
    """Test if a ValueError exception raises when passing the wrong code."""
    with pytest.raises(ValueError):
        create_concrete(fck=25, design_code='EN1995')


@pytest.mark.parametrize(
    'fck, fcm',
    [(12, 20), (35, 43), (55, 63), (90, 98), (120, 128)],
)
def test_update_attributes(fck, fcm):
    """Test update_attributes function"""
    c = create_concrete(fck=fck, design_code='mc2010')
    c.update_attributes({'fcm': fcm})
    # Test a warning is raised when a not valid key is inputted
    with pytest.warns(UserWarning):
        c.update_attributes({'not_valid_key': fcm})

    assert c.fcm is not None
    assert c.fcm == fcm


def test_not_implemented_existing_error():
    """Test existing not implemented function"""
    with pytest.raises(NotImplementedError):
        create_concrete(fck=25, design_code='mc2010', existing=True)


# Series of tests using ConcreteMC2010 class


fck_parametrized = pytest.mark.parametrize('fck', [20, 25, 30, 35, 40])


@fck_parametrized
def test_fck_getter(fck):
    """
    Test fck getter
    """
    c = ConcreteMC2010(fck)

    assert c.fck == fck


@fck_parametrized
def test_fck_setter(fck):
    """
    Test fck setter
    """
    c = ConcreteMC2010(fck)
    c.fck = fck + 5

    assert c.fck == fck + 5


def test_properties_initialized_to_none():
    """
    Test if a ConcreteMC2010 when created has the attributes set to None
    """
    c = ConcreteMC2010(fck=25)

    assert c._fcm is None
    assert c._fctm is None
    assert c._fctkmin is None
    assert c._fctkmax is None
    assert c._Gf is None


def test_reset_properties():
    """
    Test _reset_attributes function
    """
    c = ConcreteMC2010(fck=25)

    c._reset_attributes()

    assert c._fcm is None
    assert c._fctm is None
    assert c._fctkmin is None
    assert c._fctkmax is None
    assert c._Gf is None


fcm_parametrized = pytest.mark.parametrize(
    'test_input, expected',
    [(12, 20), (35, 43), (55, 63), (90, 98), (120, 128)],
)


@fcm_parametrized
def test_fcm_getter(test_input, expected):
    """Test the fcm getter."""
    c = ConcreteMC2010(fck=test_input)
    assert math.isclose(c.fcm, expected)


@fcm_parametrized
def test_fcm_setter(test_input, expected):
    """Test the fcm setter."""
    c = ConcreteMC2010(fck=test_input)
    c.fcm = expected

    assert math.isclose(c.fcm, expected)


@pytest.mark.parametrize(
    'test_input',
    [12, 35, 55, 90, 120],
)
def test_fcm_setter_exception(test_input):
    """Test the fcm setter with a wrong value."""
    c = ConcreteMC2010(fck=test_input)
    with pytest.raises(ValueError):
        c.fcm = test_input - 1


fctm_parmetrized = pytest.mark.parametrize(
    'test_input, expected',
    [
        (12, 1.6),
        (16, 1.9),
        (20, 2.2),
        (25, 2.6),
        (30, 2.9),
        (35, 3.2),
        (40, 3.5),
        (45, 3.8),
        (50, 4.1),
        (55, 4.2),
        (60, 4.4),
        (70, 4.6),
        (80, 4.8),
        (90, 5.0),
        (100, 5.2),
        (110, 5.4),
        (120, 5.6),
    ],
)


@fctm_parmetrized
def test_fctm_getter(test_input, expected):
    """Test the fctm getter function."""
    c = ConcreteMC2010(fck=test_input)
    assert math.isclose(c.fctm, expected, rel_tol=0.02)


@fctm_parmetrized
def test_fctm_setter(test_input, expected):
    """Test the fctm setter function."""
    c = ConcreteMC2010(fck=test_input)
    c.fctm = expected

    assert math.isclose(c.fctm, expected)


@pytest.mark.parametrize('test_input', [10, 15, 20, 25, 30, 35])
def test_fctm_setter_warning(test_input):
    """
    Test the fctm setter function.
    Check that a warning is raised when
    trying to set a value higher than 0.5 times fck"""
    c = ConcreteMC2010(fck=test_input)

    with pytest.warns(UserWarning):
        c.fctm = test_input * 0.5 + 1


fctkmin_parametrized = pytest.mark.parametrize(
    'test_input, expected',
    [
        (12, 1.1),
        (16, 1.3),
        (20, 1.5),
        (25, 1.8),
        (30, 2.0),
        (35, 2.2),
        (40, 2.5),
        (45, 2.7),
        (50, 2.9),
        (55, 3.0),
        (60, 3.1),
        (70, 3.2),
        (80, 3.4),
        (90, 3.5),
        (100, 3.7),
        (110, 3.8),
        (120, 3.9),
    ],
)


@fctkmin_parametrized
def test_fctkmin_getter(test_input, expected):
    """Test the fctkmin getter function."""
    c = ConcreteMC2010(fck=test_input)
    assert math.isclose(c.fctkmin, expected, rel_tol=0.031)


@fctkmin_parametrized
def test_fctkmin_setter(test_input, expected):
    """Test the fctkmin setter function."""
    c = ConcreteMC2010(fck=test_input)
    c.fctkmin = expected

    assert math.isclose(c.fctkmin, expected)


fctkmax_parmetrized = pytest.mark.parametrize(
    'test_input, expected',
    [
        (12, 2.0),
        (16, 2.5),
        (20, 2.9),
        (25, 3.3),
        (30, 3.8),
        (35, 4.2),
        (40, 4.6),
        (45, 4.8),
        (50, 5.3),
        (55, 5.5),
        (60, 5.7),
        (70, 6.0),
        (80, 6.3),
        (90, 6.6),
        (100, 6.8),
        (110, 7.0),
        (120, 7.2),
    ],
)


@fctkmax_parmetrized
def test_fctkmax_getter(test_input, expected):
    """Test the fctkmax getter function."""
    c = ConcreteMC2010(fck=test_input)
    assert math.isclose(c.fctkmax, expected, rel_tol=0.028)


@fctkmax_parmetrized
def test_fctkmax_setter(test_input, expected):
    """Test the fctkmax setter function."""
    c = ConcreteMC2010(fck=test_input)
    c.fctkmax = expected

    assert math.isclose(c.fctkmax, expected)


gf_parametrized = pytest.mark.parametrize(
    'test_input, expected',
    [
        (12, 125.172),
        (35, 143.664),
        (55, 153.888),
        (90, 166.626),
        (120, 174.832),
    ],
)


@gf_parametrized
def test_Gf_getter(test_input, expected):
    """Test the Gf getter function."""
    c = ConcreteMC2010(fck=test_input)
    assert math.isclose(c.Gf, expected, rel_tol=1e-5)


@gf_parametrized
def test_Gf_setter(test_input, expected):
    """Test the Gf setter function."""
    c = ConcreteMC2010(fck=test_input)
    c.Gf = expected

    assert math.isclose(c.Gf, expected)
