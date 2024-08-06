"""Tests for the concrete mc2010."""

import math

import pytest

from structuralcodes.codes import mc2010
from structuralcodes.materials.concrete import ConcreteMC2010, create_concrete

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
    """Test update_attributes function."""
    c = create_concrete(fck=fck, design_code='mc2010')
    c.update_attributes({'fcm': fcm})
    # Test a warning is raised when a not valid key is inputted
    with pytest.warns(UserWarning):
        c.update_attributes({'not_valid_key': fcm})

    assert c.fcm is not None
    assert c.fcm == fcm


def test_not_implemented_existing_error():
    """Test existing not implemented function."""
    with pytest.raises(NotImplementedError):
        create_concrete(fck=25, design_code='mc2010', existing=True)


# Series of tests using ConcreteMC2010 class


fck_parametrized = pytest.mark.parametrize('fck', [20, 25, 30, 35, 40])


@fck_parametrized
def test_fck_getter(fck):
    """Test fck getter."""
    c = ConcreteMC2010(fck)

    assert c.fck == fck


@fck_parametrized
def test_fck_setter(fck):
    """Test fck setter."""
    c = ConcreteMC2010(fck)
    c.fck = fck + 5

    assert c.fck == fck + 5


def test_properties_initialized_to_none():
    """Test if a ConcreteMC2010 when created has the attributes set to None."""
    c = ConcreteMC2010(fck=25)

    assert c._fcm is None
    assert c._fctm is None
    assert c._fctkmin is None
    assert c._fctkmax is None
    assert c._Gf is None
    assert c._Eci is None
    assert c._eps_c1 is None
    assert c._eps_cu1 is None
    assert c._k_sargin is None
    assert c._eps_c2 is None
    assert c._eps_cu2 is None
    assert c._n_parabolic_rectangular is None
    assert c._eps_c3 is None
    assert c._eps_cu3 is None


def test_reset_properties():
    """Test _reset_attributes function."""
    c = ConcreteMC2010(fck=25)

    _ = c.fcm
    _ = c.fctm
    _ = c.fctkmin
    _ = c.fctkmax
    _ = c.Eci
    _ = c.eps_c1
    _ = c.eps_cu1
    _ = c.k_sargin
    _ = c.eps_c2
    _ = c.eps_cu2
    _ = c.n_parabolic_rectangular
    _ = c.eps_c3
    _ = c.eps_cu3

    c._reset_attributes()

    assert c._fcm is None
    assert c._fctm is None
    assert c._fctkmin is None
    assert c._fctkmax is None
    assert c._Gf is None
    assert c._Eci is None
    assert c._eps_c1 is None
    assert c._eps_cu1 is None
    assert c._k_sargin is None
    assert c._eps_c2 is None
    assert c._eps_cu2 is None
    assert c._n_parabolic_rectangular is None
    assert c._eps_c3 is None
    assert c._eps_cu3 is None


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
    """Test the fctm setter function. Check that a warning is raised when
    trying to set a value higher than 0.5 times fck.
    """
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


def test_gamma_c():
    """Test the gamma_c property."""
    # Arrange
    concrete = ConcreteMC2010(45)

    # Assert
    assert math.isclose(concrete.gamma_c, 1.5)


@pytest.mark.parametrize(
    'fck, alpha_cc, expected',
    [
        (35, 0.85, 19.8333),
        (45, 0.85, 25.5),
        (90, 0.85, 51.0),
    ],
)
def test_fcd(fck, alpha_cc, expected):
    """Test calculating the design compressive strength."""
    # Arrange
    concrete = ConcreteMC2010(fck, alpha_cc=alpha_cc)

    # Act
    fcd = concrete.fcd()

    # Assert
    assert math.isclose(fcd, expected, rel_tol=10e-5)


@fck_parametrized
def test_eps_c1_getter(fck):
    """Test eps_c1 getter."""
    c = ConcreteMC2010(fck=fck)
    expected = mc2010.eps_c1(fck=fck)
    assert math.isclose(c.eps_c1, expected, rel_tol=1e-6)


eps_c1_parametrized = pytest.mark.parametrize(
    'test_input, expected',
    [(12, 1.8), (35, 2.25), (55, 2.5), (90, 2.8)],
)


@eps_c1_parametrized
def test_eps_c1_setter(test_input, expected):
    """Test the eps_c1 setter."""
    c = ConcreteMC2010(fck=test_input)
    c.eps_c1 = expected * 1e-3

    assert math.isclose(c.eps_c1, expected * 1e-3)


@fck_parametrized
def test_eps_cu1_getter(fck):
    """Test eps_cu1 getter."""
    c = ConcreteMC2010(fck=fck)
    expected = mc2010.eps_cu1(fck=fck)
    assert math.isclose(c.eps_cu1, expected, rel_tol=1e-6)


eps_cu1_parametrized = pytest.mark.parametrize(
    'test_input, expected',
    [(12, 3.5), (35, 3.5), (55, 3.2), (90, 2.8)],
)


@eps_cu1_parametrized
def test_eps_cu1_setter(test_input, expected):
    """Test the eps_cu1 setter."""
    c = ConcreteMC2010(fck=test_input)
    c.eps_cu1 = expected * 1e-3

    assert math.isclose(c.eps_cu1, expected * 1e-3)


@fck_parametrized
def test_k_getter(fck):
    """Test k getter."""
    c = ConcreteMC2010(fck=fck)
    expected = mc2010.k_sargin(fck=fck)
    assert math.isclose(c.k_sargin, expected, rel_tol=1e-6)


k_parametrized = pytest.mark.parametrize(
    'test_input, expected',
    [(12, 1.4), (35, 1.6), (55, 1.3), (90, 1.2)],
)


@k_parametrized
def test_k_setter(test_input, expected):
    """Test the k_sargin setter."""
    c = ConcreteMC2010(fck=test_input)
    c.k_sargin = expected

    assert math.isclose(c.k_sargin, expected)


@fck_parametrized
def test_eps_c2_getter(fck):
    """Test eps_c2 getter."""
    c = ConcreteMC2010(fck=fck)
    expected = mc2010.eps_c2(fck=fck)
    assert math.isclose(c.eps_c2, expected, rel_tol=1e-6)


eps_c2_parametrized = pytest.mark.parametrize(
    'test_input, expected',
    [(12, 2.0), (35, 2.0), (55, 2.2), (90, 2.6)],
)


@eps_c2_parametrized
def test_eps_c2_setter(test_input, expected):
    """Test the eps_c2 setter."""
    c = ConcreteMC2010(fck=test_input)
    c.eps_c2 = expected * 1e-3

    assert math.isclose(c.eps_c2, expected * 1e-3)


@fck_parametrized
def test_eps_cu2_getter(fck):
    """Test eps_cu2 getter."""
    c = ConcreteMC2010(fck=fck)
    expected = mc2010.eps_cu2(fck=fck)
    assert math.isclose(c.eps_cu2, expected, rel_tol=1e-6)


eps_cu2_parametrized = pytest.mark.parametrize(
    'test_input, expected',
    [(12, 3.5), (35, 3.5), (55, 3.1), (90, 2.6)],
)


@eps_cu2_parametrized
def test_eps_cu2_setter(test_input, expected):
    """Test the eps_cu2 setter."""
    c = ConcreteMC2010(fck=test_input)
    c.eps_cu2 = expected * 1e-3

    assert math.isclose(c.eps_cu2, expected * 1e-3)


@fck_parametrized
def test_n_getter(fck):
    """Test n getter."""
    c = ConcreteMC2010(fck=fck)
    expected = mc2010.n_parabolic_rectangular(fck=fck)
    assert math.isclose(c.n_parabolic_rectangular, expected, rel_tol=1e-6)


n_parametrized = pytest.mark.parametrize(
    'test_input, expected',
    [(12, 2), (35, 2), (55, 1.75), (90, 1.4)],
)


@n_parametrized
def test_n_setter(test_input, expected):
    """Test the n_parabolic_rettangular setter."""
    c = ConcreteMC2010(fck=test_input)
    c.n_parabolic_rectangular = expected
    assert math.isclose(c.n_parabolic_rectangular, expected)


@fck_parametrized
def test_eps_c3_getter(fck):
    """Test eps_c3 getter."""
    c = ConcreteMC2010(fck=fck)
    expected = mc2010.eps_c3(fck=fck)
    assert math.isclose(c.eps_c3, expected, rel_tol=1e-6)


eps_c3_parametrized = pytest.mark.parametrize(
    'test_input, expected',
    [(12, 1.75), (35, 1.75), (55, 1.8), (90, 2.3)],
)


@eps_c3_parametrized
def test_eps_c3_setter(test_input, expected):
    """Test the eps_c3 setter."""
    c = ConcreteMC2010(fck=test_input)
    c.eps_c3 = expected * 1e-3

    assert math.isclose(c.eps_c3, expected * 1e-3)


@fck_parametrized
def test_eps_cu3_getter(fck):
    """Test eps_cu3 getter."""
    c = ConcreteMC2010(fck=fck)
    expected = mc2010.eps_cu3(fck=fck)
    assert math.isclose(c.eps_cu3, expected, rel_tol=1e-6)


eps_cu3_parametrized = pytest.mark.parametrize(
    'test_input, expected',
    [(12, 3.5), (35, 3.5), (55, 3.1), (90, 2.6)],
)


@eps_cu3_parametrized
def test_eps_cu3_setter(test_input, expected):
    """Test the eps_cu3 setter."""
    c = ConcreteMC2010(fck=test_input)
    c.eps_cu3 = expected * 1e-3

    assert math.isclose(c.eps_cu3, expected * 1e-3)
