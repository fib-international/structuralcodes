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


fck_parametrized = pytest.mark.parametrize('fck', [20, 25, 30, 35, 40])


@fck_parametrized
def test_fck_getter(fck):
    """Test fck getter."""
    c = ConcreteEC2_2004(fck)

    assert c.fck == fck


def test_properties_initialized_to_none():
    """Test if a ConcreteEC2_2004 has the attributes set to None."""
    c = ConcreteEC2_2004(fck=25)

    assert c._fcm is None
    assert c._fctm is None
    assert c._fctk_5 is None
    assert c._fctk_95 is None
    assert c._Ecm is None
    assert c._eps_c1 is None
    assert c._eps_cu1 is None
    assert c._k_sargin is None
    assert c._eps_c2 is None
    assert c._eps_cu2 is None
    assert c._n_parabolic_rectangular is None
    assert c._eps_c3 is None
    assert c._eps_cu3 is None


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
def test_fcm_specified(test_input, expected):
    """Test specifying fcm."""
    c = ConcreteEC2_2004(fck=test_input, fcm=expected)

    assert math.isclose(c.fcm, expected)


@pytest.mark.parametrize(
    'test_input',
    [12, 35, 55, 90, 120],
)
def test_fcm_specified_exception(test_input):
    """Test the specifying fcm with a wrong value."""
    with pytest.raises(ValueError):
        ConcreteEC2_2004(fck=test_input, fcm=test_input - 1)


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
def test_fctm_specified(test_input, expected):
    """Test specifying fctm."""
    c = ConcreteEC2_2004(fck=test_input, fctm=expected)

    assert math.isclose(c.fctm, expected)


@pytest.mark.parametrize(
    'test_input',
    [12, 35, 55, 90],
)
def test_fctm_specified_warning(test_input):
    """Test specifying fctm with a wrong value."""
    with pytest.warns(UserWarning):
        ConcreteEC2_2004(fck=test_input, fctm=test_input * 0.5 + 1)


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
def test_fctk5_specified(test_input, expected):
    """Test specifying fctk_5."""
    c = ConcreteEC2_2004(fck=test_input, fctk_5=expected)

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
def test_fctk95_specified(test_input, expected):
    """Test specifying fctk_95."""
    c = ConcreteEC2_2004(fck=test_input, fctk_95=expected)

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
def test_Ecm_specified(test_input, expected):
    """Test specifying Ecm."""
    c = ConcreteEC2_2004(fck=test_input, Ecm=expected)

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
    concrete = ConcreteEC2_2004(fck=fck, alpha_cc=alpha_cc)

    # Act and assert
    assert math.isclose(
        concrete.fcd(),
        fcd,
        rel_tol=1e-4,
    )


@fck_parametrized
def test_eps_c1_getter(fck):
    """Test eps_c1 getter."""
    c = ConcreteEC2_2004(fck=fck)
    expected = ec2_2004.eps_c1(ec2_2004.fcm(fck=fck))
    assert math.isclose(c.eps_c1, expected, rel_tol=1e-6)


eps_c1_parametrized = pytest.mark.parametrize(
    'test_input, expected',
    [(12, 1.8), (35, 2.25), (55, 2.5), (90, 2.8)],
)


@eps_c1_parametrized
def test_eps_c1_specified(test_input, expected):
    """Test specifying eps_c1."""
    c = ConcreteEC2_2004(fck=test_input, eps_c1=expected * 1e-3)

    assert math.isclose(c.eps_c1, expected * 1e-3)


@pytest.mark.parametrize(
    'strain_limit',
    ('eps_c1', 'eps_cu1', 'eps_c2', 'eps_cu2', 'eps_c3', 'eps_cu3'),
)
def test_strain_limits_specified_warning(strain_limit):
    """Test specifying strain limits with a wrong value."""
    kwargs = {strain_limit: 0.15}
    with pytest.warns(UserWarning):
        ConcreteEC2_2004(fck=45, **kwargs)


@fck_parametrized
def test_eps_cu1_getter(fck):
    """Test eps_cu1 getter."""
    c = ConcreteEC2_2004(fck=fck)
    expected = ec2_2004.eps_cu1(fck=fck)
    assert math.isclose(c.eps_cu1, expected)


eps_cu1_parametrized = pytest.mark.parametrize(
    'test_input, expected',
    [(12, 3.5), (35, 3.5), (55, 3.2), (90, 2.8)],
)


@eps_cu1_parametrized
def test_eps_cu1_specified(test_input, expected):
    """Test specifying eps_cu1."""
    c = ConcreteEC2_2004(fck=test_input, eps_cu1=expected * 1e-3)

    assert math.isclose(c.eps_cu1, expected * 1e-3)


@fck_parametrized
def test_k_getter(fck):
    """Test k getter."""
    c = ConcreteEC2_2004(fck=fck)
    expected = ec2_2004.k_sargin(
        Ecm=ec2_2004.Ecm(ec2_2004.fcm(fck=fck)),
        fcm=ec2_2004.fcm(fck=fck),
        eps_c1=ec2_2004.eps_c1(ec2_2004.fcm(fck=fck)),
    )
    assert math.isclose(c.k_sargin, expected, rel_tol=1e-6)


k_parametrized = pytest.mark.parametrize(
    'test_input, expected',
    [(12, 1.4), (35, 1.6), (55, 1.3), (90, 1.2)],
)


@k_parametrized
def test_k_specified(test_input, expected):
    """Test specifying k_sargin."""
    c = ConcreteEC2_2004(fck=test_input, k_sargin=expected)

    assert math.isclose(c.k_sargin, expected)


def test_k_specified_warning():
    """Test specifying k_sargin with a wrong value."""
    with pytest.raises(ValueError):
        ConcreteEC2_2004(fck=45, k_sargin=-1.0)


@fck_parametrized
def test_eps_c2_getter(fck):
    """Test eps_c2 getter."""
    c = ConcreteEC2_2004(fck=fck)
    expected = ec2_2004.eps_c2(fck=fck)
    assert math.isclose(c.eps_c2, expected, rel_tol=1e-6)


eps_c2_parametrized = pytest.mark.parametrize(
    'test_input, expected',
    [(12, 2.0), (35, 2.0), (55, 2.2), (90, 2.6)],
)


@eps_c2_parametrized
def test_eps_c2_specified(test_input, expected):
    """Test specifying eps_c2."""
    c = ConcreteEC2_2004(fck=test_input, eps_c2=expected * 1e-3)

    assert math.isclose(c.eps_c2, expected * 1e-3)


@fck_parametrized
def test_eps_cu2_getter(fck):
    """Test eps_cu2 getter."""
    c = ConcreteEC2_2004(fck=fck)
    expected = ec2_2004.eps_cu2(fck=fck)
    assert math.isclose(c.eps_cu2, expected, rel_tol=1e-6)


eps_cu2_parametrized = pytest.mark.parametrize(
    'test_input, expected',
    [(12, 3.5), (35, 3.5), (55, 3.1), (90, 2.6)],
)


@eps_cu2_parametrized
def test_eps_cu2_specified(test_input, expected):
    """Test specifying eps_cu2."""
    c = ConcreteEC2_2004(fck=test_input, eps_cu2=expected * 1e-3)

    assert math.isclose(c.eps_cu2, expected * 1e-3)


@fck_parametrized
def test_n_getter(fck):
    """Test n getter."""
    c = ConcreteEC2_2004(fck=fck)
    expected = ec2_2004.n_parabolic_rectangular(fck=fck)
    assert math.isclose(c.n_parabolic_rectangular, expected, rel_tol=1e-6)


n_parametrized = pytest.mark.parametrize(
    'test_input, expected',
    [(12, 2), (35, 2), (55, 1.75), (90, 1.4)],
)


@n_parametrized
def test_n_specified(test_input, expected):
    """Test specifying n_parabolic_rettangular."""
    c = ConcreteEC2_2004(fck=test_input, n_parabolic_rectangular=expected)
    assert math.isclose(c.n_parabolic_rectangular, expected)


def test_n_specified_error():
    """Test specifying n_parabolic_rectangular with a wrong value."""
    with pytest.raises(ValueError):
        ConcreteEC2_2004(fck=45, n_parabolic_rectangular=-1)


def test_n_specified_warning():
    """Test specifying n_parabolic_rectangular with a wrong value."""
    with pytest.warns(UserWarning):
        ConcreteEC2_2004(fck=45, n_parabolic_rectangular=6)


@fck_parametrized
def test_eps_c3_getter(fck):
    """Test eps_c3 getter."""
    c = ConcreteEC2_2004(fck=fck)
    expected = ec2_2004.eps_c3(fck=fck)
    assert math.isclose(c.eps_c3, expected, rel_tol=1e-6)


eps_c3_parametrized = pytest.mark.parametrize(
    'test_input, expected',
    [(12, 1.75), (35, 1.75), (55, 1.8), (90, 2.3)],
)


@eps_c3_parametrized
def test_eps_c3_specified(test_input, expected):
    """Test specifying eps_c3."""
    c = ConcreteEC2_2004(fck=test_input, eps_c3=expected * 1e-3)

    assert math.isclose(c.eps_c3, expected * 1e-3)


@fck_parametrized
def test_eps_cu3_getter(fck):
    """Test eps_cu3 getter."""
    c = ConcreteEC2_2004(fck=fck)
    expected = ec2_2004.eps_cu3(fck=fck)
    assert math.isclose(c.eps_cu3, expected, rel_tol=1e-6)


eps_cu3_parametrized = pytest.mark.parametrize(
    'test_input, expected',
    [(12, 3.5), (35, 3.5), (55, 3.1), (90, 2.6)],
)


@eps_cu3_parametrized
def test_eps_cu3_specified(test_input, expected):
    """Test specifying eps_cu3."""
    c = ConcreteEC2_2004(fck=test_input, eps_cu3=expected * 1e-3)

    assert math.isclose(c.eps_cu3, expected * 1e-3)
