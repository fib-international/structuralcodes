"""Test module for testing the EC2_2023 Concrete."""

import math

import pytest

from structuralcodes.materials.concrete import ConcreteEC2_2023


@pytest.fixture
def default_concrete():
    """Returns default concrete."""
    return ConcreteEC2_2023(fck=30)


def test_initialization(default_concrete):
    """Test for initialization and attribute setting."""
    assert default_concrete.fck == 30
    assert default_concrete.name.startswith('C')
    assert default_concrete.density == 2400.0
    assert default_concrete._kE == 9500
    assert default_concrete._strength_dev_class == 'cn'
    assert default_concrete.gamma_c == 1.5


def test_invalid_strength_dev_class():
    """Test for invalid strength development class."""
    with pytest.raises(ValueError):
        ConcreteEC2_2023(fck=30, strength_dev_class='invalid_class')


def test_invalid_gamma_c():
    """Test for invalid sec. coefficient."""
    with pytest.raises(ValueError):
        ConcreteEC2_2023(fck=30, gamma_c=-1)


def test_fcm_property(default_concrete):
    """Test for fcm property."""
    expected_fcm = default_concrete.fck + 8
    assert default_concrete.fcm == expected_fcm


def test_fcm_specified_error():
    """Test for fcm setter with invalid value."""
    fck = 30
    with pytest.raises(ValueError):
        ConcreteEC2_2023(fck=fck, fcm=fck - 1)


def test_fcm_specified(default_concrete):
    """Test the setter for the fcm property."""
    expected_fcm = default_concrete.fcm + 1
    new_concrete = ConcreteEC2_2023(fck=default_concrete.fck, fcm=expected_fcm)

    assert math.isclose(new_concrete.fcm, expected_fcm)


def test_fcd_default_parameters(default_concrete):
    """Test for fcd method with default parameters."""
    expected_fcd_default = 20
    assert default_concrete.fcd() == expected_fcd_default


def test_fcd_custom_parameters(default_concrete):
    """Test for fcd method with custom parameters."""
    custom_t_ref = 15
    custom_t0 = 70
    custom_fck_ref = 30
    expected_fcd_custom = 17
    assert (
        default_concrete.fcd(
            t_ref=custom_t_ref, t0=custom_t0, fck_ref=custom_fck_ref
        )
        == expected_fcd_custom
    )


def test_fcd_invalid_parameters(default_concrete):
    """Test for fcd method with invalid parameters."""
    invalid_t_ref = -1
    with pytest.raises(ValueError):
        default_concrete.fcd(t_ref=invalid_t_ref)


def test_fctm_property(default_concrete):
    """Test for fctm property."""
    expected_fctm = 2.896
    assert math.isclose(default_concrete.fctm, expected_fctm, rel_tol=0.001)


def test_fctm_specified(default_concrete):
    """Test for custom fctm value."""
    custom_fctm = 3.5
    new_concrete = ConcreteEC2_2023(fck=default_concrete.fck, fctm=custom_fctm)
    assert new_concrete.fctm == custom_fctm


def test_fctk_5_property(default_concrete):
    """Test for fctk_5 property."""
    expected_fctk_5 = 2.028
    assert math.isclose(
        default_concrete.fctk_5, expected_fctk_5, rel_tol=0.001
    )


def test_fctk_95_property(default_concrete):
    """Test for fctk_95 property."""
    expected_fctk_95 = 3.765
    assert math.isclose(
        default_concrete.fctk_95, expected_fctk_95, rel_tol=0.001
    )


def test_Ecm_default_and_specified(default_concrete):
    """Test for Ecm property."""
    expected_Ecm = 31938.766
    assert math.isclose(default_concrete.Ecm, expected_Ecm, rel_tol=0.001)
    custom_Ecm = 35000  # Test the override of Ecm
    new_concrete = ConcreteEC2_2023(fck=default_concrete.fck, Ecm=custom_Ecm)
    assert new_concrete.Ecm == custom_Ecm


def test_fctd_default_parameter(default_concrete):
    """Test for fctd method with default parameter."""
    expected_fctd_default = 1.081
    assert math.isclose(
        default_concrete.fctd(), expected_fctd_default, rel_tol=0.001
    )


def test_fctd_custom_parameter(default_concrete):
    """Test for fctd method with custom parameter."""
    custom_t_ref = 45
    expected_fctd_custom = 0.9461
    assert math.isclose(
        default_concrete.fctd(t_ref=custom_t_ref),
        expected_fctd_custom,
        rel_tol=0.001,
    )


def test_fctd_invalid_parameter(default_concrete):
    """Test for fctd method with invalid parameter."""
    invalid_t_ref = -5
    with pytest.raises(ValueError):
        default_concrete.fctd(t_ref=invalid_t_ref)


@pytest.mark.parametrize(
    'strain_limit',
    ('eps_c1', 'eps_cu1', 'eps_c2', 'eps_cu2'),
)
def test_strain_limits_specified_warning(strain_limit):
    """Test specifying strain limits with a wrong value."""
    kwargs = {strain_limit: 0.15}
    with pytest.warns(UserWarning):
        ConcreteEC2_2023(fck=45, **kwargs)


def test_eps_c1_property(default_concrete):
    """Test for eps_c1 property."""
    expected_eps_c1 = 0.002353
    assert math.isclose(
        default_concrete.eps_c1, expected_eps_c1, rel_tol=0.001
    )


def test_eps_cu1_property(default_concrete):
    """Test for eps_cu1 property."""
    expected_eps_cu1 = 0.0035
    assert math.isclose(
        default_concrete.eps_cu1, expected_eps_cu1, rel_tol=0.001
    )


def test_k_sargin(default_concrete):
    """Test for k property for Sargin law."""
    expected = (
        1.05
        * default_concrete.Ecm
        * default_concrete.eps_c1
        / default_concrete.fcm
    )
    assert math.isclose(default_concrete.k_sargin, expected, rel_tol=0.001)


def test_k_specified_warning():
    """Test specifying k_sargin with a wrong value."""
    with pytest.raises(ValueError):
        ConcreteEC2_2023(fck=45, k_sargin=-1.0)


def test_eps_c2_property(default_concrete):
    """Test for eps_c2 property."""
    expected = 0.002
    assert math.isclose(default_concrete.eps_c2, expected)


def test_eps_cu2_property(default_concrete):
    """Test for eps_cu2 property."""
    expected = 0.0035
    assert math.isclose(default_concrete.eps_cu2, expected)


def test_setter_properties(default_concrete):
    """Test for setting properties to custom values."""
    peak_strain = 2.0e-3
    ultimate_strain = 3.5e-3
    exponent = 1.5
    fctk_5 = 2.1
    fctk_95 = 2.3
    concrete = ConcreteEC2_2023(
        fck=default_concrete.fck,
        fctk_5=fctk_5,
        fctk_95=fctk_95,
        eps_c1=peak_strain,
        eps_c2=peak_strain,
        eps_cu1=ultimate_strain,
        eps_cu2=ultimate_strain,
        n_parabolic_rectangular=exponent,
        k_sargin=exponent,
    )
    assert math.isclose(concrete.eps_c1, peak_strain)
    assert math.isclose(concrete.eps_c2, peak_strain)
    assert math.isclose(concrete.eps_cu1, ultimate_strain)
    assert math.isclose(concrete.eps_cu2, ultimate_strain)
    assert math.isclose(concrete.n_parabolic_rectangular, exponent)
    assert math.isclose(concrete.k_sargin, exponent)
    assert math.isclose(concrete.fctk_5, fctk_5)
    assert math.isclose(concrete.fctk_95, fctk_95)


def test_n_specified_error():
    """Test specifying n_parabolic_rectangular with a wrong value."""
    with pytest.raises(ValueError):
        ConcreteEC2_2023(fck=45, n_parabolic_rectangular=-1)


def test_n_specified_warning():
    """Test specifying n_parabolic_rectangular with a wrong value."""
    with pytest.warns(UserWarning):
        ConcreteEC2_2023(fck=45, n_parabolic_rectangular=6)
