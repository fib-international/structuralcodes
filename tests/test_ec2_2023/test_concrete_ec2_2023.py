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
    assert default_concrete._strength_dev_class == 'CN'
    assert default_concrete._gamma_C == 1.5


def test_invalid_strength_dev_class():
    """Test for invalid strength development class."""
    with pytest.raises(ValueError):
        ConcreteEC2_2023(fck=30, strength_dev_class='invalid_class')


def test_invalid_gamma_C():
    """Test for invalid sec. coefficient."""
    with pytest.raises(ValueError):
        ConcreteEC2_2023(fck=30, gamma_C=-1)


def test_fcm_property(default_concrete):
    """Test for fcm property."""
    expected_fcm = default_concrete.fck + 8
    assert default_concrete.fcm == expected_fcm


def test_fcm_setter_invalid_value(default_concrete):
    """Test for fcm setter with invalid value."""
    with pytest.raises(ValueError):
        default_concrete.fcm = default_concrete.fck - 1


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


def test_fctm_setter(default_concrete):
    """Test for custom fctm value."""
    custom_fctm = 3.5
    default_concrete.fctm = custom_fctm
    assert default_concrete._fctm == custom_fctm


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


def test_Ecm_property_and_setter(default_concrete):
    """Test for Ecm property and its setter."""
    expected_Ecm = 31938.766
    assert math.isclose(default_concrete.Ecm, expected_Ecm, rel_tol=0.001)
    custom_Ecm = 35000  # Test the override of Ecm
    default_concrete.Ecm = custom_Ecm
    assert default_concrete._Ecm == custom_Ecm


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


def test_sigma_c_with_valid_strain(default_concrete):
    """Test for sigma_c method with valid strain."""
    valid_strain = 0.001
    expected_sigma_c = 25.830
    assert math.isclose(
        default_concrete.sigma_c(eps_c=valid_strain),
        expected_sigma_c,
        rel_tol=0.001,
    )


def test_sigma_c_with_invalid_strain(default_concrete):
    """Test for sigma_c method with invalid strain (negative)."""
    invalid_strain = -0.001
    with pytest.raises(ValueError):
        default_concrete.sigma_c(eps_c=invalid_strain)
