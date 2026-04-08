"""Tests for concrete material properties of ACI 318-19."""

import math

import pytest

from structuralcodes.codes.aci318 import _concrete_material_properties

WC_DEFAULT = 2320.0


@pytest.mark.parametrize(
    'fc, expected',
    [
        (21, WC_DEFAULT**1.5 * 0.043 * math.sqrt(21)),
        (28, WC_DEFAULT**1.5 * 0.043 * math.sqrt(28)),
        (35, WC_DEFAULT**1.5 * 0.043 * math.sqrt(35)),
        (42, WC_DEFAULT**1.5 * 0.043 * math.sqrt(42)),
        (55, WC_DEFAULT**1.5 * 0.043 * math.sqrt(55)),
    ],
)
def test_Ec_normalweight(fc, expected):
    """Test Ec for normalweight concrete (wc=2320 kg/m3)."""
    assert math.isclose(_concrete_material_properties.Ec(fc), expected)


def test_Ec_custom_wc():
    """Test Ec with a custom unit weight."""
    fc = 28
    wc = 1800.0
    expected = wc**1.5 * 0.043 * math.sqrt(fc)
    assert math.isclose(_concrete_material_properties.Ec(fc, wc=wc), expected)


def test_Ec_invalid_fc():
    """Test Ec raises for non-positive fc."""
    with pytest.raises(ValueError):
        _concrete_material_properties.Ec(-1)


def test_Ec_invalid_wc():
    """Test Ec raises for wc outside valid range."""
    with pytest.raises(ValueError):
        _concrete_material_properties.Ec(28, wc=1000)


@pytest.mark.parametrize(
    'fc, expected',
    [
        (21, 0.62 * math.sqrt(21)),
        (28, 0.62 * math.sqrt(28)),
        (35, 0.62 * math.sqrt(35)),
        (55, 0.62 * math.sqrt(55)),
    ],
)
def test_fr(fc, expected):
    """Test modulus of rupture."""
    assert math.isclose(_concrete_material_properties.fr(fc), expected)


def test_fr_lightweight():
    """Test modulus of rupture with lightweight factor."""
    fc = 28
    lambda_s = 0.75
    expected = 0.62 * lambda_s * math.sqrt(fc)
    assert math.isclose(
        _concrete_material_properties.fr(fc, lambda_s=lambda_s),
        expected,
    )


def test_fr_invalid_fc():
    """Test fr raises for non-positive fc."""
    with pytest.raises(ValueError):
        _concrete_material_properties.fr(-1)


def test_fr_invalid_lambda():
    """Test fr raises for invalid lambda_s."""
    with pytest.raises(ValueError):
        _concrete_material_properties.fr(28, lambda_s=1.5)


@pytest.mark.parametrize(
    'fc, expected',
    [
        (21, 0.85),
        (28, 0.85),
        (35, 0.80),
        (41.5, 0.85 - 0.05 * (41.5 - 28) / 7),
        (55, 0.65),
        (69, 0.65),
    ],
)
def test_beta1(fc, expected):
    """Test Whitney stress block depth factor."""
    assert math.isclose(
        _concrete_material_properties.beta1(fc), expected, rel_tol=1e-6
    )


def test_beta1_invalid_fc():
    """Test beta1 raises for non-positive fc."""
    with pytest.raises(ValueError):
        _concrete_material_properties.beta1(-1)


def test_eps_cu():
    """Test ultimate concrete strain."""
    assert _concrete_material_properties.eps_cu() == 0.003


@pytest.mark.parametrize(
    'concrete_type, expected',
    [
        ('normalweight', 1.0),
        ('sand-lightweight', 0.85),
        ('all-lightweight', 0.75),
    ],
)
def test_lambda_factor(concrete_type, expected):
    """Test lightweight modification factor."""
    assert (
        _concrete_material_properties.lambda_factor(concrete_type) == expected
    )


def test_lambda_factor_invalid():
    """Test lambda_factor raises for unknown type."""
    with pytest.raises(ValueError):
        _concrete_material_properties.lambda_factor('unknown')


@pytest.mark.parametrize(
    'fc, expected',
    [
        (21, 0.56 * math.sqrt(21)),
        (28, 0.56 * math.sqrt(28)),
        (35, 0.56 * math.sqrt(35)),
    ],
)
def test_fct(fc, expected):
    """Test splitting tensile strength."""
    assert math.isclose(_concrete_material_properties.fct(fc), expected)


def test_alpha1():
    """Test stress block intensity factor."""
    assert _concrete_material_properties.alpha1() == 0.85
