"""Tests for concrete material properties of ACI 318-19."""

import math

import pytest

from structuralcodes.codes import aci318_19
from structuralcodes.codes.aci318_19 import _concrete_material_properties


@pytest.mark.parametrize(
    'fc, expected',
    [
        (21, 4700 * math.sqrt(21)),
        (28, 4700 * math.sqrt(28)),
        (35, 4700 * math.sqrt(35)),
        (42, 4700 * math.sqrt(42)),
        (55, 4700 * math.sqrt(55)),
    ],
)
def test_Ec_normalweight(fc, expected):
    """Test Ec for normalweight concrete."""
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
        (35, 0.85 - 0.20 / 27 * (35 - 28)),
        (41.5, 0.85 - 0.20 / 27 * (41.5 - 28)),
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


def test_eps_c0():
    """Test default peak concrete strain."""
    assert _concrete_material_properties.eps_c0() == 0.002


@pytest.mark.parametrize(
    'wc, expected',
    [
        (1600, 0.75),
        (1800, 0.0075 * 1800 / 16.01846337396014),
        (2163, 1.0),
    ],
)
def test_lambda_factor(wc, expected):
    """Test lightweight modification factor."""
    assert math.isclose(
        _concrete_material_properties.lambda_factor(wc), expected
    )


def test_lambda_factor_invalid():
    """Test lambda_factor raises for invalid density."""
    with pytest.raises(ValueError):
        _concrete_material_properties.lambda_factor(0)


def test_alpha1():
    """Test stress block intensity factor."""
    assert _concrete_material_properties.alpha1() == 0.85


def test_aci318_19_unit_conversions():
    """Test US customary conversion helpers exposed by aci318_19."""
    assert math.isclose(aci318_19.psi_to_mpa(5000), 34.473786465841806)
    assert math.isclose(aci318_19.mpa_to_psi(34.473786465841806), 5000)
    assert math.isclose(aci318_19.ksi_to_mpa(60), 413.6854375901017)
    assert math.isclose(aci318_19.mpa_to_ksi(413.6854375901017), 60)
    assert math.isclose(aci318_19.pcf_to_kg_per_m3(150), 2402.769506094021)
    assert math.isclose(aci318_19.kg_per_m3_to_pcf(2402.769506094021), 150)
    assert math.isclose(aci318_19.in_to_mm(12), 304.8)
    assert math.isclose(aci318_19.mm_to_in(304.8), 12)
    assert math.isclose(aci318_19.in2_to_mm2(1), 645.16)
    assert math.isclose(aci318_19.mm2_to_in2(645.16), 1)
    assert math.isclose(aci318_19.in4_to_mm4(1), 416231.4256)
    assert math.isclose(aci318_19.mm4_to_in4(416231.4256), 1)
    assert math.isclose(aci318_19.kip_to_n(1), 4448.2216152605)
    assert math.isclose(aci318_19.n_to_kip(4448.2216152605), 1)
    assert math.isclose(aci318_19.kip_in_to_nmm(1), 4448.2216152605 * 25.4)
    assert math.isclose(aci318_19.nmm_to_kip_in(4448.2216152605 * 25.4), 1)
