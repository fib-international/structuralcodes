"""Tests for the ConcreteACI318_19 material class."""

import math

import pytest

import structuralcodes
from structuralcodes.codes import aci318_19
from structuralcodes.materials.concrete import (
    ConcreteACI318_19,
    create_concrete,
)


@pytest.fixture(autouse=True)
def _reset_design_code():
    """Reset the global design code after each test."""
    yield
    structuralcodes.set_design_code(None)


def test_create_via_factory():
    """Test creating concrete via the factory function."""
    c = create_concrete(fck=28, design_code='aci318_19')
    assert isinstance(c, ConcreteACI318_19)


def test_default_name():
    """Test default name generation."""
    structuralcodes.set_design_code('aci318_19')
    c = create_concrete(fck=28)
    assert c.name == 'C28'


def test_fc_property():
    """Test the ACI fc property and inherited fck compatibility."""
    c = ConcreteACI318_19(fc=28)
    assert c.fc == 28
    assert c.fck == 28


def test_fc_constructor_alias():
    """Test direct ACI construction with fc."""
    c = ConcreteACI318_19(fc=28)
    assert c.fc == 28
    assert c.fck == 28


def test_fck_and_fc_constructor_alias_match():
    """Test fck and fc can both be provided when equal."""
    c = ConcreteACI318_19(fck=28, fc=28)
    assert c.fc == 28
    assert c.fck == 28


def test_fck_and_fc_constructor_alias_conflict():
    """Test fck and fc cannot conflict."""
    with pytest.raises(ValueError):
        ConcreteACI318_19(fck=28, fc=35)


def test_missing_fc_and_fck():
    """Test one compressive strength input is required."""
    with pytest.raises(ValueError):
        ConcreteACI318_19()


def test_from_psi():
    """Test creating concrete from US customary strength."""
    c = ConcreteACI318_19.from_psi(fc_psi=5000)
    expected_fc = aci318_19.psi_to_mpa(5000)
    assert c.name == 'C5000psi'
    assert math.isclose(c.fc, expected_fc)
    assert math.isclose(c.density, aci318_19.pcf_to_kg_per_m3(150))
    assert math.isclose(c.Ec, 4700 * math.sqrt(expected_fc))


def test_from_psi_customary_overrides():
    """Test customary overrides are converted to SI internals."""
    c = ConcreteACI318_19.from_psi(
        fc_psi=5000,
        density_pcf=145,
        wc_pcf=120,
        Ec_psi=4000000,
        fr_psi=500,
    )
    assert math.isclose(c.fc, aci318_19.psi_to_mpa(5000))
    assert math.isclose(c.density, aci318_19.pcf_to_kg_per_m3(145))
    assert math.isclose(c.Ec, aci318_19.psi_to_mpa(4000000))
    assert math.isclose(c.fr, aci318_19.psi_to_mpa(500))
    assert math.isclose(c.lambda_s, 0.0075 * 120)


def test_gamma_c_default():
    """Test default gamma_c is 1.0 for ACI."""
    c = ConcreteACI318_19(fc=28)
    assert c.gamma_c == 1.0


def test_Ec_property():
    """Test Ec is derived from fc."""
    c = ConcreteACI318_19(fc=28)
    expected = 4700 * math.sqrt(28)
    assert math.isclose(c.Ec, expected, rel_tol=5e-3)


def test_Ec_property_with_wc():
    """Test Ec can use ACI density-dependent expression."""
    c = ConcreteACI318_19(fc=28, wc=1800)
    expected = 1800**1.5 * 0.043 * math.sqrt(28)
    assert math.isclose(c.Ec, expected)


def test_Ec_specified():
    """Test Ec can be manually overridden."""
    c = ConcreteACI318_19(fc=28, Ec=30000)
    assert c.Ec == 30000


def test_fr_property():
    """Test fr is derived from fc."""
    c = ConcreteACI318_19(fc=28)
    expected = 0.62 * math.sqrt(28)
    assert math.isclose(c.fr, expected)


def test_fr_specified():
    """Test fr can be manually overridden."""
    c = ConcreteACI318_19(fc=28, fr=5.0)
    assert c.fr == 5.0


def test_lambda_s_property_from_wc():
    """Test lambda_s is derived from wc when not manually provided."""
    c = ConcreteACI318_19(fc=28, wc=1800)
    expected = 0.0075 * 1800 / 16.01846337396014
    assert math.isclose(c.lambda_s, expected)


def test_lambda_s_specified():
    """Test lambda_s can be manually overridden."""
    c = ConcreteACI318_19(fc=28, lambda_s=0.8)
    assert c.lambda_s == 0.8


def test_beta1_property():
    """Test beta1 is derived from fc."""
    c = ConcreteACI318_19(fc=28)
    assert c.beta1 == 0.85


def test_eps_cu_property():
    """Test eps_cu returns 0.003."""
    c = ConcreteACI318_19(fc=28)
    assert c.eps_cu == 0.003


def test_eps_c0_property():
    """Test eps_c0 returns 0.002."""
    c = ConcreteACI318_19(fc=28)
    assert c.eps_c0 == 0.002


def test_alpha1_property():
    """Test alpha1 returns 0.85."""
    c = ConcreteACI318_19(fc=28)
    assert c.alpha1 == 0.85


def test_fcd():
    """Test fcd = alpha1 * fc."""
    c = ConcreteACI318_19(fc=28)
    assert math.isclose(c.fcd(), 0.85 * 28)


def test_invalid_gamma_c():
    """Test gamma_c values other than 1.0 are rejected."""
    with pytest.raises(ValueError):
        ConcreteACI318_19(fc=28, gamma_c=1.5)


def test_invalid_fc():
    """Test ACI minimum concrete strength is enforced."""
    with pytest.raises(ValueError):
        ConcreteACI318_19(fc=16)


def test_constitutive_law_elastic():
    """Test elastic constitutive law creation."""
    c = ConcreteACI318_19(fc=28, constitutive_law='elastic')
    assert c.constitutive_law is not None


def test_constitutive_law_parabolarectangle():
    """Test parabola-rectangle constitutive law creation."""
    c = ConcreteACI318_19(fc=28, constitutive_law='parabolarectangle')
    assert c.constitutive_law is not None


def test_constitutive_law_bilinear():
    """Test bilinear compression constitutive law creation."""
    c = ConcreteACI318_19(fc=28, constitutive_law='bilinearcompression')
    assert c.constitutive_law is not None
