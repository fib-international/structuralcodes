"""Tests for the ConcreteACI318 material class."""

import math

import pytest

import structuralcodes
from structuralcodes.materials.concrete import (
    ConcreteACI318,
    create_concrete,
)


@pytest.fixture(autouse=True)
def _reset_design_code():
    """Reset the global design code after each test."""
    yield
    structuralcodes.set_design_code(None)


def test_create_via_factory():
    """Test creating concrete via the factory function."""
    c = create_concrete(fck=28, design_code='aci318')
    assert isinstance(c, ConcreteACI318)


def test_default_name():
    """Test default name generation."""
    structuralcodes.set_design_code('aci318')
    c = create_concrete(fck=28)
    assert c.name == 'C28'


def test_fc_property():
    """Test the fc alias for fck."""
    c = ConcreteACI318(fck=28)
    assert c.fc == 28
    assert c.fck == 28


def test_gamma_c_default():
    """Test default gamma_c is 1.0 for ACI."""
    c = ConcreteACI318(fck=28)
    assert c.gamma_c == 1.0


def test_Ec_property():
    """Test Ec is derived from fc."""
    c = ConcreteACI318(fck=28)
    expected = 2320**1.5 * 0.043 * math.sqrt(28)
    assert math.isclose(c.Ec, expected, rel_tol=5e-3)


def test_Ec_specified():
    """Test Ec can be manually overridden."""
    c = ConcreteACI318(fck=28, Ec=30000)
    assert c.Ec == 30000


def test_fr_property():
    """Test fr is derived from fc."""
    c = ConcreteACI318(fck=28)
    expected = 0.62 * math.sqrt(28)
    assert math.isclose(c.fr, expected)


def test_fr_specified():
    """Test fr can be manually overridden."""
    c = ConcreteACI318(fck=28, fr=5.0)
    assert c.fr == 5.0


def test_beta1_property():
    """Test beta1 is derived from fc."""
    c = ConcreteACI318(fck=28)
    assert c.beta1 == 0.85


def test_eps_cu_property():
    """Test eps_cu returns 0.003."""
    c = ConcreteACI318(fck=28)
    assert c.eps_cu == 0.003


def test_alpha1_property():
    """Test alpha1 returns 0.85."""
    c = ConcreteACI318(fck=28)
    assert c.alpha1 == 0.85


def test_fcd():
    """Test fcd = alpha1 * fc / gamma_c."""
    c = ConcreteACI318(fck=28)
    assert math.isclose(c.fcd(), 0.85 * 28)


def test_fct_property():
    """Test fct is derived from fc."""
    c = ConcreteACI318(fck=28)
    expected = 0.56 * math.sqrt(28)
    assert math.isclose(c.fct, expected)


def test_constitutive_law_elastic():
    """Test elastic constitutive law creation."""
    c = ConcreteACI318(fck=28, constitutive_law='elastic')
    assert c.constitutive_law is not None


def test_constitutive_law_parabolarectangle():
    """Test parabola-rectangle constitutive law creation."""
    c = ConcreteACI318(fck=28, constitutive_law='parabolarectangle')
    assert c.constitutive_law is not None


def test_constitutive_law_bilinear():
    """Test bilinear compression constitutive law creation."""
    c = ConcreteACI318(fck=28, constitutive_law='bilinearcompression')
    assert c.constitutive_law is not None
