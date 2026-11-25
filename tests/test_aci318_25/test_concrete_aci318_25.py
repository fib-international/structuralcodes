"""Tests for the ConcreteACI318_25 material class."""

import math

import pytest

from structuralcodes.codes import aci318_25, set_design_code
from structuralcodes.materials.concrete import (
    ConcreteACI318_25,
    create_concrete,
)


@pytest.fixture(autouse=True)
def _reset_design_code():
    """Reset global design code before and after each test."""
    set_design_code(None)
    yield
    set_design_code(None)


FCK = 27.58  # 4000 psi in MPa


class TestConstruction:
    """Tests for basic construction of ConcreteACI318_25."""

    def test_basic(self):
        """Test basic construction with fck."""
        c = ConcreteACI318_25(fck=FCK)
        assert c.fck == FCK
        assert isinstance(c, ConcreteACI318_25)

    def test_fc_alias(self):
        """Test that fc is an alias for fck."""
        c = ConcreteACI318_25(fck=FCK)
        assert c.fc == c.fck

    def test_default_name(self):
        """Test default name generation."""
        c = ConcreteACI318_25(fck=FCK)
        assert c.name == 'C28'

    def test_custom_name(self):
        """Test custom name."""
        c = ConcreteACI318_25(fck=FCK, name='Custom4000')
        assert c.name == 'Custom4000'


class TestProperties:
    """Tests for material properties."""

    def test_gamma_c(self):
        """Test that gamma_c defaults to 1.0 for ACI."""
        c = ConcreteACI318_25(fck=FCK)
        assert c.gamma_c == 1.0

    def test_fcd(self):
        """Test design compressive strength = alpha1 * f'c."""
        c = ConcreteACI318_25(fck=FCK)
        expected = 0.85 * FCK
        assert math.isclose(c.fcd(), expected, rel_tol=1e-6)

    def test_Ec_computed(self):
        """Test Ec computed from code function."""
        c = ConcreteACI318_25(fck=FCK)
        expected = aci318_25.Ec(FCK, wc=2320.0)
        assert math.isclose(c.Ec, expected, rel_tol=1e-6)

    def test_Ec_override(self):
        """Test Ec with user-specified value."""
        custom_Ec = 25000.0
        c = ConcreteACI318_25(fck=FCK, Ec=custom_Ec)
        assert math.isclose(c.Ec, custom_Ec)

    def test_fr(self):
        """Test modulus of rupture computed from code function."""
        c = ConcreteACI318_25(fck=FCK)
        expected = aci318_25.fr(FCK, lambda_s=1.0)
        assert math.isclose(c.fr, expected, rel_tol=1e-6)

    def test_fr_override(self):
        """Test fr with user-specified value."""
        custom_fr = 3.5
        c = ConcreteACI318_25(fck=FCK, fr=custom_fr)
        assert math.isclose(c.fr, custom_fr)

    def test_fct(self):
        """Test splitting tensile strength."""
        c = ConcreteACI318_25(fck=FCK)
        expected = aci318_25.fct(FCK, lambda_s=1.0)
        assert math.isclose(c.fct, expected, rel_tol=1e-6)

    def test_beta1(self):
        """Test beta1 for fck <= 28 MPa."""
        c = ConcreteACI318_25(fck=FCK)
        assert math.isclose(c.beta1, 0.85, rel_tol=1e-6)

    def test_alpha1(self):
        """Test alpha1 = 0.85."""
        c = ConcreteACI318_25(fck=FCK)
        assert math.isclose(c.alpha1, 0.85, rel_tol=1e-6)

    def test_eps_cu(self):
        """Test ultimate concrete strain = 0.003."""
        c = ConcreteACI318_25(fck=FCK)
        assert math.isclose(c.eps_cu, 0.003, rel_tol=1e-6)


class TestConstitutiveLaws:
    """Tests for constitutive law creation."""

    def test_elastic(self):
        """Test elastic constitutive law."""
        c = ConcreteACI318_25(fck=FCK, constitutive_law='elastic')
        assert c._constitutive_law is not None

    def test_parabolarectangle(self):
        """Test parabolarectangle constitutive law with correct ultimate
        strain.
        """
        c = ConcreteACI318_25(fck=FCK, constitutive_law='parabolarectangle')
        assert c._constitutive_law is not None
        # Check that the ultimate strain is 0.003
        assert math.isclose(
            c._constitutive_law.get_ultimate_strain()[0],
            -0.003,
            rel_tol=1e-6,
        )

    def test_bilinearcompression(self):
        """Test bilinearcompression constitutive law."""
        c = ConcreteACI318_25(fck=FCK, constitutive_law='bilinearcompression')
        assert c._constitutive_law is not None


class TestFactory:
    """Tests for the create_concrete factory function."""

    def test_create_via_factory(self):
        """Test creating concrete via factory with design_code string."""
        c = create_concrete(fck=FCK, design_code='aci318_25')
        assert isinstance(c, ConcreteACI318_25)
        assert c.fck == FCK

    def test_create_via_global_code(self):
        """Test creating concrete via globally set design code."""
        set_design_code('aci318_25')
        c = create_concrete(fck=FCK)
        assert isinstance(c, ConcreteACI318_25)
        assert c.fck == FCK
