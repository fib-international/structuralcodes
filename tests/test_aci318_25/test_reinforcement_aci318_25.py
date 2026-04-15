"""Tests for the ReinforcementACI318_25 material class."""

import math
import sys

import pytest

# Mock triangle module to avoid optional dependency issues
sys.modules['triangle'] = type(sys)('triangle')

from structuralcodes.codes import set_design_code  # noqa: E402
from structuralcodes.materials.constitutive_laws import (  # noqa: E402
    Elastic,
    ElasticPlastic,
)
from structuralcodes.materials.reinforcement import (  # noqa: E402
    ReinforcementACI318_25,
    create_reinforcement,
)


@pytest.fixture(autouse=True)
def _reset_design_code():
    """Reset global design code before and after each test."""
    set_design_code(None)
    yield
    set_design_code(None)


def _make_gr60(**kwargs):
    """Helper to create a Grade 60 reinforcement instance."""
    return ReinforcementACI318_25(
        fyk=420, Es=200000, ftk=550, epsuk=0.05, **kwargs
    )


class TestConstruction:
    """Tests for basic construction of ReinforcementACI318_25."""

    def test_basic(self):
        """Test basic construction with explicit parameters."""
        r = _make_gr60()
        assert r.fyk == 420
        assert r.Es == 200000
        assert r.ftk == 550
        assert r.epsuk == 0.05
        assert isinstance(r, ReinforcementACI318_25)

    def test_default_name(self):
        """Test default name generation."""
        r = _make_gr60()
        assert r.name == 'Reinforcement420'

    def test_custom_name(self):
        """Test custom name."""
        r = _make_gr60(name='Grade60')
        assert r.name == 'Grade60'

    def test_from_grade_60(self):
        """Test construction from ASTM A615 Grade 60."""
        r = ReinforcementACI318_25.from_grade('60')
        assert math.isclose(r.fyk, 420.0)
        assert math.isclose(r.ftk, 550.0)
        assert math.isclose(r.epsuk, 0.05)

    def test_from_grade_40(self):
        """Test construction from ASTM A615 Grade 40."""
        r = ReinforcementACI318_25.from_grade('40')
        assert math.isclose(r.fyk, 280.0)
        assert math.isclose(r.ftk, 420.0)

    def test_from_grade_80(self):
        """Test construction from ASTM A615 Grade 80."""
        r = ReinforcementACI318_25.from_grade('80')
        assert math.isclose(r.fyk, 550.0)
        assert math.isclose(r.ftk, 690.0)

    def test_from_grade_100(self):
        """Test construction from ASTM A615 Grade 100."""
        r = ReinforcementACI318_25.from_grade('100')
        assert math.isclose(r.fyk, 690.0)
        assert math.isclose(r.ftk, 860.0)

    def test_from_grade_invalid(self):
        """Test that an invalid grade raises ValueError."""
        with pytest.raises(ValueError):
            ReinforcementACI318_25.from_grade('999')


class TestProperties:
    """Tests for material properties."""

    def test_gamma_s_default(self):
        """Test that gamma_s defaults to 1.0 for ACI."""
        r = _make_gr60()
        assert r.gamma_s == 1.0

    def test_gamma_s_custom(self):
        """Test custom gamma_s override."""
        r = _make_gr60(gamma_s=0.9)
        assert math.isclose(r.gamma_s, 0.9)

    def test_fyd(self):
        """Test design yield strength = fyk / gamma_s = 420 / 1.0."""
        r = _make_gr60()
        assert math.isclose(r.fyd(), 420.0)

    def test_ftd(self):
        """Test design ultimate strength = ftk / gamma_s = 550 / 1.0."""
        r = _make_gr60()
        assert math.isclose(r.ftd(), 550.0)

    def test_epsud(self):
        """Test design ultimate strain = epsuk = 0.05."""
        r = _make_gr60()
        assert math.isclose(r.epsud(), 0.05)

    def test_epsyd(self):
        """Test design yield strain = fyd / Es."""
        r = _make_gr60()
        expected = 420.0 / 200000.0
        assert math.isclose(r.epsyd, expected)


class TestConstitutiveLaws:
    """Tests for constitutive law creation."""

    def test_elastic(self):
        """Test elastic constitutive law."""
        r = _make_gr60(constitutive_law='elastic')
        assert isinstance(r.constitutive_law, Elastic)
        assert math.isclose(r.constitutive_law._E, r.Es)

    def test_elasticperfectlyplastic(self):
        """Test elastic perfectly plastic constitutive law (default)."""
        r = _make_gr60(constitutive_law='elasticperfectlyplastic')
        assert isinstance(r.constitutive_law, ElasticPlastic)
        assert math.isclose(r.constitutive_law._E, r.Es)
        assert math.isclose(r.constitutive_law._fy, r.fyd())
        assert math.isclose(r.constitutive_law._Eh, 0.0)
        assert math.isclose(r.constitutive_law._eps_su, r.epsud())

    def test_elasticplastic(self):
        """Test elastic plastic constitutive law with strain hardening."""
        r = _make_gr60(constitutive_law='elasticplastic')
        assert isinstance(r.constitutive_law, ElasticPlastic)
        assert math.isclose(r.constitutive_law._E, r.Es)
        assert math.isclose(r.constitutive_law._fy, r.fyd())
        assert math.isclose(r.constitutive_law._eps_su, r.epsud())
        # Eh should be positive (strain hardening)
        assert r.constitutive_law._Eh > 0

    def test_default_constitutive_law(self):
        """Test that default constitutive law is elasticperfectlyplastic."""
        r = _make_gr60()
        assert isinstance(r.constitutive_law, ElasticPlastic)
        # elasticperfectlyplastic has Eh=0
        assert math.isclose(r.constitutive_law._Eh, 0.0)

    def test_invalid_constitutive_law(self):
        """Test that an invalid constitutive law raises ValueError."""
        with pytest.raises(ValueError):
            _make_gr60(constitutive_law='parabolarectangle')


class TestFactory:
    """Tests for the create_reinforcement factory function."""

    def test_create_via_factory(self):
        """Test creating reinforcement via factory with design_code string."""
        r = create_reinforcement(
            fyk=420,
            Es=200000,
            ftk=550,
            epsuk=0.05,
            design_code='aci318_25',
        )
        assert isinstance(r, ReinforcementACI318_25)
        assert math.isclose(r.fyk, 420.0)

    def test_create_via_global_code(self):
        """Test creating reinforcement via globally set design code."""
        set_design_code('aci318_25')
        r = create_reinforcement(fyk=420, Es=200000, ftk=550, epsuk=0.05)
        assert isinstance(r, ReinforcementACI318_25)
        assert math.isclose(r.fyk, 420.0)
