"""Tests for ACI 318-25 reinforcement material property functions (Ch. 20)."""

import math

import pytest

from structuralcodes.codes.aci318_25 import (
    _reinforcement_material_properties as rmp,
)


class TestEs:
    """Tests for the modulus of elasticity Es (Sec. 20.2.2.2)."""

    def test_value(self):
        """Es must equal 200000.0 MPa."""
        assert rmp.Es() == 200000.0


class TestFyDesign:
    """Tests for the design yield strength fy_design."""

    def test_default_phi(self):
        """With default phi=1.0, result equals fy."""
        assert math.isclose(rmp.fy_design(420.0), 420.0, rel_tol=1e-9)

    def test_phi_0_9(self):
        """With phi=0.9, result equals 0.9 * fy."""
        assert math.isclose(rmp.fy_design(420.0, phi=0.9), 378.0, rel_tol=1e-9)

    def test_invalid_fy_zero(self):
        """fy = 0 should raise ValueError."""
        with pytest.raises(ValueError):
            rmp.fy_design(0.0)

    def test_invalid_fy_negative(self):
        """Negative fy should raise ValueError."""
        with pytest.raises(ValueError):
            rmp.fy_design(-420.0)

    def test_invalid_phi_zero(self):
        """phi = 0 should raise ValueError (must be > 0)."""
        with pytest.raises(ValueError):
            rmp.fy_design(420.0, phi=0.0)

    def test_invalid_phi_above_one(self):
        """phi > 1 should raise ValueError."""
        with pytest.raises(ValueError):
            rmp.fy_design(420.0, phi=1.1)


class TestEpsyd:
    """Tests for the design yield strain epsyd (Sec. 20.2.2.2)."""

    @pytest.mark.parametrize(
        'fy, expected',
        [
            (420.0, 420.0 / 200000.0),
            (280.0, 280.0 / 200000.0),
            (550.0, 550.0 / 200000.0),
        ],
    )
    def test_epsyd_parametric(self, fy, expected):
        """Test epsyd for typical reinforcement yield strengths."""
        assert math.isclose(rmp.epsyd(fy), expected, rel_tol=1e-9)

    def test_invalid_fy_zero(self):
        """fy = 0 should raise ValueError."""
        with pytest.raises(ValueError):
            rmp.epsyd(0.0)

    def test_invalid_fy_negative(self):
        """Negative fy should raise ValueError."""
        with pytest.raises(ValueError):
            rmp.epsyd(-280.0)


class TestReinforcementGradeProps:
    """Tests for reinforcement_grade_props (Table 20.2.1.3)."""

    @pytest.mark.parametrize(
        'grade, fy, fu',
        [
            ('40', 280.0, 420.0),
            ('60', 420.0, 550.0),
            ('80', 550.0, 690.0),
            ('100', 690.0, 860.0),
        ],
    )
    def test_grade_parametric(self, grade, fy, fu):
        """Test that all four grades return correct fy and fu."""
        props = rmp.reinforcement_grade_props(grade)
        assert math.isclose(props['fy'], fy, rel_tol=1e-9)
        assert math.isclose(props['fu'], fu, rel_tol=1e-9)

    def test_invalid_grade(self):
        """Unknown grade should raise ValueError."""
        with pytest.raises(ValueError):
            rmp.reinforcement_grade_props('75')
