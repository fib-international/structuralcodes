"""Tests for ACI 318-25 one-way slab design rules (Ch. 7)."""

import math

import pytest

from structuralcodes.codes.aci318_25 import _one_way_slab as ows

# ---------------------------------------------------------------------------
# Shared test constants
# ---------------------------------------------------------------------------
SPAN = 6096.0   # mm  (20 ft)
B = 305.0       # mm  (12 in)
H = 254.0       # mm  (10 in)
FY_GR60 = 420.0  # MPa (~60 ksi)
FY_GR80 = 552.0  # MPa (~80 ksi)
D = 227.0       # mm


# ---------------------------------------------------------------------------
# TestMinThickness
# ---------------------------------------------------------------------------


class TestMinThickness:
    """Tests for min_thickness (Table 7.3.1.1)."""

    def test_simply_supported_gr60(self):
        """h = span/20 for simply-supported slab with Gr 60 rebar."""
        result = ows.min_thickness(SPAN, 'simply_supported', fy=FY_GR60)
        expected = SPAN / 20
        assert math.isclose(result, expected, rel_tol=1e-9)

    def test_one_end_continuous(self):
        """h = span/24 for one-end-continuous slab with Gr 60 rebar."""
        result = ows.min_thickness(SPAN, 'one_end_continuous', fy=FY_GR60)
        expected = SPAN / 24
        assert math.isclose(result, expected, rel_tol=1e-9)

    def test_both_ends_continuous(self):
        """h = span/28 for both-ends-continuous slab with Gr 60 rebar."""
        result = ows.min_thickness(SPAN, 'both_ends_continuous', fy=FY_GR60)
        expected = SPAN / 28
        assert math.isclose(result, expected, rel_tol=1e-9)

    def test_cantilever(self):
        """h = span/10 for cantilever slab with Gr 60 rebar."""
        result = ows.min_thickness(SPAN, 'cantilever', fy=FY_GR60)
        expected = SPAN / 10
        assert math.isclose(result, expected, rel_tol=1e-9)

    def test_fy_adjustment_gr80(self):
        """Sec. 7.3.1.1.1: Gr 80 thickness exceeds Gr 60 by fy factor."""
        h_base = SPAN / 20  # simply supported, no factor
        fy_psi = FY_GR80 * 145.038
        factor = 0.4 + fy_psi / 100000.0
        expected = h_base * factor
        result = ows.min_thickness(SPAN, 'simply_supported', fy=FY_GR80)
        assert math.isclose(result, expected, rel_tol=1e-9)

    def test_invalid_support_condition(self):
        """ValueError raised for an unknown support condition."""
        with pytest.raises(ValueError, match="Unknown support condition"):
            ows.min_thickness(SPAN, 'fixed_fixed')


# ---------------------------------------------------------------------------
# TestAsShrinkageTemperature
# ---------------------------------------------------------------------------


class TestAsShrinkageTemperature:
    """Tests for As_shrinkage_temperature (Sec. 24.4.3.2)."""

    def test_gr60_ratio(self):
        """420 MPa (~60 916 psi, just above 60 000 psi threshold).

        Falls in the >60 000 psi branch:
        ratio = max(0.0014, 0.0018 * 60 000 / fy_psi).
        """
        fy_psi = FY_GR60 * 145.038
        ratio = max(0.0014, 0.0018 * 60000.0 / fy_psi)
        expected = ratio * B * H
        result = ows.As_shrinkage_temperature(FY_GR60, B, H)
        assert math.isclose(result, expected, rel_tol=1e-9)


# ---------------------------------------------------------------------------
# TestMaxBarSpacingFlexure
# ---------------------------------------------------------------------------


class TestMaxBarSpacingFlexure:
    """Tests for max_bar_spacing_flexure (Sec. 7.7.2.3)."""

    def test_h150_governed_by_limit(self):
        """h=150 mm: 3*150=450 == 450, so result is 450."""
        assert math.isclose(ows.max_bar_spacing_flexure(150.0), 450.0, rel_tol=1e-9)

    def test_h200_governed_by_limit(self):
        """h=200 mm: 3*200=600 > 450, so result is capped at 450."""
        assert math.isclose(ows.max_bar_spacing_flexure(200.0), 450.0, rel_tol=1e-9)


# ---------------------------------------------------------------------------
# TestMaxBarSpacingShrinkage
# ---------------------------------------------------------------------------


class TestMaxBarSpacingShrinkage:
    """Tests for max_bar_spacing_shrinkage (Sec. 7.7.6.2.1)."""

    def test_h100_governed_by_limit(self):
        """h=100 mm: 5*100=500 > 450, so result is capped at 450."""
        assert math.isclose(ows.max_bar_spacing_shrinkage(100.0), 450.0, rel_tol=1e-9)

    def test_h80_governed_by_3h(self):
        """h=80 mm: 5*80=400 < 450, so result is 400."""
        assert math.isclose(ows.max_bar_spacing_shrinkage(80.0), 400.0, rel_tol=1e-9)


# ---------------------------------------------------------------------------
# TestShearCriticalSectionOffset
# ---------------------------------------------------------------------------


class TestShearCriticalSectionOffset:
    """Tests for shear_critical_section_offset (Sec. 7.4.3.2)."""

    def test_returns_d(self):
        """Critical section offset equals the effective depth d."""
        assert math.isclose(ows.shear_critical_section_offset(D), D, rel_tol=1e-9)
