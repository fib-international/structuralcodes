"""Tests for ACI 318-25 strength reduction factors (Ch. 21)."""

import math

import pytest

from structuralcodes.codes.aci318_25 import _strength_reduction as sr


class TestPhiShear:
    """Tests for phi_shear (Table 21.2.1(b))."""

    def test_value(self):
        """phi_shear must equal 0.75."""
        assert sr.phi_shear() == 0.75


class TestPhiTorsion:
    """Tests for phi_torsion (Table 21.2.1(c))."""

    def test_value(self):
        """phi_torsion must equal 0.75."""
        assert sr.phi_torsion() == 0.75


class TestPhiBearing:
    """Tests for phi_bearing (Table 21.2.1(d))."""

    def test_value(self):
        """phi_bearing must equal 0.65."""
        assert sr.phi_bearing() == 0.65


class TestPhiFlexure:
    """Tests for phi_flexure (Table 21.2.2)."""

    def test_tension_controlled_gr60(self):
        """Grade 60 (fy=420) with eps_t=0.010 is tension-controlled → 0.90."""
        assert math.isclose(sr.phi_flexure(0.010, 420.0), 0.90, rel_tol=1e-9)

    def test_tension_controlled_at_limit(self):
        """eps_t = eps_ty + 0.003 is at tension-controlled limit → 0.90."""
        eps_ty = 420.0 / 200000.0
        assert math.isclose(
            sr.phi_flexure(eps_ty + 0.003, 420.0), 0.90, rel_tol=1e-9
        )

    def test_compression_controlled_other(self):
        """eps_t = eps_ty with other transverse → compression-controlled."""
        eps_ty = 420.0 / 200000.0
        assert math.isclose(
            sr.phi_flexure(eps_ty, 420.0, transverse='other'),
            0.65,
            rel_tol=1e-9,
        )

    def test_compression_controlled_spiral(self):
        """eps_t = eps_ty with spiral transverse → compression-controlled."""
        eps_ty = 420.0 / 200000.0
        assert math.isclose(
            sr.phi_flexure(eps_ty, 420.0, transverse='spiral'),
            0.75,
            rel_tol=1e-9,
        )

    def test_transition_midpoint_other(self):
        """eps_t = eps_ty + 0.0015 (midpoint) with other → 0.775."""
        eps_ty = 420.0 / 200000.0
        expected = 0.65 + 0.25 * 0.0015 / 0.003  # = 0.775
        assert math.isclose(
            sr.phi_flexure(eps_ty + 0.0015, 420.0, transverse='other'),
            expected,
            rel_tol=1e-9,
        )

    def test_transition_midpoint_spiral(self):
        """eps_t = eps_ty + 0.0015 (midpoint) with spiral → 0.825."""
        eps_ty = 420.0 / 200000.0
        expected = 0.75 + 0.15 * 0.0015 / 0.003  # = 0.825
        assert math.isclose(
            sr.phi_flexure(eps_ty + 0.0015, 420.0, transverse='spiral'),
            expected,
            rel_tol=1e-9,
        )

    def test_gr80_tension_controlled(self):
        """Grade 80 (fy=550) with eps_t=0.010 is tension-controlled → 0.90."""
        assert math.isclose(sr.phi_flexure(0.010, 550.0), 0.90, rel_tol=1e-9)

    def test_gr80_compression_controlled(self):
        """Grade 80 (fy=550) at eps_ty with other → compression-controlled."""
        eps_ty = 550.0 / 200000.0
        assert math.isclose(
            sr.phi_flexure(eps_ty, 550.0, transverse='other'),
            0.65,
            rel_tol=1e-9,
        )

    def test_invalid_fy(self):
        """Fy <= 0 should raise ValueError."""
        with pytest.raises(ValueError):
            sr.phi_flexure(0.005, 0.0)

    def test_invalid_transverse(self):
        """Unknown transverse type should raise ValueError."""
        with pytest.raises(ValueError):
            sr.phi_flexure(0.005, 420.0, transverse='ties')


class TestSectionClassification:
    """Tests for section_classification (Table 21.2.2)."""

    def test_tension_controlled(self):
        """eps_t well above eps_ty + 0.003 → tension-controlled."""
        eps_ty = 420.0 / 200000.0
        assert (
            sr.section_classification(eps_ty + 0.005, 420.0)
            == 'tension-controlled'
        )

    def test_tension_controlled_at_limit(self):
        """eps_t = eps_ty + 0.003 → tension-controlled (boundary inclusive)."""
        eps_ty = 420.0 / 200000.0
        assert (
            sr.section_classification(eps_ty + 0.003, 420.0)
            == 'tension-controlled'
        )

    def test_compression_controlled(self):
        """eps_t = eps_ty → compression-controlled (boundary inclusive)."""
        eps_ty = 420.0 / 200000.0
        assert (
            sr.section_classification(eps_ty, 420.0)
            == 'compression-controlled'
        )

    def test_compression_controlled_below(self):
        """eps_t < eps_ty → compression-controlled."""
        eps_ty = 420.0 / 200000.0
        assert (
            sr.section_classification(eps_ty - 0.001, 420.0)
            == 'compression-controlled'
        )

    def test_transition(self):
        """eps_t between eps_ty and eps_ty + 0.003 → transition."""
        eps_ty = 420.0 / 200000.0
        assert (
            sr.section_classification(eps_ty + 0.0015, 420.0) == 'transition'
        )

    def test_invalid_fy(self):
        """Fy <= 0 should raise ValueError."""
        with pytest.raises(ValueError):
            sr.section_classification(0.005, 0.0)
