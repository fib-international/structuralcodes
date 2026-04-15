"""Tests for ACI 318-25 flexural strength functions (Ch. 22.2-22.3)."""

import math

import pytest

from structuralcodes.codes.aci318_25 import _flexure as fl

# ---------------------------------------------------------------------------
# Shared test constants
# ---------------------------------------------------------------------------
FC = 27.58      # MPa  (4000 psi)
FY = 420.0      # MPa  (~60 ksi)
B = 305.0       # mm   (12 in)
D = 227.0       # mm   (8.94 in)
H = 254.0       # mm   (10 in)
BETA1 = 0.85    # stress-block factor for fc = 4000 psi


# ---------------------------------------------------------------------------
# Equilibrium helpers
# ---------------------------------------------------------------------------


class TestStressBlockDepthSR:
    """Tests for stress_block_depth_sr."""

    def test_known_value(self):
        """Verify a = As*fy / (0.85*fc*b) for As=645."""
        As = 645.0
        expected = As * FY / (0.85 * FC * B)
        result = fl.stress_block_depth_sr(As, FY, FC, B)
        assert math.isclose(result, expected, rel_tol=1e-9)

    def test_positive(self):
        """Stress-block depth must be positive."""
        assert fl.stress_block_depth_sr(645.0, FY, FC, B) > 0.0


class TestStressBlockDepthDR:
    """Tests for stress_block_depth_dr."""

    def test_known_value(self):
        """Verify a = (As*fy - As'*fy') / (0.85*fc*b) for As=800, As'=200."""
        As, As_prime = 800.0, 200.0
        expected = (As * FY - As_prime * FY) / (0.85 * FC * B)
        result = fl.stress_block_depth_dr(As, As_prime, FY, FY, FC, B)
        assert math.isclose(result, expected, rel_tol=1e-9)

    def test_equals_sr_when_no_compression_steel(self):
        """With As'=0 the doubly-reinforced result equals the singly-reinforced."""
        As = 645.0
        dr = fl.stress_block_depth_dr(As, 0.0, FY, FY, FC, B)
        sr = fl.stress_block_depth_sr(As, FY, FC, B)
        assert math.isclose(dr, sr, rel_tol=1e-9)


class TestNeutralAxisDepth:
    """Tests for neutral_axis_depth."""

    def test_known_value(self):
        """c = a / beta1 for a=37.9, beta1=0.85."""
        a = 37.9
        beta1 = 0.85
        expected = a / beta1
        assert math.isclose(fl.neutral_axis_depth(a, beta1), expected, rel_tol=1e-9)

    def test_value_exceeds_a(self):
        """Neutral-axis depth must be >= a (since beta1 <= 1)."""
        a = 37.9
        assert fl.neutral_axis_depth(a, BETA1) >= a


class TestEpsTFromC:
    """Tests for eps_t_from_c."""

    def test_known_value(self):
        """eps_t = 0.003*(dt - c)/c for c=44.6, dt=227."""
        c, dt = 44.6, 227.0
        expected = 0.003 * (dt - c) / c
        assert math.isclose(fl.eps_t_from_c(c, dt), expected, rel_tol=1e-9)

    def test_zero_when_c_equals_dt(self):
        """When c == dt the tensile strain is zero (neutral axis at steel)."""
        assert math.isclose(fl.eps_t_from_c(D, D), 0.0, abs_tol=1e-12)


class TestEpsSPrime:
    """Tests for eps_s_prime."""

    def test_known_value(self):
        """eps_s' = 0.003*(c - d')/c for c=80, d'=40 -> 0.0015."""
        c, d_prime = 80.0, 40.0
        expected = 0.003 * (c - d_prime) / c  # = 0.0015
        assert math.isclose(fl.eps_s_prime(c, d_prime), expected, rel_tol=1e-9)
        assert math.isclose(fl.eps_s_prime(c, d_prime), 0.0015, rel_tol=1e-9)

    def test_zero_when_d_prime_equals_c(self):
        """When d' == c the compression steel sits at the neutral axis → 0."""
        assert math.isclose(fl.eps_s_prime(80.0, 80.0), 0.0, abs_tol=1e-12)


# ---------------------------------------------------------------------------
# Nominal moment strength
# ---------------------------------------------------------------------------


class TestMnSinglyReinforced:
    """Tests for Mn_singly_reinforced."""

    def test_formula(self):
        """Verify Mn = As*fy*(d - a/2) for As=645."""
        As = 645.0
        a = fl.stress_block_depth_sr(As, FY, FC, B)
        expected = As * FY * (D - a / 2.0)
        result = fl.Mn_singly_reinforced(As, FY, FC, B, D)
        assert math.isclose(result, expected, rel_tol=1e-9)

    def test_positive(self):
        """Nominal moment must be positive."""
        assert fl.Mn_singly_reinforced(645.0, FY, FC, B, D) > 0.0


class TestMnDoublyReinforced:
    """Tests for Mn_doubly_reinforced."""

    def test_formula(self):
        """Verify formula for As=800, As'=200, d'=40."""
        As, As_prime, d_prime = 800.0, 200.0, 40.0
        a = fl.stress_block_depth_dr(As, As_prime, FY, FY, FC, B)
        expected = (As * FY - As_prime * FY) * (D - a / 2.0) + As_prime * FY * (D - d_prime)
        result = fl.Mn_doubly_reinforced(As, As_prime, FY, FY, FC, B, D, d_prime)
        assert math.isclose(result, expected, rel_tol=1e-9)

    def test_exceeds_singly_reinforced(self):
        """Adding compression steel increases Mn compared to tension steel only."""
        Mn_sr = fl.Mn_singly_reinforced(800.0, FY, FC, B, D)
        Mn_dr = fl.Mn_doubly_reinforced(800.0, 200.0, FY, FY, FC, B, D, 40.0)
        assert Mn_dr > Mn_sr


# ---------------------------------------------------------------------------
# Reinforcement limits
# ---------------------------------------------------------------------------


class TestAsMinSlab:
    """Tests for As_min_slab (Sec. 7.6.1.1 -> 24.4.3.2)."""

    def test_grade60_ratio_0018(self):
        """Grade 60 (fy=420 MPa ~ 60 900 psi) -> ratio = 0.0018."""
        # 420 MPa * 145.038 psi/MPa = 60 916 psi > 60 000 -> 0.0018 * 60000/60916
        # but let's use exactly 413.7 MPa = 60 000 psi boundary for clarity;
        # 420 MPa sits > 60 000 psi, so use 413 MPa for grade-60 test.
        fy_60ksi = 413.685  # MPa corresponding to exactly 60 000 psi
        result = fl.As_min_slab(fy_60ksi, B, H)
        expected = 0.0018 * B * H
        assert math.isclose(result, expected, rel_tol=1e-6)

    def test_grade40_ratio_0020(self):
        """Grade 40 (fy ~ 276 MPa, 40 000 psi) -> ratio = 0.0020."""
        fy_40ksi = 275.79  # MPa  (40 000 psi)
        result = fl.As_min_slab(fy_40ksi, B, H)
        expected = 0.0020 * B * H
        assert math.isclose(result, expected, rel_tol=1e-6)

    def test_grade80_lower_ratio(self):
        """Grade 80 (fy=552 MPa, ~80 000 psi) -> ratio = max(0.0014, 0.0018*60000/fy_psi)."""
        fy_80ksi = 551.58  # MPa (80 000 psi)
        fy_psi = fy_80ksi * 145.038
        expected_ratio = max(0.0014, 0.0018 * 60000.0 / fy_psi)
        result = fl.As_min_slab(fy_80ksi, B, H)
        assert math.isclose(result, expected_ratio * B * H, rel_tol=1e-6)

    def test_grade80_ratio_below_0018(self):
        """Grade 80 ratio must be less than 0.0018."""
        fy_80ksi = 551.58
        ratio = fl.As_min_slab(fy_80ksi, B, H) / (B * H)
        assert ratio < 0.0018

    def test_minimum_floor(self):
        """Very high fy should not produce ratio below 0.0014."""
        fy_very_high = 700.0  # MPa
        ratio = fl.As_min_slab(fy_very_high, B, H) / (B * H)
        assert ratio >= 0.0014


class TestAsMinBeam:
    """Tests for As_min_beam (Sec. 9.6.1.2)."""

    def test_formula(self):
        """Verify As_min = max(0.25*sqrt(fc)/fy, 1.4/fy) * bw * d."""
        expected_ratio = max(0.25 * math.sqrt(FC) / FY, 1.4 / FY)
        expected = expected_ratio * B * D
        result = fl.As_min_beam(FC, FY, B, D)
        assert math.isclose(result, expected, rel_tol=1e-9)

    def test_positive(self):
        """Minimum beam steel area must be positive."""
        assert fl.As_min_beam(FC, FY, B, D) > 0.0

    def test_low_fc_governed_by_14_over_fy(self):
        """For very low fc, the 1.4/fy term should govern."""
        fc_low = 10.0  # MPa — gives 0.25*sqrt(10)/420 < 1.4/420
        ratio_concrete = 0.25 * math.sqrt(fc_low) / FY
        ratio_min = 1.4 / FY
        assert ratio_concrete < ratio_min  # confirms 1.4/fy governs
        result = fl.As_min_beam(fc_low, FY, B, D)
        assert math.isclose(result, ratio_min * B * D, rel_tol=1e-9)


class TestAsMaxCheck:
    """Tests for As_max_check."""

    def test_tension_controlled_returns_true(self):
        """eps_t = 0.010 is well into the tension-controlled zone."""
        assert fl.As_max_check(0.010, FY) is True

    def test_compression_controlled_returns_false(self):
        """eps_t = 0.002 is below the tension-controlled limit for Grade 60."""
        assert fl.As_max_check(0.002, FY) is False

    def test_at_boundary_returns_true(self):
        """eps_t exactly at fy/Es + 0.003 must return True."""
        Es = 200000.0
        eps_boundary = FY / Es + 0.003
        assert fl.As_max_check(eps_boundary, FY) is True

    def test_just_below_boundary_returns_false(self):
        """eps_t just below the tension-controlled limit must return False."""
        Es = 200000.0
        eps_just_below = FY / Es + 0.003 - 1e-10
        assert fl.As_max_check(eps_just_below, FY) is False


# ---------------------------------------------------------------------------
# Design helper
# ---------------------------------------------------------------------------


class TestAsRequired:
    """Tests for As_required."""

    def test_round_trip(self):
        """Compute As for Mu=50e6, then verify phi*Mn >= Mu."""
        Mu = 50.0e6  # N·mm
        phi = 0.9
        As = fl.As_required(Mu, phi, FY, FC, B, D)
        Mn = fl.Mn_singly_reinforced(As, FY, FC, B, D)
        assert phi * Mn >= Mu - 1.0  # allow 1 N·mm floating-point tolerance

    def test_positive_result(self):
        """Required steel area must be positive for a reasonable moment."""
        As = fl.As_required(50.0e6, 0.9, FY, FC, B, D)
        assert As > 0.0

    def test_negative_discriminant_raises(self):
        """An extremely large Mu relative to section capacity must raise ValueError."""
        with pytest.raises(ValueError, match='discriminant'):
            fl.As_required(1.0e15, 0.9, FY, FC, B, D)

    def test_phi_unity_gives_minimum_as(self):
        """phi=1.0 should produce a smaller As than phi=0.9 for the same Mu."""
        As_phi09 = fl.As_required(50.0e6, 0.9, FY, FC, B, D)
        As_phi10 = fl.As_required(50.0e6, 1.0, FY, FC, B, D)
        assert As_phi10 < As_phi09
