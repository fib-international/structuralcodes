"""Tests for ACI 318-25 one-way shear strength functions (Ch. 22.5)."""

import math

import pytest

from structuralcodes.codes.aci318_25 import _shear as sh

# ---------------------------------------------------------------------------
# Shared test constants
# ---------------------------------------------------------------------------
FC = 27.58      # MPa  (4 000 psi)
BW = 305.0      # mm   (12 in)
D = 227.0       # mm   (~8.94 in)
RHO_W = 0.009   # longitudinal reinforcement ratio


# ---------------------------------------------------------------------------
# lambda_s — size-effect factor
# ---------------------------------------------------------------------------


class TestLambdaS:
    """Tests for lambda_s (Eq. 22.5.5.1.3)."""

    def test_shallow_depth_capped_at_one(self):
        """d=227 mm is shallow enough that lambda_s should be capped at 1.0."""
        result = sh.lambda_s(D)
        # d_in = 227/25.4 = 8.937 in  -> 2/(1+8.937/10) = 2/1.8937 = 1.056 -> capped at 1.0
        assert math.isclose(result, 1.0, rel_tol=1e-9)

    def test_deep_member(self):
        """d=900 mm should give a size-effect factor less than 1.0."""
        d_deep = 900.0
        d_in = d_deep / 25.4
        expected = min(2.0 / (1.0 + d_in / 10.0), 1.0)
        result = sh.lambda_s(d_deep)
        assert math.isclose(result, expected, rel_tol=1e-9)
        assert result < 1.0

    def test_never_exceeds_one(self):
        """lambda_s must never exceed 1.0 for any positive depth."""
        for d_test in [50.0, 100.0, 200.0, 500.0, 1000.0, 2000.0]:
            assert sh.lambda_s(d_test) <= 1.0


# ---------------------------------------------------------------------------
# Vc_detailed — detailed shear strength
# ---------------------------------------------------------------------------


class TestVcDetailed:
    """Tests for Vc_detailed (Table 22.5.5.1)."""

    def test_without_min_reinforcement_uses_size_effect(self):
        """Without minimum reinforcement, lambda_s should reduce Vc."""
        sqrt_fc = min(math.sqrt(FC), 8.3)
        ls = sh.lambda_s(D)
        expected_base = 8.0 * ls * 1.0 * RHO_W ** (1.0 / 3.0) * sqrt_fc
        expected_Vc_raw = expected_base * BW * D
        # floor: lambda * sqrt_fc * bw * d
        Vc_floor = 1.0 * sqrt_fc * BW * D
        expected_Vc = max(expected_Vc_raw, Vc_floor)
        result = sh.Vc_detailed(FC, BW, D, RHO_W, Av_provided=0, Av_min=100)
        assert math.isclose(result, expected_Vc, rel_tol=1e-9)

    def test_with_min_reinforcement_no_size_effect(self):
        """With minimum reinforcement provided, lambda_s is not applied."""
        sqrt_fc = min(math.sqrt(FC), 8.3)
        expected_base = 8.0 * 1.0 * RHO_W ** (1.0 / 3.0) * sqrt_fc
        expected_Vc_raw = expected_base * BW * D
        Vc_floor = 1.0 * sqrt_fc * BW * D
        expected_Vc = max(expected_Vc_raw, Vc_floor)
        result = sh.Vc_detailed(FC, BW, D, RHO_W, Av_provided=200, Av_min=100)
        assert math.isclose(result, expected_Vc, rel_tol=1e-9)

    def test_with_min_reinforcement_exceeds_without(self):
        """Providing minimum reinforcement gives Vc >= version without."""
        Vc_no_rein = sh.Vc_detailed(FC, BW, D, RHO_W, Av_provided=0, Av_min=100)
        Vc_with_rein = sh.Vc_detailed(FC, BW, D, RHO_W, Av_provided=200, Av_min=100)
        assert Vc_with_rein >= Vc_no_rein

    def test_vc_not_negative_with_large_tension(self):
        """Large tensile force (negative Nu) must not produce negative Vc."""
        Nu_tension = -1e8  # 100 MN tension — extreme case
        Ag = BW * 300.0
        result = sh.Vc_detailed(FC, BW, D, RHO_W, Nu=Nu_tension, Ag=Ag)
        assert result >= 0.0

    def test_upper_cap_applied(self):
        """Vc must not exceed 5*lambda*sqrt(fc)*bw*d."""
        sqrt_fc = min(math.sqrt(FC), 8.3)
        cap = 5.0 * 1.0 * sqrt_fc * BW * D
        result = sh.Vc_detailed(FC, BW, D, rho_w=1.0)  # extreme rho_w to try to exceed cap
        assert result <= cap + 1e-6

    def test_compressive_axial_increases_vc(self):
        """Compressive axial load (positive Nu) should increase Vc."""
        Vc_no_axial = sh.Vc_detailed(FC, BW, D, RHO_W)
        Nu_comp = 500e3  # 500 kN compression
        Ag = BW * 300.0
        Vc_with_axial = sh.Vc_detailed(FC, BW, D, RHO_W, Nu=Nu_comp, Ag=Ag)
        assert Vc_with_axial >= Vc_no_axial


# ---------------------------------------------------------------------------
# Vc_simplified — simplified shear strength
# ---------------------------------------------------------------------------


class TestVcSimplified:
    """Tests for Vc_simplified (Table 22.5.5.1(a))."""

    def test_no_axial_formula(self):
        """With no axial load: Vc = 2*lambda*sqrt(fc)*bw*d."""
        expected = 2.0 * 1.0 * math.sqrt(FC) * BW * D
        result = sh.Vc_simplified(FC, BW, D)
        assert math.isclose(result, expected, rel_tol=1e-9)

    def test_lightweight_concrete_reduces_vc(self):
        """A lightweight factor < 1.0 should reduce Vc."""
        Vc_nw = sh.Vc_simplified(FC, BW, D)
        Vc_lw = sh.Vc_simplified(FC, BW, D, lambda_concrete=0.75)
        assert Vc_lw < Vc_nw

    def test_compressive_axial_increases_vc(self):
        """Positive Nu should increase Vc via the axial term."""
        Vc_no_axial = sh.Vc_simplified(FC, BW, D)
        Ag = BW * 300.0
        Vc_with_axial = sh.Vc_simplified(FC, BW, D, Nu=500e3, Ag=Ag)
        assert Vc_with_axial > Vc_no_axial

    def test_positive(self):
        """Vc_simplified must be positive for ordinary inputs."""
        assert sh.Vc_simplified(FC, BW, D) > 0.0


# ---------------------------------------------------------------------------
# Vs — steel shear contribution
# ---------------------------------------------------------------------------


class TestVs:
    """Tests for Vs (Eq. 22.5.8.5.3)."""

    def test_formula(self):
        """Verify Vs = Av*fyt*d/s for Av=142, fyt=420, d=227, s=150."""
        Av, fyt, d, s = 142.0, 420.0, 227.0, 150.0
        expected = Av * fyt * d / s
        result = sh.Vs(Av, fyt, d, s)
        assert math.isclose(result, expected, rel_tol=1e-9)

    def test_positive(self):
        """Vs must be positive for positive inputs."""
        assert sh.Vs(142.0, 420.0, 227.0, 150.0) > 0.0

    def test_decreases_with_larger_spacing(self):
        """Doubling the stirrup spacing should halve Vs."""
        Vs_150 = sh.Vs(142.0, 420.0, D, 150.0)
        Vs_300 = sh.Vs(142.0, 420.0, D, 300.0)
        assert math.isclose(Vs_300, Vs_150 / 2.0, rel_tol=1e-9)


# ---------------------------------------------------------------------------
# Vn — total nominal shear strength
# ---------------------------------------------------------------------------


class TestVn:
    """Tests for Vn."""

    def test_simple_sum(self):
        """Vn = Vc + Vs."""
        Vc, Vs_val = 150000.0, 90000.0
        assert math.isclose(sh.Vn(Vc, Vs_val), 240000.0, rel_tol=1e-9)

    def test_zero_vs(self):
        """With no shear steel, Vn == Vc."""
        Vc = 120000.0
        assert math.isclose(sh.Vn(Vc, 0.0), Vc, rel_tol=1e-9)


# ---------------------------------------------------------------------------
# check_cross_section
# ---------------------------------------------------------------------------


class TestCheckCrossSection:
    """Tests for check_cross_section (Eq. 22.5.1.2)."""

    def test_passes_when_vu_below_limit(self):
        """A Vu well below the limit should return True."""
        phi = 0.75
        Vc = sh.Vc_simplified(FC, BW, D)
        limit = phi * (Vc + 8.0 * math.sqrt(FC) * BW * D)
        Vu = 0.8 * limit
        assert sh.check_cross_section(Vu, phi, Vc, FC, BW, D) is True

    def test_fails_when_vu_above_limit(self):
        """A Vu above the limit should return False."""
        phi = 0.75
        Vc = sh.Vc_simplified(FC, BW, D)
        limit = phi * (Vc + 8.0 * math.sqrt(FC) * BW * D)
        Vu = 1.2 * limit
        assert sh.check_cross_section(Vu, phi, Vc, FC, BW, D) is False

    def test_at_limit_returns_true(self):
        """Vu exactly at the limit should return True (<=)."""
        phi = 0.75
        Vc = sh.Vc_simplified(FC, BW, D)
        limit = phi * (Vc + 8.0 * math.sqrt(FC) * BW * D)
        assert sh.check_cross_section(limit, phi, Vc, FC, BW, D) is True


# ---------------------------------------------------------------------------
# Av_min_per_s
# ---------------------------------------------------------------------------


class TestAvMinPerS:
    """Tests for Av_min_per_s (Sec. 9.6.3.4)."""

    def test_formula(self):
        """Verify Av/s = max(0.062*sqrt(fc), 0.35) * bw / fyt."""
        fyt = 420.0
        expected = max(0.062 * math.sqrt(FC), 0.35) * BW / fyt
        result = sh.Av_min_per_s(FC, BW, fyt)
        assert math.isclose(result, expected, rel_tol=1e-9)

    def test_positive(self):
        """Minimum Av/s must be positive."""
        assert sh.Av_min_per_s(FC, BW, 420.0) > 0.0

    def test_lower_fyt_gives_more_steel(self):
        """Lower fyt should require more area per unit length."""
        Avs_420 = sh.Av_min_per_s(FC, BW, 420.0)
        Avs_280 = sh.Av_min_per_s(FC, BW, 280.0)
        assert Avs_280 > Avs_420


# ---------------------------------------------------------------------------
# shear_reinforcement_required
# ---------------------------------------------------------------------------


class TestShearReinforcementRequired:
    """Tests for shear_reinforcement_required."""

    def test_required_when_vu_exceeds_phi_vc(self):
        """Vu > phi*Vc → reinforcement required."""
        phi_Vc = 0.75 * sh.Vc_simplified(FC, BW, D)
        Vu = phi_Vc + 1000.0  # slightly above
        assert sh.shear_reinforcement_required(Vu, phi_Vc) is True

    def test_not_required_when_vu_below_phi_vc(self):
        """Vu < phi*Vc → no reinforcement required."""
        phi_Vc = 0.75 * sh.Vc_simplified(FC, BW, D)
        Vu = phi_Vc - 1000.0  # slightly below
        assert sh.shear_reinforcement_required(Vu, phi_Vc) is False

    def test_equal_boundary_not_required(self):
        """Vu == phi*Vc is not strictly greater, so not required."""
        phi_Vc = 50000.0
        assert sh.shear_reinforcement_required(phi_Vc, phi_Vc) is False


# ---------------------------------------------------------------------------
# max_stirrup_spacing
# ---------------------------------------------------------------------------


class TestMaxStirrupSpacing:
    """Tests for max_stirrup_spacing (Sec. 9.7.6.2.2)."""

    def test_low_vs_limit_d_over_2(self):
        """For a deep member with low Vs, spacing is limited to d/2."""
        # Use a large d so that d/2 > 600 is not applicable, and d/2 < 600.
        # D = 227 mm → d/2 = 113.5 < 600 so limit is d/2 = 113.5
        Vs_low = 1.0  # essentially zero
        result = sh.max_stirrup_spacing(D, Vs_low, FC, BW)
        assert math.isclose(result, D / 2.0, rel_tol=1e-9)

    def test_low_vs_limit_capped_600(self):
        """For a very deep member, the 600 mm cap governs over d/2."""
        d_large = 1500.0  # d/2 = 750 > 600 → should be capped at 600
        Vs_low = 1.0
        result = sh.max_stirrup_spacing(d_large, Vs_low, FC, BW)
        assert math.isclose(result, 600.0, rel_tol=1e-9)

    def test_high_vs_limit_d_over_4(self):
        """For high Vs, spacing is limited to d/4 (when d/4 < 300)."""
        # D = 227 mm → d/4 = 56.75 < 300 so limit is d/4
        threshold = 4.0 * math.sqrt(FC) * BW * D
        Vs_high = threshold + 1.0  # just above threshold
        result = sh.max_stirrup_spacing(D, Vs_high, FC, BW)
        assert math.isclose(result, D / 4.0, rel_tol=1e-9)

    def test_high_vs_limit_capped_300(self):
        """For a very deep member with high Vs, the 300 mm cap governs."""
        d_large = 1500.0  # d/4 = 375 > 300 → cap at 300
        threshold = 4.0 * math.sqrt(FC) * BW * d_large
        Vs_high = threshold + 1.0
        result = sh.max_stirrup_spacing(d_large, Vs_high, FC, BW)
        assert math.isclose(result, 300.0, rel_tol=1e-9)

    def test_at_threshold_uses_lower_limit(self):
        """Vs exactly at threshold uses the less-restrictive limit (d/2, 600)."""
        threshold = 4.0 * math.sqrt(FC) * BW * D
        result = sh.max_stirrup_spacing(D, threshold, FC, BW)
        assert math.isclose(result, D / 2.0, rel_tol=1e-9)
