"""Tests for ACI 318-25 one-way shear strength functions (Ch. 22.5).

Every expected value is a hand calculation in inch-pound units, the units
in which the Code expressions are written, converted to SI with exact
factors. The reference section is bw = 12 in, d = 20 in, f'c = 4000 psi,
sqrt(f'c) = 63.246 psi, bw * d = 240 in2.
"""

import pytest

from structuralcodes.codes.aci318_25 import _shear as sh

# Exact conversion factors
IN = 25.4  # mm per in
PSI = 0.006894757293168361  # MPa per psi
KIP = 4448.2216152605  # N per kip

FC = 4000.0 * PSI  # MPa
BW = 12.0 * IN  # mm
D = 20.0 * IN  # mm
AG = 12.0 * 24.0 * IN**2  # mm2, 12 in x 24 in gross section


def kip(force_newton: float) -> float:
    """Convert N to kip."""
    return force_newton / KIP


class TestLambdaS:
    """Eq. (22.5.5.1.3) in the ACI 318-19 form: sqrt(2 / (1 + d/10)) <= 1."""

    def test_d_20_in(self):
        """Hand calculation: sqrt(2 / (1 + 20/10)) = 0.8165."""
        assert sh.lambda_s(D) == pytest.approx(0.8165, rel=1e-4)

    def test_d_40_in(self):
        """Hand calculation: sqrt(2 / (1 + 40/10)) = 0.6325."""
        assert sh.lambda_s(40.0 * IN) == pytest.approx(0.6325, rel=1e-4)

    def test_d_10_in_is_unity(self):
        """sqrt(2 / (1 + 10/10)) = 1.0."""
        assert sh.lambda_s(10.0 * IN) == pytest.approx(1.0, rel=1e-9)

    def test_shallow_member_capped_at_one(self):
        """With d = 5 in the expression gives 1.155, limited to 1.0."""
        assert sh.lambda_s(5.0 * IN) == 1.0

    def test_non_positive_depth_raises(self):
        """A depth of zero is not a section."""
        with pytest.raises(ValueError):
            sh.lambda_s(0.0)


class TestVcSimplified:
    """Table 22.5.5.1 Eq. (a): [2 lambda sqrt(f'c) + Nu/6Ag] bw d."""

    def test_no_axial(self):
        """2 x 63.246 x 240 = 30,358 lbf = 30.358 kip."""
        assert kip(sh.Vc_simplified(FC, BW, D)) == pytest.approx(
            30.358, rel=1e-4
        )

    def test_compression(self):
        """Nu = 100 kip on 288 in2: Nu/6Ag = 57.87 psi; Vc = 44.247 kip."""
        Vc = sh.Vc_simplified(FC, BW, D, Nu=100.0 * KIP, Ag=AG)
        assert kip(Vc) == pytest.approx(44.247, rel=1e-4)

    def test_tension_not_below_zero(self):
        """Nu = -300 kip: 126.5 - 173.6 psi is negative, so Vc = 0."""
        Vc = sh.Vc_simplified(FC, BW, D, Nu=-300.0 * KIP, Ag=AG)
        assert Vc == 0.0

    def test_sqrt_fc_limited_to_100_psi(self):
        """22.5.3.1: f'c = 12,000 psi; 2 x 100 x 240 = 48.0 kip."""
        Vc = sh.Vc_simplified(12000.0 * PSI, BW, D)
        assert kip(Vc) == pytest.approx(48.0, rel=1e-4)

    def test_lightweight(self):
        """Lambda = 0.75: 0.75 x 30.358 = 22.768 kip."""
        Vc = sh.Vc_simplified(FC, BW, D, lambda_concrete=0.75)
        assert kip(Vc) == pytest.approx(22.768, rel=1e-4)

    def test_axial_without_area_raises(self):
        """Nu without Ag cannot be evaluated."""
        with pytest.raises(ValueError):
            sh.Vc_simplified(FC, BW, D, Nu=100.0 * KIP)


class TestVcDetailed:
    """Table 22.5.5.1 Eq. (b) and (c), with 22.5.5.1.1 and 22.5.5.1.2."""

    def test_eq_b_with_minimum_shear_reinforcement(self):
        """8 x 0.01^(1/3) x 63.246 x 240 = 26.162 kip."""
        Vc = sh.Vc_detailed(
            FC, BW, D, rho_w=0.01, has_min_shear_reinforcement=True
        )
        assert kip(Vc) == pytest.approx(26.162, rel=1e-4)

    def test_eq_c_is_the_default(self):
        """Lambda_s = 0.8165: 0.8165 x 26.162 = 21.361 kip."""
        Vc = sh.Vc_detailed(FC, BW, D, rho_w=0.01)
        assert kip(Vc) == pytest.approx(21.361, rel=1e-4)

    def test_lower_bound_governs_lightly_reinforced(self):
        """rho_w = 0.002: Eq. (c) gives 12.492 kip.

        22.5.5.1.1 lower bound: 63.246 x 240 = 15.179 kip.
        """
        Vc = sh.Vc_detailed(FC, BW, D, rho_w=0.002)
        assert kip(Vc) == pytest.approx(15.179, rel=1e-4)

    def test_lower_bound_switched_off(self):
        """22.5.5.1.1(b): the bound does not apply, Eq. (c) stands."""
        Vc = sh.Vc_detailed(FC, BW, D, rho_w=0.002, apply_lower_bound=False)
        assert kip(Vc) == pytest.approx(12.492, rel=1e-4)

    def test_lower_bound_not_applied_in_net_tension(self):
        """22.5.5.1.1(a): Nu = -10 kip, Nu/6Ag = -5.79 psi.

        (8 x 0.8165 x 0.002^(1/3) x 63.246 - 5.787) x 240 = 11.103 kip.
        """
        Vc = sh.Vc_detailed(FC, BW, D, rho_w=0.002, Nu=-10.0 * KIP, Ag=AG)
        assert kip(Vc) == pytest.approx(11.103, rel=1e-4)

    def test_upper_limit(self):
        """rho_w = 0.02 and Nu = 400 kip, with minimum shear reinforcement.

        Nu/6Ag = 231.5 psi, limited to 0.05 f'c = 200 psi (22.5.5.1.2).
        Eq. (b) gives 80.96 kip, limited to 5 x 63.246 x 240 = 75.895 kip
        (22.5.5.1.1).
        """
        Vc = sh.Vc_detailed(
            FC,
            BW,
            D,
            rho_w=0.02,
            Nu=400.0 * KIP,
            Ag=AG,
            has_min_shear_reinforcement=True,
        )
        assert kip(Vc) == pytest.approx(75.895, rel=1e-4)

    def test_large_tension_gives_zero(self):
        """Table 22.5.5.1 note 2: Vc is not taken below zero."""
        Vc = sh.Vc_detailed(FC, BW, D, rho_w=0.01, Nu=-500.0 * KIP, Ag=AG)
        assert Vc == 0.0

    def test_invalid_reinforcement_ratio_raises(self):
        """A reinforcement ratio must be positive."""
        with pytest.raises(ValueError):
            sh.Vc_detailed(FC, BW, D, rho_w=0.0)


class TestVs:
    """Eq. (22.5.8.5.3): Vs = Av fyt d / s."""

    def test_two_leg_no_3_at_8_in(self):
        """0.22 in2 x 60 ksi x 20 in / 8 in = 33.0 kip."""
        Av = 0.22 * IN**2
        fyt = 60000.0 * PSI
        Vs = sh.Vs(Av, fyt, D, 8.0 * IN)
        assert kip(Vs) == pytest.approx(33.0, rel=1e-6)


class TestVn:
    """Vn = Vc + Vs."""

    def test_sum(self):
        """Vn is the sum of the two parts."""
        assert sh.Vn(100e3, 50e3) == 150e3


class TestCheckCrossSection:
    """Eq. (22.5.1.2): Vu <= phi (Vc + 8 sqrt(f'c) bw d)."""

    def test_limit(self):
        """Vc = 30.358 kip; 0.75 x (30.358 + 121.431) = 113.842 kip."""
        Vc = 30.358 * KIP
        assert sh.check_cross_section(113.0 * KIP, 0.75, Vc, FC, BW, D)
        assert not sh.check_cross_section(114.5 * KIP, 0.75, Vc, FC, BW, D)


class TestAvMinPerS:
    """Table 9.6.3.4 (a) and (b): the greater of 0.75 sqrt(f'c) and 50."""

    def test_50_psi_governs(self):
        """With f'c = 4000 psi, 50 x 12 / 60,000 = 0.0100 in2/in."""
        result = sh.Av_min_per_s(FC, BW, 60000.0 * PSI)
        assert result / IN == pytest.approx(0.0100, rel=1e-6)

    def test_sqrt_fc_governs(self):
        """With f'c = 6000 psi: 58.09 x 12 / 60,000 = 0.011619 in2/in."""
        result = sh.Av_min_per_s(6000.0 * PSI, BW, 60000.0 * PSI)
        assert result / IN == pytest.approx(0.011619, rel=1e-4)


class TestShearReinforcementRequired:
    """Comparison of Vu with phi Vc."""

    def test_required(self):
        """Vu above phi Vc."""
        assert sh.shear_reinforcement_required(100e3, 90e3)

    def test_not_required(self):
        """Vu at or below phi Vc."""
        assert not sh.shear_reinforcement_required(90e3, 90e3)


class TestMaxStirrupSpacing:
    """Table 9.7.6.2.2, nonprestressed beam, along the length.

    The threshold is 4 x 63.246 x 240 = 60.716 kip.
    """

    def test_low_shear_d_over_2(self):
        """Vs = 50 kip, d = 20 in: lesser of d/2 = 10 in and 24 in."""
        s = sh.max_stirrup_spacing(D, 50.0 * KIP, FC, BW)
        assert s / IN == pytest.approx(10.0, rel=1e-9)

    def test_high_shear_d_over_4(self):
        """Vs = 70 kip, d = 20 in: lesser of d/4 = 5 in and 12 in."""
        s = sh.max_stirrup_spacing(D, 70.0 * KIP, FC, BW)
        assert s / IN == pytest.approx(5.0, rel=1e-9)

    def test_deep_member_24_in(self):
        """D = 60 in, low shear: d/2 = 30 in, limited to 24 in."""
        s = sh.max_stirrup_spacing(60.0 * IN, 50.0 * KIP, FC, BW)
        assert s / IN == pytest.approx(24.0, rel=1e-9)

    def test_deep_member_12_in(self):
        """D = 60 in, threshold 182.1 kip, Vs = 200 kip: 12 in."""
        s = sh.max_stirrup_spacing(60.0 * IN, 200.0 * KIP, FC, BW)
        assert s / IN == pytest.approx(12.0, rel=1e-9)
