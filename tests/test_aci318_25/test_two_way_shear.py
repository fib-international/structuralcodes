"""Tests for ACI 318-25 two-way shear strength functions (Sec. 22.6.5).

Every expected value is a hand calculation in inch-pound units, the units
in which the Code expressions are written, converted to SI with exact
factors. For f'c = 4000 psi, sqrt(f'c) = 63.246 psi.
"""

import pytest

from structuralcodes.codes.aci318_25 import _two_way_shear as tws

# Exact conversion factors
IN = 25.4  # mm per in
PSI = 0.006894757293168361  # MPa per psi

FC = 4000.0 * PSI  # MPa


def psi(stress_mpa: float) -> float:
    """Convert MPa to psi."""
    return stress_mpa / PSI


class TestAlphaS:
    """Sec. 22.6.5.3: 40 interior, 30 edge, 20 corner."""

    @pytest.mark.parametrize(
        'location, expected',
        [('interior', 40.0), ('edge', 30.0), ('corner', 20.0)],
    )
    def test_values(self, location, expected):
        """The three values of the provision."""
        assert tws.alpha_s(location) == expected

    def test_unknown_location_raises(self):
        """Any other name is rejected."""
        with pytest.raises(ValueError):
            tws.alpha_s('re-entrant')


class TestVcTwoWay:
    """Table 22.6.5.2: the least of expressions (a), (b) and (c)."""

    def test_expression_a_governs(self):
        """Interior 24 in square column, d = 20 in, bo = 4 x 44 = 176 in.

        (a) 4; (b) 2 + 4/1 = 6; (c) 2 + 40 x 20 / 176 = 6.55.
        lambda_s = sqrt(2 / (1 + 20/10)) = 0.8165.
        vc = 4 x 0.8165 x 63.246 = 206.56 psi.
        """
        vc = tws.vc_two_way(
            FC, d=20.0 * IN, bo=176.0 * IN, beta=1.0, alpha_s=40.0
        )
        assert psi(vc) == pytest.approx(206.56, rel=1e-4)

    def test_expression_b_governs(self):
        """Interior 12 in x 48 in column, d = 20 in, bo = 200 in.

        (a) 4; (b) 2 + 4/4 = 3; (c) 2 + 40 x 20 / 200 = 6.
        vc = 3 x 0.8165 x 63.246 = 154.92 psi.
        """
        vc = tws.vc_two_way(
            FC, d=20.0 * IN, bo=200.0 * IN, beta=4.0, alpha_s=40.0
        )
        assert psi(vc) == pytest.approx(154.92, rel=1e-4)

    def test_expression_c_governs(self):
        """Corner section, d = 8 in, bo = 100 in, alpha_s = 20.

        (a) 4; (b) 6; (c) 2 + 20 x 8 / 100 = 3.6. With d less than 10 in,
        lambda_s = 1. vc = 3.6 x 63.246 = 227.68 psi.
        """
        vc = tws.vc_two_way(
            FC, d=8.0 * IN, bo=100.0 * IN, beta=1.0, alpha_s=20.0
        )
        assert psi(vc) == pytest.approx(227.68, rel=1e-4)

    def test_lightweight(self):
        """As above with lambda = 0.75: 0.75 x 227.68 = 170.76 psi."""
        vc = tws.vc_two_way(
            FC,
            d=8.0 * IN,
            bo=100.0 * IN,
            beta=1.0,
            alpha_s=20.0,
            lambda_concrete=0.75,
        )
        assert psi(vc) == pytest.approx(170.76, rel=1e-4)

    def test_sqrt_fc_limited_to_100_psi(self):
        """Sec. 22.6.3.1 with f'c = 12,000 psi and d = 8 in.

        Interior square column, bo = 80 in: (c) 2 + 40 x 8 / 80 = 6, so (a)
        governs. vc = 4 x 100 = 400 psi.
        """
        vc = tws.vc_two_way(
            12000.0 * PSI, d=8.0 * IN, bo=80.0 * IN, beta=1.0, alpha_s=40.0
        )
        assert psi(vc) == pytest.approx(400.0, rel=1e-9)

    def test_beta_below_one_raises(self):
        """Beta is the long side over the short side."""
        with pytest.raises(ValueError):
            tws.vc_two_way(
                FC, d=8.0 * IN, bo=80.0 * IN, beta=0.5, alpha_s=40.0
            )

    def test_zero_perimeter_raises(self):
        """A section with no perimeter has no strength to report."""
        with pytest.raises(ValueError):
            tws.vc_two_way(FC, d=8.0 * IN, bo=0.0, beta=1.0, alpha_s=40.0)
