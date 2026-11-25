"""Tests for ACI 318-25 two-way shear with moment transfer (Sec. 8.4.4.2).

Expected values are hand calculations in inch-pound units converted to SI.
The interior case is a 24 in square column with d = 20 in: the critical
section is a 44 in square and bo d = 176 x 20 = 3520 in2. The edge case is
the same column flush with a slab edge at its +x face: two legs of 34 in at
y = +/-22 in running from x = -22 in to x = 12 in, and a side of 44 in at
x = -22 in.
"""

import pytest

from structuralcodes.codes.aci318_25 import _moment_transfer as mt
from structuralcodes.codes.aci318_25._critical_section import (
    critical_section_rectangular,
)

IN = 25.4  # mm per in
PSI = 0.006894757293168361  # MPa per psi
LBF = 4.4482216152605  # N per lbf

C = 24.0 * IN
D = 20.0 * IN


@pytest.fixture
def interior():
    """Critical section of the interior column."""
    return critical_section_rectangular(C, C, D)


@pytest.fixture
def edge():
    """Critical section of the edge column, free edge at the +x face."""
    return critical_section_rectangular(C, C, D, edge_distances={'+x': 0.0})


class TestGamma:
    """Eq. (8.4.2.2.2) and Eq. (8.4.4.2.2)."""

    def test_gamma_f_square(self):
        """With b1 = b2: 1 / (1 + 2/3) = 0.600."""
        assert mt.gamma_f(44.0, 44.0) == pytest.approx(0.600, rel=1e-9)

    def test_gamma_f_rectangular(self):
        """With b1 = 32, b2 = 68: 1 / (1 + (2/3) sqrt(32/68)) = 0.6862."""
        assert mt.gamma_f(32.0, 68.0) == pytest.approx(0.6862, rel=1e-4)

    def test_gamma_v(self):
        """With b1 = b2: 1 - 0.600 = 0.400."""
        assert mt.gamma_v(44.0, 44.0) == pytest.approx(0.400, rel=1e-9)

    def test_non_positive_dimension_raises(self):
        """Both dimensions of the section must be positive."""
        with pytest.raises(ValueError):
            mt.gamma_f(0.0, 44.0)


class TestCentroid:
    """Centroid of the critical section."""

    def test_interior(self, interior):
        """A closed symmetric section is centred on the column."""
        x, y = mt.critical_section_centroid(interior)
        assert x == pytest.approx(0.0, abs=1e-9)
        assert y == pytest.approx(0.0, abs=1e-9)

    def test_edge(self, edge):
        """Legs centred at x = -5 in, side at x = -22 in.

        x = (2 x 34 x (-5) + 44 x (-22)) / 112 = -11.679 in.
        """
        x, y = mt.critical_section_centroid(edge)
        assert x / IN == pytest.approx(-11.679, rel=1e-4)
        assert y == pytest.approx(0.0, abs=1e-9)


class TestJc:
    """Property Jc of the critical section, R8.4.4.2.3."""

    def test_interior_square(self, interior):
        """Commentary expression with c1 + d = c2 + d = 44 in, d = 20 in.

        20 x 44^3 / 6 + 44 x 20^3 / 6 + 20 x 44 x 44^2 / 2
        = 283,947 + 58,667 + 851,840 = 1,194,453 in4.
        """
        Jc = mt.Jc(interior, D, 'x')
        assert Jc / IN**4 == pytest.approx(1194453.0, rel=1e-6)

    def test_interior_rectangle(self):
        """12 in x 48 in column, arms along x: b1 = 32 in, b2 = 68 in.

        20 x 32^3 / 6 + 32 x 20^3 / 6 + 20 x 68 x 32^2 / 2 = 848,213 in4.
        """
        section = critical_section_rectangular(12.0 * IN, 48.0 * IN, D)
        Jc = mt.Jc(section, D, 'x')
        assert Jc / IN**4 == pytest.approx(848213.0, rel=1e-6)

    def test_edge_arms_perpendicular_to_edge(self, edge):
        """Parallel axis theorem about the centroid at x = -11.679 in.

        Each leg: 20 x 34^3 / 12 + 34 x 20^3 / 12 + 680 x 6.679^2.
        Side: 880 x 10.321^2. Total 330,755 in4.
        """
        Jc = mt.Jc(edge, D, 'x')
        assert Jc / IN**4 == pytest.approx(330755.0, rel=1e-5)

    def test_edge_arms_parallel_to_edge(self, edge):
        """Arms along y, centroid at y = 0.

        Legs: 2 x 680 x 22^2. Side: 20 x 44^3 / 12 + 44 x 20^3 / 12.
        Total 829,547 in4.
        """
        Jc = mt.Jc(edge, D, 'y')
        assert Jc / IN**4 == pytest.approx(829547.0, rel=1e-5)

    def test_unknown_direction_raises(self, interior):
        """Directions are 'x' and 'y'."""
        with pytest.raises(ValueError):
            mt.Jc(interior, D, 'z')


class TestVuMax:
    """Maximum factored shear stress, 8.4.4.2.1 to 8.4.4.2.3."""

    def test_shear_only(self, interior):
        """Vu = 300 kip: 300,000 / 3520 = 85.23 psi."""
        vu = mt.vu_max(300e3 * LBF, interior, D)
        assert vu / PSI == pytest.approx(85.23, rel=1e-4)

    def test_interior_with_moment(self, interior):
        """Vu = 300 kip, Msc = 200 kip-ft, gamma_v = 0.4.

        0.4 x 2,400,000 x 22 / 1,194,453 = 17.68 psi.
        vu = 85.23 + 17.68 = 102.91 psi.
        """
        moment = 0.4 * 200e3 * 12.0 * LBF * IN
        vu = mt.vu_max(300e3 * LBF, interior, D, gamma_v_Msc_x=moment)
        assert vu / PSI == pytest.approx(102.91, rel=1e-4)

    def test_interior_sign_of_moment_does_not_matter(self, interior):
        """A symmetric section gives the same maximum for either sign."""
        moment = 0.4 * 200e3 * 12.0 * LBF * IN
        vu = mt.vu_max(300e3 * LBF, interior, D, gamma_v_Msc_x=-moment)
        assert vu / PSI == pytest.approx(102.91, rel=1e-4)

    def test_interior_two_directions(self, interior):
        """The same moment in both directions: 85.23 + 2 x 17.68 = 120.59."""
        moment = 0.4 * 200e3 * 12.0 * LBF * IN
        vu = mt.vu_max(
            300e3 * LBF,
            interior,
            D,
            gamma_v_Msc_x=moment,
            gamma_v_Msc_y=moment,
        )
        assert vu / PSI == pytest.approx(120.59, rel=1e-4)

    def test_edge_moment_raising_stress_at_the_free_edge(self, edge):
        """Vu = 100 kip, gamma_v Msc = +50 kip-ft, arms along x.

        vuv = 100,000 / (112 x 20) = 44.64 psi. The ends of the legs at
        the free edge are 12 + 11.679 = 23.679 in from the centroid:
        600,000 x 23.679 / 330,755 = 42.95 psi. vu = 87.60 psi.
        """
        moment = 50e3 * 12.0 * LBF * IN
        vu = mt.vu_max(100e3 * LBF, edge, D, gamma_v_Msc_x=moment)
        assert vu / PSI == pytest.approx(87.60, rel=1e-4)

    def test_edge_moment_raising_stress_at_the_inner_side(self, edge):
        """As above with the moment reversed.

        The inner side is 22 - 11.679 = 10.321 in from the centroid:
        600,000 x 10.321 / 330,755 = 18.72 psi. vu = 63.37 psi.
        """
        moment = -50e3 * 12.0 * LBF * IN
        vu = mt.vu_max(100e3 * LBF, edge, D, gamma_v_Msc_x=moment)
        assert vu / PSI == pytest.approx(63.37, rel=1e-4)
