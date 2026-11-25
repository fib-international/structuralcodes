"""Tests for ACI 318-19 critical sections for two-way shear (Sec. 22.6.4).

Expected perimeters are hand calculations. Dimensions are worked in inches
and converted to mm. The reference case is a 24 in square column with
d = 20 in, so the critical section is a 44 in square, 22 in from the
column centroid on each side, and bo = 4 x 44 = 176 in.
"""

import pytest
from shapely.geometry import Polygon

from structuralcodes.codes.aci318_19 import _critical_section as cs

IN = 25.4  # mm per in

C = 24.0 * IN
D = 20.0 * IN
H = 24.0 * IN


def rect(x0: float, x1: float, y0: float, y1: float) -> Polygon:
    """Rectangle with corners given in inches."""
    return Polygon(
        [
            (x0 * IN, y0 * IN),
            (x1 * IN, y0 * IN),
            (x1 * IN, y1 * IN),
            (x0 * IN, y1 * IN),
        ]
    )


def bo_in(section) -> float:
    """Perimeter of a section in inches."""
    return section.length / IN


class TestWithoutOpenings:
    """Sec. 22.6.4.1 and 22.6.4.1.1: straight sides at d/2."""

    def test_interior_square(self):
        """4 x (24 + 20) = 176 in, four sides."""
        section = cs.critical_section_rectangular(C, C, D)
        assert bo_in(section) == pytest.approx(176.0, rel=1e-9)
        assert len(section.geoms) == 4

    def test_interior_rectangle(self):
        """12 in x 48 in column: 2 x (32 + 68) = 200 in."""
        section = cs.critical_section_rectangular(12.0 * IN, 48.0 * IN, D)
        assert bo_in(section) == pytest.approx(200.0, rel=1e-9)

    def test_edge_column_flush_with_edge(self):
        """Free edge at the +x column face.

        Two legs of 24 + 10 = 34 in and one side of 44 in: 112 in.
        """
        section = cs.critical_section_rectangular(
            C, C, D, edge_distances={'+x': 0.0}
        )
        assert bo_in(section) == pytest.approx(112.0, rel=1e-9)
        assert len(section.geoms) == 3

    def test_edge_column_set_back_from_edge(self):
        """Free edge 6 in beyond the +x face: legs of 40 in, 124 in."""
        section = cs.critical_section_rectangular(
            C, C, D, edge_distances={'+x': 6.0 * IN}
        )
        assert bo_in(section) == pytest.approx(124.0, rel=1e-9)

    def test_corner_column_flush(self):
        """Free edges at the +x and -y faces: 34 + 34 = 68 in."""
        section = cs.critical_section_rectangular(
            C, C, D, edge_distances={'+x': 0.0, '-y': 0.0}
        )
        assert bo_in(section) == pytest.approx(68.0, rel=1e-9)
        assert len(section.geoms) == 2

    def test_unknown_face_raises(self):
        """Faces are named '+x', '-x', '+y' and '-y'."""
        with pytest.raises(ValueError):
            cs.critical_section_rectangular(
                C, C, D, edge_distances={'north': 0.0}
            )

    def test_negative_edge_distance_raises(self):
        """The slab edge cannot be inside the column."""
        with pytest.raises(ValueError):
            cs.critical_section_rectangular(
                C, C, D, edge_distances={'+x': -1.0}
            )

    def test_four_free_edges_raises(self):
        """No side remains."""
        with pytest.raises(ValueError):
            cs.critical_section_rectangular(
                C,
                C,
                D,
                edge_distances={'+x': 0.0, '-x': 0.0, '+y': 0.0, '-y': 0.0},
            )


class TestOpenings:
    """Sec. 22.6.4.3: openings closer than 4h, here 4 x 24 = 96 in."""

    def test_opening_on_the_column_axis(self):
        """Opening 10 in square, 18 in clear of the +x face.

        Its near corners are at (30, +/-5). The tangent lines from the
        centroid cut the side at x = 22 at y = +/-22 x 5 / 30 = +/-3.667,
        so 7.333 in is ineffective: bo = 176 - 7.333 = 168.667 in.
        """
        section = cs.critical_section_rectangular(
            C, C, D, openings=[rect(30, 40, -5, 5)], h=H
        )
        assert bo_in(section) == pytest.approx(168.667, rel=1e-5)
        assert len(section.geoms) == 5

    def test_opening_beyond_4h_is_ignored(self):
        """Opening 100 in clear of the face: no reduction."""
        section = cs.critical_section_rectangular(
            C, C, D, openings=[rect(112, 122, -5, 5)], h=H
        )
        assert bo_in(section) == pytest.approx(176.0, rel=1e-9)

    def test_opening_at_exactly_4h_is_ignored(self):
        """The provision is for openings closer than 4h."""
        section = cs.critical_section_rectangular(
            C, C, D, openings=[rect(108, 118, -5, 5)], h=H
        )
        assert bo_in(section) == pytest.approx(176.0, rel=1e-9)

    def test_opening_on_the_diagonal(self):
        """Opening 10 in square with corners (30, 30) and (40, 40).

        The tangent lines pass through (40, 30) and (30, 40). The first
        meets x = 22 at y = 22 x 30 / 40 = 16.5, so 22 - 16.5 = 5.5 in of
        the +x side is ineffective, and likewise 5.5 in of the +y side:
        bo = 176 - 11 = 165 in.
        """
        section = cs.critical_section_rectangular(
            C, C, D, openings=[rect(30, 40, 30, 40)], h=H
        )
        assert bo_in(section) == pytest.approx(165.0, rel=1e-9)

    def test_two_openings(self):
        """The axis opening mirrored on -x: 176 - 2 x 7.333 = 161.333."""
        section = cs.critical_section_rectangular(
            C,
            C,
            D,
            openings=[rect(30, 40, -5, 5), rect(-40, -30, -5, 5)],
            h=H,
        )
        assert bo_in(section) == pytest.approx(161.333, rel=1e-5)

    def test_opening_given_as_coordinates(self):
        """An opening may be a sequence of (x, y) points."""
        coords = [
            (30 * IN, -5 * IN),
            (40 * IN, -5 * IN),
            (40 * IN, 5 * IN),
            (30 * IN, 5 * IN),
        ]
        section = cs.critical_section_rectangular(
            C, C, D, openings=[coords], h=H
        )
        assert bo_in(section) == pytest.approx(168.667, rel=1e-5)

    def test_opening_with_free_edge(self):
        """Edge column flush at +x, axis opening on the -x side.

        112 - 7.333 = 104.667 in.
        """
        section = cs.critical_section_rectangular(
            C,
            C,
            D,
            edge_distances={'+x': 0.0},
            openings=[rect(-40, -30, -5, 5)],
            h=H,
        )
        assert bo_in(section) == pytest.approx(104.667, rel=1e-5)

    def test_openings_without_h_raises(self):
        """The 4h rule needs the slab thickness."""
        with pytest.raises(ValueError):
            cs.critical_section_rectangular(
                C, C, D, openings=[rect(30, 40, -5, 5)]
            )

    def test_opening_over_the_centroid_raises(self):
        """No tangent lines exist from a centroid inside the opening."""
        with pytest.raises(ValueError):
            cs.critical_section_rectangular(
                C, C, D, openings=[rect(-5, 5, -5, 5)], h=H
            )
