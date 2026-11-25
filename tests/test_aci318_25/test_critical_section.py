"""Tests for ACI 318-25 critical sections for two-way shear (Sec. 22.6.4).

The provisions are worded as in ACI 318-19 and the function delegates to
the ACI 318-19 one, which carries the full set of hand calculations. These
tests confirm the ACI 318-25 entry point on two of them.
"""

import pytest
from shapely.geometry import Polygon

from structuralcodes.codes import aci318_25

IN = 25.4  # mm per in


def test_interior_square():
    """24 in square column, d = 20 in: 4 x 44 = 176 in."""
    section = aci318_25.critical_section_rectangular(
        24.0 * IN, 24.0 * IN, 20.0 * IN
    )
    assert section.length / IN == pytest.approx(176.0, rel=1e-9)


def test_opening_on_the_column_axis():
    """Opening with near corners at (30, +/-5) in, h = 24 in.

    The tangent lines cut the side at x = 22 in at y = +/-3.667 in:
    bo = 176 - 7.333 = 168.667 in.
    """
    opening = Polygon(
        [
            (30 * IN, -5 * IN),
            (40 * IN, -5 * IN),
            (40 * IN, 5 * IN),
            (30 * IN, 5 * IN),
        ]
    )
    section = aci318_25.critical_section_rectangular(
        24.0 * IN, 24.0 * IN, 20.0 * IN, openings=[opening], h=24.0 * IN
    )
    assert section.length / IN == pytest.approx(168.667, rel=1e-5)
