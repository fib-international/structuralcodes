"""Tests for the Marin integration."""

import math

import numpy as np

from structuralcodes.geometry import PointGeometry
from structuralcodes.materials.concrete import ConcreteMC2010


def test_point_geometry():
    """Test creating a PointGeometry object."""

    # Create two points with default naming (uses global counter)
    # TODO: to change with reinforcement as soon as we have its implementation
    C25 = ConcreteMC2010(25)
    for i in range(2):
        p = PointGeometry(np.array([2, 3]), 12, C25)
        assert p.name == f"Geometry_{i}"
        assert math.isclose(p.diameter, 12)
        assert math.isclose(p.point.coords[0][0], 2)
        assert math.isclose(p.point.coords[0][1], 3)
    # Create two points with custom label for filtering
    for i in range(2):
        p = PointGeometry(np.array([2, 3]), 12, C25, group_label='Bottom')
        assert p.name == f"Geometry_{i+2}"
        assert p.group_label == 'Bottom'
        assert math.isclose(p.diameter, 12)
        assert math.isclose(p.point.coords[0][0], 2)
        assert math.isclose(p.point.coords[0][1], 3)
