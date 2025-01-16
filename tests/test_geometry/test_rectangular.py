"""Tests for the Geometry."""

import math

import numpy as np
import pytest

from structuralcodes.geometry import (
    RectangularGeometry,
    add_reinforcement_line,
)
from structuralcodes.materials.concrete import ConcreteMC2010
from structuralcodes.materials.constitutive_laws import Elastic
from structuralcodes.materials.reinforcement import ReinforcementMC2010


# Test create a rectangular geometry
@pytest.mark.parametrize(
    'w, h',
    [(200, 400), (200, 200), (400, 200), (0.2, 0.4), (0.4, 0.4), (0.4, 0.2)],
)
def test_create_rectangular_geometry(w, h):
    """Test creating a RectangularGeometry."""
    mat = Elastic(300000)
    rect = RectangularGeometry(w, h, mat)

    assert math.isclose(rect.width, w)
    assert math.isclose(rect.height, h)

    area = rect.area
    expected_area = w * h
    assert math.isclose(area, expected_area, rel_tol=0.05)


# Test wrong input
@pytest.mark.parametrize(
    'wrong_width, wrong_height',
    [
        (-100, -200),
        (0, -200),
        (0, 200),
        (-100, 200),
        (0, -200),
        (100, -200),
        (100, 0),
    ],
)
def test_create_rectangular_geometry_exception(wrong_width, wrong_height):
    """Test raising exception when inputing wrong value."""
    mat = Elastic(30000)
    with pytest.raises(ValueError):
        RectangularGeometry(wrong_width, wrong_height, mat)


# Test providing origin
@pytest.mark.parametrize(
    'origin',
    [(0, 0), (30, 600), (1.0,), None],
)
def test_rectangle_with_origin(origin):
    """Test creating a rectangle with an origin."""
    # Arrange
    mat = Elastic(30000)
    height = 600
    width = 200

    # Act and assert
    if origin is not None and len(origin) == 2:
        geom = RectangularGeometry(
            width=width, height=height, material=mat, origin=origin
        )
        assert np.allclose(
            (geom.polygon.centroid.xy[0][0], geom.polygon.centroid.xy[1][0]),
            origin,
        )
    elif origin is None:
        geom = RectangularGeometry(
            width=width, height=height, material=mat, origin=origin
        )
        assert np.allclose(
            (geom.polygon.centroid.xy[0][0], geom.polygon.centroid.xy[1][0]),
            (0, 0),
        )
    else:
        with pytest.raises(ValueError):
            geom = RectangularGeometry(
                width=width, height=height, material=mat, origin=origin
            )


# Test adding reinforcement in circular pattern
# Whole circle
@pytest.mark.parametrize('w, h, c, n', [(200, 400, 40, 2), (400, 400, 50, 4)])
def test_reinforced_rectangular_geometry(w, h, c, n):
    """Test adding reinforcemnet to a RectangularGeometry."""
    # Create materials
    concrete = ConcreteMC2010(40)
    steel = ReinforcementMC2010(450, 200000, 450, 0.065)

    # create the circular geometry
    rect = RectangularGeometry(w, h, concrete)

    assert rect.concrete

    # add reinforcement in the circular pattern
    diam = 16
    rc = add_reinforcement_line(
        geo=rect,
        coords_i=(-w / 2 + c, -h / 2 + c),
        coords_j=(w / 2 - c, -h / 2 + c),
        diameter=diam,
        material=steel,
        n=n,
    )
    rc = add_reinforcement_line(
        geo=rc,
        coords_i=(-w / 2 + c, h / 2 - c),
        coords_j=(w / 2 - c, h / 2 - c),
        diameter=diam,
        material=steel,
        n=n,
    )
    assert len(rc.point_geometries) == 2 * n
