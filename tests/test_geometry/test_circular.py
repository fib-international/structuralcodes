"""Tests for the Geometry."""

import math

import numpy as np
import pytest

from structuralcodes.geometry import CircularGeometry, add_reinforcement_circle
from structuralcodes.materials.concrete import ConcreteMC2010
from structuralcodes.materials.constitutive_laws import Elastic
from structuralcodes.materials.reinforcement import ReinforcementMC2010


# Test create a circular geometry
@pytest.mark.parametrize(
    'diameter, n_points',
    [(100, 20), (200, 20), (300, 20), (300, 30), (400, 40), (500, 24)],
)
def test_create_circular_geometry(diameter, n_points):
    """Test creating a CircularGeometry."""
    mat = Elastic(300000)
    circle = CircularGeometry(diameter, mat, n_points)

    assert len(circle.polygon.exterior.coords) == n_points + 2

    assert math.isclose(circle.diameter, diameter)
    assert math.isclose(circle.radius, diameter / 2.0)

    area = circle.area
    expected_area = np.pi * diameter**2 / 4.0
    assert math.isclose(area, expected_area, rel_tol=0.05)


# Test wrong input
@pytest.mark.parametrize('wrong_diameter', [-100, 0])
def test_create_circular_geometry_exception(wrong_diameter):
    """Test raising exception when inputing wrong value."""
    mat = Elastic(300000)
    with pytest.raises(ValueError):
        CircularGeometry(wrong_diameter, mat)


# Test providing origin
@pytest.mark.parametrize(
    'origin',
    [(0, 0), (30, 600), (1.0,)],
)
def test_circle_with_origin(origin):
    """Test creating a circle with an origin."""
    # Arrange
    mat = Elastic(300000)
    diameter = 100

    # Act and assert
    if origin is not None and len(origin) == 2:
        geom = CircularGeometry(diameter=diameter, material=mat, origin=origin)
        assert np.allclose(
            (geom.polygon.centroid.xy[0][0], geom.polygon.centroid.xy[1][0]),
            origin,
        )
    elif origin is None:
        geom = CircularGeometry(diameter=diameter, material=mat, origin=origin)
        assert np.allclose(
            (geom.polygon.centroid.xy[0][0], geom.polygon.centroid.xy[1][0]),
            (0, 0),
        )
    else:
        with pytest.raises(ValueError):
            geom = CircularGeometry(
                diameter=diameter, material=mat, origin=origin
            )


# Test adding reinforcement in circular pattern
# Whole circle
@pytest.mark.parametrize('diameter, cover, n', [(200, 40, 8), (600, 50, 16)])
def test_reinforced_circular_geometry(diameter, cover, n):
    """Test adding reinforcemnet to a CircularGeometry."""
    # Create materials
    concrete = ConcreteMC2010(40)
    steel = ReinforcementMC2010(450, 200000, 450, 0.065)

    # create the circular geometry
    circle = CircularGeometry(diameter, concrete)

    assert circle.concrete

    # add reinforcement in the circular pattern
    rc = add_reinforcement_circle(
        geo=circle,
        center=(0, 0),
        radius=circle.radius - cover,
        diameter=16,
        material=steel,
        n=n,
    )

    assert len(rc.point_geometries) == n


# Test adding reinforcement in circular pattern
# circular arch
@pytest.mark.parametrize('diameter, cover, n', [(200, 40, 8), (600, 50, 16)])
def test_reinforced_circular_geometry_arch(diameter, cover, n):
    """Test adding reinforcemnet to a CircularGeometry on circle arches."""
    # Create materials
    concrete = ConcreteMC2010(40)
    steel = ReinforcementMC2010(450, 200000, 450, 0.065)

    # create the circular geometry
    circle = CircularGeometry(diameter, concrete)

    assert circle.concrete

    # add reinforcement in the circular arch pattern
    rc = add_reinforcement_circle(
        geo=circle,
        center=(0, 0),
        radius=circle.radius - cover,
        diameter=16,
        material=steel,
        n=int(n / 2) + 1,
        stop_angle=np.pi,
    )

    assert len(rc.point_geometries) == int(n / 2) + 1

    # add further reinforcement in a circular arch pattern

    rc = add_reinforcement_circle(
        geo=rc,
        center=(0, 0),
        radius=circle.radius - cover,
        diameter=16,
        material=steel,
        n=int(n / 2) + 1,
        start_angle=np.pi,
        stop_angle=np.pi * 2,
        first=False,
        last=False,
    )

    # These two operations should give same result as reinforcement in whole
    # circumference. Let's check this

    assert len(rc.point_geometries) == n


# Test adding reinforcement in circular pattern
# circular arch
@pytest.mark.parametrize('diameter, cover, n', [(200, 40, 8), (600, 50, 16)])
def test_reinforced_circular_geometry_exception(diameter, cover, n):
    """Test excpetion when adding reinforcemnet to a CircularGeometry."""
    # Create materials
    concrete = ConcreteMC2010(40)
    steel = ReinforcementMC2010(450, 200000, 450, 0.065)

    # create the circular geometry
    circle = CircularGeometry(diameter, concrete)

    assert circle.concrete

    # add reinforcement in the circular arch pattern
    with pytest.raises(ValueError):
        _ = add_reinforcement_circle(
            geo=circle,
            center=(0, 0),
            radius=circle.radius - cover,
            diameter=16,
            material=steel,
            n=int(n / 2) + 1,
            stop_angle=-np.pi,
        )


# add reinforcement providing spacing
def test_reinforcement_circular_spacing():
    """Test adding reinforcement in circular pattern providing spacing."""
    # Create materials
    concrete = ConcreteMC2010(40)
    steel = ReinforcementMC2010(450, 200000, 450, 0.065)

    # create the circular geometry
    circle = CircularGeometry(300, concrete)

    assert circle.concrete

    # add reinforcement in the circular arch pattern
    rc = add_reinforcement_circle(
        geo=circle,
        center=(0, 0),
        radius=120,
        diameter=16,
        material=steel,
        s=50,
        start_angle=0,
        stop_angle=np.pi,
    )

    # This should create 8 bars
    assert len(rc.point_geometries) == 8


# add reinforcement providing spacing and number
def test_reinforcement_circular_spacing_number():
    """Test adding reinforcement in circular pattern providing spacing and
    number.
    """
    # Create materials
    concrete = ConcreteMC2010(40)
    steel = ReinforcementMC2010(450, 200000, 450, 0.065)

    # create the circular geometry
    circle = CircularGeometry(300, concrete)

    assert circle.concrete

    # add reinforcement in the circular arch pattern
    rc = add_reinforcement_circle(
        geo=circle,
        center=(0, 0),
        radius=120,
        diameter=16,
        material=steel,
        s=50,
        n=6,
        start_angle=0,
        stop_angle=np.pi,
    )

    # This should create 6 bars as requested
    assert len(rc.point_geometries) == 6

    # The first bar should be in position 103.589, 60.574
    assert math.isclose(
        rc.point_geometries[0].point.coords[0][0], 103.58960753977033
    )
    assert math.isclose(
        rc.point_geometries[0].point.coords[0][1], 60.57386573231362
    )


# add reinforcement providing spacing and number
def test_reinforcement_circular_spacing_number_error():
    """Test raising exception when adding reinforcement in circular pattern
    providing spacing and number.
    """
    # Create materials
    concrete = ConcreteMC2010(40)
    steel = ReinforcementMC2010(450, 200000, 450, 0.065)

    # create the circular geometry
    circle = CircularGeometry(300, concrete)

    assert circle.concrete

    # add reinforcement in the circular arch pattern
    with pytest.raises(ValueError):
        _ = add_reinforcement_circle(
            geo=circle,
            center=(0, 0),
            radius=120,
            diameter=16,
            material=steel,
            s=50,
            n=10,
            start_angle=0,
            stop_angle=np.pi,
        )
