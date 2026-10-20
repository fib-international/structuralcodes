"""Tests for adding Geometry objects."""

import typing as t

import pytest
from shapely import Polygon

from structuralcodes.geometry import (
    CompoundGeometry,
    PointGeometry,
    SurfaceGeometry,
)
from structuralcodes.materials.basic import ElasticMaterial


@pytest.fixture
def two_rectangles() -> t.Tuple[SurfaceGeometry]:
    """Fixture for two rectangles."""
    concrete = ElasticMaterial(E=30000, density=2500)
    width_web = 250
    width_flange = 1000
    height_web = 650
    height_flange = 150
    polygon_web = Polygon(
        (
            (-width_web / 2, 0),
            (width_web / 2, 0),
            (width_web / 2, height_web),
            (-width_web / 2, height_web),
        )
    )
    polygon_flange = Polygon(
        (
            (-width_flange / 2, height_web),
            (width_flange / 2, height_web),
            (width_flange / 2, height_web + height_flange),
            (-width_flange / 2, height_web + height_flange),
        ),
    )

    return SurfaceGeometry(
        poly=polygon_web, material=concrete
    ), SurfaceGeometry(poly=polygon_flange, material=concrete)


@pytest.fixture
def two_points() -> t.Tuple[PointGeometry]:
    """Fixture for two points."""
    steel = ElasticMaterial(E=200000, density=7850)
    return PointGeometry(
        point=(0, 50), diameter=16, material=steel
    ), PointGeometry(point=(100, 50), diameter=16, material=steel)


def test_add_surface_geometries(two_rectangles):
    """Test adding two surface geometries."""
    geom_web, geom_flange = two_rectangles

    compound = geom_web + geom_flange

    assert isinstance(compound, CompoundGeometry)
    assert geom_web, geom_flange in compound.geometries


def test_add_surface_and_point_geometries(two_rectangles, two_points):
    """Test adding a surface geometry and a point geometry."""
    surface, _ = two_rectangles
    point, _ = two_points

    # Add surface and point
    compound = surface + point

    assert isinstance(compound, CompoundGeometry)
    assert surface in compound.geometries
    assert point in compound.point_geometries
    assert len(compound.geometries) == 1
    assert len(compound.point_geometries) == 1

    # Add point and surface
    compound = point + surface

    assert isinstance(compound, CompoundGeometry)
    assert surface in compound.geometries
    assert point in compound.point_geometries
    assert len(compound.geometries) == 1
    assert len(compound.point_geometries) == 1


def test_add_point_geometries(two_points):
    """Test adding two point geometries."""
    point_1, point_2 = two_points

    compound = point_1 + point_2

    assert isinstance(compound, CompoundGeometry)
    assert point_1, point_2 in compound.point_geometries


def test_add_to_compound_geometries(two_rectangles, two_points):
    """Test adding to compound geometries."""
    geom_web, geom_flange = two_rectangles
    point_1, _ = two_points

    compound_1 = geom_web + geom_flange
    compound_2 = geom_web + geom_flange
    compound_3 = compound_1 + compound_2

    assert isinstance(compound_1, CompoundGeometry)
    assert isinstance(compound_2, CompoundGeometry)
    assert isinstance(compound_3, CompoundGeometry)
    assert len(compound_1.geometries) == 2
    assert len(compound_2.geometries) == 2
    assert len(compound_3.geometries) == 4
    assert len(compound_3.point_geometries) == 0

    # Add a point and a surface
    compound_4 = compound_3 + point_1
    compound_4 += geom_web
    assert isinstance(compound_4, CompoundGeometry)
    assert len(compound_4.geometries) == 5
    assert len(compound_4.point_geometries) == 1

    # Add a compound to a point
    compound_5 = point_1 + compound_4
    assert isinstance(compound_5, CompoundGeometry)
    assert len(compound_5.geometries) == 5
    assert len(compound_5.point_geometries) == 2

    # Add a compound to a surface
    compound_6 = geom_web + compound_5
    assert isinstance(compound_6, CompoundGeometry)
    assert len(compound_6.geometries) == 6
    assert len(compound_6.point_geometries) == 2
