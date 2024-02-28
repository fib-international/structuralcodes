"""Tests for the Geometry."""

import math

import numpy as np

import pytest

from structuralcodes.geometry import (
    Geometry,
    PointGeometry,
    SurfaceGeometry,
    create_line_point_angle,
)
from structuralcodes.sections._reinforcement import (
    add_reinforcement,
    add_reinforcement_line,
)
from structuralcodes.materials.concrete import ConcreteMC2010
from structuralcodes.materials.constitutive_laws import (
    ElasticPlastic,
    ParabolaRectangle,
)
from shapely import Polygon, MultiLineString, LineString
from shapely.testing import assert_geometries_equal


# Test PointGeometry
def test_point_geometry():
    """Test creating a PointGeometry object."""

    Geometry.section_counter = 0
    # Create a consitutive law to use
    steel = ElasticPlastic(210000, 450)

    # Create two points with default naming (uses global counter)
    for i in range(2):
        p = PointGeometry(np.array([2, 3]), 12, steel)
        assert p.name == f"Geometry_{i}"
        assert math.isclose(p.diameter, 12)
        assert math.isclose(p.point.coords[0][0], 2)
        assert math.isclose(p.point.coords[0][1], 3)
        assert math.isclose(p.x, 2)
        assert math.isclose(p.y, 3)
        assert math.isclose(p.area, 6**2 * np.pi)
        assert p._repr_svg_() == (
            '<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w'
            '3.org/1999/xlink" width="100.0" height="100.0" viewBox="1.0 2.0 '
            '2.0 2.0" preserveAspectRatio="xMinYMin meet"><g transform="matrix'
            '(1,0,0,-1,0,6.0)"><circle cx="2.0" cy="3.0" r="0.06" stroke="#555'
            '555" stroke-width="0.02" fill="#66cc99" opacity="0.6" /></g></svg>'
        )

    # Create two points with custom label for filtering
    for i in range(2):
        p = PointGeometry(np.array([2, 3]), 12, steel, group_label='Bottom')
        assert p.name == f"Geometry_{i+2}"
        assert p.group_label == 'Bottom'
        assert math.isclose(p.diameter, 12)
        assert math.isclose(p.point.coords[0][0], 2)
        assert math.isclose(p.point.coords[0][1], 3)

    # Create two points with custom name
    for i in range(2):
        p = PointGeometry(np.array([2, 3]), 12, steel, name='Rebar')
        assert p.name == f"Rebar"
        assert math.isclose(p.diameter, 12)
        assert math.isclose(p.point.coords[0][0], 2)
        assert math.isclose(p.point.coords[0][1], 3)

    # Pass less than two coordinates
    with pytest.raises(ValueError) as excinfo:
        p = PointGeometry(np.array([2]), 12, steel)
    assert str(excinfo.value) == 'Two coordinates are needed'

    # Pass more than two coordinates: check that a warning is thrown
    with pytest.warns(UserWarning) as record:
        p = PointGeometry(np.array([2, 3, 4]), 12, steel)
    assert len(record) == 1
    msg = f'Two coordinates are needed. 3'
    msg += ' coords provided. The extra entries will be'
    msg += ' discarded'
    assert str(record[0].message) == msg

    # Pass something else than a material or constitutive law
    with pytest.raises(TypeError) as excinfo:
        p = PointGeometry(np.array([2, 3]), 12, 'steel')

    # Trick for now since we don't hav a steel material
    C25 = ConcreteMC2010(25)
    p = PointGeometry(np.array([2, 3]), 12, C25)
    assert isinstance(p.material, ParabolaRectangle) == True


# Test Surface Geometry
def test_surface_geometry():
    """Test creating a SurfaceGeometry object."""

    # Create a material to use
    C25 = ConcreteMC2010(25)

    # Create a constitutive law to use
    C25_const = ParabolaRectangle(25)

    # Create a rectangular geometry
    poly = Polygon(((0, 0), (200, 0), (200, 400), (0, 400)))
    for mat in (C25, C25_const):
        geo = SurfaceGeometry(poly, mat)
        assert isinstance(geo.material, ParabolaRectangle) == True
        assert geo.area == 200 * 400
        assert geo.centroid[0] == 100
        assert geo.centroid[1] == 200
        geo_t = geo.translate(-100, -200)
        assert isinstance(geo_t.material, ParabolaRectangle) == True
        assert geo_t.area == 200 * 400
        assert geo_t.centroid[0] == 0
        assert geo_t.centroid[1] == 0

    # Test splitting of polygon
    ab, bl = geo_t.split(line=((0, 0), 0))
    assert len(ab) == 1
    assert len(bl) == 1
    assert_geometries_equal(
        ab[0],
        Polygon(((-100, 0), (-100, 200), (100, 200), (100, 0), (-100, 0))),
    )
    assert_geometries_equal(
        bl[0],
        Polygon(((100, 0), (100, -200), (-100, -200), (-100, 0), (100, 0))),
    )
    # Not intersecting: all above
    ab, bl = geo.split(line=((0, -10), 0))
    assert len(ab) == 1
    assert len(bl) == 0

    # Splitting with two lines
    line1 = LineString([(-120, -50), (120, -50)])
    line2 = LineString([(-120, -60), (120, -60)])
    lines = MultiLineString((line1, line2))
    result = geo_t.split_two_lines(lines)
    assert_geometries_equal(
        result,
        Polygon(
            ((100, -60), (-100, -60), (-100, -50), (100, -50), (100, -60))
        ),
    )

    # Alternatively one can create the Lines
    line1 = create_line_point_angle((0, 0), 0, geo_t.polygon.bounds)
    line2 = create_line_point_angle((0, 5), 0, geo_t.polygon.bounds)
    result = geo_t.split_two_lines((line1, line2))
    assert_geometries_equal(
        result, Polygon(((100, 0), (-100, 0), (-100, 5), (100, 5), (100, 0)))
    )

    # Rotate the geometry
    geo_r = geo_t.rotate(np.pi / 2)
    assert_geometries_equal(
        geo_r.polygon,
        Polygon(
            ((200, -100), (200, 100), (-200, 100), (-200, -100), (200, -100))
        ),
    )


def test_compound_geometry():
    """Test creating a SurfaceGeometry object."""

    # Create a material to use
    C25 = ConcreteMC2010(25)
    steel = ElasticPlastic(210000, 450)

    # Create a rectangular geometry
    poly = Polygon(((0, 0), (200, 0), (200, 400), (0, 400)))
    # Create a Surface Geometry
    geo = SurfaceGeometry(poly, C25)
    # Add single reinforcement
    geo = add_reinforcement(geo, (40, 360), 20, steel)
    geo = add_reinforcement(geo, (160, 360), 20, steel)
    # Add a line of reinforcement
    geo2 = add_reinforcement_line(geo, (40, 40), (160, 40), 20, steel, n=4)
    assert len(geo2.geometries) == 1
    assert len(geo2.point_geometries) == 6

    # it is possible to add also a line specifying the spacing
    geo2 = add_reinforcement_line(geo, (40, 40), (160, 40), 20, steel, s=30)
    assert len(geo2.geometries) == 1
    assert len(geo2.point_geometries) == 7

    # it is possible to add also a line specifying the spacing and the number
    geo2 = add_reinforcement_line(
        geo, (40, 40), (160, 40), 20, steel, n=3, s=30
    )
    assert len(geo2.geometries) == 1
    assert len(geo2.point_geometries) == 5

    # it is possible to add also to skip the first and/or last bar
    geo2 = add_reinforcement_line(
        geo, (40, 40), (160, 40), 20, steel, n=3, s=30, first=False, last=False
    )
    assert len(geo2.geometries) == 1
    assert len(geo2.point_geometries) == 3

    geo = add_reinforcement_line(geo, (40, 40), (160, 40), 20, steel, n=4)
    assert len(geo.geometries) == 1
    assert len(geo.point_geometries) == 6

    # Translate the whole CompoundGeometry
    geo_t = geo.translate(-100, -200)
    assert math.isclose(geo_t.area, 200 * 400)
    assert math.isclose(geo_t.geometries[0].area, 200 * 400)
    assert math.isclose(geo_t.geometries[0].centroid[0], 0)
    assert math.isclose(geo_t.geometries[0].centroid[1], 0)

    # Rotate the whole CompountGeometry
    geo_r = geo_t.rotate(np.pi / 2)
    assert math.isclose(geo_r.point_geometries[0].x, -160)
    assert math.isclose(geo_r.point_geometries[0].y, -60)
