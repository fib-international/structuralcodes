"""Tests for the Geometry."""

import math

import numpy as np
import pytest
from shapely import LineString, MultiLineString, MultiPolygon, Polygon
from shapely.affinity import translate
from shapely.testing import assert_geometries_equal

from structuralcodes.geometry import (
    CompoundGeometry,
    Geometry,
    PointGeometry,
    SurfaceGeometry,
    add_reinforcement,
    add_reinforcement_line,
    create_line_point_angle,
)
from structuralcodes.materials.concrete import ConcreteMC2010
from structuralcodes.materials.constitutive_laws import (
    Elastic,
    ElasticPlastic,
    ParabolaRectangle,
)
from structuralcodes.materials.reinforcement import ReinforcementMC2010


# Test create line
@pytest.mark.parametrize(
    'xp, yp, th, xb1, yb1, xb2, yb2, x1, y1, x2, y2',
    [
        (0, 0, 0, -10, -10, 10, 10, -10 - 1e-3, 0, 10 + 1e-3, 0),
        (0, 0.5, 0, -1.0, -1.0, 1.0, 1.0, -1 - 1e-3, 0.5, 1 + 1e-3, 0.5),
        (
            0,
            0,
            np.pi / 4,
            -1.0,
            -1.0,
            1.0,
            1.0,
            -1 - 1e-3,
            -1 - 1e-3,
            1 + 1e-3,
            1 + 1e-3,
        ),
        (0, 0, np.pi / 2, -1.0, -1.0, 1.0, 1.0, 0, -1 - 1e-3, 0, 1 + 1e-3),
    ],
)
def test_create_line_point_angle(
    xp, yp, th, xb1, yb1, xb2, yb2, x1, y1, x2, y2
):
    """Test creating a line from point and angle."""
    line = create_line_point_angle(
        point=(xp, yp), theta=th, bbox=(xb1, yb1, xb2, yb2)
    )
    assert_geometries_equal(line, LineString([(x1, y1), (x2, y2)]))


# Test PointGeometry
def test_point_geometry():
    """Test creating a PointGeometry object."""
    Geometry.section_counter = 0
    # Create a consitutive law to use
    steel = ElasticPlastic(210000, 450)

    # Create two points with default naming (uses global counter)
    for i in range(2):
        p = PointGeometry(np.array([2, 3]), 12, steel)
        assert p.name == f'Geometry_{i}'
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
            '555" stroke-width="0.02" fill="#66cc99" opacity="0.6" /></g>'
            '</svg>'
        )

    # Create two points with custom label for filtering
    for i in range(2):
        p = PointGeometry(np.array([2, 3]), 12, steel, group_label='Bottom')
        assert p.name == f'Geometry_{i+2}'
        assert p.group_label == 'Bottom'
        assert math.isclose(p.diameter, 12)
        assert math.isclose(p.point.coords[0][0], 2)
        assert math.isclose(p.point.coords[0][1], 3)

    # Create two points with custom name
    for i in range(2):
        p = PointGeometry(np.array([2, 3]), 12, steel, name='Rebar')
        assert p.name == 'Rebar'
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
    msg = 'Two coordinates are needed. 3'
    msg += ' coords provided. The extra entries will be'
    msg += ' discarded'
    assert str(record[0].message) == msg

    # Pass something else than a material or constitutive law
    with pytest.raises(TypeError) as excinfo:
        p = PointGeometry(np.array([2, 3]), 12, 'steel')

    # Trick for now since we don't hav a steel material
    C25 = ConcreteMC2010(25)
    p = PointGeometry(np.array([2, 3]), 12, C25)
    assert isinstance(p.material, ParabolaRectangle)


# Test Surface Geometry
def test_surface_geometry():  # noqa: PLR0915
    """Test creating a SurfaceGeometry object."""
    # Create a material to use
    C25 = ConcreteMC2010(25)

    # Create a constitutive law to use
    C25_const = ParabolaRectangle(25)

    # Create a rectangular geometry
    poly = Polygon(((0, 0), (200, 0), (200, 400), (0, 400)))
    for mat in (C25, C25_const):
        geo = SurfaceGeometry(poly, mat)
        assert isinstance(geo.material, ParabolaRectangle)
        assert geo.area == 200 * 400
        assert geo.centroid[0] == 100
        assert geo.centroid[1] == 200
        geo_t = geo.translate(-100, -200)
        assert isinstance(geo_t.material, ParabolaRectangle)
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
    # Not interecting: all below
    ab, bl = geo.split(line=((0, 410), 0))
    assert len(ab) == 0
    assert len(bl) == 1

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

    # Check assertion if passed more than two lines
    with pytest.raises(RuntimeError) as excinfo:
        geo_t.split_two_lines((line1, line2, line1))
    assert str(excinfo.value) == 'Two lines must be input'

    # Rotate the geometry
    geo_r = geo_t.rotate(np.pi / 2)
    assert_geometries_equal(
        geo_r.polygon,
        Polygon(
            ((200, -100), (200, 100), (-200, 100), (-200, -100), (200, -100))
        ),
    )

    # Pass something else than a polygon
    with pytest.raises(TypeError) as excinfo:
        SurfaceGeometry(poly=[0, 0, 1, 1], material=C25)
    assert (
        str(excinfo.value)
        == f'poly need to be a valid shapely.geometry.Polygon object. \
                {repr([0, 0, 1, 1])}'
    )
    # Pass something else than a polygon
    with pytest.raises(TypeError) as excinfo:
        SurfaceGeometry(poly=poly, material=1)
    assert (
        str(excinfo.value)
        == f'mat should be a valid structuralcodes.base.Material \
                or structuralcodes.base.ConstitutiveLaw object. \
                {repr(1)}'
    )

    # add two surfaces
    poly1 = Polygon(((-100, -200), (100, -200), (100, 200), (-100, 200)))
    poly2 = Polygon(((-500, 200), (500, 200), (500, 300), (-500, 300)))
    geo1 = SurfaceGeometry(poly=poly1, material=C25)
    geo2 = SurfaceGeometry(poly=poly2, material=C25)
    geo_add = geo1 + geo2
    assert isinstance(geo_add, CompoundGeometry)
    assert_geometries_equal(geo_add.geometries[0].polygon, geo1.polygon)
    assert_geometries_equal(geo_add.geometries[1].polygon, geo2.polygon)

    # check svg representation
    assert geo1._repr_svg_() == (
        '<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w'
        '3.org/1999/xlink" width="232.0" height="300" viewBox="-116.0 -216.0'
        ' 232.0 432.0" preserveAspectRatio="xMinYMin meet"><g transform='
        '"matrix(1,0,0,-1,0,0.0)"><path fill-rule="evenodd" fill="#66cc99" '
        'stroke="#555555" stroke-width="2.88" opacity="0.6" d="M -100.0,-200.0'
        ' L 100.0,-200.0 L 100.0,200.0 L -100.0,200.0 L -100.0,-200.0 z" />'
        '</g></svg>'
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

    # Create a CompoundGeometry with a MultiPolygon
    web = Polygon([(-100, -200), (100, -200), (100, 200), (-100, 200)])
    coords = [(-500, 200), (500, 200), (500, 300), (-500, 300)]
    flange = Polygon(coords)
    multi_pol = MultiPolygon([web, flange])
    # Create a compound geometry with the same material
    geo = CompoundGeometry(multi_pol, C25)
    assert_geometries_equal(geo.geometries[0].polygon, web)
    assert_geometries_equal(geo.geometries[1].polygon, flange)
    # Create a compound geometry with different materials
    geo = CompoundGeometry(multi_pol, [C25, steel])
    assert_geometries_equal(geo.geometries[0].polygon, web)
    assert_geometries_equal(geo.geometries[1].polygon, flange)
    assert isinstance(geo.geometries[0].material, ParabolaRectangle)
    assert isinstance(geo.geometries[1].material, ElasticPlastic)
    # check error is raised when the number of materials is incorrect
    with pytest.raises(ValueError) as excinfo:
        CompoundGeometry(multi_pol, [C25, C25, C25])
    assert (
        str(excinfo.value)
        == 'geometries and materials should have the same length'
    )

    # check svg representation
    assert geo._repr_svg_() == (
        '<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w'
        '3.org/1999/xlink" width="300" height="300" viewBox="-540.0 -240.0'
        ' 1080.0 580.0" preserveAspectRatio="xMinYMin meet"><g transform='
        '"matrix(1,0,0,-1,0,100.0)"><g><path fill-rule="evenodd" '
        'fill="#ff3333" '
        'stroke="#555555" stroke-width="7.2" opacity="0.6" d="M -100.0,-200.0'
        ' L 100.0,-200.0 L 100.0,200.0 L -100.0,200.0 L -100.0,-200.0 z" />'
        '<path fill-rule="evenodd" fill="#ff3333" stroke="#555555" '
        'stroke-width'
        '="7.2" opacity="0.6" d="M -500.0,200.0 L 500.0,200.0 L 500.0,300.0 L '
        '-500.0,300.0 L -500.0,200.0 z" /></g></g></svg>'
    )


def test_add_geometries():
    """Test addition for different geometries."""
    polys = []
    polys.append(
        Polygon([(-150, -300), (0, -300), (0, -289.3), (-150, -289.3)])
    )
    mat = ElasticPlastic(E=206000, fy=300)
    geo1 = SurfaceGeometry(polys[-1], mat)

    polys.append(Polygon([(-150, 289.3), (0, 289.3), (0, 300), (-150, 300)]))
    geo2 = SurfaceGeometry(polys[-1], mat)

    # add two surface geometries
    geo = geo1 + geo2
    assert isinstance(geo, CompoundGeometry)

    # add a further geometry
    polys.append(
        Polygon(
            [
                (-78.55, -289.3),
                (-71.45, -289.3),
                (-71.45, 289.3),
                (-78.55, 289.3),
            ]
        )
    )

    geo3 = SurfaceGeometry(polys[-1], mat)
    # add one compound geometry with a surface geometry
    geo += geo3
    assert isinstance(geo, CompoundGeometry)

    for i in range(3):
        polys.append(translate(polys[i], 150, 0))

    geo_t = geo.translate(dx=150, dy=0)

    # add one compound geometry with another compoud geometry
    geo += geo_t
    assert isinstance(geo, CompoundGeometry)

    # check the sum was fine
    assert len(polys) == len(geo.geometries)
    for i, p in enumerate(geo.geometries):
        assert_geometries_equal(polys[i], p.polygon)


def test_sub_geometries():
    """Test subtraction between geometries."""
    mat = ElasticPlastic(E=206000, fy=300)
    # Create the exptected polygons
    poly_1 = Polygon(
        shell=[(-100, -200), (100, -200), (100, 200), (-100, 200)],
        holes=[[(-50, -150), (50, -150), (50, -50), (-50, -50)]],
    )
    poly_2 = Polygon(
        shell=[(-100, -200), (100, -200), (100, 200), (-100, 200)],
        holes=[
            [(-50, -150), (50, -150), (50, -50), (-50, -50)],
            [(-50, 50), (50, 50), (50, 150), (-50, 150)],
        ],
    )
    # Create a rectangle
    geo_rect = SurfaceGeometry(
        Polygon([(-100, -200), (100, -200), (100, 200), (-100, 200)]), mat
    )

    # Create one hole
    geo_hole = SurfaceGeometry(
        Polygon([(-50, -150), (50, -150), (50, -50), (-50, -50)]), mat
    )

    # Create a compound geometry with two holes
    geo_holes = geo_hole + geo_hole.translate(dx=0, dy=200)

    # Subtract from surface one hole
    assert_geometries_equal(
        poly_1, (geo_rect - geo_hole).polygon, normalize=True
    )

    # Subtract from surface two holes (a compund)
    assert_geometries_equal(
        poly_2, (geo_rect - geo_holes).polygon, normalize=True
    )

    # Create a compound geometry
    geo_rects = geo_rect.translate(-100, 0) + geo_rect.translate(100, 0)
    geo_holes = geo_holes.translate(-100, 0) + geo_holes.translate(100, 0)

    # Subtract from Compound another compound
    geo_sub = geo_rects - geo_holes
    assert len(geo_sub.geometries) == 2
    assert_geometries_equal(
        geo_sub.geometries[0].polygon,
        translate(poly_2, xoff=-100),
        normalize=True,
    )
    assert_geometries_equal(
        geo_sub.geometries[1].polygon,
        translate(poly_2, xoff=100),
        normalize=True,
    )


@pytest.mark.parametrize(
    'w, h',
    [
        (100, 300),
        (200, 500),
        (400, 200),
        (500, 300),
    ],
)
def test_extents_calculation(w, h):
    """Test extents calculation for SurfaceGeometry and CompoundGeometry."""
    mat = Elastic(E=206000)
    # Create a rectangle
    geo_rect = SurfaceGeometry(
        Polygon(
            [
                (-w / 2, -h / 2),
                (w / 2, -h / 2),
                (w / 2, h / 2),
                (-w / 2, h / 2),
            ]
        ),
        mat,
    )
    x_min, x_max, y_min, y_max = geo_rect.calculate_extents()
    assert math.isclose(x_min, -w / 2, rel_tol=1e-5)
    assert math.isclose(x_max, w / 2, rel_tol=1e-5)
    assert math.isclose(y_min, -h / 2, rel_tol=1e-5)
    assert math.isclose(y_max, h / 2, rel_tol=1e-5)

    geo_rect = geo_rect.translate(w / 2, h / 2)
    x_min, x_max, y_min, y_max = geo_rect.calculate_extents()
    assert math.isclose(x_min, 0, abs_tol=1e-5)
    assert math.isclose(x_max, w, rel_tol=1e-5)
    assert math.isclose(y_min, 0, abs_tol=1e-5)
    assert math.isclose(y_max, h, rel_tol=1e-5)

    comp_geo = CompoundGeometry([geo_rect])
    x_min, x_max, y_min, y_max = comp_geo.calculate_extents()
    assert math.isclose(x_min, 0, abs_tol=1e-5)
    assert math.isclose(x_max, w, rel_tol=1e-5)
    assert math.isclose(y_min, 0, abs_tol=1e-5)
    assert math.isclose(y_max, h, rel_tol=1e-5)


@pytest.mark.parametrize(
    'w, h, c',
    [
        (100, 300, 40),
        (200, 500, 40),
        (400, 200, 40),
        (500, 300, 40),
    ],
)
def test_property_reinforced_concrete(w, h, c):
    """Test property reinforced_concrete."""
    mat = Elastic(E=206000)
    # Create a rectangle
    geo_rect = SurfaceGeometry(
        Polygon(
            [
                (-w / 2, -h / 2),
                (w / 2, -h / 2),
                (w / 2, h / 2),
                (-w / 2, h / 2),
            ]
        ),
        mat,
    )
    steel = Elastic(E=206000)
    geo_rc = add_reinforcement_line(
        geo_rect,
        (-w / 2 + c, -h / 2 + c),
        (w / 2 + c, -h / 2 + c),
        20,
        steel,
        n=4,
    )
    assert not geo_rc.reinforced_concrete

    geo_rect = SurfaceGeometry(
        poly=Polygon(
            [
                (-w / 2, -h / 2),
                (w / 2, -h / 2),
                (w / 2, h / 2),
                (-w / 2, h / 2),
            ]
        ),
        material=mat,
        concrete=True,
    )
    geo_rc = add_reinforcement_line(
        geo_rect,
        (-w / 2 + c, -h / 2 + c),
        (w / 2 + c, -h / 2 + c),
        20,
        steel,
        n=4,
    )
    assert geo_rc.reinforced_concrete

    mat = ConcreteMC2010(fck=35)
    geo_rect = SurfaceGeometry(
        poly=Polygon(
            [
                (-w / 2, -h / 2),
                (w / 2, -h / 2),
                (w / 2, h / 2),
                (-w / 2, h / 2),
            ]
        ),
        material=mat,
    )
    geo_rc = add_reinforcement_line(
        geo_rect,
        (-w / 2 + c, -h / 2 + c),
        (w / 2 + c, -h / 2 + c),
        20,
        steel,
        n=4,
    )
    assert geo_rc.reinforced_concrete


def test_reinforcement_group_label_one():
    """Test adding one reinforcement with a group label to a section."""
    # Arrange
    group_label = 'bottom-reinf'
    width = 200
    height = 500
    cover = 50
    diameter = 16

    concrete = ConcreteMC2010(fck=35)
    reinforcement = ReinforcementMC2010(
        fyk=500, ftk=520, epsuk=0.05, Es=200000
    )

    geometry = SurfaceGeometry(
        poly=Polygon(
            (
                (-width / 2, -height / 2),
                (width / 2, -height / 2),
                (width / 2, height / 2),
                (-width / 2, height / 2),
            )
        ),
        material=concrete,
    )

    # Act
    geometry = add_reinforcement(
        geometry,
        (0, -height / 2 + cover),
        diameter,
        reinforcement,
        group_label=group_label,
    )
    point_geometries = geometry.point_geometries

    # Assert
    assert len(point_geometries) == 1
    assert point_geometries[0].group_label == group_label


def test_reinforcement_group_label_line():
    """Test adding reinforcement with a group label along a line."""
    # Arrange
    group_label = 'bottom-reinf'
    width = 200
    height = 500
    cover = 50
    diameter = 16
    n_bars = 3

    concrete = ConcreteMC2010(fck=35)
    reinforcement = ReinforcementMC2010(
        fyk=500, ftk=520, epsuk=0.05, Es=200000
    )

    geometry = SurfaceGeometry(
        poly=Polygon(
            (
                (-width / 2, -height / 2),
                (width / 2, -height / 2),
                (width / 2, height / 2),
                (-width / 2, height / 2),
            )
        ),
        material=concrete,
    )

    # Act
    geometry = add_reinforcement_line(
        geometry,
        (-width / 2 + cover, -height / 2 + cover),
        (width / 2 - cover, -height / 2 + cover),
        diameter,
        reinforcement,
        n_bars,
        group_label=group_label,
    )
    point_geometries = geometry.point_geometries

    # Assert
    assert len(point_geometries) == n_bars
    for point in point_geometries:
        assert point.group_label == group_label


def test_surface_geometry_name_group_label():
    """Test the name and group label attribute of a SurfaceGeometry."""
    # Arrange
    width = 200
    height = 500
    concrete = ConcreteMC2010(fck=35)
    name = 'concrete_geometry'
    group_label = 'concrete'

    # Act
    geometry = SurfaceGeometry(
        poly=Polygon(
            (
                (-width / 2, -height / 2),
                (width / 2, -height / 2),
                (width / 2, height / 2),
                (-width / 2, height / 2),
            )
        ),
        material=concrete,
        name=name,
        group_label=group_label,
    )

    # Assert
    assert geometry.name == name
    assert geometry.group_label == group_label
