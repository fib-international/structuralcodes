"""Tests for the Geometry."""

import math

import pytest
from shapely import Point, Polygon, set_precision
from shapely.affinity import scale, translate
from shapely.ops import unary_union

from structuralcodes.geometry import CompoundGeometry, SurfaceGeometry
from structuralcodes.materials.concrete import ConcreteMC2010
from structuralcodes.materials.constitutive_laws import (
    Elastic,
    ElasticPlastic,
    Sargin,
    UserDefined,
)
from structuralcodes.materials.reinforcement import ReinforcementMC2010
from structuralcodes.sections._generic import GenericSection
from structuralcodes.sections._reinforcement import (
    add_reinforcement,
    add_reinforcement_line,
)


# Test rectangular section
def test_rectangular_section():
    """Test rectangular section."""
    # Create materials to use
    concrete = ConcreteMC2010(25)
    steel = ReinforcementMC2010(fyk=450, Es=210000, ftk=450, epsuk=0.0675)

    # The section
    poly = Polygon(((0, 0), (200, 0), (200, 400), (0, 400)))
    geo = SurfaceGeometry(poly, concrete)
    geo = add_reinforcement_line(geo, (40, 40), (160, 40), 16, steel, n=4)
    geo = add_reinforcement_line(geo, (40, 360), (160, 360), 16, steel, n=4)
    geo = geo.translate(-100, -200)
    assert geo.geometries[0].centroid[0] == 0
    assert geo.geometries[0].centroid[1] == 0

    # Create the section (default Marin integrator)
    sec = GenericSection(geo)
    assert sec.name == 'GenericSection'

    assert math.isclose(sec.gross_properties.area, 200 * 400)

    # Compute bending strength
    res_marin = sec.section_analyzer.calculate_bending_strength(theta=0, n=0)

    # Compute moment curvature
    res_mc_marin = sec.section_analyzer.calculate_moment_curvature(
        theta=0, n=0
    )

    # Use fiber integration
    sec = GenericSection(geo, integrator='Fiber', mesh_size=0.0001)
    assert math.isclose(sec.gross_properties.area, 200 * 400)

    # Compute bending strength
    res_fiber = sec.section_analyzer.calculate_bending_strength(theta=0, n=0)

    assert math.isclose(res_marin.m_y, res_fiber.m_y, rel_tol=1e-2)

    # Compute moment curvature
    res_mc_fiber = sec.section_analyzer.calculate_moment_curvature(
        theta=0, n=0
    )

    assert math.isclose(
        res_mc_marin.m_y[-1], res_mc_fiber.m_y[-1], rel_tol=1e-2
    )


# Test rectangular section with Sargin Model
def test_rectangular_section_Sargin():
    """Test rectangular section."""
    # Create materials to use
    concrete = ConcreteMC2010(25)
    steel = ReinforcementMC2010(fyk=450, Es=210000, ftk=450, epsuk=0.0675)
    # Set a different constitutive law respect to default Parabola-Rectangle
    # Here we use Sargin law (MC2010 eq 5.1-26) with parameters taken from
    # MC2010 table 5.1-8
    concrete.constitutive_law = Sargin(
        fc=-35, eps_c1=-2.3e-3, eps_cu1=-3.5e-3, k=1.92
    )

    # The section
    poly = Polygon(((0, 0), (200, 0), (200, 400), (0, 400)))
    geo = SurfaceGeometry(poly, concrete)
    geo = add_reinforcement_line(geo, (40, 40), (160, 40), 16, steel, n=4)
    geo = add_reinforcement_line(geo, (40, 360), (160, 360), 16, steel, n=4)
    geo = geo.translate(-100, -200)
    assert geo.geometries[0].centroid[0] == 0
    assert geo.geometries[0].centroid[1] == 0

    # Create the section (default Marin integrator)
    sec = GenericSection(geo)

    assert math.isclose(sec.gross_properties.area, 200 * 400)

    # Compute bending strength
    res_marin = sec.section_analyzer.calculate_bending_strength(theta=0, n=0)

    # Use fiber integration
    sec = GenericSection(geo, integrator='Fiber', mesh_size=0.0001)
    assert math.isclose(sec.gross_properties.area, 200 * 400)

    # Compute bending strength
    res_fiber = sec.section_analyzer.calculate_bending_strength(theta=0, n=0)

    assert math.isclose(res_marin.m_y, res_fiber.m_y, rel_tol=1e-2)


# Test holed section
def test_holed_section():
    """Test a section with a hole."""
    # Create materials to use
    concrete = ConcreteMC2010(25)
    steel = ReinforcementMC2010(fyk=450, Es=210000, ftk=450, epsuk=0.0675)

    # The section
    poly = Polygon(
        ((0, 0), (500, 0), (500, 1000), (0, 1000)),
        [((100, 100), (400, 100), (400, 900), (100, 900))],
    )
    geo = SurfaceGeometry(poly, concrete)

    # Add reinforcement
    geo = add_reinforcement_line(geo, (50, 50), (450, 50), 28, steel, n=4)
    geo = add_reinforcement_line(geo, (50, 950), (450, 950), 28, steel, n=4)
    # Translate the section
    geo = geo.translate(-250, -500)
    assert geo.geometries[0].centroid[0] == 0
    assert geo.geometries[0].centroid[1] == 0

    # Create the section (default Marin integrator)
    sec = GenericSection(geo)
    assert sec.name == 'GenericSection'

    assert math.isclose(sec.gross_properties.area, 260000)

    # Compute bending strength
    res_marin = sec.section_analyzer.calculate_bending_strength(theta=0, n=0)

    assert math.isclose(abs(res_marin.m_y * 1e-6), 1012, rel_tol=1e-2)

    # Use fiber integration
    sec = GenericSection(geo, integrator='Fiber', mesh_size=0.0001)
    assert math.isclose(sec.gross_properties.area, 260000)

    # Compute bending strength
    res_fiber = sec.section_analyzer.calculate_bending_strength(theta=0, n=0)

    assert math.isclose(res_marin.m_y, res_fiber.m_y, rel_tol=1e-2)


# Test U section
def test_u_section():
    """Test a section with a U shape."""
    # Create materials to use
    concrete = ConcreteMC2010(25)
    steel = ReinforcementMC2010(fyk=450, Es=210000, ftk=450, epsuk=0.0675)

    # The section
    # bottom flange
    poly = Polygon(
        ((0, 0), (500, 0), (500, 100), (0, 100)),
    )
    poly = poly.union(Polygon(((0, 100), (100, 100), (100, 1000), (0, 1000))))
    poly = poly.union(
        Polygon(((400, 100), (500, 100), (500, 1000), (400, 1000)))
    )
    geo = SurfaceGeometry(poly, concrete)

    # Add reinforcement
    geo = add_reinforcement_line(geo, (50, 50), (450, 50), 28, steel, n=4)
    geo = add_reinforcement(geo, (50, 950), 28, steel)
    geo = add_reinforcement(geo, (450, 950), 28, steel)
    # Translate the section
    geo = geo.translate(-250, -441.30434782608694)
    assert geo.geometries[0].centroid[0] == 0
    assert geo.geometries[0].centroid[1] == 0

    # Create the section (default Marin integrator)
    sec = GenericSection(geo)
    assert sec.name == 'GenericSection'

    assert math.isclose(sec.gross_properties.area, 230000)

    # Compute bending strength
    res_marin = sec.section_analyzer.calculate_bending_strength(theta=0, n=0)

    assert math.isclose(abs(res_marin.m_y * 1e-6), 993.3, rel_tol=1e-2)

    # Use fiber integration
    sec = GenericSection(geo, integrator='Fiber', mesh_size=0.0001)
    assert math.isclose(sec.gross_properties.area, 230000)

    # Compute bending strength
    res_fiber = sec.section_analyzer.calculate_bending_strength(theta=0, n=0)

    assert math.isclose(res_marin.m_y, res_fiber.m_y, rel_tol=1e-2)


# Test steel I section
def _iyy_IPE(h: float, b: float, tw: float, tf: float, r: float) -> float:
    return (
        1 / 12 * b * h**3
        - 1 / 12 * (b - tw) * (h - 2 * tf) ** 3
        + 4 * (r**4 / 12.0 + r**2 * (h / 2 - tf - r / 2) ** 2)
        - 4
        * (
            (math.pi / 16 - 4 / (9 * math.pi)) * r**4
            + math.pi
            * r**2
            / 4
            * (h / 2 - tf - r * (3 * math.pi - 4) / (3 * math.pi)) ** 2
        )
    )


def _izz_IPE(h: float, b: float, tw: float, tf: float, r: float) -> float:
    return (
        2 / 12 * tf * b**3
        + 1 / 12 * (h - 2 * tf) * tw**3
        + 4 * (r**4 / 12.0 + r**2 * ((tw + r) / 2) ** 2)
        - 4
        * (
            (math.pi / 16 - 4 / (9 * math.pi)) * r**4
            + math.pi
            * r**2
            / 4
            * (tw / 2 + r * (3 * math.pi - 4) / (3 * math.pi)) ** 2
        )
    )


def _wyel_IPE(h: float, b: float, tw: float, tf: float, r: float) -> float:
    return _iyy_IPE(h, b, tw, tf, r) / (h / 2.0)


def _wzel_IPE(h: float, b: float, tw: float, tf: float, r: float) -> float:
    return _izz_IPE(h, b, tw, tf, r) / (b / 2.0)


def _wypl_IPE(h: float, b: float, tw: float, tf: float, r: float) -> float:
    return 2 * (
        tw * h**2 / 8.0
        + tf * (b - tw) * (h - tf) / 2.0
        + 2.0 * r**2 * (h / 2.0 - tf - r / 2.0)
        - 2
        * math.pi
        * r**2
        / 4.0
        * (h / 2 - tf - r * (3 * math.pi - 4) / (3 * math.pi))
    )


def _wzpl_IPE(h: float, b: float, tw: float, tf: float, r: float) -> float:
    return 2 * (
        2 * tf * b**2 / 8.0
        + tw**2 / 8.0 * (h - 2 * tf)
        + 2 * r**2 * (tw + r) / 2.0
        - 2
        * math.pi
        * r**2
        / 4
        * (tw / 2 + r * (3 * math.pi - 4) / (3 * math.pi))
    )


def _create_I_section(h: float, b: float, tw: float, tf: float, r: float):
    # create the geometry of section
    # top flange
    top_flange = Polygon(
        [
            (-b / 2, -h / 2),
            (b / 2, -h / 2),
            (b / 2, -h / 2 + tf),
            (-b / 2, -h / 2 + tf),
        ]
    )
    # bottom flange
    bottom_flange = translate(top_flange, xoff=0, yoff=h - tf)
    web = Polygon(
        [
            (-tw / 2, -h / 2 + tf),
            (tw / 2, -h / 2 + tf),
            (tw / 2, h / 2 - tf),
            (-tw / 2, h / 2 - tf),
        ]
    )
    # fillets
    p = Point([tw / 2 + r, -h / 2 + tf + r]).buffer(r)
    s = Polygon(
        [
            (tw / 2, -h / 2 + tf),
            (tw / 2 + r, -h / 2 + tf),
            (tw / 2 + r, -h / 2 + tf + r),
            (tw / 2, -h / 2 + tf + r),
        ]
    )
    fillet = s.difference(p)
    p = Point([-tw / 2 - r, -h / 2 + tf + r]).buffer(r)
    s = Polygon(
        [
            (-tw / 2 - r, -h / 2 + tf),
            (-tw / 2, -h / 2 + tf),
            (-tw / 2, -h / 2 + tf + r),
            (-tw / 2 - r, -h / 2 + tf + r),
        ]
    )
    fillet = s.difference(p).union(fillet)
    fillet = translate(
        scale(fillet, 1, -1), xoff=0, yoff=h - 2 * tf - r
    ).union(fillet)
    # Create the geometry
    geometries = [
        set_precision(geometry, grid_size=1e-13)
        for geometry in [fillet, top_flange, bottom_flange, web]
    ]
    return unary_union(geometries)


@pytest.mark.parametrize(
    'h, b, tw, tf, r',
    [
        (80, 46, 3.8, 5.2, 5),
        (100, 55, 4.1, 5.7, 7),
        (120, 64, 4.4, 6.3, 7),
        (140, 73, 4.7, 6.9, 7),
        (160, 82, 5, 7.4, 9),
        (180, 91, 5.3, 8, 9),
        (200, 100, 5.6, 8.5, 12),
        (220, 110, 5.9, 9.2, 12),
        (240, 120, 6.2, 9.8, 15),
        (270, 135, 6.6, 10.2, 15),
        (300, 150, 7.1, 10.7, 15),
        (330, 160, 7.5, 11.5, 18),
        (360, 170, 8, 12.7, 18),
        (400, 180, 8.6, 13.5, 21),
        (450, 190, 9.4, 14.6, 21),
        (500, 200, 10.2, 16, 21),
        (550, 210, 11.1, 17.2, 21),
        (600, 220, 12, 19, 24),
    ],
)
def test_Isection_elastic_fiber(h, b, tw, tf, r):
    """Test Steel I section elastic strength."""
    Es = 206000
    fy = 355
    eps_su = fy / Es
    steel = ElasticPlastic(E=Es, fy=fy, eps_su=eps_su)

    # Create geometry
    geo = CompoundGeometry(
        [SurfaceGeometry(_create_I_section(h, b, tw, tf, r), steel)]
    )

    # Compute expected values
    wy_el = _wyel_IPE(h, b, tw, tf, r)
    wz_el = _wzel_IPE(h, b, tw, tf, r)
    my_expected = wy_el * fy * 1e-6
    mz_expected = wz_el * fy * 1e-6
    # Create the section with fiber
    sec = GenericSection(geo, integrator='Fiber', mesh_size=0.001)
    results = sec.section_analyzer.calculate_bending_strength(theta=0, n=0)
    assert math.isclose(-results.m_y * 1e-6, my_expected, rel_tol=1e-3)
    results = sec.section_analyzer.calculate_bending_strength(
        theta=math.pi / 2, n=0
    )
    assert math.isclose(-results.m_z * 1e-6, mz_expected, rel_tol=1e-3)


@pytest.mark.parametrize(
    'h, b, tw, tf, r',
    [
        (80, 46, 3.8, 5.2, 5),
        (100, 55, 4.1, 5.7, 7),
        (120, 64, 4.4, 6.3, 7),
        (140, 73, 4.7, 6.9, 7),
        (160, 82, 5, 7.4, 9),
        (180, 91, 5.3, 8, 9),
        (200, 100, 5.6, 8.5, 12),
        (220, 110, 5.9, 9.2, 12),
        (240, 120, 6.2, 9.8, 15),
        (270, 135, 6.6, 10.2, 15),
        (300, 150, 7.1, 10.7, 15),
        (330, 160, 7.5, 11.5, 18),
        (360, 170, 8, 12.7, 18),
        (400, 180, 8.6, 13.5, 21),
        (450, 190, 9.4, 14.6, 21),
        (500, 200, 10.2, 16, 21),
        (550, 210, 11.1, 17.2, 21),
        (600, 220, 12, 19, 24),
    ],
)
def test_Isection_elastic_marin(h, b, tw, tf, r):
    """Test Steel I section elastic strength."""
    Es = 206000
    fy = 355
    eps_su = fy / Es
    steel = ElasticPlastic(E=Es, fy=fy, eps_su=eps_su)

    # Create geometry
    geo = CompoundGeometry(
        [SurfaceGeometry(_create_I_section(h, b, tw, tf, r), steel)]
    )

    # Compute expected values
    wy_el = _wyel_IPE(h, b, tw, tf, r)
    wz_el = _wzel_IPE(h, b, tw, tf, r)
    my_expected = wy_el * fy * 1e-6
    mz_expected = wz_el * fy * 1e-6
    # Create the section with Marin integrator
    sec = GenericSection(geo)
    results = sec.section_analyzer.calculate_bending_strength(theta=0, n=0)
    assert math.isclose(-results.m_y * 1e-6, my_expected, rel_tol=1e-3)
    results = sec.section_analyzer.calculate_bending_strength(
        theta=math.pi / 2, n=0
    )
    assert math.isclose(-results.m_z * 1e-6, mz_expected, rel_tol=1e-3)


@pytest.mark.parametrize(
    'h, b, tw, tf, r',
    [
        (80, 46, 3.8, 5.2, 5),
        (100, 55, 4.1, 5.7, 7),
        (120, 64, 4.4, 6.3, 7),
        (140, 73, 4.7, 6.9, 7),
        (160, 82, 5, 7.4, 9),
        (180, 91, 5.3, 8, 9),
        (200, 100, 5.6, 8.5, 12),
        (220, 110, 5.9, 9.2, 12),
        (240, 120, 6.2, 9.8, 15),
        (270, 135, 6.6, 10.2, 15),
        (300, 150, 7.1, 10.7, 15),
        (330, 160, 7.5, 11.5, 18),
        (360, 170, 8, 12.7, 18),
        (400, 180, 8.6, 13.5, 21),
        (450, 190, 9.4, 14.6, 21),
        (500, 200, 10.2, 16, 21),
        (550, 210, 11.1, 17.2, 21),
        (600, 220, 12, 19, 24),
    ],
)
def test_Isection_plastic_fiber(h, b, tw, tf, r):
    """Test Steel I section elastic strength."""
    Es = 206000
    fy = 355
    eps_su = 0.15
    steel = ElasticPlastic(E=Es, fy=fy, eps_su=eps_su)

    # Create geometry
    geo = CompoundGeometry(
        [SurfaceGeometry(_create_I_section(h, b, tw, tf, r), steel)]
    )

    # Compute expected values
    wy_pl = _wypl_IPE(h, b, tw, tf, r)
    wz_pl = _wzpl_IPE(h, b, tw, tf, r)
    my_expected = wy_pl * fy * 1e-6
    mz_expected = wz_pl * fy * 1e-6
    # Create the section with fiber
    sec = GenericSection(geo, integrator='Fiber', mesh_size=0.001)
    results = sec.section_analyzer.calculate_bending_strength(theta=0, n=0)
    assert math.isclose(-results.m_y * 1e-6, my_expected, rel_tol=1e-2)
    results = sec.section_analyzer.calculate_bending_strength(
        theta=math.pi / 2, n=0
    )
    assert math.isclose(-results.m_z * 1e-6, mz_expected, rel_tol=1e-2)


@pytest.mark.parametrize(
    'h, b, tw, tf, r',
    [
        (80, 46, 3.8, 5.2, 5),
        (100, 55, 4.1, 5.7, 7),
        (120, 64, 4.4, 6.3, 7),
        (140, 73, 4.7, 6.9, 7),
        (160, 82, 5, 7.4, 9),
        (180, 91, 5.3, 8, 9),
        (200, 100, 5.6, 8.5, 12),
        (220, 110, 5.9, 9.2, 12),
        (240, 120, 6.2, 9.8, 15),
        (270, 135, 6.6, 10.2, 15),
        (300, 150, 7.1, 10.7, 15),
        (330, 160, 7.5, 11.5, 18),
        (360, 170, 8, 12.7, 18),
        (400, 180, 8.6, 13.5, 21),
        (450, 190, 9.4, 14.6, 21),
        (500, 200, 10.2, 16, 21),
        (550, 210, 11.1, 17.2, 21),
        (600, 220, 12, 19, 24),
    ],
)
def test_Isection_plastic_marin(h, b, tw, tf, r):
    """Test Steel I section elastic strength."""
    Es = 206000
    fy = 355
    eps_su = 0.15
    steel = ElasticPlastic(E=Es, fy=fy, eps_su=eps_su)

    # Create geometry
    geo = CompoundGeometry(
        [SurfaceGeometry(_create_I_section(h, b, tw, tf, r), steel)]
    )

    # Compute expected values
    wy_pl = _wypl_IPE(h, b, tw, tf, r)
    wz_pl = _wzpl_IPE(h, b, tw, tf, r)
    my_expected = wy_pl * fy * 1e-6
    mz_expected = wz_pl * fy * 1e-6
    # Create the section with marin integrator
    sec = GenericSection(geo)
    results = sec.section_analyzer.calculate_bending_strength(theta=0, n=0)
    assert math.isclose(-results.m_y * 1e-6, my_expected, rel_tol=1e-2)
    results = sec.section_analyzer.calculate_bending_strength(
        theta=math.pi / 2, n=0
    )
    assert math.isclose(-results.m_z * 1e-6, mz_expected, rel_tol=1e-2)


@pytest.mark.parametrize(
    'h, b, tw, tf, r',
    [
        (80, 46, 3.8, 5.2, 5),
        (100, 55, 4.1, 5.7, 7),
        (120, 64, 4.4, 6.3, 7),
        (140, 73, 4.7, 6.9, 7),
        (160, 82, 5, 7.4, 9),
        (180, 91, 5.3, 8, 9),
        (200, 100, 5.6, 8.5, 12),
        (220, 110, 5.9, 9.2, 12),
        (240, 120, 6.2, 9.8, 15),
        (270, 135, 6.6, 10.2, 15),
        (300, 150, 7.1, 10.7, 15),
        (330, 160, 7.5, 11.5, 18),
        (360, 170, 8, 12.7, 18),
        (400, 180, 8.6, 13.5, 21),
        (450, 190, 9.4, 14.6, 21),
        (500, 200, 10.2, 16, 21),
        (550, 210, 11.1, 17.2, 21),
        (600, 220, 12, 19, 24),
    ],
)
def test_Isection_elastic_material_marin(h, b, tw, tf, r):
    """Test Steel I section elastic strength."""
    Es = 206000
    fy = 355
    steel = Elastic(E=Es)
    steel.set_ultimate_strain(fy / Es)
    # Create geometry
    geo = CompoundGeometry(
        [SurfaceGeometry(_create_I_section(h, b, tw, tf, r), steel)]
    )

    # Compute expected values
    wy_el = _wyel_IPE(h, b, tw, tf, r)
    wz_el = _wzel_IPE(h, b, tw, tf, r)
    my_expected = wy_el * fy * 1e-6
    mz_expected = wz_el * fy * 1e-6
    # Create the section with marin integrator
    sec = GenericSection(geo)
    results = sec.section_analyzer.calculate_bending_strength(theta=0, n=0)
    assert math.isclose(-results.m_y * 1e-6, my_expected, rel_tol=1e-3)
    results = sec.section_analyzer.calculate_bending_strength(
        theta=math.pi / 2, n=0
    )
    assert math.isclose(-results.m_z * 1e-6, mz_expected, rel_tol=1e-3)


@pytest.mark.parametrize(
    'h, b, tw, tf, r',
    [
        (80, 46, 3.8, 5.2, 5),
        (100, 55, 4.1, 5.7, 7),
        (120, 64, 4.4, 6.3, 7),
        (140, 73, 4.7, 6.9, 7),
        (160, 82, 5, 7.4, 9),
        (180, 91, 5.3, 8, 9),
        (200, 100, 5.6, 8.5, 12),
        (220, 110, 5.9, 9.2, 12),
        (240, 120, 6.2, 9.8, 15),
        (270, 135, 6.6, 10.2, 15),
        (300, 150, 7.1, 10.7, 15),
        (330, 160, 7.5, 11.5, 18),
        (360, 170, 8, 12.7, 18),
        (400, 180, 8.6, 13.5, 21),
        (450, 190, 9.4, 14.6, 21),
        (500, 200, 10.2, 16, 21),
        (550, 210, 11.1, 17.2, 21),
        (600, 220, 12, 19, 24),
    ],
)
def test_Isection_user_material_marin(h, b, tw, tf, r):
    """Test Steel I section elastic strength."""
    Es = 206000
    fy = 355
    eps_su = 7e-2
    steel = UserDefined(
        x=[-eps_su, -fy / Es, 0, fy / Es, eps_su], y=[-fy, -fy, 0, fy, fy]
    )
    # Create geometry
    geo = CompoundGeometry(
        [SurfaceGeometry(_create_I_section(h, b, tw, tf, r), steel)]
    )

    # Compute expected values
    wy_pl = _wypl_IPE(h, b, tw, tf, r)
    wz_pl = _wzpl_IPE(h, b, tw, tf, r)
    my_expected = wy_pl * fy * 1e-6
    mz_expected = wz_pl * fy * 1e-6
    # Create the section with fiber
    sec = GenericSection(geo)
    results = sec.section_analyzer.calculate_bending_strength(theta=0, n=0)
    assert math.isclose(-results.m_y * 1e-6, my_expected, rel_tol=1e-3)
    results = sec.section_analyzer.calculate_bending_strength(
        theta=math.pi / 2, n=0
    )
    assert math.isclose(-results.m_z * 1e-6, mz_expected, rel_tol=1e-2)
