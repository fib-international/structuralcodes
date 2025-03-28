"""Tests for the Generic Section."""

import math

import numpy as np
import pytest
from shapely import Polygon

from structuralcodes.codes.ec2_2004 import reinforcement_duct_props
from structuralcodes.geometry import (
    CircularGeometry,
    RectangularGeometry,
    SurfaceGeometry,
    add_reinforcement,
    add_reinforcement_circle,
    add_reinforcement_line,
)
from structuralcodes.materials.concrete import ConcreteEC2_2004, ConcreteMC2010
from structuralcodes.materials.constitutive_laws import Elastic, Sargin
from structuralcodes.materials.reinforcement import (
    ReinforcementEC2_2004,
    ReinforcementMC2010,
)
from structuralcodes.sections import GenericSection


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

    # Compute max / min axial load
    n_max_marin = sec.section_calculator.n_max
    n_min_marin = sec.section_calculator.n_min

    # Compute bending strength
    res_marin = sec.section_calculator.calculate_bending_strength(theta=0, n=0)

    # Use integrate_strain_response
    N, My, Mz = sec.section_calculator.integrate_strain_profile(
        (res_marin.eps_a, res_marin.chi_y, res_marin.chi_z)
    )

    assert math.isclose(N, res_marin.n)
    assert math.isclose(My, res_marin.m_y)
    assert math.isclose(Mz, res_marin.m_z)

    # Compute moment curvature
    res_mc_marin = sec.section_calculator.calculate_moment_curvature(
        theta=0, n=0
    )

    # Use fiber integration
    sec = GenericSection(geo, integrator='Fiber', mesh_size=0.0001)
    assert math.isclose(sec.gross_properties.area, 200 * 400)

    # Compute bending strength
    res_fiber = sec.section_calculator.calculate_bending_strength(theta=0, n=0)

    assert math.isclose(res_marin.m_y, res_fiber.m_y, rel_tol=2e-3)

    # Compute moment curvature
    res_mc_fiber = sec.section_calculator.calculate_moment_curvature(
        theta=0, n=0
    )

    assert math.isclose(
        res_mc_marin.m_y[-1], res_mc_fiber.m_y[-1], rel_tol=2e-3
    )

    n_min_fiber = sec.section_calculator.n_min
    n_max_fiber = sec.section_calculator.n_max

    # Calculate moment-curvature for a given array of curvatures
    res_mc_fiber_same_curvature = (
        sec.section_calculator.calculate_moment_curvature(
            theta=0, n=0, chi=res_mc_fiber.chi_y
        )
    )
    assert math.isclose(
        res_mc_fiber.m_y[-1], res_mc_fiber_same_curvature.m_y[-1]
    )

    # Calculate moment-curvature for a given array of curvatures, but with
    # significant axial compression. This should raise a ValueError since we
    # cannot find equilibrium.
    with pytest.raises(ValueError):
        sec.section_calculator.calculate_moment_curvature(
            theta=0, n=0.95 * n_min_fiber, chi=res_mc_fiber.chi_y
        )

    # check if axial limit forces are the same for marin and fiber
    assert math.isclose(n_min_marin, n_min_fiber, rel_tol=2e-2)
    assert math.isclose(n_max_marin, n_max_fiber, rel_tol=2e-2)

    # confirm assertion if too much axial load
    with pytest.raises(ValueError):
        sec.section_calculator.calculate_bending_strength(
            theta=0, n=n_min_marin * 1.5
        )

    with pytest.raises(ValueError):
        sec.section_calculator.calculate_bending_strength(
            theta=0, n=n_max_marin * 1.5
        )
    with pytest.raises(ValueError):
        sec.section_calculator.calculate_moment_curvature(
            theta=0, n=n_min_marin * 1.5
        )

    with pytest.raises(ValueError):
        sec.section_calculator.calculate_moment_curvature(
            theta=0, n=n_max_marin * 1.5
        )


# Test rectangular section tangent stiffness
@pytest.mark.parametrize(
    'b, h, E',
    [
        (250, 500, 30000),
        (300, 600, 10000),
        (400, 200, 20000),
        (600, 350, 30000),
    ],
)
@pytest.mark.parametrize('integrator', ['fiber', 'marin'])
def test_rectangular_section_tangent_stiffness(b, h, E, integrator):
    """Test stiffness matrix of elastic rectangular section."""
    # Create materials to use
    elastic = Elastic(E)

    # The section
    poly = Polygon(
        ((-b / 2, -h / 2), (b / 2, -h / 2), (b / 2, h / 2), (-b / 2, h / 2))
    )
    geo = SurfaceGeometry(poly, elastic)

    assert geo.polygon.centroid.coords[0][0] == 0
    assert geo.polygon.centroid.coords[0][1] == 0

    # Create the section with fiber integrator
    sec = GenericSection(geo, integrator=integrator, mesh_size=0.0001)
    assert sec.name == 'GenericSection'

    assert math.isclose(sec.gross_properties.area, b * h)

    # compute stiffness matrix
    stiffness = sec.section_calculator.integrate_strain_profile(
        [0, 0, 0], 'modulus'
    )
    assert stiffness.shape == (3, 3)
    stiffness /= E

    # compute stiffness matrix analitically
    # The first moment and product moment should be zero in this case
    stiffness_expected = np.zeros((3, 3))
    stiffness_expected[0, 0] = b * h
    stiffness_expected[1, 1] = 1 / 12.0 * b * h**3
    stiffness_expected[2, 2] = 1 / 12.0 * h * b**3
    assert np.allclose(
        stiffness,
        stiffness_expected,
        atol=stiffness.max() * 1e-6,
        rtol=1e-2,
    )


# Test rectangular RC section tangent stiffness initial
@pytest.mark.parametrize(
    'b, h, c, diameter, n_bottom, n_top',
    [
        (250, 500, 40, 20, 3, 2),
    ],
)
@pytest.mark.parametrize('integrator', ['fiber', 'marin'])
def test_rectangular_rc_section_initial_tangent_stiffness(
    b, h, c, diameter, n_bottom, n_top, integrator
):
    """Test stiffness matrix of RC rectangular section."""
    # Create materials to use
    concrete = ConcreteEC2_2004(30)
    Es = 200000
    steel = ReinforcementEC2_2004(
        450, Es=Es, ftk=450, epsuk=0.075, gamma_s=1.15, gamma_eps=0.9
    )

    # The section
    poly = Polygon(
        ((-b / 2, -h / 2), (b / 2, -h / 2), (b / 2, h / 2), (-b / 2, h / 2))
    )
    geo = SurfaceGeometry(poly, concrete)

    assert geo.polygon.centroid.coords[0][0] == 0
    assert geo.polygon.centroid.coords[0][1] == 0

    # Add reinforcement
    geo = add_reinforcement_line(
        geo,
        coords_i=(-b / 2 + c, -h / 2 + c),
        coords_j=(b / 2 - c, -h / 2 + c),
        diameter=diameter,
        material=steel,
        n=n_bottom,
    )
    geo = add_reinforcement_line(
        geo,
        coords_i=(-b / 2 + c, h / 2 - c),
        coords_j=(b / 2 - c, h / 2 - c),
        diameter=diameter,
        material=steel,
        n=n_top,
    )

    # Manual computations
    Ec = concrete.constitutive_law.get_tangent(0)
    EA_gross = Ec * b * h
    EIyy_gross = Ec * 1 / 12.0 * b * h**3
    EIzz_gross = Ec * 1 / 12.0 * h * b**3
    for p in geo.point_geometries:
        y, z = p.point.coords.xy
        As = p.area
        EA_gross += Es * As
        EIyy_gross += Es * As * z[0] ** 2
        EIzz_gross += Es * As * y[0] ** 2

    # Create the section with fiber integrator
    sec = GenericSection(geo, integrator=integrator, mesh_size=0.0001)

    # compute initial stiffness matrix (gross)
    stiffness = sec.section_calculator.integrate_strain_profile(
        [0, 0, 0], 'modulus'
    )
    assert stiffness.shape == (3, 3)

    assert math.isclose(EA_gross, stiffness[0, 0], rel_tol=1e-3)
    assert math.isclose(EIyy_gross, stiffness[1, 1], rel_tol=1e-3)
    assert math.isclose(EIzz_gross, stiffness[2, 2], rel_tol=1e-3)


# Test rectangular RC section tangent stiffness
@pytest.mark.parametrize(
    'b, h, c, diameter, n_bottom, n_top',
    [
        (250, 500, 40, 20, 3, 2),
    ],
)
@pytest.mark.parametrize('integrator', ['fiber', 'marin'])
def test_rectangular_rc_section_tangent_stiffness(
    b, h, c, diameter, n_bottom, n_top, integrator
):
    """Test stiffness matrix of RC rectangular section."""
    # Create materials to use
    concrete = ConcreteEC2_2004(30)
    Es = 200000
    steel = ReinforcementEC2_2004(
        450, Es=Es, ftk=450, epsuk=0.075, gamma_s=1.15, gamma_eps=0.9
    )

    # The section
    poly = Polygon(
        ((-b / 2, -h / 2), (b / 2, -h / 2), (b / 2, h / 2), (-b / 2, h / 2))
    )
    geo = SurfaceGeometry(poly, concrete)

    assert geo.polygon.centroid.coords[0][0] == 0
    assert geo.polygon.centroid.coords[0][1] == 0

    # Add reinforcement
    geo = add_reinforcement_line(
        geo,
        coords_i=(-b / 2 + c, -h / 2 + c),
        coords_j=(b / 2 - c, -h / 2 + c),
        diameter=diameter,
        material=steel,
        n=n_bottom,
    )
    geo = add_reinforcement_line(
        geo,
        coords_i=(-b / 2 + c, h / 2 - c),
        coords_j=(b / 2 - c, h / 2 - c),
        diameter=diameter,
        material=steel,
        n=n_top,
    )

    # Manual computations
    Ec = concrete.constitutive_law.get_tangent(0)
    As = diameter**2 / 4.0 * np.pi

    # Position of neutral axis
    n = Es / Ec
    As_bottom = n_bottom * As
    As_top = n_top * As
    d = h - c
    x = n * (As_bottom + As_top) / b
    x *= (
        -1
        + (
            1
            + (2 * b * (As_bottom * d + As_top * c))
            / (n * (As_bottom + As_top) ** 2)
        )
        ** 0.5
    )

    EA = Ec * b * x
    EIyy = Ec * 1 / 12.0 * b * x**3 + Ec * b * x * (x / 2) ** 2
    EIzz = Ec * 1 / 12.0 * x * b**3
    for p in geo.point_geometries:
        y, z = p.point.coords.xy
        As = p.area
        EA += Es * As
        EIyy += Es * As * (z[0] - h / 2 + x) ** 2
        EIzz += Es * As * y[0] ** 2
    EIyy += EA * (h / 2 - x) ** 2

    # Create the section with fiber integrator
    sec = GenericSection(geo, integrator=integrator, mesh_size=0.0001)

    # compute tangent stiffness matrix (cracked)
    # for the given position of n.a. and a very small curvature
    chi_y = -1e-9
    eps_a = -chi_y * (h / 2 - x)

    stiffness = sec.section_calculator.integrate_strain_profile(
        [eps_a, chi_y, 0], 'modulus'
    )

    assert stiffness.shape == (3, 3)

    assert math.isclose(EA, stiffness[0, 0], rel_tol=1e-3)
    assert math.isclose(EIyy, stiffness[1, 1], rel_tol=1e-3)
    assert math.isclose(EIzz, stiffness[2, 2], rel_tol=1e-3)

    # For comparison, let's compare this with a elastic section of only
    # reacting concrete:
    elastic = Elastic(Ec)
    geo = SurfaceGeometry(
        Polygon(
            (
                (-b / 2, h / 2 - x),
                (b / 2, h / 2 - x),
                (b / 2, h / 2),
                (-b / 2, h / 2),
            )
        ),
        elastic,
    )
    # Add reinforcement
    geo = add_reinforcement_line(
        geo,
        coords_i=(-b / 2 + c, -h / 2 + c),
        coords_j=(b / 2 - c, -h / 2 + c),
        diameter=diameter,
        material=steel,
        n=n_bottom,
    )
    geo = add_reinforcement_line(
        geo,
        coords_i=(-b / 2 + c, h / 2 - c),
        coords_j=(b / 2 - c, h / 2 - c),
        diameter=diameter,
        material=steel,
        n=n_top,
    )

    # create this effective elastic section
    sec = GenericSection(geo, integrator=integrator, mesh_size=0.0001)

    stiffness2 = sec.section_calculator.integrate_strain_profile(
        [0, 0, 0], 'modulus'
    )

    assert np.allclose(
        stiffness,
        stiffness2,
        atol=stiffness.max() * 1e-3,
        rtol=1e-2,
    )


# Test rectangular section tangent stiffness
# respect to bottm left point
@pytest.mark.parametrize(
    'b, h, E',
    [
        (250, 500, 30000),
        (300, 600, 10000),
        (400, 200, 20000),
        (600, 350, 30000),
    ],
)
@pytest.mark.parametrize('integrator', ['fiber', 'marin'])
def test_rectangular_section_tangent_stiffness_translated(b, h, E, integrator):
    """Test stiffness matrix of elastic rectangular section."""
    # Create materials to use
    elastic = Elastic(E)

    # The section
    poly = Polygon(((0, 0), (b, 0), (b, h), (0, h)))
    geo = SurfaceGeometry(poly, elastic)

    assert geo.polygon.centroid.coords[0][0] == b / 2
    assert geo.polygon.centroid.coords[0][1] == h / 2

    # Create the section with fiber integrator
    sec = GenericSection(geo, integrator=integrator, mesh_size=0.0001)
    assert sec.name == 'GenericSection'

    assert math.isclose(sec.gross_properties.area, b * h)

    # compute stiffness matrix
    stiffness, _ = (
        sec.section_calculator.integrator.integrate_strain_response_on_geometry(
            sec.geometry,
            [0, 0, 0],
            integrate='modulus',
            mesh_size=sec.section_calculator.mesh_size,
        )
    )
    assert stiffness.shape == (3, 3)
    stiffness /= E

    # compute stiffness matrix analitically
    stiffness_expected = np.zeros((3, 3))
    stiffness_expected[0, 0] = b * h
    stiffness_expected[1, 1] = 1 / 3.0 * b * h**3
    stiffness_expected[2, 2] = 1 / 3.0 * h * b**3
    stiffness_expected[0, 1] = stiffness_expected[1, 0] = b * h**2 / 2.0
    stiffness_expected[0, 2] = stiffness_expected[2, 0] = -h * b**2 / 2.0
    stiffness_expected[1, 2] = stiffness_expected[2, 1] = -(h**2) * b**2 / 4.0
    assert np.allclose(
        stiffness,
        stiffness_expected,
        rtol=1e-2,
    )


@pytest.mark.parametrize('fck', [25, 35, 55, 65])
@pytest.mark.parametrize('fyk', [450, 500, 550])
@pytest.mark.parametrize(
    'ductility_class',
    ['a', 'b', 'c'],
)
def test_rectangular_section_parabola_rectangle(fck, fyk, ductility_class):
    """Test a rectangular section with different concretes and n."""
    # crete the materials to use
    concrete = ConcreteEC2_2004(fck=fck)
    props = reinforcement_duct_props(fyk=fyk, ductility_class=ductility_class)

    steel = ReinforcementEC2_2004(
        fyk=fyk, Es=200000, ftk=props['ftk'], epsuk=props['epsuk']
    )

    # The section
    poly = Polygon(((0, 0), (200, 0), (200, 400), (0, 400)))
    geo = SurfaceGeometry(poly, concrete)
    geo = add_reinforcement_line(geo, (40, 40), (160, 40), 16, steel, n=4)
    geo = add_reinforcement_line(geo, (40, 360), (160, 360), 16, steel, n=4)
    geo = geo.translate(-100, -200)

    # Create the section with fiber integrator
    sec_fiber = GenericSection(geo, integrator='fiber', mesh_size=0.001)

    # Compute bending strength My-
    res_fiber = sec_fiber.section_calculator.calculate_bending_strength()

    # Create the section with default marin integrator
    sec_marin = GenericSection(geo)

    # Compute bending strength My-
    res_marin = sec_marin.section_calculator.calculate_bending_strength()

    assert math.isclose(res_fiber.m_y, res_marin.m_y, rel_tol=1e-3)


def test_rectangular_section_mn_domain():
    """Test rectangular section interaction domain."""
    # Create materials to use
    concrete = ConcreteMC2010(25)
    steel = ReinforcementMC2010(fyk=450, Es=210000, ftk=450, epsuk=0.0675)

    # The section
    poly = Polygon(((0, 0), (200, 0), (200, 400), (0, 400)))
    geo = SurfaceGeometry(poly, concrete)
    geo = add_reinforcement_line(geo, (40, 40), (160, 40), 16, steel, n=4)
    geo = add_reinforcement_line(geo, (40, 360), (160, 360), 16, steel, n=4)
    geo = geo.translate(-100, -200)

    # Create the section (default Marin integrator)
    sec_marin = GenericSection(geo)

    # Compute MN domain
    mn_res_marin = (
        sec_marin.section_calculator.calculate_nm_interaction_domain(theta=0)
    )

    # Use fiber integration
    sec_fiber = GenericSection(geo, integrator='Fiber', mesh_size=0.0001)

    # compute MN domain
    # Fast version
    mn_res_fiber = (
        sec_fiber.section_calculator.calculate_nm_interaction_domain(theta=0)
    )

    assert math.isclose(
        mn_res_marin.m_y.flat[np.abs(mn_res_marin.m_y).argmax()],
        mn_res_fiber.m_y.flat[np.abs(mn_res_fiber.m_y).argmax()],
        rel_tol=2e-2,
    )


def test_rectangular_section_mm_domain():
    """Test rectangular section MM interaction domain."""
    # Create materials to use
    concrete = ConcreteMC2010(25)
    steel = ReinforcementMC2010(fyk=450, Es=210000, ftk=450, epsuk=0.075)

    # The section
    poly = Polygon(((0, 0), (200, 0), (200, 400), (0, 400)))
    geo = SurfaceGeometry(poly, concrete)
    geo = add_reinforcement_line(geo, (40, 40), (160, 40), 16, steel, n=4)
    geo = add_reinforcement_line(geo, (40, 360), (160, 360), 16, steel, n=4)
    geo = geo.translate(-100, -200)

    # Create the section (default Marin integrator)
    sec_marin = GenericSection(geo)
    # Compute MN domain
    mm_res_marin = (
        sec_marin.section_calculator.calculate_mm_interaction_domain(n=0)
    )

    # Use fiber integration
    sec_fiber = GenericSection(geo, integrator='Fiber', mesh_size=0.0001)
    # compute MN domain
    mm_res_fiber = (
        sec_fiber.section_calculator.calculate_mm_interaction_domain(n=0)
    )

    assert math.isclose(
        mm_res_marin.m_y.max(),
        mm_res_fiber.m_y.max(),
        rel_tol=2e-2,
    )
    assert math.isclose(
        mm_res_marin.m_z.max(),
        mm_res_fiber.m_z.max(),
        rel_tol=2e-2,
    )


def test_rectangular_section_nmm_domain():
    """Test rectangular section NMM interaction domain."""
    # Create materials to use
    concrete = ConcreteMC2010(25)
    steel = ReinforcementMC2010(fyk=450, Es=210000, ftk=450, epsuk=0.075)

    # The section
    poly = Polygon(((0, 0), (200, 0), (200, 400), (0, 400)))
    geo = SurfaceGeometry(poly, concrete)
    geo = add_reinforcement_line(geo, (40, 40), (160, 40), 16, steel, n=4)
    geo = add_reinforcement_line(geo, (40, 360), (160, 360), 16, steel, n=4)
    geo = geo.translate(-100, -200)

    # Create the section (default Marin integrator)
    sec_marin = GenericSection(geo)
    # Compute MN domain
    mm_res_marin = (
        sec_marin.section_calculator.calculate_nmm_interaction_domain()
    )

    # Use fiber integration
    sec_fiber = GenericSection(geo, integrator='Fiber', mesh_size=0.0001)
    # compute MN domain
    mm_res_fiber = (
        sec_fiber.section_calculator.calculate_nmm_interaction_domain()
    )

    assert math.isclose(
        mm_res_fiber.forces.max(),
        mm_res_fiber.forces.max(),
        rel_tol=2e-2,
    )
    assert math.isclose(
        mm_res_marin.forces.max(),
        mm_res_fiber.forces.max(),
        rel_tol=2e-2,
    )


# Test rectangular section with Sargin Model
def test_rectangular_section_Sargin():
    """Test rectangular section."""
    # Set a different constitutive law respect to default Parabola-Rectangle
    # Here we use Sargin law (MC2010 eq 5.1-26) with parameters taken from
    # MC2010 table 5.1-8
    constitutive_law = Sargin(fc=-35, eps_c1=-2.3e-3, eps_cu1=-3.5e-3, k=1.92)
    # Create materials to use
    concrete = ConcreteMC2010(25, constitutive_law=constitutive_law)
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

    assert math.isclose(sec.gross_properties.area, 200 * 400)

    # Compute bending strength
    res_marin = sec.section_calculator.calculate_bending_strength(theta=0, n=0)

    # Use fiber integration
    sec = GenericSection(geo, integrator='Fiber', mesh_size=0.0001)
    assert math.isclose(sec.gross_properties.area, 200 * 400)

    # Compute bending strength
    res_fiber = sec.section_calculator.calculate_bending_strength(theta=0, n=0)

    assert math.isclose(res_marin.m_y, res_fiber.m_y, rel_tol=1e-2)


# Test holed section
def test_holed_section():
    """Test a section with a hole."""
    # Create materials to use
    concrete = ConcreteMC2010(25, gamma_c=1.0, alpha_cc=1.0)
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

    # Compute bending strength with Marin
    res_marin = sec.section_calculator.calculate_bending_strength(theta=0, n=0)

    # Use fiber integration
    sec = GenericSection(geo, integrator='Fiber', mesh_size=0.0001)
    assert math.isclose(sec.gross_properties.area, 260000)

    # Compute bending strength
    res_fiber = sec.section_calculator.calculate_bending_strength(theta=0, n=0)

    assert math.isclose(res_marin.m_y, res_fiber.m_y, rel_tol=1e-3)


# Test U section
def test_u_section():
    """Test a section with a U shape."""
    # Create materials to use
    concrete = ConcreteMC2010(25, gamma_c=1.0)
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
    res_marin = sec.section_calculator.calculate_bending_strength(theta=0, n=0)

    # Use fiber integration
    sec = GenericSection(geo, integrator='Fiber', mesh_size=0.0001)
    assert math.isclose(sec.gross_properties.area, 230000)

    # Compute bending strength
    res_fiber = sec.section_calculator.calculate_bending_strength(theta=0, n=0)

    assert math.isclose(res_marin.m_y, res_fiber.m_y, rel_tol=1e-2)


def test_refined_moment_curvature():
    """Test overriding defaults when calculating moment-curvature."""
    # Create materials to use
    concrete = ConcreteMC2010(25)
    steel = ReinforcementMC2010(fyk=450, Es=210000, ftk=450, epsuk=0.0675)

    # The section
    poly = Polygon(((0, 0), (200, 0), (200, 400), (0, 400)))
    geo = SurfaceGeometry(poly, concrete)
    geo = add_reinforcement_line(geo, (40, 40), (160, 40), 16, steel, n=4)
    geo = add_reinforcement_line(geo, (40, 360), (160, 360), 16, steel, n=4)
    geo = geo.translate(-100, -200)

    # Create the section
    sec = GenericSection(geo, integrator='fiber')

    # Calculate default moment-curvature relation
    res_default = sec.section_calculator.calculate_moment_curvature()

    # Get the yield curvature. Since we are using default values, this is found
    # at place 9 in the chi_y array
    chi_yield = res_default.chi_y[9]

    # Calculate moment-curvature relation with more points
    res_more_points = sec.section_calculator.calculate_moment_curvature(
        num_pre_yield=15, num_post_yield=33
    )

    # Calculate moment-curvature relation when providing a first curvature
    # larger than the yield curvature
    res_large_chi_first = sec.section_calculator.calculate_moment_curvature(
        chi_first=1.5 * chi_yield
    )

    # Assert
    assert len(res_more_points.chi_y) == 48
    assert math.isclose(res_large_chi_first.chi_y[0], chi_yield / 10)


def test_refined_mn_domain():
    """Test overriding defaults when calculating moment-curvature."""
    # Create materials to use
    concrete = ConcreteMC2010(25)
    steel = ReinforcementMC2010(fyk=450, Es=210000, ftk=450, epsuk=0.0675)

    # The section
    poly = Polygon(((0, 0), (200, 0), (200, 400), (0, 400)))
    geo = SurfaceGeometry(poly, concrete)
    geo = add_reinforcement_line(geo, (40, 40), (160, 40), 16, steel, n=4)
    geo = add_reinforcement_line(geo, (40, 360), (160, 360), 16, steel, n=4)
    geo = geo.translate(-100, -200)

    # Create the section
    sec = GenericSection(geo, integrator='fiber')

    # Calculate default moment-curvature relation
    res_default = sec.section_calculator.calculate_nm_interaction_domain()

    # Specify more discretized domain
    res_refined = sec.section_calculator.calculate_nm_interaction_domain(
        num=100
    )

    # Specify detailed discretization for each field
    res_detailed = sec.section_calculator.calculate_nm_interaction_domain(
        num_1=2,
        num_2=10,
        num_3=40,
        num_4=40,
        num_5=10,
        num_6=5,
        type_4='geometric',
    )

    # Assertion
    assert len(res_default.strains) == 35
    assert len(res_refined.strains) >= 0.9 * 100
    assert len(res_detailed.strains) == 107


@pytest.mark.parametrize(
    'theta',
    [
        0,
        np.pi / 2,
        np.pi / 4,
    ],
)
def test_rectangular_section_biaxial_moment(theta):
    """Test for rectangular section under biaxial moment."""
    # Create materials to use
    concrete = ConcreteMC2010(25)
    steel = ReinforcementMC2010(fyk=450, Es=210000, ftk=450, epsuk=0.0675)

    # The section
    poly = Polygon(((0, 0), (200, 0), (200, 400), (0, 400)))
    geo = SurfaceGeometry(poly, concrete)
    geo = add_reinforcement_line(geo, (40, 40), (160, 40), 16, steel, n=4)
    geo = add_reinforcement_line(geo, (40, 360), (160, 360), 16, steel, n=4)
    geo = geo.translate(-100, -200)

    # Create the section
    sec = GenericSection(geo, integrator='fiber')

    # Calculate default moment-curvature relation
    res = sec.section_calculator.calculate_bending_strength(theta=theta)

    theta_inverse = math.atan2(res.chi_z, res.chi_y) + np.pi

    theta_inverse = theta_inverse % (2 * np.pi)

    assert math.isclose(theta, theta_inverse, rel_tol=1e-3)


n = [(x) for x in np.linspace(-100e3, 100e3, 6)]
my = [(x) for x in np.linspace(-100e3, 100e3, 6)]
mz = [(x) for x in np.linspace(-100e3, 100e3, 6)]


@pytest.mark.parametrize('n', n)
@pytest.mark.parametrize('my', my)
@pytest.mark.parametrize('mz', mz)
@pytest.mark.parametrize(
    'Ec, b, h',
    [
        (30000, 400, 400),
        (30000, 200, 400),
        (20000, 200, 400),
        (20000, 400, 400),
        (20000, 400, 200),
    ],
)
def test_strain_plane_calculation_elastic_Nmm(n, my, mz, Ec, b, h):
    """Test deeply calculate strain profile method.

    Elastic materials, test many load combinations.
    """
    # Create materials to use
    concrete = Elastic(Ec)

    # Create the section
    geom = SurfaceGeometry(
        Polygon(
            (
                (-b / 2, -h / 2),
                (b / 2, -h / 2),
                (b / 2, h / 2),
                (-b / 2, h / 2),
            )
        ),
        concrete,
    )

    # Fiber
    section = GenericSection(geom, integrator='fiber', mesh_size=0.001)
    strain_fiber = section.section_calculator.calculate_strain_profile(
        n=n, my=my, mz=mz, tol=1e-7
    )

    # Marin
    section = GenericSection(geom)
    strain_marin = section.section_calculator.calculate_strain_profile(
        n=n, my=my, mz=mz, tol=1e-7
    )

    # Evaluate expected response
    # Ec = concrete.constitutive_law.get_tangent(0)
    EA = Ec * b * h
    EIyy = Ec * 1 / 12.0 * b * h**3
    EIzz = Ec * 1 / 12.0 * h * b**3

    eps_0 = n / EA
    chi_y = my / EIyy
    chi_z = mz / EIzz

    expected_strain = np.array([eps_0, chi_y, chi_z])

    # Compare
    assert np.allclose(
        np.array(strain_fiber), expected_strain, rtol=1e-2, atol=1e-8
    )
    assert np.allclose(
        np.array(strain_marin), expected_strain, rtol=1e-2, atol=1e-8
    )


n = [(x) for x in np.linspace(-100, 100, 6)]
my = [(x) for x in np.linspace(-100, 100, 6)]
mz = [(x) for x in np.linspace(-100, 100, 6)]


@pytest.mark.parametrize('n', n)
@pytest.mark.parametrize('my', my)
@pytest.mark.parametrize('mz', mz)
@pytest.mark.parametrize(
    'Ec, b, h',
    [
        (30000e3, 0.4, 0.4),
        (30000e3, 0.2, 0.4),
        (20000e3, 0.2, 0.4),
        (20000e3, 0.4, 0.4),
        (20000e3, 0.4, 0.2),
    ],
)
def test_strain_plane_calculation_elastic_kNm(n, my, mz, Ec, b, h):
    """Test deeply calculate strain profile method.

    Elastic materials, test many load combinations.

    Units in kN, m and kPa
    """
    # Create materials to use
    concrete = Elastic(Ec)

    # Create the section
    geom = SurfaceGeometry(
        Polygon(
            (
                (-b / 2, -h / 2),
                (b / 2, -h / 2),
                (b / 2, h / 2),
                (-b / 2, h / 2),
            )
        ),
        concrete,
    )

    # Fiber
    section = GenericSection(geom, integrator='fiber', mesh_size=0.001)
    strain_fiber = section.section_calculator.calculate_strain_profile(
        n=n, my=my, mz=mz, tol=1e-7
    )

    # Marin
    section = GenericSection(geom)
    strain_marin = section.section_calculator.calculate_strain_profile(
        n=n, my=my, mz=mz, tol=1e-7
    )

    # Evaluate expected response
    # Ec = concrete.constitutive_law.get_tangent(0)
    EA = Ec * b * h
    EIyy = Ec * 1 / 12.0 * b * h**3
    EIzz = Ec * 1 / 12.0 * h * b**3

    eps_0 = n / EA
    chi_y = my / EIyy
    chi_z = mz / EIzz

    expected_strain = np.array([eps_0, chi_y, chi_z])

    # Compare
    assert np.allclose(
        np.array(strain_fiber), expected_strain, rtol=5e-2, atol=1e-5
    )
    assert np.allclose(
        np.array(strain_marin), expected_strain, rtol=5e-2, atol=1e-5
    )


@pytest.mark.parametrize(
    'n, my, mz, fck, b, h',
    [
        (0.0, -50e6, 0, 30, 200, 400),
        (-100e3, -60e6, 0, 40, 250, 400),
        (-100e3, 50e6, 0, 40, 250, 400),
        (50e3, -30e6, 0, 30, 200, 400),
    ],
)
def test_strain_plane_calculation_rectangular_rc(n, my, mz, fck, b, h):
    """Test calculate strain profile method."""
    # Create materials to use
    concrete = ConcreteMC2010(fck)
    steel = ReinforcementMC2010(fyk=450, Es=200000, ftk=450, epsuk=0.075)

    # The section
    poly = Polygon(
        ((-b / 2, -h / 2), (b / 2, -h / 2), (b / 2, h / 2), (-b / 2, h / 2))
    )
    geo = SurfaceGeometry(poly, concrete)

    # Add reinforcement
    c = 40
    diameter = 16
    n_bottom = 3
    n_top = 2
    geo = add_reinforcement_line(
        geo,
        coords_i=(-b / 2 + c, -h / 2 + c),
        coords_j=(b / 2 - c, -h / 2 + c),
        diameter=diameter,
        material=steel,
        n=n_bottom,
    )
    geo = add_reinforcement_line(
        geo,
        coords_i=(-b / 2 + c, h / 2 - c),
        coords_j=(b / 2 - c, h / 2 - c),
        diameter=diameter,
        material=steel,
        n=n_top,
    )

    # Check with given loads that both marin and fiber gives same result
    section = GenericSection(geo)
    strain_marin = section.section_calculator.calculate_strain_profile(
        n, my, mz, tol=1e-7
    )
    strain_marin = np.array(strain_marin)

    section = GenericSection(geo, integrator='fiber', mesh_size=0.0001)

    strain_fiber = section.section_calculator.calculate_strain_profile(
        n, my, mz, tol=1e-7
    )
    strain_fiber = np.array(strain_fiber)

    # check that initial tangent gives the same solution at the end
    strain_fiber_initial = section.section_calculator.calculate_strain_profile(
        n, my, mz, tol=1e-7, initial=True, max_iter=80
    )
    strain_fiber_initial = np.array(strain_fiber_initial)

    assert np.allclose(strain_marin, strain_fiber, rtol=2e-2, atol=1e-6)
    assert np.allclose(
        strain_fiber, strain_fiber_initial, rtol=2e-2, atol=1e-6
    )


@pytest.mark.parametrize(
    'n, my, mz, fck, b, h',
    [
        (-100e3, -100e6, 0, 40, 250, 400),
    ],
)
def test_strain_plane_calculation_rectangular_rc_high_load(
    n, my, mz, fck, b, h
):
    """Test calculate strain profile method.

    check a case with high load (non convergence).
    """
    # Create materials to use
    concrete = ConcreteMC2010(fck)
    steel = ReinforcementMC2010(fyk=450, Es=200000, ftk=450, epsuk=0.075)

    # The section
    poly = Polygon(
        ((-b / 2, -h / 2), (b / 2, -h / 2), (b / 2, h / 2), (-b / 2, h / 2))
    )
    geo = SurfaceGeometry(poly, concrete)

    # Add reinforcement
    c = 40
    diameter = 16
    n_bottom = 3
    n_top = 2
    geo = add_reinforcement_line(
        geo,
        coords_i=(-b / 2 + c, -h / 2 + c),
        coords_j=(b / 2 - c, -h / 2 + c),
        diameter=diameter,
        material=steel,
        n=n_bottom,
    )
    geo = add_reinforcement_line(
        geo,
        coords_i=(-b / 2 + c, h / 2 - c),
        coords_j=(b / 2 - c, h / 2 - c),
        diameter=diameter,
        material=steel,
        n=n_top,
    )

    # Check that with given loads we don't reach convergence
    section = GenericSection(geo)
    with pytest.raises(
        StopIteration, match='Maximum number of iterations reached'
    ):
        section.section_calculator.calculate_strain_profile(
            n, my, mz, tol=1e-7
        )


@pytest.mark.parametrize(
    'w, h, fck, fyk, ftk, Es, eps_su, nbars, diameter, cover',
    [(250, 1000, 75.0, 500.0, 550, 200000, 7.0e-2, 3, 25, 50)],
)
@pytest.mark.parametrize('n', [(x) for x in np.linspace(-3000e3, 200e3, 10)])
def test_moment_curvature_slender(
    w, h, fck, fyk, ftk, Es, eps_su, nbars, diameter, cover, n
):
    """Test moment-curvature for slender section.

    This test did not pass with the current implementation of the
    moment-curvature calculation. Fixed with PR #224.
    """
    # Create the materials
    concrete = ConcreteMC2010(fck=fck)
    reinforcement = ReinforcementMC2010(fyk=fyk, Es=Es, ftk=ftk, epsuk=eps_su)

    # Create the geometry
    geometry = RectangularGeometry(width=w, height=h, material=concrete)

    # Add reinforcement
    geometry = add_reinforcement_line(
        geo=geometry,
        coords_i=(-w / 2 + cover, -h / 2 + cover),
        coords_j=(w / 2 - cover, -h / 2 + cover),
        diameter=diameter,
        material=reinforcement,
        n=nbars,
    )

    # Create fiber section
    section_fiber = GenericSection(geometry=geometry, integrator='fiber')

    # Create marin section
    section_marin = GenericSection(
        geometry=geometry,
    )

    # Calculation
    if (
        n >= section_fiber.section_calculator.n_max
        or n <= section_fiber.section_calculator.n_min
    ):
        return
    results_fiber = (
        section_fiber.section_calculator.calculate_moment_curvature(n=n)
    )

    results_marin = (
        section_marin.section_calculator.calculate_moment_curvature(n=n)
    )

    # check
    # Do not check last couple of points because with a not so fine fiber
    # discretization some differences can arise.
    # Also first point for some values of n has shows a diference of about
    # 7-8% (with this discretization). That is why the rtol is set to 10%
    assert np.allclose(
        results_fiber.m_y[:-2], results_marin.m_y[:-2], rtol=1e-1, atol=1e-1
    )
    assert np.allclose(
        results_fiber.chi_y[:-2],
        results_marin.chi_y[:-2],
        rtol=1e-1,
        atol=1e-4,
    )


def test_moment_curvature_large_circular_section():
    """Test moment-curvature for large hollow circular section.

    This test did not pass with the current implementation of the
    moment-curvature calculation. Fixed with PR #224.
    """
    # Geometry
    outer_diameter = 10000
    shaft_thickness = 500

    # Reinforcement
    diameter_reinf = 20
    cc_reinf = 200
    cover = 50

    # Materials
    fck = 45
    gamma_c = 1.5
    alpha_cc = 0.85
    fyk = 500
    Es = 200000

    # Modelling
    n_points_circle = 20

    # Materials
    concrete = ConcreteEC2_2004(fck=fck, gamma_c=gamma_c, alpha_cc=alpha_cc)
    reinforcement = ReinforcementEC2_2004(
        fyk=fyk,
        Es=Es,
        **reinforcement_duct_props(fyk=fyk, ductility_class='C'),
    )

    # Establish geometry
    shaft = CircularGeometry(
        diameter=outer_diameter,
        material=concrete,
        n_points=n_points_circle,
    ) - CircularGeometry(
        diameter=(outer_diameter - 2 * shaft_thickness),
        material=concrete,
        n_points=n_points_circle,
    )

    for radius in (
        0.5 * outer_diameter - cover,
        0.5 * outer_diameter - shaft_thickness + cover,
    ):
        shaft = add_reinforcement_circle(
            shaft,
            center=(0, 0),
            radius=radius,
            diameter=diameter_reinf,
            material=reinforcement,
            s=cc_reinf,
        )

    # Section
    section = GenericSection(geometry=shaft, integrator='fiber')

    # Moment-curvature
    result = section.section_calculator.calculate_moment_curvature(n=-1e6)

    # Compare
    exp_result_chi = np.array(
        [
            -1e-08,
            -3.7599839056175666e-08,
            -6.519967811235132e-08,
            -9.279951716852699e-08,
            -1.2039935622470266e-07,
            -1.4799919528087832e-07,
            -1.75599034337054e-07,
            -2.0319887339322965e-07,
            -2.307987124494053e-07,
            -2.5839855150558097e-07,
            -8.015207515413767e-07,
            -1.3446429515771726e-06,
            -1.8877651516129684e-06,
            -2.4308873516487644e-06,
            -2.97400955168456e-06,
            -3.517131751720356e-06,
            -4.060253951756152e-06,
            -4.603376151791947e-06,
            -5.146498351827743e-06,
            -5.689620551863539e-06,
        ]
    )

    exp_result_m = np.array(
        [
            -8.61365389e09,
            -2.27251866e10,
            -3.67794409e10,
            -5.08062579e10,
            -6.48065933e10,
            -7.87797566e10,
            -9.27250247e10,
            -1.06641639e11,
            -1.20528802e11,
            -1.34385677e11,
            -1.83835739e11,
            -1.91419082e11,
            -1.95481133e11,
            -1.98304138e11,
            -2.00740666e11,
            -2.02971492e11,
            -2.04990133e11,
            -2.06880606e11,
            -2.08682290e11,
            -2.10425206e11,
        ]
    )

    assert np.allclose(result.chi_y, exp_result_chi, rtol=1e-5, atol=1e-6)
    assert np.allclose(result.m_y, exp_result_m, rtol=1e-5, atol=1e-2)
