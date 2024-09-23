"""Tests for the Generic Section."""

import math

import numpy as np
import pytest
from shapely import Polygon

from structuralcodes.codes.ec2_2004 import reinforcement_duct_props
from structuralcodes.geometry import SurfaceGeometry
from structuralcodes.materials.concrete import ConcreteEC2_2004, ConcreteMC2010
from structuralcodes.materials.constitutive_laws import Sargin
from structuralcodes.materials.reinforcement import (
    ReinforcementEC2_2004,
    ReinforcementMC2010,
)
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
    res_marin = sec.section_calculator.calculate_bending_strength(theta=0, n=0)

    # Compute moment curvature
    res_mc_marin = sec.section_calculator.calculate_moment_curvature(
        theta=0, n=0
    )

    n_min_marin = sec.section_calculator.n_min
    n_max_marin = sec.section_calculator.n_max

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


def test_calculate_strain_profile():
    # region CREATE THE SECTION
    concrete = ConcreteEC2_2004(fck=25)
    reinforcemnet = ReinforcementEC2_2004(
        fyk=500,
        Es=200000,
        density=7850,
        ftk=500,
        epsuk=0.07,
    )
    poly = Polygon(((0, 0), (350, 0), (350, 500), (0, 500)))
    geo = SurfaceGeometry(poly, concrete)
    geo = add_reinforcement_line(
        geo, (50, 450), (300, 450), 12, reinforcemnet, n=3
    )
    geo = add_reinforcement_line(
        geo, (50, 50), (300, 50), 20, reinforcemnet, n=6
    )
    sec = GenericSection(
        geo,
    )
    sec.geometry = sec.geometry.translate(-175, -250)
    # endregion

    n = -100 * 1e3
    my = 50 * 1e6
    mz = -40 * 1e6
    strain_profile = sec.section_calculator.calculate_strain_profile(n, my, mz)
    assert strain_profile is not None
    res = sec.section_calculator.integrate_strain_profile(strain_profile)
    assert math.isclose(n, res[0], rel_tol=1e-3)
    assert math.isclose(my, res[1], rel_tol=1e-3)
    assert math.isclose(mz, res[2], rel_tol=1e-3)

    # confirm assertion if too much load
    with pytest.raises(ValueError):
        n = 0 * 1e3
        my = 500 * 1e6
        mz = 0 * 1e6
        sec.section_calculator.calculate_strain_profile(n, my, mz)
