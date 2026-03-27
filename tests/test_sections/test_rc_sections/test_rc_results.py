"""Testing section results for RC cross sections."""

import math

import numpy as np
import pytest

from structuralcodes.core import _section_results as s_res
from structuralcodes.geometry import (
    RectangularGeometry,
    add_reinforcement_line,
)
from structuralcodes.materials.concrete import ConcreteMC2010
from structuralcodes.materials.reinforcement import ReinforcementMC2010
from structuralcodes.sections import BeamSection


@pytest.fixture
def simple_rc_section():
    """Create a simple RC section with known properties."""
    # Square section centered at (0, 0)
    width = 400
    cover = 40
    diameter_reinf = 20

    concrete = ConcreteMC2010(fck=30)
    reinforcement = ReinforcementMC2010(
        fyk=500, Es=200000, ftk=500, epsuk=0.075
    )

    geo = RectangularGeometry(
        width, width, material=concrete, name='concrete', group_label='simple'
    )

    geo = add_reinforcement_line(
        geo,
        coords_i=(-width / 2 + cover, -width / 2 + cover),
        coords_j=(width / 2 - cover, -width / 2 + cover),
        diameter=diameter_reinf,
        n=4,
        material=reinforcement,
        group_label='bottom',
    )
    geo = add_reinforcement_line(
        geo,
        coords_i=(-width / 2 + cover, width / 2 - cover),
        coords_j=(width / 2 - cover, width / 2 - cover),
        diameter=diameter_reinf,
        n=4,
        material=reinforcement,
        group_label='top',
    )

    return BeamSection(geo, integrator='fiber')


## Test utility functions developed as non public


def test_matching_geometries_intersection_case_sensitive(simple_rc_section):
    """Test _matching_geometries with intersection and case sensitivity."""
    section = simple_rc_section

    # Case sensitive: "con*" matches "concrete" but "CON*" does not
    surfaces, points = s_res._matching_geometries(
        section.geometry, name='con*'
    )
    assert len(surfaces) == 1
    assert surfaces[0].name == 'concrete'
    assert len(points) == 0

    surfaces, _ = s_res._matching_geometries(section.geometry, name='CON*')
    assert len(surfaces) == 0

    # Case insensitive: "CON*" matches "concrete"
    surfaces, _ = s_res._matching_geometries(
        section.geometry, name='CON*', case_sensitive=False
    )
    assert len(surfaces) == 1
    assert surfaces[0].name == 'concrete'


def test_strain_from_kinematics_scalar_and_vectorized():
    """Test _strain_from_kinematics with both scalar and vector inputs."""
    # scalar
    eps = s_res._strain_from_kinematics(
        eps_a=1.0, chi_y=2.0, chi_z=3.0, y=4.0, z=5.0
    )
    assert eps == pytest.approx(1.0 - 3.0 * 4.0 + 2.0 * 5.0)

    # vectorized
    eps_a = np.array([0.0, 1.0, 2.0])
    chi_y = np.array([1.0, 1.0, 1.0])
    chi_z = np.array([0.0, 2.0, 4.0])
    out = s_res._strain_from_kinematics(eps_a, chi_y, chi_z, y=2.0, z=3.0)
    expected = eps_a - chi_z * 2.0 + chi_y * 3.0
    np.testing.assert_allclose(out, expected)


def test_point_inside_surface_geometry(simple_rc_section):
    """Test _point_inside_geometry_surface."""
    section = simple_rc_section
    surface = section.geometry.geometries[0]

    # Point inside
    assert s_res._point_inside_geometry_surface(surface, y=0, z=0)

    # Point on boundary (should be considered inside)
    assert s_res._point_inside_geometry_surface(surface, y=200, z=0)

    # Point outside
    assert not s_res._point_inside_geometry_surface(surface, y=300, z=0)


def test_point_matches_geometry_point_isclose(simple_rc_section):
    """Test _point_matches_geometry_point."""
    section = simple_rc_section
    pg = section.geometry.point_geometries[0]

    y, z = pg.x, pg.y

    assert (
        s_res._point_matches_geometry_point(pg, y + 1e-12, z - 1e-12) is True
    )
    assert s_res._point_matches_geometry_point(pg, y + 0.5, z + 0.5) is False


def test_get_point_response_section_none():
    """Test that get point response returns None if section is None."""
    out = s_res._get_point_response(
        section=None,
        eps_a=0.0,
        chi_y=0.0,
        chi_z=0.0,
        y=0.0,
        z=0.0,
        response_type='strain',
    )
    assert out is None


def test_get_point_response(simple_rc_section):
    """Test that get point response returns the first matching geometry."""
    section = simple_rc_section

    # We know that the point (-160, -160) is in the bottom reinforcement line
    y, z = -160.0, -160.0

    # Get strain at this point
    strain = s_res._get_point_response(
        section=section,
        eps_a=0.001,
        chi_y=0.01,
        chi_z=0.02,
        y=y,
        z=z,
        response_type='strain',
    )
    assert isinstance(strain, float)
    assert strain == pytest.approx(0.001 - 0.02 * y + 0.01 * z)

    # Get stress at this point
    stress = s_res._get_point_response(
        section=section,
        eps_a=0.001,
        chi_y=0.01,
        chi_z=0.02,
        y=y,
        z=z,
        response_type='stress',
    )
    assert isinstance(stress, (int, float))
    # Because it should find only the surface geometry first, which is
    # concrete, and the strain is positive so giving zero stress.
    assert stress == 0

    # Get stress at the point for all geometries
    stresses = s_res._get_point_response(
        section=section,
        eps_a=0.001,
        chi_y=0.01,
        chi_z=0.02,
        y=y,
        z=z,
        response_type='stress',
        all_results=True,
    )
    # Now we should get a dictionary
    assert isinstance(stresses, dict)
    # with length 2 (one for each geometry)
    assert len(stresses) == 2
    assert stresses.get(('concrete', 'simple'), None) is not None
    # The second geometry should be the bottom reinforcement
    bottom_found = False
    for key in stresses:
        if key[1] == 'bottom':
            bottom_found = True
            break
    assert bottom_found is True

    # no hits => returns None even if all_results=True
    res2 = s_res._get_point_response(
        section=section,
        eps_a=0.001,
        chi_y=0.01,
        chi_z=0.02,
        y=999.0,
        z=999.0,
        response_type='strain',
        all_results=True,
    )
    assert res2 is None


def test_bending_strength_result_point(simple_rc_section):
    """Test that bending strength result is computed correctly."""
    section = simple_rc_section

    res = section.section_calculator.calculate_bending_strength()

    # We know that the point (-160, -160) is in the bottom reinforcement line
    y1, z1 = -160.0, -160.0
    stress_1 = s_res._get_point_response(
        section=section,
        eps_a=res.eps_a,
        chi_y=res.chi_y,
        chi_z=res.chi_z,
        y=y1,
        z=z1,
        response_type='stress',
        group_label='bottom',
    )
    assert isinstance(stress_1, (int, float))

    stress_2 = res.get_point_stress(y=y1, z=z1, group_label='bottom')
    assert isinstance(stress_2, (int, float))

    assert np.isclose(stress_1, stress_2)

    # We know that the point (160, 160) is in the top reinforcement line
    y2, z2 = 160.0, 160.0
    stress1 = s_res._get_point_response(
        section=section,
        eps_a=res.eps_a,
        chi_y=res.chi_y,
        chi_z=res.chi_z,
        y=y2,
        z=z2,
        response_type='stress',
        group_label='top',
    )
    assert isinstance(stress1, (int, float))

    stress2 = res.get_point_stress(y=y2, z=z2, group_label='top')
    assert isinstance(stress2, (int, float))

    assert np.isclose(stress1, stress2)


def test_moment_curvature_result_point(simple_rc_section):
    """Test that moment curvature result is computed correctly."""
    section = simple_rc_section

    res = section.section_calculator.calculate_moment_curvature()

    # We know that the point (-160, -160) is in the bottom reinforcement line
    y1, z1 = -160.0, -160.0
    stress = res.get_point_stress(y=y1, z=z1, group_label='bottom')
    assert isinstance(stress, np.ndarray)

    assert len(stress) == len(res.eps_a)

    # Check that at the end the stress on the bar is the yield strength of
    # reinforcement.
    assert np.isclose(
        stress[-1], section.geometry.point_geometries[0].material.fyd()
    )

    # We know that the point (160, 160) is in the top reinforcement line
    y2, z2 = 160.0, 160.0
    stress = res.get_point_stress(y=y2, z=z2, group_label='bottom')
    strain = res.get_point_strain(y=y2, z=z2, group_label='bottom')
    assert stress is None
    assert strain is None
    stress = res.get_point_stress(y=y2, z=z2, group_label='top')
    assert isinstance(stress, np.ndarray)

    assert len(stress) == len(res.eps_a)


def test_moment_curvature_result(simple_rc_section):
    """Test that moment curvature result is computed correctly."""
    section = simple_rc_section

    res = section.section_calculator.calculate_moment_curvature()

    stresses_bottom = res.get_point_stress(
        y=-160.0, z=-160.0, group_label='bottom'
    )
    stresses_top = res.get_point_stress(y=160.0, z=160.0, group_label='top')

    assert len(res.detailed_result.surface_data['stress']) == 1000
    assert len(res.detailed_result.point_data['stress']) == 8

    for i in range(len(res.chi_y)):
        stress_bottom = res.detailed_result.point_data['stress'][0]
        assert np.isclose(stress_bottom, stresses_bottom[i])

        stress_top = res.detailed_result.point_data['stress'][-1]
        assert np.isclose(stress_top, stresses_top[i])

        res.next_step()

    res.previous_step()
    stress_bottom = res.detailed_result.point_data['stress'][0]
    assert np.isclose(stress_bottom, stresses_bottom[-2])

    res.set_step(10)
    stress_bottom = res.detailed_result.point_data['stress'][0]
    assert np.isclose(stress_bottom, stresses_bottom[10])


def test_bending_strength_result(simple_rc_section):
    """Test that bending strength result is computed correctly."""
    section = simple_rc_section

    res = section.section_calculator.calculate_bending_strength()

    stress_bottom = res.get_point_stress(
        y=-160.0, z=-160.0, group_label='bottom'
    )
    stress_top = res.get_point_stress(y=160.0, z=160.0, group_label='top')

    assert len(res.detailed_result.surface_data['stress']) == 1000
    assert len(res.detailed_result.point_data['stress']) == 8

    assert np.isclose(
        stress_bottom, res.detailed_result.point_data['stress'][0]
    )
    assert np.isclose(
        stress_bottom, res.detailed_result.point_data['material'][0].fyd()
    )

    assert np.isclose(stress_top, res.detailed_result.point_data['stress'][-1])


@pytest.mark.parametrize(
    'max_iter, tol, initial',
    [
        (10, 1e-6, False),
        (20, 1e-8, False),
        (100, 1e-6, True),
    ],
)
def test_calc_strain_profile_results(
    simple_rc_section, max_iter, tol, initial
):
    """Test that the calculate_strain_profile result."""
    section = simple_rc_section

    # External forces
    n = 0
    my = -175 * 1e6
    mz = 0
    f_ext = np.array([n, my, mz])

    # Get the bending strength result
    strain_res = section.section_calculator.calculate_strain_profile(
        *f_ext, max_iter=max_iter, tol=tol, initial=initial
    )

    # Check from result that the residual = f_ext - f_int
    f_int = np.array([strain_res.n, strain_res.m_y, strain_res.m_z])
    assert np.allclose(
        strain_res.residual, f_ext - f_int, rtol=1e-7, atol=1e-7
    )
    assert math.isclose(
        strain_res.residual_norm, np.linalg.norm(f_ext - f_int), abs_tol=1e-7
    )

    # Check the the correct data are stored from convergence
    assert math.isclose(strain_res.tolerance, tol)
    assert strain_res.max_iter == max_iter
    assert strain_res.used_initial_tangent == initial

    # check that length of residual_history is the same as iterations
    assert len(strain_res.residual_history) == strain_res.iterations

    assert len(strain_res.response_history) == strain_res.iterations
    assert len(strain_res.strain_history) == strain_res.iterations
    assert len(strain_res.residual_norm_history) == strain_res.iterations

    # Check that length of delta_strain is iterations - 1
    assert len(strain_res.delta_strain_history) == strain_res.iterations - 1
    assert (
        len(strain_res.delta_strain_norm_history) == strain_res.iterations - 1
    )

    # Check point stress and strain at the bottom reinforcement line
    stress_bottom = strain_res.get_point_stress(
        y=-160.0, z=-160.0, group_label='bottom'
    )
    strain_bottom = strain_res.get_point_strain(
        y=-160.0, z=-160.0, group_label='bottom'
    )

    assert np.isclose(
        stress_bottom, strain_res.detailed_result.point_data['stress'][0]
    )
    assert np.isclose(
        stress_bottom,
        strain_res.detailed_result.point_data['material'][0].fyd(),
    )
    assert np.isclose(
        strain_bottom, strain_res.detailed_result.point_data['strain'][0]
    )
