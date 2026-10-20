"""Test for results objects."""

import math

import numpy as np
import pytest

from structuralcodes.geometry import RectangularGeometry
from structuralcodes.materials.basic import ElasticMaterial
from structuralcodes.sections import BeamSection

width = 200
height = 400
E = 10000
eps_u = 0.01


@pytest.fixture
def elastic_rec_section() -> BeamSection:
    """Create a simple rectangular section with known properties."""
    # Square section centered at (0, 0)
    mat = ElasticMaterial(E=E, density=700, ultimate_strain=eps_u)

    geo = RectangularGeometry(width, height, mat)

    # It should not be concrete
    assert not geo.concrete

    return BeamSection(geo, integrator='marin')


# Test the string representation of gross properties
def test_gross_properties_text(elastic_rec_section):
    """A test for gross properties result string representation."""
    sec = elastic_rec_section

    str_res = sec.gross_properties.__str__()

    # Check the string starts with the expected title
    assert str_res.startswith('Section Properties:')

    # Check that there are 30+1 rows of properties
    assert len(str_res.split('\n')) == 31

    # Check that last item is principal axis angle
    last_entry = str_res.split('\n')[-2]
    assert last_entry.startswith('Principal axis angle (theta):')

    str_fmt = sec.gross_properties.__format__('.2f')
    # Check the string starts with the expected title
    assert str_fmt.startswith('Section Properties:')

    # Check that there are 30+1 rows of properties
    assert len(str_fmt.split('\n')) == 31

    # Check that last item is principal axis angle
    last_entry = str_fmt.split('\n')[-2]
    assert last_entry.startswith('Principal axis angle (theta):')

    last_val = last_entry.split(':')[-1].strip()
    assert last_val == '0.00'

    # Check the value of Iyy as an example meeting the format required
    Iyy = sec.gross_properties.iyy
    Iyy_str = format(Iyy, '.2f')
    for entry in str_fmt.split('\n'):
        if entry.startswith('Second moment (Iyy):'):
            val = entry.split(':')[-1].strip()
            assert val == Iyy_str
            break


def test_detailed_result_none_point(elastic_rec_section):
    """Test that point data in detailed result is empty."""
    sec = elastic_rec_section

    mult = sec.section_calculator.calculate_bending_strength()

    mult.create_detailed_result()

    # Now detailed result should not be None anymore
    assert mult.detailed_result

    # But point data should be None
    assert not mult.detailed_result.point_data


def test_detailed_result_bending_strength(elastic_rec_section):
    """Test detailed result for beinding strength."""
    sec = elastic_rec_section

    Iyy = sec.gross_properties.iyy

    mult = sec.section_calculator.calculate_bending_strength()

    # Check detailed_result data correspond with results
    assert mult.n == mult.detailed_result.n
    assert mult.m_y == mult.detailed_result.m_y
    assert mult.m_z == mult.detailed_result.m_z
    assert np.allclose(
        np.array([mult.eps_a, mult.chi_y, mult.chi_z]),
        mult.detailed_result.strain,
    )

    # Check the ultimate moment result
    fy = eps_u * E
    # - due to sign convention
    m_y_expected = -width * height**2 / 6 * fy

    assert math.isclose(m_y_expected, mult.m_y)

    # Check stress and strain in top and bottom fiber
    stress_top_expected = m_y_expected / Iyy * height / 2
    strain_top_expected = stress_top_expected / E

    stress_bottom_expected = m_y_expected / Iyy * (-height / 2)
    strain_bottom_expected = stress_bottom_expected / E

    stress_top = mult.get_point_stress(0, height / 2)
    stress_bottom = mult.get_point_stress(0, -height / 2)

    strain_top = mult.get_point_strain(0, height / 2)
    strain_bottom = mult.get_point_strain(0, -height / 2)

    assert math.isclose(stress_top_expected, stress_top)
    assert math.isclose(stress_bottom_expected, stress_bottom)
    assert math.isclose(strain_top_expected, strain_top)
    assert math.isclose(strain_bottom_expected, strain_bottom)


def test_detailed_result_moment_curvature(elastic_rec_section):
    """Test detailed result for moment curvature."""
    sec = elastic_rec_section

    mcurv = sec.section_calculator.calculate_moment_curvature()

    # Check detailed_result data correspond with results
    assert mcurv.n == mcurv.detailed_result.n

    for i in range(len(mcurv.chi_y)):
        mcurv.set_step(i)
        assert mcurv.chi_y[i] == mcurv.detailed_result.chi_y
        assert mcurv.chi_z[i] == mcurv.detailed_result.chi_z

        assert np.allclose(
            np.array([mcurv.eps_a[i], mcurv.chi_y[i], mcurv.chi_z[i]]),
            mcurv.detailed_result.strain,
        )


@pytest.mark.parametrize(
    'n, my',
    [
        (0, -30e6),  # pure bending
        (-10e3, -30e6),  # bending with compression
        (10e3, -30e6),  # bending with tension
    ],
)
def test_detailed_result_plane_strain(elastic_rec_section, n, my):
    """Test detailed result for plane strain calculation."""
    sec = elastic_rec_section

    Iyy = sec.gross_properties.iyy
    A = sec.gross_properties.area

    strain_res = sec.section_calculator.calculate_strain_profile(
        n=n, my=my, mz=0
    )

    assert math.isclose(strain_res.n_ext, n)
    assert math.isclose(strain_res.m_y_ext, my)
    assert math.isclose(strain_res.m_z_ext, 0)

    # Check detailed_result data correspond with results
    assert strain_res.n == strain_res.detailed_result.n
    assert strain_res.m_y == strain_res.detailed_result.m_y
    assert strain_res.m_z == strain_res.detailed_result.m_z

    stress_top_expected = n / A + my / Iyy * height / 2
    strain_top_expected = stress_top_expected / E

    stress_bottom_expected = n / A - my / Iyy * height / 2
    strain_bottom_expected = stress_bottom_expected / E

    stress_top = strain_res.get_point_stress(0, height / 2)
    stress_bottom = strain_res.get_point_stress(0, -height / 2)

    strain_top = strain_res.get_point_strain(0, height / 2)
    strain_bottom = strain_res.get_point_strain(0, -height / 2)

    assert math.isclose(stress_top_expected, stress_top)
    assert math.isclose(stress_bottom_expected, stress_bottom)
    assert math.isclose(strain_top_expected, strain_top)
    assert math.isclose(strain_bottom_expected, strain_bottom)


@pytest.mark.parametrize(
    'eps_a, chi_y',
    [
        (0, -1e-7),  # pure bending
        (-1e-4, -1e-7),  # bending with compression
        (1e-4, -2e-7),  # bending with tension
    ],
)
def test_detailed_result_integrate_plane_strain(
    elastic_rec_section, eps_a, chi_y
):
    """Test detailed result for integration of stresses."""
    sec = elastic_rec_section

    int_res = sec.section_calculator.integrate_strain_profile(
        strain=np.array([eps_a, chi_y, 0])
    )

    # Check detailed_result data correspond with results
    assert int_res.n == int_res.detailed_result.n
    assert int_res.m_y == int_res.detailed_result.m_y
    assert int_res.m_z == int_res.detailed_result.m_z

    strain_top_expected = eps_a + height / 2 * chi_y
    stress_top_expected = strain_top_expected * E

    strain_bottom_expected = eps_a - height / 2 * chi_y
    stress_bottom_expected = strain_bottom_expected * E

    stress_top = int_res.get_point_stress(0, height / 2)
    stress_bottom = int_res.get_point_stress(0, -height / 2)

    strain_top = int_res.get_point_strain(0, height / 2)
    strain_bottom = int_res.get_point_strain(0, -height / 2)

    assert math.isclose(stress_top_expected, stress_top)
    assert math.isclose(stress_bottom_expected, stress_bottom)
    assert math.isclose(strain_top_expected, strain_top)
    assert math.isclose(strain_bottom_expected, strain_bottom)
