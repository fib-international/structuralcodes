"""Tests for test_rc_utils."""

import pytest
from shapely import Polygon

from structuralcodes import set_design_code
from structuralcodes.geometry import SurfaceGeometry, add_reinforcement
from structuralcodes.materials.basic import ElasticMaterial
from structuralcodes.materials.concrete import create_concrete
from structuralcodes.materials.reinforcement import create_reinforcement
from structuralcodes.sections import GenericSection
from structuralcodes.sections._rc_utils import (
    calculate_elastic_cracked_properties,
)


def test_calculate_elastic_cracked_properties_comparison():
    """Test calculating cracked properties."""
    # region create materials
    set_design_code('ec2_2004')
    concrete = create_concrete(fck=45)
    reinforcement = create_reinforcement(
        fyk=500, Es=210000, ftk=550, epsuk=0.07
    )
    # endregion

    # region create section1
    width = 250
    height = 500
    polygon = Polygon(
        [
            (-width / 2, -height / 2),
            (width / 2, -height / 2),
            (width / 2, height / 2),
            (-width / 2, height / 2),
        ]
    )
    geometry = SurfaceGeometry(poly=polygon, material=concrete)
    diameter_reinf = 25
    cover = 50
    geometry = add_reinforcement(
        geometry,
        (
            -width / 2 + cover + diameter_reinf / 2,
            -height / 2 + cover + diameter_reinf / 2,
        ),
        diameter_reinf,
        reinforcement,
    )
    geometry = add_reinforcement(
        geometry,
        (
            width / 2 - cover - diameter_reinf / 2,
            -height / 2 + cover + diameter_reinf / 2,
        ),
        diameter_reinf,
        reinforcement,
    )

    # Create section
    section1 = GenericSection(geometry)
    # endregion

    # region Create cracked section to compare (section2)
    poly_cracked = polygon = Polygon(
        [
            (-width / 2, 119.98),
            (width / 2, 119.98),
            (width / 2, height / 2),
            (-width / 2, height / 2),
        ]
    )
    geometry_cracked = SurfaceGeometry(poly=poly_cracked, material=concrete)
    geometry_cracked = add_reinforcement(
        geometry_cracked,
        (
            -width / 2 + cover + diameter_reinf / 2,
            -height / 2 + cover + diameter_reinf / 2,
        ),
        diameter_reinf,
        reinforcement,
    )
    geometry_cracked = add_reinforcement(
        geometry_cracked,
        (
            width / 2 - cover - diameter_reinf / 2,
            -height / 2 + cover + diameter_reinf / 2,
        ),
        diameter_reinf,
        reinforcement,
    )

    section2 = GenericSection(geometry_cracked)
    # end region

    # Calculate cracked properties for both sections
    cracked_prop1 = calculate_elastic_cracked_properties(section1, theta=0)
    cracked_prop2 = section2.gross_properties

    # Compare specific properties between the two objects
    assert cracked_prop1.area == pytest.approx(cracked_prop2.area, rel=1e-3), (
        'Areas do not match'
    )
    assert cracked_prop1.cy == pytest.approx(cracked_prop2.cy, rel=1e-3), (
        'cy do not match'
    )
    assert cracked_prop1.cz == pytest.approx(cracked_prop2.cz, rel=1e-3), (
        'cz do not match'
    )
    assert cracked_prop1.e_i11 == pytest.approx(
        cracked_prop2.e_i11, rel=1e-3
    ), 'cze_i11 do not match'
    assert cracked_prop1.e_i22 == pytest.approx(
        cracked_prop2.e_i22, rel=1e-3
    ), 'e_i22 do not match'
    assert cracked_prop1.sy == pytest.approx(cracked_prop2.sy, rel=1e-3), (
        'sy do not match'
    )
    assert cracked_prop1.sz == pytest.approx(cracked_prop2.sz, rel=1e-3), (
        'sz do not match'
    )

    # check if works when geometry is required
    result = calculate_elastic_cracked_properties(
        section1, 0, return_cracked_section=True
    )
    assert isinstance(result, tuple)
    assert len(result) == 2

    # check if it works with multiples geometries
    geometry_dummy = SurfaceGeometry(poly=polygon, material=concrete)
    geometry += geometry_dummy
    section1 = GenericSection(geometry, integrator='fiber')
    result = calculate_elastic_cracked_properties(section1, 0)
    # check if it is not reinforced concrete
    section1.geometry._reinforced_concrete = False
    result = calculate_elastic_cracked_properties(section1)
    assert result is None


def test_calculate_elastic_cracked_properties_return_geometry():
    """Test calculating cracked properties when returning also geometry."""
    # region create materials
    set_design_code('ec2_2004')
    concrete = create_concrete(fck=45)
    reinforcement = create_reinforcement(
        fyk=500, Es=210000, ftk=550, epsuk=0.07
    )
    # endregion

    # region create section1
    width = 250
    height = 500
    polygon = Polygon(
        [
            (-width / 2, -height / 2),
            (width / 2, -height / 2),
            (width / 2, height / 2),
            (-width / 2, height / 2),
        ]
    )
    geometry = SurfaceGeometry(poly=polygon, material=concrete)
    diameter_reinf = 25
    cover = 50
    geometry = add_reinforcement(
        geometry,
        (
            -width / 2 + cover + diameter_reinf / 2,
            -height / 2 + cover + diameter_reinf / 2,
        ),
        diameter_reinf,
        reinforcement,
    )
    geometry = add_reinforcement(
        geometry,
        (
            width / 2 - cover - diameter_reinf / 2,
            -height / 2 + cover + diameter_reinf / 2,
        ),
        diameter_reinf,
        reinforcement,
    )

    # Create section
    section1 = GenericSection(geometry)
    # endregion

    # check if works when geometry is required
    result = calculate_elastic_cracked_properties(
        section1, 0, return_cracked_section=True
    )
    assert isinstance(result, tuple)
    assert len(result) == 2

    # check if it works with multiples geometries
    geometry_dummy = SurfaceGeometry(poly=polygon, material=concrete)
    geometry += geometry_dummy
    section1 = GenericSection(geometry, integrator='fiber')
    result = calculate_elastic_cracked_properties(section1, 0)
    # check if it is not reinforced concrete
    section1.geometry._reinforced_concrete = False
    result = calculate_elastic_cracked_properties(section1)
    assert result is None


def test_calculate_elastic_cracked_properties_no_reinforced_concrete():
    """Test calculating cracked properties when the section is not RC."""
    # Create an elastic material
    elastic = ElasticMaterial(E=30000, density=2400)

    # Create the geometry and section
    width = 250
    height = 500
    polygon = Polygon(
        [
            (-width / 2, -height / 2),
            (width / 2, -height / 2),
            (width / 2, height / 2),
            (-width / 2, height / 2),
        ]
    )
    geometry = SurfaceGeometry(poly=polygon, material=elastic)
    section1 = GenericSection(geometry, integrator='fiber')

    assert not section1.geometry.reinforced_concrete
    # Calculate cracked properties
    result = calculate_elastic_cracked_properties(section1)
    assert result is None
