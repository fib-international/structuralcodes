"""Tests for L profiles."""

import json
import math
from pathlib import Path

import pytest
from shapely.testing import assert_geometries_equal

from structuralcodes.geometry import (
    SurfaceGeometry,
)
from structuralcodes.geometry.profiles import (
    L,
)
from structuralcodes.materials.basic import (
    ElasticMaterial,
)
from structuralcodes.sections._generic import GenericSection


def load_L_profiles_data():
    """Load L profiles data from l.json file."""
    json_file = Path(__file__).parent / 'l.json'
    with open(json_file, 'r') as f:
        profiles_data = []
        for line in f:
            profile_data = json.loads(line.strip())
            profile_name = profile_data['ProfileName']
            h_mm = profile_data['h_mm']
            b_mm = profile_data['b_mm']
            t_mm = profile_data['t_mm']
            r1_mm = profile_data['r1_mm']
            r2_mm = profile_data['r2_mm']
            A_cm2 = profile_data['A_cm2']
            Iy_cm4 = profile_data['Iy_cm4']
            Iz_cm4 = profile_data['Iz_cm4']
            Wely_cm3 = profile_data['Wely_cm3']
            Welz_cm3 = profile_data['Welz_cm3']
            zs_cm = profile_data['zs_cm']
            ys_cm = profile_data['ys_cm']

            profiles_data.append(
                (
                    profile_name,
                    h_mm,
                    b_mm,
                    t_mm,
                    r1_mm,
                    r2_mm,
                    A_cm2,
                    Iy_cm4,
                    Iz_cm4,
                    Wely_cm3,
                    Welz_cm3,
                    zs_cm,
                    ys_cm,
                )
            )
    return profiles_data


@pytest.mark.parametrize(
    (
        'profile_name, h_mm, b_mm, '
        't_mm, r1_mm, r2_mm, '
        'A_cm2, Iy_cm4, Iz_cm4, '
        'Wely_cm3, Welz_cm3, zs_cm, ys_cm'
    ),
    load_L_profiles_data(),
)
def test_L_geometric_data(
    profile_name,
    h_mm,
    b_mm,
    t_mm,
    r1_mm,
    r2_mm,
    A_cm2,
    Iy_cm4,
    Iz_cm4,
    Wely_cm3,
    Welz_cm3,
    zs_cm,
    ys_cm,
):
    """Test L profile geometric data matches JSON data."""
    # Delete not used variables to avoid linting errors
    del (
        A_cm2,
        Iy_cm4,
        Iz_cm4,
        Wely_cm3,
        Welz_cm3,
        zs_cm,
        ys_cm,
    )

    # Create the profile
    profile = L(profile_name)

    # Check geometric properties
    assert math.isclose(profile.h, h_mm, rel_tol=1e-3)
    assert math.isclose(profile.b, b_mm, rel_tol=1e-3)
    assert math.isclose(profile.t, t_mm, rel_tol=1e-3)
    assert math.isclose(profile.r1, r1_mm, rel_tol=1e-3)
    assert math.isclose(profile.r2, r2_mm, rel_tol=1e-3)


@pytest.mark.parametrize(
    (
        'profile_name, h_mm, b_mm, '
        't_mm, r1_mm, r2_mm, '
        'A_cm2, Iy_cm4, Iz_cm4, '
        'Wely_cm3, Welz_cm3, zs_cm, ys_cm'
    ),
    load_L_profiles_data(),
)
def test_L_massprop_data(
    profile_name,
    h_mm,
    b_mm,
    t_mm,
    r1_mm,
    r2_mm,
    A_cm2,
    Iy_cm4,
    Iz_cm4,
    Wely_cm3,
    Welz_cm3,
    zs_cm,
    ys_cm,
):
    """Test L profile mass property matches JSON data."""
    # Delete not used variables to avoid linting errors
    del (
        h_mm,
        b_mm,
        t_mm,
        r1_mm,
        r2_mm,
        zs_cm,
        ys_cm,
    )
    # Create the profile
    profile = L(profile_name)

    # Check mass properties
    assert math.isclose(profile.A, A_cm2 * 1e2, rel_tol=2e-2)
    assert math.isclose(profile.Iy, Iy_cm4 * 1e4, rel_tol=1.5e-2)
    assert math.isclose(profile.Iz, Iz_cm4 * 1e4, rel_tol=1e-2)
    assert math.isclose(profile.Wely, Wely_cm3 * 1e3, rel_tol=2e-2)
    assert math.isclose(profile.Welz, Welz_cm3 * 1e3, rel_tol=2.5e-2)


@pytest.mark.parametrize(
    (
        'profile_name, h_mm, b_mm, '
        't_mm, r1_mm, r2_mm, '
        'A_cm2, Iy_cm4, Iz_cm4, '
        'Wely_cm3, Welz_cm3, zs_cm, ys_cm'
    ),
    load_L_profiles_data(),
)
def test_L_polygon(
    profile_name,
    h_mm,
    b_mm,
    t_mm,
    r1_mm,
    r2_mm,
    A_cm2,
    Iy_cm4,
    Iz_cm4,
    Wely_cm3,
    Welz_cm3,
    zs_cm,
    ys_cm,
):
    """Test L profile polygon."""
    # Delete not used variables to avoid linting errors
    del (
        h_mm,
        b_mm,
        t_mm,
        r1_mm,
        r2_mm,
        A_cm2,
        Iy_cm4,
        Iz_cm4,
        Wely_cm3,
        Welz_cm3,
        zs_cm,
        ys_cm,
    )

    # Ceate the profile
    profile = L(profile_name)

    # Get the polygon from the object
    poly1 = profile.polygon
    # Get the polygon from class method passing the full name
    poly2 = L.get_polygon(profile_name)
    assert_geometries_equal(poly1, poly2)


@pytest.mark.parametrize('invalid_name', ['L40x10', '40x40x10', 'L400'])
def test_L_invalid_name(
    invalid_name,
):
    """Test L profile with invalid name raises ValueError."""
    # Invalid name for constructor
    with pytest.raises(ValueError):
        L(invalid_name)
    # Invalid name for class method
    with pytest.raises(ValueError):
        L.get_polygon(invalid_name)


@pytest.mark.parametrize(
    (
        'profile_name, h_mm, b_mm, '
        't_mm, r1_mm, r2_mm, '
        'A_cm2, Iy_cm4, Iz_cm4, '
        'Wely_cm3, Welz_cm3, zs_cm, ys_cm'
    ),
    load_L_profiles_data(),
)
def test_L_massprops_polygon(
    profile_name,
    h_mm,
    b_mm,
    t_mm,
    r1_mm,
    r2_mm,
    A_cm2,
    Iy_cm4,
    Iz_cm4,
    Wely_cm3,
    Welz_cm3,
    zs_cm,
    ys_cm,
):
    """Test L profile massprops from polygon."""
    # Delete not used variables to avoid linting errors
    del (
        h_mm,
        b_mm,
        t_mm,
        r1_mm,
        r2_mm,
    )

    # Material properties
    Es = 200000
    steel = ElasticMaterial(E=Es, density=7850)

    # Create geometry
    geo = SurfaceGeometry(L.get_polygon(profile_name), material=steel)
    # Create the section
    sec = GenericSection(geo)

    # Get values from section
    a_section = sec.gross_properties.area
    iy_section = sec.gross_properties.iyy_c
    iz_section = sec.gross_properties.izz_c

    bbox = geo.polygon.bounds
    zmin, zmax = bbox[1], bbox[3]
    ymin, ymax = bbox[0], bbox[2]
    zs_section = abs(zmin)
    ys_section = abs(ymin)
    wy_section = iy_section / zmax
    wz_section = iz_section / ymax

    # Assert
    assert math.isclose(zs_cm * 10, zs_section, rel_tol=1e-2)
    assert math.isclose(ys_cm * 10, ys_section, rel_tol=1e-2)
    assert math.isclose(A_cm2 * 1e2, a_section, rel_tol=2e-2)
    assert math.isclose(Iy_cm4 * 1e4, iy_section, rel_tol=1.5e-2)
    assert math.isclose(Iz_cm4 * 1e4, iz_section, rel_tol=1e-2)
    assert math.isclose(Wely_cm3 * 1e3, wy_section, rel_tol=2e-2)
    assert math.isclose(Welz_cm3 * 1e3, wz_section, rel_tol=2.5e-2)
