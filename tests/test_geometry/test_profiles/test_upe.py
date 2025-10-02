"""Tests for UPE profiles."""

import json
import math
from pathlib import Path

import pytest
from shapely.testing import assert_geometries_equal

from structuralcodes.geometry import (
    SurfaceGeometry,
)
from structuralcodes.geometry.profiles import (
    UPE,
)
from structuralcodes.materials.basic import (
    ElasticMaterial,
)
from structuralcodes.sections._generic import GenericSection


def load_upe_profiles_data():
    """Load UPE profiles data from upe.json file."""
    json_file = Path(__file__).parent / 'upe.json'
    with open(json_file, 'r') as f:
        profiles_data = []
        for line in f:
            profile_data = json.loads(line.strip())
            profile_name = profile_data['ProfileName']
            h_mm = profile_data['h_mm']
            b_mm = profile_data['b_mm']
            tw_mm = profile_data['tw_mm']
            tf_mm = profile_data['tf_mm']
            r_mm = profile_data['r_mm']
            A_cm2 = profile_data['A_cm2']
            Iy_cm4 = profile_data['Iy_cm4']
            Iz_cm4 = profile_data['Iz_cm4']
            Wely_cm3 = profile_data['Wely_cm3']
            Welz_cm3 = profile_data['Welz_cm3']
            ys_cm = profile_data['ys_cm']

            profiles_data.append(
                (
                    profile_name,
                    h_mm,
                    b_mm,
                    tw_mm,
                    tf_mm,
                    r_mm,
                    A_cm2,
                    Iy_cm4,
                    Iz_cm4,
                    Wely_cm3,
                    Welz_cm3,
                    ys_cm,
                )
            )
    return profiles_data


@pytest.mark.parametrize(
    (
        'profile_name, h_mm, b_mm, '
        'tw_mm, tf_mm, r_mm, '
        'A_cm2, Iy_cm4, Iz_cm4, '
        'Wely_cm3, Welz_cm3, ys_cm'
    ),
    load_upe_profiles_data(),
)
def test_upe_geometric_data(
    profile_name,
    h_mm,
    b_mm,
    tw_mm,
    tf_mm,
    r_mm,
    A_cm2,
    Iy_cm4,
    Iz_cm4,
    Wely_cm3,
    Welz_cm3,
    ys_cm,
):
    """Test UPE profile geometric data matches JSON data."""
    # Delete not used variables to avoid linting errors
    del (
        A_cm2,
        Iy_cm4,
        Iz_cm4,
        Wely_cm3,
        Welz_cm3,
        ys_cm,
    )

    # Create the profile
    profile = UPE(profile_name)

    # Check geometric properties
    assert math.isclose(profile.h, h_mm, rel_tol=1e-3)
    assert math.isclose(profile.b, b_mm, rel_tol=1e-3)
    assert math.isclose(profile.tw, tw_mm, rel_tol=1e-3)
    assert math.isclose(profile.tf, tf_mm, rel_tol=1e-3)
    assert math.isclose(profile.r, r_mm, rel_tol=1e-3)


@pytest.mark.parametrize(
    (
        'profile_name, h_mm, b_mm, '
        'tw_mm, tf_mm, r_mm, '
        'A_cm2, Iy_cm4, Iz_cm4, '
        'Wely_cm3, Welz_cm3, ys_cm'
    ),
    load_upe_profiles_data(),
)
def test_upe_massprop_data(
    profile_name,
    h_mm,
    b_mm,
    tw_mm,
    tf_mm,
    r_mm,
    A_cm2,
    Iy_cm4,
    Iz_cm4,
    Wely_cm3,
    Welz_cm3,
    ys_cm,
):
    """Test UPE profile mass property matches JSON data."""
    # Delete not used variables to avoid linting errors
    del (
        h_mm,
        b_mm,
        tw_mm,
        tf_mm,
        r_mm,
        ys_cm,
    )
    # Create the profile
    profile = UPE(profile_name)

    # Check mass properties
    assert math.isclose(profile.A, A_cm2 * 1e2, rel_tol=2e-2)
    assert math.isclose(profile.Iy, Iy_cm4 * 1e4, rel_tol=1e-3)
    assert math.isclose(profile.Iz, Iz_cm4 * 1e4, rel_tol=1e-2)
    assert math.isclose(profile.Wely, Wely_cm3 * 1e3, rel_tol=1e-3)
    assert math.isclose(profile.Welz, Welz_cm3 * 1e3, rel_tol=1e-2)


@pytest.mark.parametrize(
    (
        'profile_name, h_mm, b_mm, '
        'tw_mm, tf_mm, r_mm, '
        'A_cm2, Iy_cm4, Iz_cm4, '
        'Wely_cm3, Welz_cm3, ys_cm'
    ),
    load_upe_profiles_data(),
)
def test_upe_polygon(
    profile_name,
    h_mm,
    b_mm,
    tw_mm,
    tf_mm,
    r_mm,
    A_cm2,
    Iy_cm4,
    Iz_cm4,
    Wely_cm3,
    Welz_cm3,
    ys_cm,
):
    """Test UPE profile polygon."""
    # Delete not used variables to avoid linting errors
    del (
        h_mm,
        b_mm,
        tw_mm,
        tf_mm,
        r_mm,
        A_cm2,
        Iy_cm4,
        Iz_cm4,
        Wely_cm3,
        Welz_cm3,
        ys_cm,
    )

    # Ceate the profile
    profile = UPE(profile_name)

    # Get the polygon from the object
    poly1 = profile.polygon
    # Get the polygon from class method passing the full name
    poly2 = UPE.get_polygon(profile_name)
    assert_geometries_equal(poly1, poly2)


@pytest.mark.parametrize(
    (
        'profile_name, h_mm, b_mm, '
        'tw_mm, tf_mm, r_mm, '
        'A_cm2, Iy_cm4, Iz_cm4, '
        'Wely_cm3, Welz_cm3, ys_cm'
    ),
    load_upe_profiles_data(),
)
def test_upe_reduced_names(
    profile_name,
    h_mm,
    b_mm,
    tw_mm,
    tf_mm,
    r_mm,
    A_cm2,
    Iy_cm4,
    Iz_cm4,
    Wely_cm3,
    Welz_cm3,
    ys_cm,
):
    """Test UPE profile with short name."""
    # Delete not used variables to avoid linting errors
    del (
        h_mm,
        b_mm,
        tw_mm,
        tf_mm,
        r_mm,
        A_cm2,
        Iy_cm4,
        Iz_cm4,
        Wely_cm3,
        Welz_cm3,
        ys_cm,
    )

    # Ceate the profile
    profile = UPE(profile_name)
    profile_short = UPE(int(profile_name[3:]))

    assert_geometries_equal(profile.polygon, profile_short.polygon)

    # Use the class method
    assert_geometries_equal(
        UPE.get_polygon(profile_name), UPE.get_polygon(int(profile_name[3:]))
    )


@pytest.mark.parametrize(
    'invalid_name', ['UPE1000', 'UPE318', '224', 1400.0, -10]
)
def test_upe_invalid_name(
    invalid_name,
):
    """Test UPE profile with invalid name raises ValueError."""
    # Invalid name for constructor
    with pytest.raises(ValueError):
        UPE(invalid_name)
    # Invalid name for class method
    with pytest.raises(ValueError):
        UPE.get_polygon(invalid_name)


@pytest.mark.parametrize(
    (
        'profile_name, h_mm, b_mm, '
        'tw_mm, tf_mm, r_mm, '
        'A_cm2, Iy_cm4, Iz_cm4, '
        'Wely_cm3, Welz_cm3, ys_cm'
    ),
    load_upe_profiles_data(),
)
def test_upe_massprops_polygon(
    profile_name,
    h_mm,
    b_mm,
    tw_mm,
    tf_mm,
    r_mm,
    A_cm2,
    Iy_cm4,
    Iz_cm4,
    Wely_cm3,
    Welz_cm3,
    ys_cm,
):
    """Test UPE profile massprops from polygon."""
    # Delete not used variables to avoid linting errors
    del (
        tw_mm,
        tf_mm,
        r_mm,
    )

    # Material properties
    Es = 200000
    steel = ElasticMaterial(E=Es, density=7850)

    # Create geometry
    geo = SurfaceGeometry(UPE.get_polygon(profile_name), material=steel)
    # Create the section
    sec = GenericSection(geo)

    # Get values from section
    a_section = sec.gross_properties.area
    iy_section = sec.gross_properties.iyy_c
    iz_section = sec.gross_properties.izz_c
    wy_section = iy_section / (h_mm / 2)
    wz_section = iz_section / (b_mm - ys_cm * 10)

    # Assert
    assert math.isclose(A_cm2 * 1e2, a_section, rel_tol=2e-2)
    assert math.isclose(Iy_cm4 * 1e4, iy_section, rel_tol=1e-3)
    assert math.isclose(Iz_cm4 * 1e4, iz_section, rel_tol=1e-2)
    assert math.isclose(Wely_cm3 * 1e3, wy_section, rel_tol=1e-3)
    assert math.isclose(Welz_cm3 * 1e3, wz_section, rel_tol=1e-2)
