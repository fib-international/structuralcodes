"""Test for the EC2_2004 module."""

import math

import numpy as np
import pytest
from shapely import Polygon

from structuralcodes.geometry import SurfaceGeometry, add_reinforcement
from structuralcodes.materials.concrete import ConcreteEC2_2004
from structuralcodes.materials.reinforcement import ReinforcementEC2_2004
from structuralcodes.sections._rc_shear._EC2_2004 import (
    BeamSection,
    ShearReinforcement,
    max_area_shear_reinf,
    required_shear_reinf,
    shearcap_rectangular_section,
    shearcap_rectangular_uncracked_prestressed,
    shearcap_reinf_rectangular_section,
)


@pytest.mark.parametrize(
    'diameter, s, material, expected',
    [
        (10, 200, ReinforcementEC2_2004(500, 200000, 510, 0.06), 10),
        (12, 200, ReinforcementEC2_2004(500, 200000, 510, 0.06), 12),
        (14, 200, ReinforcementEC2_2004(500, 200000, 510, 0.06), 14),
    ],
)
def test_ShearReinforcement_diameter(diameter, s, material, expected):
    """Test diameter property of ShearReninforcement class."""
    assert math.isclose(
        ShearReinforcement(diameter, s, material).diameter,
        expected,
        rel_tol=10e-5,
    )


@pytest.mark.parametrize(
    'diameter, s, material, n, expected',
    [
        (
            10,
            200,
            ReinforcementEC2_2004(500, 200000, 510, 0.06),
            2,
            2 * 10**2 / 4 * np.pi,
        ),
        (
            12,
            200,
            ReinforcementEC2_2004(500, 200000, 510, 0.06),
            2,
            2 * 12**2 / 4 * np.pi,
        ),
        (
            14,
            200,
            ReinforcementEC2_2004(500, 200000, 510, 0.06),
            3,
            3 * 14**2 / 4 * np.pi,
        ),
    ],
)
def test_ShearReinforcement_Asw(diameter, s, material, n, expected):
    """Test Asw property of ShearReninforcement class."""
    assert math.isclose(
        ShearReinforcement(diameter, s, material, n=n).Asw,
        expected,
        rel_tol=10e-5,
    )


@pytest.mark.parametrize(
    'diameter, s, material, n, expected',
    [
        (
            10,
            200,
            ReinforcementEC2_2004(500, 200000, 510, 0.06),
            2,
            2 * 10**2 / 4 * np.pi / 200,
        ),
        (
            10,
            250,
            ReinforcementEC2_2004(500, 200000, 510, 0.06),
            3,
            3 * 10**2 / 4 * np.pi / 250,
        ),
        (
            12,
            200,
            ReinforcementEC2_2004(500, 200000, 510, 0.06),
            2,
            2 * 12**2 / 4 * np.pi / 200,
        ),
    ],
)
def test_ShearReinforcement_Asw_s(diameter, s, material, n, expected):
    """Test Asw_s property of ShearReninforcement class."""
    assert math.isclose(
        ShearReinforcement(diameter, s, material, n=n).Asw_s,
        expected,
        rel_tol=10e-5,
    )


@pytest.mark.parametrize(
    'diameter, s, material, expected',
    [
        (10, 200, ReinforcementEC2_2004(500, 200000, 510, 0.06), 200),
        (10, 220, ReinforcementEC2_2004(500, 200000, 510, 0.06), 220),
        (10, 240, ReinforcementEC2_2004(500, 200000, 510, 0.06), 240),
    ],
)
def test_ShearReinforcement_s(diameter, s, material, expected):
    """Test s property of ShearReinforcement class."""
    assert math.isclose(
        ShearReinforcement(diameter, s, material).s,
        expected,
        rel_tol=10e-5,
    )


@pytest.mark.parametrize(
    'diameter, s, material, n, expected',
    [
        (10, 200, ReinforcementEC2_2004(500, 200000, 510, 0.06), 2, 2),
        (10, 220, ReinforcementEC2_2004(500, 200000, 510, 0.06), 3, 3),
        (10, 240, ReinforcementEC2_2004(500, 200000, 510, 0.06), 4, 4),
    ],
)
def test_ShearReinforcement_n(diameter, s, material, n, expected):
    """Test n property of ShearReinforcement class."""
    assert math.isclose(
        ShearReinforcement(diameter, s, material, n=n).n,
        expected,
        rel_tol=10e-5,
    )


@pytest.mark.parametrize(
    'diameter, s, material, alpha, expected',
    [
        (10, 200, ReinforcementEC2_2004(500, 200000, 510, 0.06), 90, 90),
        (10, 220, ReinforcementEC2_2004(500, 200000, 510, 0.06), 45, 45),
        (10, 240, ReinforcementEC2_2004(500, 200000, 510, 0.06), 60, 60),
    ],
)
def test_ShearReinforcement_alpha(diameter, s, material, alpha, expected):
    """Test alpha property of ShearReinforcement class."""
    assert math.isclose(
        ShearReinforcement(diameter, s, material, alpha=alpha).alpha,
        expected,
        rel_tol=10e-5,
    )


@pytest.mark.parametrize(
    'fcd, fck, bw, s, fywd, NEd, Ac, alpha, expected',
    [
        (20 / 1.5, 20, 100, 200, 435, 100e3, 100 * 400, 90, 201),
        (20 / 1.5, 20, 100, 200, 435, 100e3, 100 * 400, 45, 284),
    ],
)
def test_max_area_shear_reinf(fcd, fck, bw, s, fywd, NEd, Ac, alpha, expected):
    """Test the max_area_shear_reinf function."""
    # Test values are reused from test_ec2_2004_shear.py
    concrete = ConcreteEC2_2004(
        fck=fck,
        fcd=fcd,
    )
    polygon = Polygon(
        [
            (-bw / 2, -Ac / bw / 2),
            (bw / 2, -Ac / bw / 2),
            (bw / 2, Ac / bw / 2),
            (-bw / 2, Ac / bw / 2),
        ]
    )
    geometry = SurfaceGeometry(
        poly=polygon,
        material=concrete,
    )
    section = BeamSection(geometry)
    shear_material = ReinforcementEC2_2004(
        # Es, ftk and ftk are not used in calculations, but are requred as
        # input
        fyk=fywd * 1.15,
        Es=200000,
        ftk=510,
        epsuk=0.06,
        gamma_s=1.15,
    )
    shear_reinf = ShearReinforcement(
        diameter=10,
        material=shear_material,
        s=s,
        alpha=alpha,
    )
    assert math.isclose(
        max_area_shear_reinf(section, shear_reinf, NEd),
        expected,
        rel_tol=0.01,
    )


@pytest.mark.parametrize(
    'VEd, z, theta, fywd, expected',
    [
        (100e3, 300, 45, 500 / 1.15, 0.76666),
        (150e3, 300, 45, 500 / 1.15, 1.14999),
        (200e3, 350, 30, 500 / 1.15, 0.75880),
        (250e3, 350, 45, 500 / 1.15, 1.64285),
        (120e3, 400, 45, 500 / 1.15, 0.68999),
        (180e3, 350, 35, 500 / 1.15, 0.82824),
    ],
)
def test_required_shear_reinf(VEd, z, theta, fywd, expected):
    """Test the required_shear_reinf function."""
    # Concrete properties are not used in the calculation
    concrete = ConcreteEC2_2004(35)
    polygon = Polygon(
        # Divide z by 0.9 to obtain correct value for d
        [
            (0, 0),
            (300, 0),
            (300, z / 0.9),
            (0, z / 0.9),
        ]
    )
    geometry = SurfaceGeometry(
        poly=polygon,
        material=concrete,
    )
    # Add some reinforcement at the bottom of cross section so d can be
    # obtained. Longitudinal reinforcement is not included in the calculation
    material = ReinforcementEC2_2004(
        fyk=fywd * 1.15,
        Es=200000,
        ftk=510,
        epsuk=0.06,
        gamma_s=1.15,
    )
    geometry = add_reinforcement(
        geometry,
        (0, 0),
        10,
        material,
    )
    section = BeamSection(geometry)
    shear_reinf = required_shear_reinf(section, material, VEd, theta)
    assert math.isclose(
        shear_reinf.Asw_s,
        expected,
        rel_tol=0.01,
    )


@pytest.mark.parametrize(
    'fck, d, Asl, bw, NEd, Ac, k1, gamma_c, expected',
    [
        (20, 250, 200, 100, 5e4, 30000, 0.15, 1.5, 20538),
        (37.5, 450, 1000, 500, 1e8, 250000, 0.15, 1.5, 283334),
        (20, 250, 0, 100, 5e4, 30000, 0.15, 1.5, 16425),
        (37.5, 450, 0, 500, 1e8, 250000, 0.15, 1.5, 272475),
    ],
)
def test_shearcap_rectangular_section_VRdc(
    fck,
    d,
    Asl,
    bw,
    NEd,
    Ac,
    k1,
    gamma_c,
    expected,
):
    """Test VRdc in the shearcap_rectangular_section functoin."""
    concrete = ConcreteEC2_2004(fck=fck, gamma_c=gamma_c)
    reinforcement = ReinforcementEC2_2004(
        fyk=500,
        Es=200000,
        ftk=510,
        epsuk=0.06,
    )
    c = Ac / bw - d
    polygon = Polygon(
        [
            (-bw / 2, -d / 2 - c),
            (bw / 2, -d / 2 - c),
            (bw / 2, d / 2),
            (-bw / 2, d / 2),
        ]
    )
    geometry = SurfaceGeometry(
        poly=polygon,
        material=concrete,
    )
    geometry = add_reinforcement(
        geometry, (0, -d / 2), np.sqrt(Asl * 4 / np.pi), reinforcement
    )
    section = BeamSection(geometry)
    result = shearcap_rectangular_section(section, NEd, k1)
    assert math.isclose(result[0], expected, rel_tol=0.01)


@pytest.mark.parametrize(
    'bw, d, fck, expected',
    [
        (100, 250, 20, 91770.0),
        (100, 250, 37.5, 159375.0),
        (500, 450, 37.5, 1434375.0),
    ],
)
def test_shearcap_rectangular_section_VEdmax_unreinf(
    bw,
    d,
    fck,
    expected,
):
    """Test VEdmax_unreinf in the shearcap_rectangular_section function."""
    concrete = ConcreteEC2_2004(fck)
    reinforcement = ReinforcementEC2_2004(
        fyk=500,
        Es=200000,
        ftk=510,
        epsuk=0.06,
    )
    polygon = Polygon(
        [
            (-bw / 2, -d / 2),
            (bw / 2, -d / 2),
            (bw / 2, d / 2),
            (-bw / 2, d / 2),
        ]
    )
    geometry = SurfaceGeometry(polygon, concrete)
    # Longitudinal reinforcement not used in calculation, but is required to
    # run the driver function
    geometry = add_reinforcement(geometry, (0, -d / 2), 10, reinforcement)
    section = BeamSection(geometry)
    result = shearcap_rectangular_section(section)
    assert math.isclose(result[1], expected, rel_tol=0.01)


@pytest.mark.parametrize(
    'bw, fctd, NEd, Ac, L_x, L_pt2, expected',
    [
        (300, 1.2, 1.5e6, 150e3, 0, 2, 120000),
        (300, 1.2, 1.5e6, 150e3, 0.5, 2, 210713),
        (300, 1.2, 1.5e6, 150e3, 1.0, 2, 272764),
        (300, 1.2, 1.5e6, 150e3, 2.0, 2, 366606),
        (300, 1.2, 1.5e6, 150e3, None, 2, 366606),
        (300, 1.2, 1.5e6, 150e3, None, None, 366606),
    ],
)
def test_shearcap_rectangular_uncracked_prestressed(
    bw,
    fctd,
    NEd,
    Ac,
    L_x,
    L_pt2,
    expected,
):
    """Test the shearcap_rectangular_uncracked_prestressed function."""
    h = Ac / bw
    polygon = Polygon(
        [
            (-bw / 2, -h / 2),
            (bw / 2, -h / 2),
            (bw / 2, h / 2),
            (-bw / 2, h / 2),
        ]
    )
    concrete = ConcreteEC2_2004(25, fctk_5=fctd * 1.5)
    geometry = SurfaceGeometry(polygon, concrete)
    section = BeamSection(geometry)
    assert math.isclose(
        shearcap_rectangular_uncracked_prestressed(
            section,
            NEd,
            L_x=L_x,
            L_pt2=L_pt2,
        ),
        expected,
        rel_tol=0.01,
    )


@pytest.mark.parametrize(
    'diam, s, z, theta, fyk, alpha, gamma_s, expected',
    [
        (16, 200, 4880, 45, 400, 90, 1.15, 3413e3),
        (20, 200, 4880, 45, 400, 90, 1.15, 5332e3),
        # Check if ValueError is raised.
        (20, 200, 1757, 100, 400, 90, 1.15, 1),
        (16, 200, 4880, 45, 400, 45, 1.15, 4825e3),
        (20, 200, 4880, 45, 400, 45, 1.15, 7537e3),
        # Check if ValueError is raised.
        (20, 200, 1757, 100, 400, 45, 1.15, 1),
    ],
)
def test_shearcap_reinf_rectangular_section_VRds(
    diam, s, z, theta, fyk, alpha, gamma_s, expected
):
    """Test VRds in the shearcap_reinf_rectangular_section function."""
    reinforcement = ReinforcementEC2_2004(fyk, 200000, 510, 0.06, gamma_s)
    shear_reinf = ShearReinforcement(diam, s, reinforcement, 2, alpha=alpha)
    d = z / 0.9
    concrete = ConcreteEC2_2004(35)
    polygon = Polygon(
        [
            (-250, -d / 2),
            (250, -d / 2),
            (250, d / 2),
            (-250, -d / 2),
        ]
    )
    geometry = SurfaceGeometry(polygon, concrete)
    geometry = add_reinforcement(geometry, (0, -d / 2), 10, reinforcement)
    section = BeamSection(geometry)
    try:
        result = shearcap_reinf_rectangular_section(
            section, shear_reinf, theta=theta
        )
        assert math.isclose(result[0], expected, rel_tol=0.01)
    except ValueError:
        with pytest.raises(ValueError) as exc_info:
            result = shearcap_reinf_rectangular_section(
                section, shear_reinf, theta=theta
            )
            assert math.isclose(result[0], expected, rel_tol=0.01)
        assert str(exc_info.value).startswith(
            'Wrong value for theta is chosen.'
        )


@pytest.mark.parametrize(
    (
        'bw, z, fck, theta, NEd, Ac, gamma_c, alpha, alpha_cc, limit_fyd, '
        'expected'
    ),
    [
        (100, 300, 20, 45, 100e3, 100 * 400, 1.5, 90, 1.0, False, 131100),
        (100, 300, 20, 21.8, 100e3, 100 * 400, 1.5, 90, 1.0, False, 90409.12),
        (100, 300, 20, 45, 100e3, 100 * 400, 1.5, 45, 1.0, False, 262200),
        (100, 300, 20, 21.8, 100e3, 100 * 400, 1.5, 45, 1.0, False, 126570.19),
        (100, 300, 70, 45, 100e3, 100 * 400, 1.5, 90, 1.0, False, 318600),
        (100, 300, 70, 21.8, 100e3, 100 * 400, 1.5, 90, 1.0, False, 219712.79),
        (100, 300, 70, 45, 100e3, 100 * 400, 1.5, 45, 1.0, False, 637200),
        (100, 300, 70, 21.8, 100e3, 100 * 400, 1.5, 45, 1.0, False, 307591.63),
        (100, 300, 20, 45, 100e3, 100 * 400, 1.5, 90, 1.0, True, 142500),
        (100, 300, 20, 21.8, 100e3, 100 * 400, 1.5, 90, 1.0, True, 98270.78),
        (100, 300, 20, 45, 100e3, 100 * 400, 1.5, 45, 1.0, True, 285000),
        (100, 300, 20, 21.8, 100e3, 100 * 400, 1.5, 45, 1.0, True, 137576.29),
        (100, 300, 70, 45, 100e3, 100 * 400, 1.5, 90, 1.0, True, 405625),
        (100, 300, 70, 21.8, 100e3, 100 * 400, 1.5, 90, 1.0, True, 279726.93),
        (100, 300, 70, 45, 100e3, 100 * 400, 1.5, 45, 1.0, True, 811250),
        (100, 300, 70, 21.8, 100e3, 100 * 400, 1.5, 45, 1.0, True, 391609.72),
    ],
)
def test_shearcap_reinf_rectangular_section_VRdmax(
    bw, z, fck, theta, NEd, Ac, gamma_c, alpha, alpha_cc, limit_fyd, expected
):
    """Test VRdmax in the shearcap_reinf_rectangular_section function."""
    concrete = ConcreteEC2_2004(fck, alpha_cc=alpha_cc, gamma_c=gamma_c)
    reinforcement = ReinforcementEC2_2004(500, 200000, 510, 0.06)
    h = Ac / bw
    d = z / 0.9
    polygon = Polygon(
        [
            (-bw / 2, -h / 2),
            (bw / 2, -h / 2),
            (bw / 2, h / 2),
            (-bw / 2, h / 2),
        ]
    )
    geometry = SurfaceGeometry(polygon, concrete)
    geometry = add_reinforcement(geometry, (0, h / 2 - d), 10, reinforcement)
    section = BeamSection(geometry)
    shear_reinf = ShearReinforcement(10, 200, reinforcement, alpha=alpha)
    assert math.isclose(
        shearcap_reinf_rectangular_section(
            section, shear_reinf, NEd, theta=theta, limit_fyd=limit_fyd
        )[1],
        expected,
        rel_tol=0.01,
    )
