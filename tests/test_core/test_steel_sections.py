"""Tests for the Steel Sections."""

import math

import pytest

from structuralcodes.geometry import (
    HE,
    IPE,
    IPN,
    UB,
    UBP,
    UC,
    UPN,
    CompoundGeometry,
    SurfaceGeometry,
)
from structuralcodes.materials.constitutive_laws import (
    Elastic,
    ElasticPlastic,
    UserDefined,
)
from structuralcodes.sections._generic import GenericSection
from structuralcodes.sections.section_integrators._marin_integration import (
    marin_integration,
)


# Test steel I section
def _iyy_I_beam(h: float, b: float, tw: float, tf: float, r: float) -> float:
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


def _izz_I_beam(h: float, b: float, tw: float, tf: float, r: float) -> float:
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


def _wyel_I_beam(h: float, b: float, tw: float, tf: float, r: float) -> float:
    return _iyy_I_beam(h, b, tw, tf, r) / (h / 2.0)


def _wzel_I_beam(h: float, b: float, tw: float, tf: float, r: float) -> float:
    return _izz_I_beam(h, b, tw, tf, r) / (b / 2.0)


def _wypl_I_beam(h: float, b: float, tw: float, tf: float, r: float) -> float:
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


def _wzpl_I_beam(h: float, b: float, tw: float, tf: float, r: float) -> float:
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


@pytest.mark.parametrize(
    'cls, invalid_name',
    [
        (IPE, 'IPE95'),
        (IPE, 'az132'),
        (HE, 'HEA3000'),
        (HE, 'HEB3000'),
        (HE, 'HEM3000'),
        (UB, 'UB203'),
        (UC, 'UC152'),
        (UBP, 'UBP203'),
        (IPN, 'IPN125'),
        (UPN, 'UPN123'),
    ],
)
def test_names_invalid(cls, invalid_name):
    """Test Steel I section elastic strength."""
    # Invalid name for class method
    with pytest.raises(ValueError):
        cls.get_polygon(name=invalid_name)
    # Invalid name for constructure
    with pytest.raises(ValueError):
        cls(invalid_name).polygon


@pytest.mark.parametrize(
    'cls, name',
    [
        (IPE, 'IPE80'),
        (IPE, 'IPE100'),
        (IPE, 'IPE120'),
        (IPE, 'IPE140'),
        (IPE, 'IPE160'),
        (IPE, 'IPE180'),
        (IPE, 'IPE200'),
        (IPE, 'IPE220'),
        (IPE, 'IPE240'),
        (IPE, 'IPE270'),
        (IPE, 'IPE300'),
        (IPE, 'IPE330'),
        (IPE, 'IPE360'),
        (IPE, 'IPE400'),
        (IPE, 'IPE450'),
        (IPE, 'IPE500'),
        (IPE, 'IPE550'),
        (IPE, 'IPE600'),
        (HE, 'HEA300'),
        (HE, 'HEB300'),
        (HE, 'HEM300'),
        (HE, 'HEA400'),
        (HE, 'HEB400'),
        (HE, 'HEM400'),
        (UB, 'UB203x102x23'),
        (UB, 'UB203x133x30'),
        (UC, 'UC152x152x44'),
        (UC, 'UC203x203x52'),
        (UBP, 'UBP203x203x54'),
        (UBP, 'UBP254x254x71'),
        (IPN, 'IPN80'),
        (IPN, 'IPN300'),
        (IPN, 'IPN320'),
        (IPN, 'IPN340'),
        (IPN, 'IPN360'),
        (IPN, 'IPN380'),
        (IPN, 'IPN400'),
        (UPN, 'UPN120'),
        (UPN, 'UPN140'),
        (UPN, 'UPN160'),
        (UPN, 'UPN180'),
        (UPN, 'UPN200'),
        (UPN, 'UPN300'),
    ],
)
def test_section_get_polygon(cls, name):
    """Test Steel I section elastic strength."""
    Es = 206000
    fy = 355
    eps_su = fy / Es
    steel = ElasticPlastic(E=Es, fy=fy, eps_su=eps_su)

    # Create geometry
    geo = CompoundGeometry([SurfaceGeometry(cls.get_polygon(name), steel)])
    sec = GenericSection(geo)

    # Compute expected values
    xy = sec.geometry.geometries[0].polygon.exterior.coords.xy

    A = marin_integration(xy[0], xy[1], 0, 0)

    A_sec = sec.gross_properties.area
    assert math.isclose(A_sec, A)


@pytest.mark.parametrize(
    'cls, name, A, Iy, Wely, Wply, iy, Iz, Welz, Wplz, iz',
    [
        (
            IPE,
            'IPE80',
            7.6,
            80.13,
            20.03,
            23.21,
            3.2,
            8.489,
            3.69,
            5.817,
            1.05,
        ),
        (IPE, 'IPE100', 10.3, 171, 34.2, 39.4, 4, 15.91, 5.788, 9.145, 1.24),
        (
            IPE,
            'IPE120',
            13.2,
            317.7,
            52.95,
            60.72,
            4.9,
            27.66,
            8.646,
            13.58,
            1.45,
        ),
        (
            IPE,
            'IPE140',
            16.4,
            541.2,
            77.31,
            88.34,
            5.7,
            44.91,
            12.3,
            19.24,
            1.65,
        ),
        (
            IPE,
            'IPE160',
            20.1,
            869.2,
            108.6,
            123.8,
            6.5,
            68.31,
            16.66,
            26.09,
            1.84,
        ),
        (
            IPE,
            'IPE180',
            23.9,
            1316,
            146.3,
            166.4,
            7.4,
            100.8,
            22.16,
            34.59,
            2.05,
        ),
        (
            IPE,
            'IPE200',
            28.5,
            1943,
            194.3,
            220.6,
            8.2,
            142.3,
            28.47,
            44.61,
            2.24,
        ),
        (
            IPE,
            'IPE220',
            33.4,
            2771,
            251.9,
            285.4,
            9.1,
            204.8,
            37.25,
            58.11,
            2.48,
        ),
        (
            IPE,
            'IPE240',
            39.1,
            3891,
            324.3,
            366.6,
            9.9,
            283.6,
            47.27,
            73.92,
            2.69,
        ),
        (IPE, 'IPE270', 45.9, 5789, 428.8, 483.9, 11.2, 419.8, 62.2, 96.95, 3),
        (IPE, 'IPE300', 53.8, 8356, 557, 628.3, 12.4, 603.7, 80.5, 125.2, 3.3),
        (
            IPE,
            'IPE330',
            62.6,
            11770,
            713.1,
            804.3,
            13.7,
            788.1,
            98.51,
            153.6,
            3.5,
        ),
        (
            IPE,
            'IPE360',
            72.7,
            16260,
            903.6,
            1019,
            14.9,
            1043,
            122.7,
            191,
            3.79,
        ),
        (
            IPE,
            'IPE400',
            84.5,
            23120,
            1156,
            1307,
            16.5,
            1317,
            146.4,
            229,
            3.95,
        ),
        (
            IPE,
            'IPE450',
            98.8,
            33740,
            1499,
            1701,
            18.4,
            1675,
            176.4,
            276.3,
            4.1,
        ),
        (
            IPE,
            'IPE500',
            115.5,
            48190,
            1927,
            2194,
            20.4,
            2141,
            214.1,
            335.8,
            4.3,
        ),
        (IPE, 'IPE550', 134.4, 67110, 2440, 2787, 22.3, 2667, 254, 400.5, 4.4),
        (IPE, 'IPE600', 156, 92080, 3069, 3512, 24.2, 3387, 307.9, 485.6, 4.6),
        (
            HE,
            'HEA300',
            112.5,
            18260,
            1260,
            1383,
            12.74,
            6310,
            420.6,
            641.2,
            7.49,
        ),
        (
            UB,
            'UB203x133x25',
            31.97,
            2340,
            230.3,
            257.7,
            8.56,
            307.6,
            46.19,
            70.94,
            3.1,
        ),
        (
            UC,
            'UC152x152x37',
            47.11,
            2210,
            273.2,
            308.8,
            6.85,
            706.2,
            91.48,
            139.6,
            3.87,
        ),
        (IPN, 'IPN200', 33.4, 2140, 214, 250, 8, 117, 26, 43.5, 1.87),
        (UPN, 'UPN180', 28, 1350, 150, 179, 6.95, 114, 22.4, 42.9, 2.02),
    ],
)
def test_section_properties(
    cls, name, A, Iy, Wely, Wply, iy, Iz, Welz, Wplz, iz
):
    """Test Steel section comparing section properties with tables."""
    # Create polygon representing geometyr of profile
    i_beam = cls(name)

    # Compare wth stored values with a 2% tolerance
    assert math.isclose(A * 1e2, i_beam.A, rel_tol=2e-2)
    assert math.isclose(Iy * 1e4, i_beam.Iy, rel_tol=2e-2)
    assert math.isclose(Iz * 1e4, i_beam.Iz, rel_tol=2e-2)
    assert math.isclose(Wely * 1e3, i_beam.Wely, rel_tol=2e-2)
    assert math.isclose(Welz * 1e3, i_beam.Welz, rel_tol=2e-2)
    assert math.isclose(Wply * 1e3, i_beam.Wply, rel_tol=2e-2)
    assert math.isclose(Wplz * 1e3, i_beam.Wplz, rel_tol=2e-2)
    assert math.isclose(iy * 10, i_beam.iy, rel_tol=2e-2)
    assert math.isclose(iz * 10, i_beam.iz, rel_tol=2e-2)


@pytest.mark.parametrize(
    'cls, name',
    [
        (IPE, 'IPE80'),
        (IPE, 'IPE100'),
        (IPE, 'IPE120'),
        (IPE, 'IPE140'),
        (IPE, 'IPE160'),
        (IPE, 'IPE180'),
        (IPE, 'IPE200'),
        (IPE, 'IPE220'),
        (IPE, 'IPE240'),
        (IPE, 'IPE270'),
        (IPE, 'IPE300'),
        (IPE, 'IPE330'),
        (IPE, 'IPE360'),
        (IPE, 'IPE400'),
        (IPE, 'IPE450'),
        (IPE, 'IPE500'),
        (IPE, 'IPE550'),
        (IPE, 'IPE600'),
        (HE, 'HEA300'),
        (HE, 'HEB300'),
        (HE, 'HEM300'),
        (HE, 'HEA400'),
        (HE, 'HEB400'),
        (HE, 'HEM400'),
        (UB, 'UB203x102x23'),
        (UB, 'UB203x133x30'),
        (UC, 'UC152x152x44'),
        (UC, 'UC203x203x52'),
        (UBP, 'UBP203x203x54'),
        (UBP, 'UBP254x254x71'),
    ],
)
def test_Isection_elastic_fiber(cls, name):
    """Test Steel I section elastic strength."""
    Es = 206000
    fy = 355
    eps_su = fy / Es
    steel = ElasticPlastic(E=Es, fy=fy, eps_su=eps_su)

    # Create geometry
    i_beam = cls(name)
    geo = CompoundGeometry([SurfaceGeometry(i_beam.polygon, steel)])

    # Compute expected values
    wy_el = _wyel_I_beam(i_beam.h, i_beam.b, i_beam.tw, i_beam.tf, i_beam.r)
    wz_el = _wzel_I_beam(i_beam.h, i_beam.b, i_beam.tw, i_beam.tf, i_beam.r)
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
    'cls, name',
    [
        (IPE, 'IPE80'),
        (IPE, 'IPE100'),
        (IPE, 'IPE120'),
        (IPE, 'IPE140'),
        (IPE, 'IPE160'),
        (IPE, 'IPE180'),
        (IPE, 'IPE200'),
        (IPE, 'IPE220'),
        (IPE, 'IPE240'),
        (IPE, 'IPE270'),
        (IPE, 'IPE300'),
        (IPE, 'IPE330'),
        (IPE, 'IPE360'),
        (IPE, 'IPE400'),
        (IPE, 'IPE450'),
        (IPE, 'IPE500'),
        (IPE, 'IPE550'),
        (IPE, 'IPE600'),
        (HE, 'HEA300'),
        (HE, 'HEB300'),
        (HE, 'HEM300'),
        (HE, 'HEA400'),
        (HE, 'HEB400'),
        (HE, 'HEM400'),
        (UB, 'UB203x102x23'),
        (UB, 'UB203x133x30'),
        (UC, 'UC152x152x44'),
        (UC, 'UC203x203x52'),
        (UBP, 'UBP203x203x54'),
        (UBP, 'UBP254x254x71'),
    ],
)
def test_Isection_elastic_marin(cls, name):
    """Test Steel I section elastic strength."""
    Es = 206000
    fy = 355
    eps_su = fy / Es
    steel = ElasticPlastic(E=Es, fy=fy, eps_su=eps_su)

    # Create geometry
    i_beam = cls(name)
    geo = CompoundGeometry([SurfaceGeometry(i_beam.polygon, steel)])

    # Compute expected values
    wy_el = _wyel_I_beam(i_beam.h, i_beam.b, i_beam.tw, i_beam.tf, i_beam.r)
    wz_el = _wzel_I_beam(i_beam.h, i_beam.b, i_beam.tw, i_beam.tf, i_beam.r)
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
    'cls, name',
    [
        (IPE, 'IPE80'),
        (IPE, 'IPE100'),
        (IPE, 'IPE120'),
        (IPE, 'IPE140'),
        (IPE, 'IPE160'),
        (IPE, 'IPE180'),
        (IPE, 'IPE200'),
        (IPE, 'IPE220'),
        (IPE, 'IPE240'),
        (IPE, 'IPE270'),
        (IPE, 'IPE300'),
        (IPE, 'IPE330'),
        (IPE, 'IPE360'),
        (IPE, 'IPE400'),
        (IPE, 'IPE450'),
        (IPE, 'IPE500'),
        (IPE, 'IPE550'),
        (IPE, 'IPE600'),
        (HE, 'HEA300'),
        (HE, 'HEB300'),
        (HE, 'HEM300'),
        (HE, 'HEA400'),
        (HE, 'HEB400'),
        (HE, 'HEM400'),
        (UB, 'UB203x102x23'),
        (UB, 'UB203x133x30'),
        (UC, 'UC152x152x44'),
        (UC, 'UC203x203x52'),
        (UBP, 'UBP203x203x54'),
        (UBP, 'UBP254x254x71'),
    ],
)
def test_Isection_plastic_fiber(cls, name):
    """Test Steel I section elastic strength."""
    Es = 206000
    fy = 355
    eps_su = 0.15
    steel = ElasticPlastic(E=Es, fy=fy, eps_su=eps_su)

    # Create geometry
    i_beam = cls(name)
    geo = CompoundGeometry([SurfaceGeometry(i_beam.polygon, steel)])

    # Compute expected values
    wy_pl = _wypl_I_beam(i_beam.h, i_beam.b, i_beam.tw, i_beam.tf, i_beam.r)
    wz_pl = _wzpl_I_beam(i_beam.h, i_beam.b, i_beam.tw, i_beam.tf, i_beam.r)
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
    'cls, name',
    [
        (IPE, 'IPE80'),
        (IPE, 'IPE100'),
        (IPE, 'IPE120'),
        (IPE, 'IPE140'),
        (IPE, 'IPE160'),
        (IPE, 'IPE180'),
        (IPE, 'IPE200'),
        (IPE, 'IPE220'),
        (IPE, 'IPE240'),
        (IPE, 'IPE270'),
        (IPE, 'IPE300'),
        (IPE, 'IPE330'),
        (IPE, 'IPE360'),
        (IPE, 'IPE400'),
        (IPE, 'IPE450'),
        (IPE, 'IPE500'),
        (IPE, 'IPE550'),
        (IPE, 'IPE600'),
        (HE, 'HEA300'),
        (HE, 'HEB300'),
        (HE, 'HEM300'),
        (HE, 'HEA400'),
        (HE, 'HEB400'),
        (HE, 'HEM400'),
        (UB, 'UB203x102x23'),
        (UB, 'UB203x133x30'),
        (UC, 'UC152x152x44'),
        (UC, 'UC203x203x52'),
        (UBP, 'UBP203x203x54'),
        (UBP, 'UBP254x254x71'),
    ],
)
def test_Isection_plastic_marin(cls, name):
    """Test Steel I section elastic strength."""
    Es = 206000
    fy = 355
    eps_su = 0.15
    steel = ElasticPlastic(E=Es, fy=fy, eps_su=eps_su)

    # Create geometry
    i_beam = cls(name)
    geo = CompoundGeometry([SurfaceGeometry(i_beam.polygon, steel)])

    # Compute expected values
    wy_pl = _wypl_I_beam(i_beam.h, i_beam.b, i_beam.tw, i_beam.tf, i_beam.r)
    wz_pl = _wzpl_I_beam(i_beam.h, i_beam.b, i_beam.tw, i_beam.tf, i_beam.r)
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
    'cls, name',
    [
        (IPE, 'IPE80'),
        (IPE, 'IPE100'),
        (IPE, 'IPE120'),
        (IPE, 'IPE140'),
        (IPE, 'IPE160'),
        (IPE, 'IPE180'),
        (IPE, 'IPE200'),
        (IPE, 'IPE220'),
        (IPE, 'IPE240'),
        (IPE, 'IPE270'),
        (IPE, 'IPE300'),
        (IPE, 'IPE330'),
        (IPE, 'IPE360'),
        (IPE, 'IPE400'),
        (IPE, 'IPE450'),
        (IPE, 'IPE500'),
        (IPE, 'IPE550'),
        (IPE, 'IPE600'),
        (HE, 'HEA300'),
        (HE, 'HEB300'),
        (HE, 'HEM300'),
        (HE, 'HEA400'),
        (HE, 'HEB400'),
        (HE, 'HEM400'),
        (UB, 'UB203x102x23'),
        (UB, 'UB203x133x30'),
        (UC, 'UC152x152x44'),
        (UC, 'UC203x203x52'),
        (UBP, 'UBP203x203x54'),
        (UBP, 'UBP254x254x71'),
    ],
)
def test_Isection_elastic_material_marin(cls, name):
    """Test Steel I section elastic strength."""
    Es = 206000
    fy = 355
    steel = Elastic(E=Es)
    steel.set_ultimate_strain(fy / Es)
    # Create geometry
    i_beam = cls(name)
    geo = CompoundGeometry([SurfaceGeometry(i_beam.polygon, steel)])

    # Compute expected values
    wy_el = _wyel_I_beam(i_beam.h, i_beam.b, i_beam.tw, i_beam.tf, i_beam.r)
    wz_el = _wzel_I_beam(i_beam.h, i_beam.b, i_beam.tw, i_beam.tf, i_beam.r)
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
    'cls, name',
    [
        (IPE, 'IPE80'),
        (IPE, 'IPE100'),
        (IPE, 'IPE120'),
        (IPE, 'IPE140'),
        (IPE, 'IPE160'),
        (IPE, 'IPE180'),
        (IPE, 'IPE200'),
        (IPE, 'IPE220'),
        (IPE, 'IPE240'),
        (IPE, 'IPE270'),
        (IPE, 'IPE300'),
        (IPE, 'IPE330'),
        (IPE, 'IPE360'),
        (IPE, 'IPE400'),
        (IPE, 'IPE450'),
        (IPE, 'IPE500'),
        (IPE, 'IPE550'),
        (IPE, 'IPE600'),
        (HE, 'HEA300'),
        (HE, 'HEB300'),
        (HE, 'HEM300'),
        (HE, 'HEA400'),
        (HE, 'HEB400'),
        (HE, 'HEM400'),
        (UB, 'UB203x102x23'),
        (UB, 'UB203x133x30'),
        (UC, 'UC152x152x44'),
        (UC, 'UC203x203x52'),
        (UBP, 'UBP203x203x54'),
        (UBP, 'UBP254x254x71'),
    ],
)
def test_Isection_user_material_marin(cls, name):
    """Test Steel I section elastic strength."""
    Es = 206000
    fy = 355
    eps_su = 7e-2
    steel = UserDefined(
        x=[-eps_su, -fy / Es, 0, fy / Es, eps_su], y=[-fy, -fy, 0, fy, fy]
    )
    # Create geometry
    i_beam = cls(name)
    geo = CompoundGeometry([SurfaceGeometry(i_beam.polygon, steel)])
    # Compute expected values
    wy_pl = _wypl_I_beam(i_beam.h, i_beam.b, i_beam.tw, i_beam.tf, i_beam.r)
    wz_pl = _wzpl_I_beam(i_beam.h, i_beam.b, i_beam.tw, i_beam.tf, i_beam.r)
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


@pytest.mark.parametrize(
    'cls, name, Wyel, Wzel, Wypl, Wzpl',
    [
        (IPE, 'IPE80', 20.03e3, 3.69e3, 23.21e3, 5.817e3),
        (IPE, 'IPE100', 34.2e3, 5.788e3, 39.4e3, 9.145e3),
        (IPE, 'IPE120', 52.95e3, 8.646e3, 60.72e3, 13.58e3),
        (IPE, 'IPE140', 77.31e3, 12.3e3, 88.34e3, 19.24e3),
        (IPE, 'IPE160', 108.6e3, 16.66e3, 123.8e3, 26.09e3),
        (IPE, 'IPE180', 146.3e3, 22.16e3, 166.4e3, 34.59e3),
        (IPE, 'IPE200', 194.3e3, 28.47e3, 220.6e3, 44.61e3),
        (IPE, 'IPE220', 251.9e3, 37.25e3, 285.4e3, 58.11e3),
        (IPE, 'IPE240', 324.3e3, 47.27e3, 366.6e3, 73.92e3),
        (IPE, 'IPE270', 428.8e3, 62.2e3, 483.9e3, 96.95e3),
        (IPE, 'IPE300', 557e3, 80.5e3, 628.3e3, 125.2e3),
        (IPE, 'IPE330', 713.1e3, 98.51e3, 804.3e3, 153.6e3),
        (IPE, 'IPE360', 903.6e3, 122.7e3, 1019e3, 191e3),
        (IPE, 'IPE400', 1156e3, 146.4e3, 1307e3, 229e3),
        (IPE, 'IPE450', 1499e3, 176.4e3, 1701e3, 276.3e3),
        (IPE, 'IPE500', 1927e3, 214.1e3, 2194e3, 335.8e3),
        (IPE, 'IPE550', 2440e3, 254e3, 2787e3, 400.5e3),
        (IPE, 'IPE600', 3069e3, 307.9e3, 3512e3, 485.6e3),
        (HE, 'HEA100', 72.75e3, 26.76e3, 83.01e3, 41.14e3),
        (HE, 'HEA120', 106.3e3, 38.48e3, 119.5e3, 58.85e3),
        (HE, 'HEA140', 155.3e3, 55.61e3, 173.4e3, 84.84e3),
        (HE, 'HEA160', 220.1e3, 76.94e3, 245.1e3, 117.6e3),
        (HE, 'HEA180', 293.6e3, 102.7e3, 324.8e3, 156.4e3),
        (HE, 'HEA200', 388.6e3, 133.5e3, 429.4e3, 203.8e3),
        (HE, 'HEA220', 515.2e3, 177.6e3, 568.4e3, 270.5e3),
        (HE, 'HEA240', 675e3, 230.7e3, 744.6e3, 351.6e3),
        (HE, 'HEA260', 836.3e3, 282.1e3, 919.7e3, 430.1e3),
        (HE, 'HEA280', 1013e3, 340.2e3, 1112e3, 518.1e3),
        (HE, 'HEA300', 1259e3, 420.6e3, 1383e3, 641.1e3),
        (HE, 'HEA320', 1479e3, 465.6e3, 1628e3, 709.7e3),
        (HE, 'HEA340', 1678e3, 495.7e3, 1850e3, 755.9e3),
        (HE, 'HEA360', 1890e3, 525.7e3, 2088e3, 802.2e3),
        (HE, 'HEA400', 2311e3, 570.9e3, 2561e3, 872.8e3),
        (HE, 'HEA450', 2896e3, 631e3, 3215e3, 965.5e3),
        (HE, 'HEA500', 3549e3, 691.1e3, 3948e3, 1058e3),
        (HE, 'HEA550', 4145e3, 721.2e3, 4621e3, 1106e3),
        (HE, 'HEA600', 4787e3, 751.4e3, 5350e3, 1156e3),
        (HE, 'HEA650', 5474e3, 781.5e3, 6136e3, 1204e3),
        (HE, 'HEA700', 6240e3, 811.9e3, 7031e3, 1256e3),
        (HE, 'HEA800', 7682e3, 842.5e3, 8699e3, 1312e3),
        (HE, 'HEA900', 9484e3, 903.1e3, 10810e3, 1414e3),
        (HE, 'HEA1000', 11180e3, 933.6e3, 12820e3, 1469e3),
        (HE, 'HEB100', 89.9e3, 33.45e3, 104.2e3, 51.42e3),
        (HE, 'HEB120', 144e3, 52.92e3, 165.2e3, 80.96e3),
        (HE, 'HEB140', 215.6e3, 78.52e3, 245.4e3, 119.7e3),
        (HE, 'HEB160', 311.4e3, 111.1e3, 353.9e3, 169.9e3),
        (HE, 'HEB180', 425.6e3, 151.4e3, 481.4e3, 231e3),
        (HE, 'HEB200', 569.6e3, 200.3e3, 642.5e3, 305.8e3),
        (HE, 'HEB220', 735.5e3, 258.4e3, 827e3, 393.8e3),
        (HE, 'HEB240', 938.2e3, 326.8e3, 1053e3, 498.4e3),
        (HE, 'HEB260', 1147e3, 394.9e3, 1282e3, 602.2e3),
        (HE, 'HEB280', 1376e3, 471e3, 1534e3, 717.5e3),
        (HE, 'HEB300', 1677e3, 570.8e3, 1868e3, 870.1e3),
        (HE, 'HEB320', 1926e3, 615.9e3, 2149e3, 939e3),
        (HE, 'HEB340', 2156e3, 645.9e3, 2408e3, 985.7e3),
        (HE, 'HEB360', 2399e3, 676e3, 2682e3, 1032e3),
        (HE, 'HEB400', 2884e3, 721.2e3, 3231e3, 1104e3),
        (HE, 'HEB450', 3550e3, 781.4e3, 3982e3, 1197e3),
        (HE, 'HEB500', 4287e3, 841.6e3, 4815e3, 1292e3),
        (HE, 'HEB550', 4971e3, 871.8e3, 5591e3, 1341e3),
        (HE, 'HEB600', 5701e3, 902e3, 6425e3, 1391e3),
        (HE, 'HEB650', 6480e3, 932.2e3, 7319e3, 1441e3),
        (HE, 'HEB700', 7339e3, 962.7e3, 8327e3, 1495e3),
        (HE, 'HEB800', 8977e3, 993.5e3, 10220e3, 1553e3),
        (HE, 'HEB900', 10970e3, 1054e3, 12580e3, 1658e3),
        (HE, 'HEB1000', 12890e3, 1085e3, 14850e3, 1716e3),
        (HE, 'HEM100', 190.4e3, 75.31e3, 235.8e3, 116.3e3),
        (HE, 'HEM120', 288.2e3, 111.5e3, 350.6e3, 171.6e3),
        (HE, 'HEM140', 411.4e3, 156.7e3, 493.8e3, 240.5e3),
        (HE, 'HEM160', 566.4e3, 211.8e3, 674.5e3, 325.4e3),
        (HE, 'HEM180', 748.3e3, 277.4e3, 883.4e3, 425.1e3),
        (HE, 'HEM200', 967.4e3, 354.4e3, 1135e3, 543.2e3),
        (HE, 'HEM220', 1217e3, 443.5e3, 1419e3, 678.5e3),
        (HE, 'HEM240', 1799e3, 657.4e3, 2116e3, 1005e3),
        (HE, 'HEM260', 2159e3, 779.7e3, 2523e3, 1192e3),
        (HE, 'HEM280', 2551e3, 914e3, 2965e3, 1396e3),
        (HE, 'HEM300', 3482e3, 1251e3, 4077e3, 1913e3),
        (HE, 'HEM320', 3795e3, 1275e3, 4435e3, 1950e3),
        (HE, 'HEM340', 4051e3, 1275e3, 4717e3, 1952e3),
        (HE, 'HEM360', 4297e3, 1267e3, 4989e3, 1942e3),
        (HE, 'HEM400', 4820e3, 1260e3, 5571e3, 1934e3),
        (HE, 'HEM450', 5501e3, 1259e3, 6331e3, 1939e3),
        (HE, 'HEM500', 6180e3, 1251e3, 7094e3, 1932e3),
        (HE, 'HEM550', 6922e3, 1252e3, 7932e3, 1937e3),
        (HE, 'HEM600', 7659e3, 1244e3, 8772e3, 1930e3),
        (HE, 'HEM650', 8433e3, 1244e3, 9656e3, 1935e3),
        (HE, 'HEM700', 9197e3, 1236e3, 10530e3, 1928e3),
        (HE, 'HEM800', 10870e3, 1229e3, 12480e3, 1930e3),
        (HE, 'HEM900', 12530e3, 1221e3, 14440e3, 1928e3),
        (HE, 'HEM1000', 14330e3, 1222e3, 16560e3, 1939e3),
        (IPN, 'IPN80', 19.5e3, 3e3, 22.8e3, 5e3),
        (IPN, 'IPN100', 34.2e3, 4.88e3, 39.8e3, 8.1e3),
        (IPN, 'IPN120', 54.7e3, 7.41e3, 63.6e3, 12.4e3),
        (IPN, 'IPN140', 81.9e3, 10.7e3, 95.4e3, 17.9e3),
        (IPN, 'IPN160', 117e3, 14.8e3, 136e3, 24.9e3),
        (IPN, 'IPN180', 161e3, 19.8e3, 187e3, 33.2e3),
        (IPN, 'IPN200', 214e3, 26e3, 250e3, 43.5e3),
        (IPN, 'IPN220', 278e3, 33.1e3, 324e3, 55.7e3),
        (IPN, 'IPN240', 354e3, 41.7e3, 412e3, 70e3),
        (IPN, 'IPN260', 442e3, 51e3, 514e3, 85.9e3),
        (IPN, 'IPN280', 542e3, 61.2e3, 632e3, 103e3),
        (IPN, 'IPN300', 653e3, 72.2e3, 762e3, 121e3),
        (IPN, 'IPN320', 782e3, 84.7e3, 914e3, 143e3),
        (IPN, 'IPN340', 923e3, 98.4e3, 1080e3, 166e3),
        (IPN, 'IPN360', 1090e3, 114e3, 1276e3, 194e3),
        (IPN, 'IPN380', 1260e3, 131e3, 1482e3, 221e3),
        (IPN, 'IPN400', 1460e3, 149e3, 1714e3, 253e3),
        (IPN, 'IPN450', 2040e3, 203e3, 2400e3, 345e3),
        (IPN, 'IPN500', 2750e3, 268e3, 3240e3, 456e3),
        (IPN, 'IPN550', 3610e3, 349e3, 4240e3, 592e3),
        (IPN, 'IPN600', 4630e3, 434e3, 5452e3, 752e3),
        (UPN, 'UPN50', 10.6e3, 3.75e3, 13.10e3, 6.78e3),
        (UPN, 'UPN65', 17.7e3, 5.07e3, 21.7e3, 9.38e3),
        (UPN, 'UPN80', 26.5e3, 6.36e3, 32.3e3, 11.9e3),
        (UPN, 'UPN100', 41.2e3, 8.49e3, 49e3, 16.2e3),
        (UPN, 'UPN120', 60.7e3, 11.1e3, 72.6e3, 21.2e3),
        (UPN, 'UPN140', 86.4e3, 14.8e3, 103e3, 28.3e3),
        (UPN, 'UPN160', 116e3, 18.3e3, 138e3, 35.2e3),
        (UPN, 'UPN180', 150e3, 22.4e3, 179e3, 42.9e3),
        (UPN, 'UPN200', 191e3, 27e3, 228e3, 51.8e3),
        (UPN, 'UPN220', 245e3, 33.6e3, 292e3, 64.1e3),
        (UPN, 'UPN240', 300e3, 39.6e3, 358e3, 75.7e3),
        (UPN, 'UPN260', 371e3, 47.7e3, 442e3, 91.6e3),
        (UPN, 'UPN280', 448e3, 57.2e3, 532e3, 109e3),
        (UPN, 'UPN300', 535e3, 67.8e3, 632e3, 130e3),
        (UPN, 'UPN320', 679e3, 80.6e3, 826e3, 152e3),
        (UPN, 'UPN350', 734e3, 75e3, 918e3, 143e3),
        (UPN, 'UPN380', 829e3, 78.7e3, 1010e3, 148e3),
        (UPN, 'UPN400', 1020e3, 102e3, 1240e3, 190e3),
    ],
)
def test_profiles(cls, name, Wyel, Wzel, Wypl, Wzpl):
    """Test bending strength elastic and plastic."""
    Es = 206000
    fy = 355
    eps_su = 7e-2
    steel = UserDefined(
        x=[-eps_su, -fy / Es, 0, fy / Es, eps_su], y=[-fy, -fy, 0, fy, fy]
    )
    # Create geometry
    beam = cls(name)
    geo = CompoundGeometry([SurfaceGeometry(beam.polygon, steel)])

    # Compute expected values
    mye_expected = Wyel * fy * 1e-6
    mze_expected = Wzel * fy * 1e-6
    myp_expected = Wypl * fy * 1e-6
    mzp_expected = Wzpl * fy * 1e-6
    # Create the section with fiber
    sec = GenericSection(geo)
    # Elastic strength
    steel.set_ultimate_strain(fy / Es)
    results = sec.section_analyzer.calculate_bending_strength(theta=0, n=0)
    assert math.isclose(-results.m_y * 1e-6, mye_expected, rel_tol=2.5e-2)
    results = sec.section_analyzer.calculate_bending_strength(
        theta=math.pi / 2, n=0
    )
    assert math.isclose(-results.m_z * 1e-6, mze_expected, rel_tol=2.5e-2)
    # Plastic strength
    steel.set_ultimate_strain(0.07)
    results = sec.section_analyzer.calculate_bending_strength(theta=0, n=0)
    assert math.isclose(-results.m_y * 1e-6, myp_expected, rel_tol=2.5e-2)
    results = sec.section_analyzer.calculate_bending_strength(
        theta=math.pi / 2, n=0
    )
    assert math.isclose(-results.m_z * 1e-6, mzp_expected, rel_tol=2.5e-2)
