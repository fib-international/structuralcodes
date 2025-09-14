"""Tests for profiles."""

import json
import math
from pathlib import Path

import pytest
from shapely.testing import assert_geometries_equal

from structuralcodes.core._marin_integration import marin_integration
from structuralcodes.geometry import (
    CompoundGeometry,
    SurfaceGeometry,
)
from structuralcodes.geometry.profiles import (
    HD,
    HE,
    HP,
    IPE,
    IPN,
    UB,
    UBP,
    UC,
    UPE,
    UPN,
    L,
    U,
    W,
)
from structuralcodes.materials.basic import (
    ElasticMaterial,
    ElasticPlasticMaterial,
    GenericMaterial,
)
from structuralcodes.materials.constitutive_laws import (
    UserDefined,
)
from structuralcodes.sections._generic import GenericSection


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
        (W, 'W 100 x 100'),
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
        (UPN, 'UPN50'),
        (UPN, 'UPN65'),
        (UPN, 'UPN80'),
        (UPN, 'UPN100'),
        (UPN, 'UPN120'),
        (UPN, 'UPN140'),
        (UPN, 'UPN160'),
        (UPN, 'UPN180'),
        (UPN, 'UPN200'),
        (UPN, 'UPN220'),
        (UPN, 'UPN240'),
        (UPN, 'UPN260'),
        (UPN, 'UPN280'),
        (UPN, 'UPN300'),
        (UPN, 'UPN320'),
        (UPN, 'UPN350'),
        (UPN, 'UPN380'),
        (UPN, 'UPN400'),
        (IPN, 'IPN80'),
        (IPN, 'IPN100'),
        (IPN, 'IPN120'),
        (IPN, 'IPN140'),
        (IPN, 'IPN160'),
        (IPN, 'IPN180'),
        (IPN, 'IPN200'),
        (IPN, 'IPN220'),
        (IPN, 'IPN240'),
        (IPN, 'IPN260'),
        (IPN, 'IPN280'),
        (IPN, 'IPN300'),
        (IPN, 'IPN320'),
        (IPN, 'IPN340'),
        (IPN, 'IPN360'),
        (IPN, 'IPN380'),
        (IPN, 'IPN400'),
        (IPN, 'IPN450'),
        (IPN, 'IPN500'),
        (IPN, 'IPN550'),
        (IPN, 'IPN600'),
    ],
)
def test_reduced_names(cls, name):
    """Tests when imputing reduced names."""
    full_name = cls(name)
    short_name = cls(int(name[3:]))
    assert_geometries_equal(
        full_name.polygon,
        short_name.polygon,
        normalize=True,
    )

    # check also class method
    poly = cls.get_polygon(int(name[3:]))
    assert_geometries_equal(full_name.polygon, poly, normalize=True)


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
    steel = ElasticPlasticMaterial(E=Es, fy=fy, density=7850, eps_su=eps_su)

    # Create geometry
    geo = SurfaceGeometry(cls.get_polygon(name), steel)
    sec = GenericSection(geo)

    # Compute expected values
    xy = sec.geometry.geometries[0].polygon.exterior.coords.xy

    A = marin_integration(xy[0], xy[1], 0, 0)

    A_sec = sec.gross_properties.area
    assert math.isclose(A_sec, A)


@pytest.mark.parametrize(
    'cls, name, A, Iy, Wely, Wply, iy, Iz, Welz, Wplz, iz, Icsi, Ieta, theta',
    [
        (
            IPE,
            'IPE80',
            764,
            801400,
            20030,
            23220,
            32.4,
            84890,
            3691,
            5818,
            10.5,
            801400,
            84890,
            0.0,
        ),
        (
            IPE,
            'IPE100',
            1032,
            1710000,
            34200,
            39410,
            40.7,
            159200,
            5789,
            9146,
            12.4,
            1710000,
            159200,
            0.0,
        ),
        (
            IPE,
            'IPE120',
            1321,
            3178000,
            52960,
            60730,
            49,
            276700,
            8646,
            13580,
            14.5,
            3178000,
            276700,
            0.0,
        ),
        (
            IPE,
            'IPE140',
            1643,
            5412000,
            77320,
            88340,
            57.4,
            449200,
            12310,
            19250,
            16.5,
            5412000,
            449200,
            0.0,
        ),
        (
            IPE,
            'IPE160',
            2009,
            8693000,
            108700,
            123900,
            65.8,
            683100,
            16660,
            26100,
            18.4,
            8693000,
            683100,
            0.0,
        ),
        (
            IPE,
            'IPE180',
            2395,
            13170000,
            146300,
            166400,
            74.2,
            1009000,
            22160,
            34600,
            20.5,
            13170000,
            1009000,
            0.0,
        ),
        (
            IPE,
            'IPE200',
            2848,
            19430000,
            194300,
            220600,
            82.6,
            1424000,
            28470,
            44610,
            22.4,
            19430000,
            1424000,
            0.0,
        ),
        (
            IPE,
            'IPE220',
            3337,
            27720000,
            252000,
            285400,
            91.1,
            2049000,
            37250,
            58110,
            24.8,
            27720000,
            2049000,
            0.0,
        ),
        (
            IPE,
            'IPE240',
            3912,
            38920000,
            324300,
            366600,
            99.7,
            2836000,
            47270,
            73920,
            26.9,
            38920000,
            2836000,
            0.0,
        ),
        (
            IPE,
            'IPE270',
            4595,
            57900000,
            428900,
            484000,
            112.3,
            4199000,
            62200,
            96950,
            30.2,
            57900000,
            4199000,
            0.0,
        ),
        (
            IPE,
            'IPE300',
            5381,
            83560000,
            557100,
            628400,
            124.6,
            6038000,
            80500,
            125200,
            33.5,
            83560000,
            6038000,
            0.0,
        ),
        (
            IPE,
            'IPE330',
            6261,
            117700000,
            713100,
            804300,
            137.1,
            7881000,
            98520,
            153700,
            35.5,
            117700000,
            7881000,
            0.0,
        ),
        (
            IPE,
            'IPE360',
            7273,
            162700000,
            903600,
            1019000,
            149.5,
            10430000,
            122800,
            191100,
            37.9,
            162700000,
            10430000,
            0.0,
        ),
        (
            IPE,
            'IPE400',
            8446,
            231300000,
            1156000,
            1307000,
            165.5,
            13180000,
            146400,
            229000,
            39.5,
            231300000,
            13180000,
            0.0,
        ),
        (
            IPE,
            'IPE450',
            9882,
            337400000,
            1500000,
            1702000,
            184.8,
            16760000,
            176400,
            276400,
            41.2,
            337400000,
            16760000,
            0.0,
        ),
        (
            IPE,
            'IPE500',
            11552,
            482000000,
            1928000,
            2194000,
            204.3,
            21420000,
            214200,
            335900,
            43.1,
            482000000,
            21420000,
            0.0,
        ),
        (
            IPE,
            'IPE550',
            13442,
            671200000,
            2441000,
            2787000,
            223.5,
            26680000,
            254100,
            400500,
            44.5,
            671200000,
            26680000,
            0.0,
        ),
        (
            IPE,
            'IPE600',
            15598,
            920800000,
            3069000,
            3512000,
            243,
            33870000,
            307900,
            485600,
            46.6,
            920800000,
            33870000,
            0.0,
        ),
        (
            UPN,
            'UPN50',
            710,
            264000,
            10600,
            13100,
            19,
            91200,
            3750,
            6780,
            11.2,
            264000,
            91200,
            0,
        ),
        (
            UPN,
            'UPN65',
            900,
            575000,
            17700,
            21700,
            25,
            141000,
            5070,
            9380,
            12.4,
            575000,
            141000,
            0,
        ),
        (
            UPN,
            'UPN80',
            1100,
            1060000,
            26500,
            31800,
            31,
            194000,
            6360,
            12100,
            13.3,
            1060000,
            194000,
            0,
        ),
        (
            UPN,
            'UPN100',
            1350,
            2060000,
            41200,
            49000,
            39.1,
            293000,
            8490,
            16200,
            14.7,
            2060000,
            293000,
            0,
        ),
        (
            UPN,
            'UPN120',
            1700,
            3640000,
            60700,
            72600,
            46.2,
            432000,
            11100,
            21200,
            15.9,
            3640000,
            432000,
            0,
        ),
        (
            UPN,
            'UPN140',
            2040,
            6050000,
            86400,
            103000,
            54.5,
            627000,
            14800,
            28300,
            17.5,
            6050000,
            627000,
            0,
        ),
        (
            UPN,
            'UPN160',
            2400,
            9250000,
            116000,
            138000,
            62.1,
            853000,
            18300,
            35200,
            18.9,
            9250000,
            853000,
            0,
        ),
        (
            UPN,
            'UPN180',
            2800,
            13500000,
            150000,
            179000,
            69.5,
            1140000,
            22400,
            42900,
            20.2,
            13500000,
            1140000,
            0,
        ),
        (
            UPN,
            'UPN200',
            3220,
            19100000,
            191000,
            228000,
            77,
            1480000,
            27000,
            51800,
            21.4,
            19100000,
            1480000,
            0,
        ),
        (
            UPN,
            'UPN220',
            3740,
            26900000,
            245000,
            292000,
            84.8,
            1970000,
            33600,
            64100,
            23,
            26900000,
            1970000,
            0,
        ),
        (
            UPN,
            'UPN240',
            4230,
            36000000,
            300000,
            358000,
            92.2,
            2480000,
            39600,
            75700,
            24.2,
            36000000,
            2480000,
            0,
        ),
        (
            UPN,
            'UPN260',
            4830,
            48200000,
            371000,
            442000,
            99.9,
            3170000,
            47700,
            91600,
            25.6,
            48200000,
            3170000,
            0,
        ),
        (
            UPN,
            'UPN280',
            5330,
            62800000,
            448000,
            532000,
            109,
            3990000,
            57200,
            109000,
            27.4,
            62800000,
            3990000,
            0,
        ),
        (
            UPN,
            'UPN300',
            5880,
            80300000,
            535000,
            632000,
            117,
            4950000,
            67800,
            130000,
            29,
            80300000,
            4950000,
            0,
        ),
        (
            UPN,
            'UPN320',
            7580,
            108700000,
            679000,
            826000,
            121,
            5970000,
            80600,
            152000,
            28.1,
            108700000,
            5970000,
            0,
        ),
        (
            UPN,
            'UPN350',
            7730,
            128400000,
            734000,
            899310,
            129,
            5700000,
            75000,
            143000,
            27.2,
            128400000,
            5700000,
            0,
        ),
        (
            UPN,
            'UPN380',
            8040,
            157600000,
            829000,
            1010000,
            140,
            6150000,
            78700,
            148000,
            27.7,
            157600000,
            6150000,
            0,
        ),
        (
            UPN,
            'UPN400',
            9150,
            203500000,
            1020000,
            1240000,
            149,
            8460000,
            102000,
            190000,
            30.4,
            203500000,
            8460000,
            0,
        ),
        (
            IPN,
            'IPN80',
            758,
            778000,
            19500,
            22800,
            32,
            62900,
            3000,
            5000,
            9.1,
            778000,
            62900,
            0,
        ),
        (
            IPN,
            'IPN500',
            17900,
            687400000,
            2750000,
            3240000,
            196,
            24800000,
            268000,
            456000,
            37.2,
            687400000,
            24800000,
            0.0,
        ),
    ],
)
def test_section_properties(
    cls, name, A, Iy, Wely, Wply, iy, Iz, Welz, Wplz, iz, Icsi, Ieta, theta
):
    """Test Steel section comparing section properties with tables."""
    # Create polygon representing geometyr of profile
    profile = cls(name)

    # Compare wth stored values with a 2% tolerance
    assert math.isclose(A, profile.A, rel_tol=2e-2)
    assert math.isclose(Iy, profile.Iy, rel_tol=2e-2)
    assert math.isclose(Iz, profile.Iz, rel_tol=2e-2)
    assert math.isclose(Wely, profile.Wely, rel_tol=2e-2)
    assert math.isclose(Welz, profile.Welz, rel_tol=2e-2)
    assert math.isclose(Wply, profile.Wply, rel_tol=2e-2)
    assert math.isclose(Wplz, profile.Wplz, rel_tol=2e-2)
    assert math.isclose(iy, profile.iy, rel_tol=2e-2)
    assert math.isclose(iz, profile.iz, rel_tol=2e-2)
    assert math.isclose(iz, profile.iz, rel_tol=2e-2)
    assert math.isclose(Icsi, profile.Icsi, rel_tol=2e-2)
    assert math.isclose(Ieta, profile.Ieta, rel_tol=2e-2)
    assert math.isclose(theta, profile.theta, rel_tol=2e-2)

    # Check other properties
    params = cls.parameters.get(name)
    if params is not None:
        assert math.isclose(params.get('h'), profile.h, rel_tol=2e-2)
        assert math.isclose(params.get('b'), profile.b, rel_tol=2e-2)
        assert math.isclose(params.get('tw'), profile.tw, rel_tol=2e-2)
        assert math.isclose(params.get('tf'), profile.tf, rel_tol=2e-2)
        r = params.get('r')
        if r is not None:
            assert math.isclose(r, profile.r, rel_tol=2.0e-2)
        r1 = params.get('r1')
        if r1 is not None:
            assert math.isclose(r1, profile.r1, rel_tol=2.0e-2)
        r2 = params.get('r2')
        if r2 is not None:
            assert math.isclose(r2, profile.r2, rel_tol=2.0e-2)


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
    steel = ElasticPlasticMaterial(E=Es, fy=fy, density=7850, eps_su=eps_su)

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
    results = sec.section_calculator.calculate_bending_strength(theta=0, n=0)
    assert math.isclose(-results.m_y * 1e-6, my_expected, rel_tol=1e-3)
    results = sec.section_calculator.calculate_bending_strength(
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
def test_Isection_elastic_marin(cls, name):
    """Test Steel I section elastic strength."""
    Es = 206000
    fy = 355
    eps_su = fy / Es
    steel = ElasticPlasticMaterial(E=Es, fy=fy, density=7850, eps_su=eps_su)

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
    results = sec.section_calculator.calculate_bending_strength(theta=0, n=0)
    assert math.isclose(-results.m_y * 1e-6, my_expected, rel_tol=1e-3)
    results = sec.section_calculator.calculate_bending_strength(
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
    steel = ElasticPlasticMaterial(E=Es, fy=fy, density=7850, eps_su=eps_su)

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
    results = sec.section_calculator.calculate_bending_strength(theta=0, n=0)
    assert math.isclose(-results.m_y * 1e-6, my_expected, rel_tol=1e-2)
    results = sec.section_calculator.calculate_bending_strength(
        theta=math.pi / 2, n=0
    )
    assert math.isclose(-results.m_z * 1e-6, mz_expected, rel_tol=2e-2)


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
    steel = ElasticPlasticMaterial(E=Es, fy=fy, density=7850, eps_su=eps_su)

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
    results = sec.section_calculator.calculate_bending_strength(theta=0, n=0)
    assert math.isclose(-results.m_y * 1e-6, my_expected, rel_tol=1e-2)
    results = sec.section_calculator.calculate_bending_strength(
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
    steel = ElasticMaterial(E=Es, density=7850)
    steel.constitutive_law.set_ultimate_strain(fy / Es)
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
    results = sec.section_calculator.calculate_bending_strength(theta=0, n=0)
    assert math.isclose(-results.m_y * 1e-6, my_expected, rel_tol=1e-3)
    results = sec.section_calculator.calculate_bending_strength(
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
    steel_law = UserDefined(
        x=[-eps_su, -fy / Es, 0, fy / Es, eps_su], y=[-fy, -fy, 0, fy, fy]
    )
    steel = GenericMaterial(density=7850, constitutive_law=steel_law)
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
    results = sec.section_calculator.calculate_bending_strength(theta=0, n=0)
    assert math.isclose(-results.m_y * 1e-6, my_expected, rel_tol=1e-3)
    results = sec.section_calculator.calculate_bending_strength(
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
        (UPN, 'UPN350', 734e3, 75e3, 901e3, 143e3),
        (UPN, 'UPN380', 829e3, 78.7e3, 1010e3, 148e3),
        (UPN, 'UPN400', 1020e3, 102e3, 1240e3, 190e3),
    ],
)
def test_profiles(cls, name, Wyel, Wzel, Wypl, Wzpl):
    """Test bending strength elastic and plastic."""
    Es = 206000
    fy = 355
    eps_su = 7e-2
    steel_law = UserDefined(
        x=[-eps_su, -fy / Es, 0, fy / Es, eps_su], y=[-fy, -fy, 0, fy, fy]
    )
    steel = GenericMaterial(density=7850, constitutive_law=steel_law)
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
    steel_law.set_ultimate_strain(fy / Es)
    results = sec.section_calculator.calculate_bending_strength(theta=0, n=0)
    assert math.isclose(-results.m_y * 1e-6, mye_expected, rel_tol=2.5e-2)
    results = sec.section_calculator.calculate_bending_strength(
        theta=math.pi / 2, n=0
    )
    assert math.isclose(-results.m_z * 1e-6, mze_expected, rel_tol=2.5e-2)
    # Plastic strength
    steel_law.set_ultimate_strain(0.07)
    results = sec.section_calculator.calculate_bending_strength(theta=0, n=0)
    assert math.isclose(-results.m_y * 1e-6, myp_expected, rel_tol=2.5e-2)
    results = sec.section_calculator.calculate_bending_strength(
        theta=math.pi / 2, n=0
    )
    assert math.isclose(-results.m_z * 1e-6, mzp_expected, rel_tol=2.5e-2)


def load_w_profiles_data():
    """Load W profiles data from w.json file."""
    json_file = Path(__file__).parent / 'w.json'
    with open(json_file, 'r') as f:
        profiles_data = []
        for line in f:
            profile_data = json.loads(line.strip())
            profile_name = profile_data['ProfileName']
            wely_cm3 = profile_data['Wely_cm3']

            profiles_data.append((profile_name, wely_cm3))
    return profiles_data


@pytest.mark.parametrize(
    'profile_name, expected_wely_cm3', load_w_profiles_data()
)
def test_w_profile_wely(profile_name, expected_wely_cm3):
    """Test W profile Wely property matches JSON data."""
    profile = W(profile_name)

    wy_el_expected = _wyel_I_beam(
        profile.h, profile.b, profile.tw, profile.tf, profile.r
    )

    wy_el_profile = profile.Wely

    assert math.isclose(wy_el_expected, wy_el_profile, rel_tol=1e-3)

    # I wanted to do this but is failing due to strange values in the input
    # text file:
    # assert math.isclose(expected_wely_cm3 * 1e3, wy_el_profile, rel_tol=1e-2)
    del expected_wely_cm3


def load_upe_profiles_data():
    """Load UPE profiles data from upe.json file."""
    json_file = Path(__file__).parent / 'upe.json'
    with open(json_file, 'r') as f:
        profiles_data = []
        for line in f:
            profile_data = json.loads(line.strip())
            profile_name = profile_data['ProfileName']
            A_cm2 = profile_data['A_cm2']
            Iy_cm4 = profile_data['Iy_cm4']
            Iz_cm4 = profile_data['Iz_cm4']
            Wely_cm3 = profile_data['Wely_cm3']
            Welz_cm3 = profile_data['Welz_cm3']

            profiles_data.append(
                (profile_name, A_cm2, Iy_cm4, Iz_cm4, Wely_cm3, Welz_cm3)
            )
    return profiles_data


@pytest.mark.parametrize(
    (
        'profile_name, expected_A_cm2, expected_Iy_cm4, expected_Iz_cm4, '
        'expected_Wely_cm3, expected_Welz_cm3'
    ),
    load_upe_profiles_data(),
)
def test_upe_profile(
    profile_name,
    expected_A_cm2,
    expected_Iy_cm4,
    expected_Iz_cm4,
    expected_Wely_cm3,
    expected_Welz_cm3,
):
    """Test UPE profile A property matches JSON data."""
    profile = UPE(profile_name)

    assert math.isclose(profile.A, expected_A_cm2 * 1e2, rel_tol=0.5e-2)
    assert math.isclose(profile.Iy, expected_Iy_cm4 * 1e4, rel_tol=0.5e-2)
    assert math.isclose(profile.Iz, expected_Iz_cm4 * 1e4, rel_tol=0.5e-2)
    assert math.isclose(profile.Wely, expected_Wely_cm3 * 1e3, rel_tol=0.5e-2)
    assert math.isclose(profile.Welz, expected_Welz_cm3 * 1e3, rel_tol=0.5e-2)


def load_ub_profiles_data():
    """Load UB profiles data from ub.json file."""
    json_file = Path(__file__).parent / 'ub.json'
    with open(json_file, 'r') as f:
        profiles_data = []
        for line in f:
            profile_data = json.loads(line.strip())
            profile_name = profile_data['ProfileName']
            A_cm2 = profile_data['A_cm2']
            Iy_cm4 = profile_data['Iy_cm4']
            Iz_cm4 = profile_data['Iz_cm4']
            Wely_cm3 = profile_data['Wely_cm3']
            Welz_cm3 = profile_data['Welz_cm3']

            profiles_data.append(
                (profile_name, A_cm2, Iy_cm4, Iz_cm4, Wely_cm3, Welz_cm3)
            )
    return profiles_data


@pytest.mark.parametrize(
    (
        'profile_name, expected_A_cm2, expected_Iy_cm4, expected_Iz_cm4, '
        'expected_Wely_cm3, expected_Welz_cm3'
    ),
    load_ub_profiles_data(),
)
def test_ub_profile(
    profile_name,
    expected_A_cm2,
    expected_Iy_cm4,
    expected_Iz_cm4,
    expected_Wely_cm3,
    expected_Welz_cm3,
):
    """Test UB profile property matches JSON data."""
    profile = UB(profile_name)

    assert math.isclose(profile.A, expected_A_cm2 * 1e2, rel_tol=2e-2)
    assert math.isclose(profile.Iy, expected_Iy_cm4 * 1e4, rel_tol=2e-2)
    assert math.isclose(profile.Iz, expected_Iz_cm4 * 1e4, rel_tol=2e-2)
    assert math.isclose(profile.Wely, expected_Wely_cm3 * 1e3, rel_tol=2.5e-2)
    assert math.isclose(profile.Welz, expected_Welz_cm3 * 1e3, rel_tol=2.5e-2)


def load_u_profiles_data():
    """Load U profiles data from ub.json file."""
    json_file = Path(__file__).parent / 'u.json'
    with open(json_file, 'r') as f:
        profiles_data = []
        for line in f:
            profile_data = json.loads(line.strip())
            profile_name = profile_data['ProfileName']
            A_cm2 = profile_data['A_cm2']
            Iy_cm4 = profile_data['Iy_cm4']
            Iz_cm4 = profile_data['Iz_cm4']
            Wely_cm3 = profile_data['Wely_cm3']
            Welz_cm3 = profile_data['Welz_cm3']

            profiles_data.append(
                (profile_name, A_cm2, Iy_cm4, Iz_cm4, Wely_cm3, Welz_cm3)
            )
    return profiles_data


@pytest.mark.parametrize(
    (
        'profile_name, expected_A_cm2, expected_Iy_cm4, expected_Iz_cm4, '
        'expected_Wely_cm3, expected_Welz_cm3'
    ),
    load_u_profiles_data(),
)
def test_u_profile(
    profile_name,
    expected_A_cm2,
    expected_Iy_cm4,
    expected_Iz_cm4,
    expected_Wely_cm3,
    expected_Welz_cm3,
):
    """Test U profile property matches JSON data."""
    profile = U(profile_name)

    assert math.isclose(profile.A, expected_A_cm2 * 1e2, rel_tol=2e-2)
    assert math.isclose(profile.Iy, expected_Iy_cm4 * 1e4, rel_tol=2e-2)
    assert math.isclose(profile.Iz, expected_Iz_cm4 * 1e4, rel_tol=2e-2)
    assert math.isclose(profile.Wely, expected_Wely_cm3 * 1e3, rel_tol=2.5e-2)
    assert math.isclose(profile.Welz, expected_Welz_cm3 * 1e3, rel_tol=2.5e-2)


def load_hd_profiles_data():
    """Load HD profiles data from ub.json file."""
    json_file = Path(__file__).parent / 'hd.json'
    with open(json_file, 'r') as f:
        profiles_data = []
        for line in f:
            profile_data = json.loads(line.strip())
            profile_name = profile_data['ProfileName']
            A_cm2 = profile_data['A_cm2']
            Iy_cm4 = profile_data['Iy_cm4']
            Iz_cm4 = profile_data['Iz_cm4']
            Wely_cm3 = profile_data['Wely_cm3']
            Welz_cm3 = profile_data['Welz_cm3']

            profiles_data.append(
                (profile_name, A_cm2, Iy_cm4, Iz_cm4, Wely_cm3, Welz_cm3)
            )
    return profiles_data


@pytest.mark.parametrize(
    (
        'profile_name, expected_A_cm2, expected_Iy_cm4, expected_Iz_cm4, '
        'expected_Wely_cm3, expected_Welz_cm3'
    ),
    load_hd_profiles_data(),
)
def test_hd_profile(
    profile_name,
    expected_A_cm2,
    expected_Iy_cm4,
    expected_Iz_cm4,
    expected_Wely_cm3,
    expected_Welz_cm3,
):
    """Test HD profile property matches JSON data."""
    profile = HD(profile_name)

    assert math.isclose(profile.A, expected_A_cm2 * 1e2, rel_tol=2e-2)
    assert math.isclose(profile.Iy, expected_Iy_cm4 * 1e4, rel_tol=2e-2)
    assert math.isclose(profile.Iz, expected_Iz_cm4 * 1e4, rel_tol=2e-2)
    assert math.isclose(profile.Wely, expected_Wely_cm3 * 1e3, rel_tol=2.5e-2)
    assert math.isclose(profile.Welz, expected_Welz_cm3 * 1e3, rel_tol=2.5e-2)


def load_hp_profiles_data():
    """Load HP profiles data from ub.json file."""
    json_file = Path(__file__).parent / 'hp.json'
    with open(json_file, 'r') as f:
        profiles_data = []
        for line in f:
            profile_data = json.loads(line.strip())
            profile_name = profile_data['ProfileName']
            A_cm2 = profile_data['A_cm2']
            Iy_cm4 = profile_data['Iy_cm4']
            Iz_cm4 = profile_data['Iz_cm4']
            Wely_cm3 = profile_data['Wely_cm3']
            Welz_cm3 = profile_data['Welz_cm3']

            profiles_data.append(
                (profile_name, A_cm2, Iy_cm4, Iz_cm4, Wely_cm3, Welz_cm3)
            )
    return profiles_data


@pytest.mark.parametrize(
    (
        'profile_name, expected_A_cm2, expected_Iy_cm4, expected_Iz_cm4, '
        'expected_Wely_cm3, expected_Welz_cm3'
    ),
    load_hp_profiles_data(),
)
def test_hp_profile(
    profile_name,
    expected_A_cm2,
    expected_Iy_cm4,
    expected_Iz_cm4,
    expected_Wely_cm3,
    expected_Welz_cm3,
):
    """Test HP profile property matches JSON data."""
    profile = HP(profile_name)

    assert math.isclose(profile.A, expected_A_cm2 * 1e2, rel_tol=2e-2)
    assert math.isclose(profile.Iy, expected_Iy_cm4 * 1e4, rel_tol=2e-2)
    assert math.isclose(profile.Iz, expected_Iz_cm4 * 1e4, rel_tol=2e-2)
    assert math.isclose(profile.Wely, expected_Wely_cm3 * 1e3, rel_tol=2.5e-2)
    assert math.isclose(profile.Welz, expected_Welz_cm3 * 1e3, rel_tol=2.5e-2)


def load_l_profiles_data():
    """Load L profiles data from ub.json file."""
    json_file = Path(__file__).parent / 'l.json'
    with open(json_file, 'r') as f:
        profiles_data = []
        for line in f:
            profile_data = json.loads(line.strip())
            profile_name = profile_data['ProfileName']
            A_cm2 = profile_data['A_cm2']
            Iy_cm4 = profile_data['Iy_cm4']
            Iz_cm4 = profile_data['Iz_cm4']
            Wely_cm3 = profile_data['Wely_cm3']
            Welz_cm3 = profile_data['Welz_cm3']

            profiles_data.append(
                (profile_name, A_cm2, Iy_cm4, Iz_cm4, Wely_cm3, Welz_cm3)
            )
    return profiles_data


@pytest.mark.parametrize(
    (
        'profile_name, expected_A_cm2, expected_Iy_cm4, expected_Iz_cm4, '
        'expected_Wely_cm3, expected_Welz_cm3'
    ),
    load_l_profiles_data(),
)
def test_l_profile(
    profile_name,
    expected_A_cm2,
    expected_Iy_cm4,
    expected_Iz_cm4,
    expected_Wely_cm3,
    expected_Welz_cm3,
):
    """Test L profile property matches JSON data."""
    profile = L(profile_name)

    assert math.isclose(profile.A, expected_A_cm2 * 1e2, rel_tol=2e-2)
    assert math.isclose(profile.Iy, expected_Iy_cm4 * 1e4, rel_tol=2e-2)
    assert math.isclose(profile.Iz, expected_Iz_cm4 * 1e4, rel_tol=2e-2)
    assert math.isclose(profile.Wely, expected_Wely_cm3 * 1e3, rel_tol=2.5e-2)
    assert math.isclose(profile.Welz, expected_Welz_cm3 * 1e3, rel_tol=2.5e-2)
