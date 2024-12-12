"""Tests for the marin coefficients from the constitutive_laws module."""

import math

import numpy as np
import pytest

from structuralcodes.materials.constitutive_laws import Elastic, ElasticPlastic

eps0 = (
    [(x) for x in np.linspace(0, 0.01, 5)]
    + [(x) for x in np.linspace(0.02, 0.08, 3)]
    + [(x) for x in np.linspace(-0.001, -0.005, 5)]
    + [(x) for x in np.linspace(-0.01, -0.09, 3)]
)


@pytest.mark.parametrize('eps_0', eps0)
@pytest.mark.parametrize(
    'E',
    [(20000), (30000), (210000), (2100000)],
)
def test_marin_elastic_uniform_strain(eps_0, E):
    """Test the marin coefficients for a strain plane of pure tension or
    compression for Elastic constitutive law.
    """
    law = Elastic(E)

    # Check marin coefficients for stress integration
    expected = (None, (E * eps_0, 0.0))

    strain, coeffs = law.__marin__([eps_0, 0])

    assert strain == expected[0]
    assert len(coeffs[0]) == 2
    assert math.isclose(coeffs[0][0], expected[1][0])
    assert math.isclose(coeffs[0][1], expected[1][1])

    # Check marin coefficients for modulus integration
    expected = (None, (E,))

    strain, coeffs = law.__marin_tangent__([eps_0, 0])

    assert strain == expected[0]
    assert len(coeffs[0]) == 1
    assert math.isclose(coeffs[0][0], expected[1][0])


@pytest.mark.parametrize(
    'eps_max', [(0.08), (0.01), (0.005), (0.004), (0.003), (0.002), (0.001)]
)
@pytest.mark.parametrize(
    'E, h',
    [(20000, 400), (30000, 400), (210000, 400), (2100000, 0.4)],
)
def test_marin_elastic_bending_strain(eps_max, E, h):
    """Test the marin coefficients for a strain plane of uniaxial bending
    for Elastic constitutive law.
    """
    law = Elastic(E)

    # Check marin coefficients for stress integration
    kappa_y = -eps_max * 2 / h
    expected = (None, (0.0, E * kappa_y))

    strain, coeffs = law.__marin__([0, kappa_y])

    assert strain == expected[0]
    assert len(coeffs[0]) == 2
    assert math.isclose(coeffs[0][0], expected[1][0])
    assert math.isclose(coeffs[0][1], expected[1][1])

    # Check marin coefficients for modulus integration
    expected = (None, (E,))

    strain, coeffs = law.__marin_tangent__([0, kappa_y])

    assert strain == expected[0]
    assert len(coeffs[0]) == 1
    assert math.isclose(coeffs[0][0], expected[1][0])


@pytest.mark.parametrize('eps_0', eps0)
@pytest.mark.parametrize(
    'eps_min', [(-0.01), (-0.005), (-0.004), (-0.003), (-0.002), (-0.001)]
)
@pytest.mark.parametrize(
    'E, h',
    [(20000, 400), (30000, 400), (210000, 400), (2100000, 0.4)],
)
def test_marin_elastic_axial_bending_strain(eps_0, eps_min, E, h):
    """Test the marin coefficients for a strain plane of uniaxial bending
    with tension / compression for Elastic constitutive law.
    """
    law = Elastic(E)

    # Check marin coefficients for stress integration
    kappa_y = (eps_min - eps_0) * 2 / h
    expected = (None, (E * eps_0, E * kappa_y))

    strain, coeffs = law.__marin__([eps_0, kappa_y])

    assert strain == expected[0]
    assert len(coeffs[0]) == 2
    assert math.isclose(coeffs[0][0], expected[1][0])
    assert math.isclose(coeffs[0][1], expected[1][1])

    # Check marin coefficients for modulus integration
    expected = (None, (E,))

    strain, coeffs = law.__marin_tangent__([eps_0, kappa_y])

    assert strain == expected[0]
    assert len(coeffs[0]) == 1
    assert math.isclose(coeffs[0][0], expected[1][0])


@pytest.mark.parametrize('eps_0', eps0)
@pytest.mark.parametrize(
    'E, fy, Eh, eps_su',
    [
        (200000, 450, 0.0, None),
        (200000, 500, 2000.0, None),
        (200000, 500, 0.0, 0.0675),
        (200000000, 500000, 0.0, 0.0675),
        (200000000, 500000, 0.0, None),
    ],
)
def test_marin_elasticplastic_uniform_strain(eps_0, E, fy, Eh, eps_su):
    """Test the marin coefficients for a strain plane of pure tension or
    compression for ElasticPlastic constitutive law.
    """
    law = ElasticPlastic(E, fy, Eh, eps_su)

    eps_y = law._eps_sy
    _, eps_su_p = law.get_ultimate_strain()

    delta_sigma = fy * (1 - Eh / E)

    # Check marin coefficients for stress integration
    sign = 1 if eps_0 >= 0 else -1
    if abs(eps_0) <= eps_y:
        expected = (None, (E * eps_0, 0.0))
    elif abs(eps_0) <= eps_su_p:
        expected = (None, (Eh * eps_0 + sign * delta_sigma, 0.0))
    else:
        expected = (None, (0.0, 0.0))

    strain, coeffs = law.__marin__([eps_0, 0])

    assert strain == expected[0]
    assert len(coeffs[0]) == 2
    assert math.isclose(coeffs[0][0], expected[1][0])
    assert math.isclose(coeffs[0][1], expected[1][1])

    # Check marin coefficients for modulus integration
    if abs(eps_0) <= eps_y:
        expected = (None, (E,))
    elif abs(eps_0) <= eps_su_p:
        expected = (None, (Eh,))
    else:
        expected = (None, (0.0,))

    strain, coeffs = law.__marin_tangent__([eps_0, 0])

    assert strain == expected[0]
    assert len(coeffs[0]) == 1
    assert math.isclose(coeffs[0][0], expected[1][0])


@pytest.mark.parametrize('eps_0', eps0)
@pytest.mark.parametrize(
    'eps_min', [(-0.01), (-0.005), (-0.004), (-0.003), (-0.002), (-0.001)]
)
@pytest.mark.parametrize(
    'E, fy, Eh, eps_su, h',
    [
        (200000, 450, 0.0, None, 400),
        (200000, 500, 2000.0, None, 400),
        (200000, 500, 0.0, 0.0675, 300),
        (200000000, 500000, 0.0, 0.0675, 0.4),
        (200000000, 500000, 0.0, None, 0.3),
    ],
)
def test_marin_elasticplastic_axial_bending_strain(
    eps_0, eps_min, E, fy, Eh, eps_su, h
):
    """Test the marin coefficients for a strain plane of uniaxial bending
    with tension / compression for ElasticPlastic constitutive law.
    """
    law = ElasticPlastic(E, fy, Eh, eps_su)

    eps_y = law._eps_sy
    eps_su_n, eps_su_p = law.get_ultimate_strain()

    delta_sigma = fy * (1 - Eh / E)

    # Check marin coefficients for stress integration
    kappa_y = (eps_min - eps_0) * 2 / h

    # We are not checking here case with uniform compression / tension, skip
    if kappa_y == 0:
        return

    expected = (
        [(eps_su_n, -eps_y), (-eps_y, eps_y), (eps_y, eps_su_p)],
        [
            (Eh * eps_0 - delta_sigma, Eh * kappa_y),
            (E * eps_0, E * kappa_y),
            (Eh * eps_0 + delta_sigma, Eh * kappa_y),
        ],
    )

    strain, coeffs = law.__marin__([eps_0, kappa_y])

    assert len(strain) == len(expected[0])
    assert len(expected[1]) == len(coeffs)

    for coeff, expect in zip(coeffs, expected[1]):
        assert len(coeff) == len(expect)
        assert math.isclose(coeff[0], expect[0])
        assert math.isclose(coeff[1], expect[1])

    # Check marin coefficients for modulus integration
    expected = (
        [(eps_su_n, -eps_y), (-eps_y, eps_y), (eps_y, eps_su_p)],
        [
            (Eh,),
            (E,),
            (Eh,),
        ],
    )

    strain, coeffs = law.__marin_tangent__([eps_0, kappa_y])

    assert len(strain) == len(expected[0])
    assert len(expected[1]) == len(coeffs)

    for coeff, expect in zip(coeffs, expected[1]):
        assert len(coeff) == len(expect)
        assert math.isclose(coeff[0], expect[0])
