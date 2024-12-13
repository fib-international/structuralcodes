"""Tests for the marin coefficients from the constitutive_laws module."""

import math

import numpy as np
import pytest

from structuralcodes.materials.constitutive_laws import (
    BilinearCompression,
    Elastic,
)

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
    'fc, eps_c, eps_cu',
    [
        (-30, -0.002, -0.0035),
        (-20, -0.0015, -0.004),
        (-35, -0.002, -0.003),
        (-20000, -0.002, -0.0025),
        (-30000, -0.0025, -0.0035),
    ],
)
def test_marin_bilinearcompression_uniform_strain(eps_0, fc, eps_c, eps_cu):
    """Test the marin coefficients for a strain plane of pure tension or
    compression for BilinearCompression constitutive law.
    """
    law = BilinearCompression(fc, eps_c, eps_cu)

    # Check marin coefficients for stress integration
    if eps_0 > 0:
        expected = (None, (0.0,))
    elif eps_0 > eps_c:
        expected = (None, (fc / eps_c * eps_0, 0.0))
    elif eps_0 >= eps_cu:
        expected = (None, (fc,))
    else:
        expected = (None, (0.0,))

    strain, coeffs = law.__marin__([eps_0, 0])

    assert strain == expected[0]
    assert len(coeffs[0]) == len(expected[1])
    for i in range(len(coeffs[0])):
        assert math.isclose(coeffs[0][i], expected[1][i])

    # Check marin coefficients for modulus integration
    if eps_0 > 0:
        expected = (None, (0.0,))
    elif eps_0 > eps_c:
        expected = (None, (fc / eps_c,))
    else:
        expected = (None, (0.0,))

    strain, coeffs = law.__marin_tangent__([eps_0, 0])

    assert strain == expected[0]
    assert len(coeffs[0]) == len(expected[1])
    for i in range(len(coeffs[0])):
        assert math.isclose(coeffs[0][i], expected[1][i])


@pytest.mark.parametrize('eps_0', eps0)
@pytest.mark.parametrize(
    'eps_min', [(-0.01), (-0.005), (-0.004), (-0.003), (-0.002), (-0.001)]
)
@pytest.mark.parametrize(
    'fc, eps_c, eps_cu, h',
    [
        (-30, -0.002, -0.0035, 400),
        (-20, -0.0015, -0.004, 400),
        (-35, -0.002, -0.003, 350),
        (-20000, -0.002, -0.0025, 0.4),
        (-30000, -0.0025, -0.0035, 0.3),
    ],
)
def test_marin_bilinearcompression_axial_bending_strain(
    eps_0, eps_min, fc, eps_c, eps_cu, h
):
    """Test the marin coefficients for a strain plane of uniaxial bending
    with tension / compression for BilinearCompression constitutive law.
    """
    law = BilinearCompression(fc, eps_c, eps_cu)

    # Check marin coefficients for stress integration
    kappa_y = (eps_min - eps_0) * 2 / h

    # We are not checking here case with uniform compression / tension, skip
    if kappa_y == 0:
        return

    expected = (
        [(eps_c, 0), (eps_cu, eps_c)],
        [
            (fc / eps_c * eps_0, fc / eps_c * kappa_y),
            (fc,),
        ],
    )

    strain, coeffs = law.__marin__([eps_0, kappa_y])

    assert len(strain) == len(expected[0])
    assert len(expected[1]) == len(coeffs)

    for coeff, expect in zip(coeffs, expected[1]):
        assert len(coeff) == len(expect)
        for i in range(len(coeff)):
            assert math.isclose(coeff[i], expect[i])

    # Check marin coefficients for modulus integration
    expected = (
        [(eps_c, 0.0), (eps_cu, eps_c)],
        [
            (fc / eps_c,),
            (0.0,),
        ],
    )

    strain, coeffs = law.__marin_tangent__([eps_0, kappa_y])

    assert len(strain) == len(expected[0])
    assert len(expected[1]) == len(coeffs)

    for coeff, expect in zip(coeffs, expected[1]):
        assert len(coeff) == len(expect)
        for i in range(len(coeff)):
            assert math.isclose(coeff[i], expect[i])
