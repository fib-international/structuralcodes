"""Tests for the marin coefficients from the constitutive_laws module."""

import math

import numpy as np
import pytest

from structuralcodes.materials.constitutive_laws import (
    BilinearCompression,
    Elastic,
    ElasticPlastic,
    InitStrain,
    ParabolaRectangle,
    UserDefined,
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
        expected = (None, (0.0,))

    strain, coeffs = law.__marin__([eps_0, 0])

    assert strain == expected[0]
    assert len(coeffs[0]) == len(expected[1])
    for i in range(len(coeffs[0])):
        assert math.isclose(coeffs[0][i], expected[1][i])

    # Check marin coefficients for modulus integration
    if abs(eps_0) <= eps_y:
        expected = (None, (E,))
    elif abs(eps_0) <= eps_su_p:
        expected = (None, (Eh,))
    else:
        expected = (None, (0.0,))

    strain, coeffs = law.__marin_tangent__([eps_0, 0])

    assert strain == expected[0]
    assert len(coeffs[0]) == len(expected[1])
    for i in range(len(coeffs[0])):
        assert math.isclose(coeffs[0][i], expected[1][i])


@pytest.mark.parametrize('eps_0', eps0)
@pytest.mark.parametrize('init_strain', [0.0, 0.001])
@pytest.mark.parametrize(
    'E, fy, Eh, eps_su',
    [
        (200000, 450, 0.0, 0.0675),
        (200000, 500, 2000.0, 0.0675),
        (200000, 500, 0.0, 0.0675),
    ],
)
def test_marin_init_strain_uniform_strain(
    eps_0, E, fy, Eh, eps_su, init_strain
):
    """Test the marin coefficients for a strain plane of pure tension or
    compression for InitialStrain constitutive law.
    """
    base_law = ElasticPlastic(E, fy, Eh, eps_su)
    law = InitStrain(base_law, init_strain)

    eps_y = base_law._eps_sy
    _, eps_su_p = law.get_ultimate_strain()

    delta_sigma = fy * (1 - Eh / E)

    # Check marin coefficients for stress integration
    sign = 1 if eps_0 >= 0 else -1
    if abs(eps_0 + init_strain) <= eps_y:
        expected = (None, (E * (eps_0 + init_strain), 0.0))
    elif abs(eps_0 + init_strain) <= eps_su_p:
        expected = (
            None,
            (Eh * (eps_0 + init_strain) + sign * delta_sigma, 0.0),
        )
    else:
        expected = (None, (0.0,))

    strain, coeffs = law.__marin__([eps_0, 0])

    assert strain == expected[0]
    assert len(coeffs[0]) == len(expected[1])
    for i in range(len(coeffs[0])):
        assert math.isclose(coeffs[0][i], expected[1][i])

    # Check marin coefficients for modulus integration
    if abs(eps_0 + init_strain) <= eps_y:
        expected = (None, (E,))
    elif abs(eps_0 + init_strain) <= eps_su_p:
        expected = (None, (Eh,))
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
        for i in range(len(coeff)):
            assert math.isclose(coeff[i], expect[i])

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
        for i in range(len(coeff)):
            assert math.isclose(coeff[i], expect[i])


@pytest.mark.parametrize('eps_0', eps0)
@pytest.mark.parametrize(
    'fc, eps_c0, eps_cu',
    [
        (-30, -0.002, -0.0035),
        (-20, -0.0015, -0.004),
        (-35, -0.002, -0.003),
        (-20000, -0.002, -0.0025),
        (-30000, -0.0025, -0.0035),
    ],
)
def test_marin_parabolarectangle_uniform_strain(eps_0, fc, eps_c0, eps_cu):
    """Test the marin coefficients for a strain plane of pure tension or
    compression for ParabolaRectangle constitutive law.
    """
    law = ParabolaRectangle(fc, eps_c0, eps_cu)

    # Check marin coefficients for stress integration
    if eps_0 > 0:
        expected = (None, (0.0,))
    elif eps_0 > eps_c0:
        expected = (
            None,
            (fc / eps_c0 * eps_0 * (2 - eps_0 / eps_c0), 0.0, 0.0),
        )
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
    elif eps_0 > eps_c0:
        expected = (None, (2 * fc / eps_c0 * (1 - eps_0 / eps_c0), 0.0))
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
    'fc, eps_c0, eps_cu, h',
    [
        (-30, -0.002, -0.0035, 400),
        (-20, -0.0015, -0.004, 400),
        (-35, -0.002, -0.003, 350),
        (-20000, -0.002, -0.0025, 0.4),
        (-30000, -0.0025, -0.0035, 0.3),
    ],
)
def test_marin_parabolarectangle_axial_bending_strain(
    eps_0, eps_min, fc, eps_c0, eps_cu, h
):
    """Test the marin coefficients for a strain plane of uniaxial bending
    with tension / compression for ParabolaRectangle constitutive law.
    """
    law = ParabolaRectangle(fc, eps_c0, eps_cu)

    # Check marin coefficients for stress integration
    kappa_y = (eps_min - eps_0) * 2 / h

    # We are not checking here case with uniform compression / tension, skip
    if kappa_y == 0:
        return

    expected = (
        [(eps_c0, 0), (eps_cu, eps_c0)],
        [
            (
                fc / eps_c0 * eps_0 * (2 - eps_0 / eps_c0),
                2 * fc * kappa_y / eps_c0 * (1 - eps_0 / eps_c0),
                -fc * kappa_y**2 / eps_c0**2,
            ),
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
        [(eps_c0, 0), (eps_cu, eps_c0)],
        [
            (
                2 * fc / eps_c0 * (1 - eps_0 / eps_c0),
                -2 * fc * kappa_y / eps_c0**2,
            ),
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


@pytest.mark.parametrize('eps_0', eps0)
@pytest.mark.parametrize(
    'x, y',
    [
        ([0, 0.002, 0.01], [0, 400, 420]),
        ([0, 0.002, 0.07], [0, 400000, 500000]),
    ],
)
def test_marin_puserdefined_uniform_strain(eps_0, x, y):
    """Test the marin coefficients for a strain plane of pure tension or
    compression for UserDefined constitutive law.
    """
    law = UserDefined(x, y)

    # Get the ordered numpy arrays
    x = law._x
    y = law._y

    # Check marin coefficients for stress integration
    # Find the index where x_target would fit in the ordered array
    eps_0 = law.preprocess_strains_with_limits(eps_0)
    found = False
    for i in range(len(x) - 1):
        if x[i] <= eps_0 <= x[i + 1]:
            stiffness = (y[i + 1] - y[i]) / (x[i + 1] - x[i])
            a0 = stiffness * (eps_0 - x[i]) + y[i]
            expected = (None, (a0, 0.0))
            found = True
            break
    if not found:
        expected = (None, (0.0,))

    strain, coeffs = law.__marin__([eps_0, 0])

    assert strain == expected[0]
    assert len(coeffs[0]) == len(expected[1])
    for i in range(len(coeffs[0])):
        assert math.isclose(coeffs[0][i], expected[1][i])

    # Check marin coefficients for modulus integration
    found = False
    for i in range(len(x) - 1):
        if x[i] <= eps_0 <= x[i + 1]:
            stiffness = (y[i + 1] - y[i]) / (x[i + 1] - x[i])
            expected = (None, (stiffness,))
            found = True
            break
    if not found:
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
    'x, y, h',
    [
        ([0, 0.002, 0.01], [0, 400, 420], 400),
        ([0, 0.002, 0.07], [0, 400000, 500000], 0.3),
    ],
)
def test_marin_userdefined_axial_bending_strain(eps_0, eps_min, x, y, h):
    """Test the marin coefficients for a strain plane of uniaxial bending
    with tension / compression for UserDefined constitutive law.
    """
    law = UserDefined(x, y)

    # Get the ordered numpy arrays
    x = law._x
    y = law._y

    # Check marin coefficients for stress integration
    kappa_y = (eps_min - eps_0) * 2 / h

    # We are not checking here case with uniform compression / tension, skip
    if kappa_y == 0:
        return

    strains = []
    coeffs = []
    for i in range(len(x) - 1):
        # For each branch of the linear piecewise function
        stiffness = (y[i + 1] - y[i]) / (x[i + 1] - x[i])
        strains.append((x[i], x[i + 1]))
        a0 = stiffness * (eps_0 - x[i]) + y[i]
        a1 = stiffness * kappa_y
        coeffs.append((a0, a1))
    expected = (strains, coeffs)

    strain, coeffs = law.__marin__([eps_0, kappa_y])

    assert len(strain) == len(expected[0])
    assert len(expected[1]) == len(coeffs)

    for coeff, expect in zip(coeffs, expected[1]):
        assert len(coeff) == len(expect)
        for i in range(len(coeff)):
            assert math.isclose(coeff[i], expect[i])

    # Check marin coefficients for modulus integration
    strains = []
    coeffs = []
    for i in range(len(x) - 1):
        # For each branch of the linear piecewise function
        stiffness = (y[i + 1] - y[i]) / (x[i + 1] - x[i])
        strains.append((x[i], x[i + 1]))
        coeffs.append((stiffness,))
    expected = (strains, coeffs)

    strain, coeffs = law.__marin_tangent__([eps_0, kappa_y])

    assert len(strain) == len(expected[0])
    assert len(expected[1]) == len(coeffs)

    for coeff, expect in zip(coeffs, expected[1]):
        assert len(coeff) == len(expect)
        for i in range(len(coeff)):
            assert math.isclose(coeff[i], expect[i])
