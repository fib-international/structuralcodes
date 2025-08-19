"""Tests for the constitutive_laws module."""

import math

import numpy as np
import pytest
from numpy.testing import assert_allclose

from structuralcodes.materials.constitutive_laws import (
    BilinearCompression,
    Elastic,
    ElasticPlastic,
    InitialStrain,
    ParabolaRectangle,
    Popovics,
    Sargin,
    UserDefined,
)


@pytest.mark.parametrize(
    'E, strain, expected',
    [
        (210000, 0.003, 210000 * 0.003),
        (200000, 0.002, 200000 * 0.002),
        (200000, -0.002, -200000 * 0.002),
        (200000, -0.01, -200000 * 0.01),
    ],
)
def test_elastic_floats(E, strain, expected):
    """Test the elastic material."""
    assert math.isclose(Elastic(E).get_stress(strain), expected)
    assert math.isclose(Elastic(E).get_tangent(strain), E)
    assert math.isclose(Elastic(E).get_ultimate_strain()[0], -100)
    assert math.isclose(Elastic(E).get_ultimate_strain()[1], 100)


@pytest.mark.parametrize(
    'E, eps_su',
    [
        (210000, 0.07),
        (210000, 0.03),
        (200000, 0.002),
    ],
)
def test_elastic_set_ultimate_strain_float(E, eps_su):
    """Test elastic material set ultimate strain with a float."""
    material = Elastic(E)
    material.set_ultimate_strain(eps_su)
    assert math.isclose(material.get_ultimate_strain()[0], -eps_su)
    assert math.isclose(material.get_ultimate_strain()[1], eps_su)


@pytest.mark.parametrize(
    'E, eps_su',
    [
        (210000, (0.07, -0.07)),
        (210000, (0.03, -0.02)),
        (200000, (0.002, -0.001)),
    ],
)
def test_elastic_set_ultimate_strain_tuple(E, eps_su):
    """Test elastic material set ultimate strain with a tuple."""
    material = Elastic(E)
    material.set_ultimate_strain(eps_su=eps_su)
    assert math.isclose(material.get_ultimate_strain()[0], min(eps_su))
    assert math.isclose(material.get_ultimate_strain()[1], max(eps_su))


@pytest.mark.parametrize(
    'E, eps_su',
    [
        (210000, (0.07, 0.02)),
        (210000, (-0.03, -0.02)),
        (200000, (-0.002, -0.001, 0.003)),
        (200000, (-0.002,)),
        (210000, [0.07, -0.07]),
    ],
)
def test_elastic_set_ultimate_strain_tuple_error(E, eps_su):
    """Test elastic material set ultimate strain raise error."""
    material = Elastic(E)
    with pytest.raises(ValueError):
        material.set_ultimate_strain(eps_su)


def test_elastic_numpy():
    """Test the elastic material with numpy input."""
    E = 200000
    strain = np.linspace(0, 0.01, 100)
    sig_expected = E * strain
    sig = Elastic(E).get_stress(strain)
    assert np.allclose(sig, sig_expected)


@pytest.mark.parametrize(
    'E, fy, strain, expected',
    [
        (210000, 410, 0.001, 210.0),
        (210000, 410, 0.003, 410.0),
        (210000, 410, 0.010, 410.0),
        (200000, 450, 0.002, 400.0),
        (200000, 450, 0.004, 450.0),
        (200000, 450, 0.010, 450.0),
    ],
)
def test_elasticplastic_get_stress(E, fy, strain, expected):
    """Test the elasticPlastic material."""
    assert math.isclose(ElasticPlastic(E, fy).get_stress(strain), expected)


def test_elasticplastic_input_correct():
    """Test invalid input values for ElasticPlastic."""
    with pytest.raises(ValueError) as excinfo:
        ElasticPlastic(-210000, 450)
    assert str(excinfo.value) == 'Elastic modulus E must be greater than zero'


@pytest.mark.parametrize(
    'E, fy, strain, expected',
    [
        (210000, 410, 0.001, 210000),
        (210000, 410, 0.003, 0.0),
        (210000, 410, 0.010, 0.0),
        (200000, 450, 0.002, 200000),
        (200000, 450, 0.004, 0.0),
        (200000, 450, 0.010, 0.0),
    ],
)
def test_elasticplastic_get_tangent(E, fy, strain, expected):
    """Test the elasticPlastic material."""
    assert math.isclose(ElasticPlastic(E, fy).get_tangent(strain), expected)


@pytest.mark.parametrize(
    'fc, eps_0, eps_u, strain, stress, tangent',
    [
        (-30.0, -0.002, -0.0035, 0.001, 0.0, 0.0),
        (-30.0, -0.002, -0.0035, -0.002, -30, 0.0),
        (-30.0, -0.002, -0.0035, -0.003, -30, 0.0),
        (-30.0, -0.002, -0.0035, -0.004, 0, 0.0),
        (-30.0, -0.002, -0.0035, -0.001, -22.5, 15000.0),
        (-30.0, -0.002, -0.0035, 0.0, 0.0, 30000.0),
        (-45.0, -0.002, -0.004, 0.001, 0.0, 0.0),
        (-45.0, -0.002, -0.004, -0.002, -45, 0.0),
        (-45.0, -0.002, -0.004, -0.0035, -45, 0.0),
        (-45.0, -0.002, -0.004, -0.0045, 0.0, 0.0),
        (-45.0, -0.002, -0.004, -0.001, -33.75, 22500.0),
    ],
)
def test_parabola_rectangle_floats(fc, eps_0, eps_u, strain, stress, tangent):
    """Test the parabola-rectangle material."""
    mat = ParabolaRectangle(fc, eps_0, eps_u)
    assert math.isclose(mat.get_stress(strain), stress)
    assert math.isclose(mat.get_tangent(strain), tangent)


@pytest.mark.parametrize(
    'x, y, flag, strain, expected',
    [
        ([0, 0.002, 0.005], [0, 20, 23], 0, 0.001, 10),
        ([0, 0.002, 0.005], [0, 20, 23], 0, 0.002, 20),
        ([0, 0.002, 0.005], [0, 20, 23], 0, 0.003, 21),
        ([0, 0.002, 0.005], [0, 20, 23], 0, 0.004, 22),
        ([0, 0.002, 0.005], [0, 20, 23], 0, 0.005, 23),
        ([0, 0.002, 0.005], [0, 20, 23], 0, 0.006, 0),
        ([0, 0.002, 0.005], [0, 20, 23], 0, -0.001, -10),
        ([0, 0.002, 0.005], [0, 20, 23], 0, -0.002, -20),
        ([0, 0.002, 0.005], [0, 20, 23], 0, -0.003, -21),
        ([0, 0.002, 0.005], [0, 20, 23], 0, -0.004, -22),
        ([0, 0.002, 0.005], [0, 20, 23], 0, -0.005, -23),
        ([0, 0.002, 0.005], [0, 20, 23], 0, -0.006, 0),
        ([0, 0.002, 0.005], [0, 20, 23], 1, -0.006, -23),
        ([0, 0.002, 0.005], [0, 20, 23], 1, 0.006, 23),
        ([0, 0.002, 0.005], [0, 20, 23], 2, -0.006, -24),
        ([0, 0.002, 0.005], [0, 20, 23], 2, 0.006, 24),
        ([0, 0.002, 0.005], [0, 20, 23], 3, -0.006, -27.6),
        ([0, 0.002, 0.005], [0, 20, 23], 3, 0.006, 27.6),
        (
            [-0.0035, -0.002, 0, 0.001, 0.005],
            [-35, -35, 0, 3, 0],
            0,
            -0.005,
            0.0,
        ),
        (
            [-0.0035, -0.002, 0, 0.001, 0.005],
            [-35, -35, 0, 3, 0],
            1,
            -0.005,
            -35,
        ),
        (
            [-0.0035, -0.002, 0, 0.001, 0.005],
            [-35, -35, 0, 3, 0],
            2,
            -0.005,
            -35,
        ),
        (
            [-0.0035, -0.002, 0, 0.001, 0.005],
            [-35, -35, 0, 3, 0],
            3,
            -0.005,
            -50,
        ),
        (
            [-0.0035, -0.002, 0, 0.001, 0.005],
            [-35, -35, 0, 3, 0],
            0,
            -0.003,
            -35,
        ),
        (
            [-0.0035, -0.002, 0, 0.001, 0.005],
            [-35, -35, 0, 3, 0],
            0,
            -0.001,
            -17.5,
        ),
        ([-0.0035, -0.002, 0, 0.001, 0.005], [-35, -35, 0, 3, 0], 0, 0.0, 0.0),
        (
            [-0.0035, -0.002, 0, 0.001, 0.005],
            [-35, -35, 0, 3, 0],
            0,
            0.0005,
            1.5,
        ),
        (
            [-0.0035, -0.002, 0, 0.001, 0.005],
            [-35, -35, 0, 3, 0],
            0,
            0.001,
            3.0,
        ),
        (
            [-0.0035, -0.002, 0, 0.001, 0.005],
            [-35, -35, 0, 3, 0],
            0,
            0.003,
            1.5,
        ),
        (
            [-0.0035, -0.002, 0, 0.001, 0.005],
            [-35, -35, 0, 3, 0],
            0,
            0.006,
            0.0,
        ),
        (
            [-0.0035, -0.002, 0, 0.001, 0.005],
            [-35, -35, 0, 3, 0],
            1,
            0.006,
            0.0,
        ),
        (
            [-0.0035, -0.002, 0, 0.001, 0.005],
            [-35, -35, 0, 3, 0],
            2,
            0.006,
            -0.75,
        ),
        (
            [-0.0035, -0.002, 0, 0.001, 0.005],
            [-35, -35, 0, 3, 0],
            3,
            0.006,
            0.0,
        ),
    ],
)
def test_user_defined_floats(x, y, flag, strain, expected):
    """Test the parabola-rectangle material."""
    assert math.isclose(
        UserDefined(x, y, flag=flag).get_stress(strain), expected
    )
    xmin = x[0] if x[0] < 0 else -x[-1]
    xmax = x[-1]
    assert math.isclose(UserDefined(x, y).get_ultimate_strain()[0], xmin)
    assert math.isclose(UserDefined(x, y).get_ultimate_strain()[1], xmax)

    with pytest.raises(ValueError) as excinfo:
        UserDefined(x[:-1], y)
    assert str(excinfo.value) == 'The two arrays should have the same length'
    with pytest.raises(ValueError) as excinfo:
        UserDefined(x, y, flag=5)
    assert str(excinfo.value) == 'Flag can assume values 0, 1, 2 or 3.'


@pytest.mark.parametrize(
    'E, fy, eps_su',
    [
        (210000, 350, 0.07),
        (210000, 250, 0.03),
        (200000, 350, 0.002),
    ],
)
def test_userdefined_set_ultimate_strain_float(E, fy, eps_su):
    """Test UserDefined material set ultimate strain with a float."""
    x = [-3 * fy / E, 0, 3 * fy / E]
    y = [-fy, 0, fy]
    material = UserDefined(x=x, y=y)
    material.set_ultimate_strain(eps_su)
    assert math.isclose(material.get_ultimate_strain()[0], -eps_su)
    assert math.isclose(material.get_ultimate_strain()[1], eps_su)


@pytest.mark.parametrize(
    'E, fy, eps_su',
    [
        (210000, 350, (0.07, -0.07)),
        (210000, 250, (0.03, -0.02)),
        (200000, 350, (0.002, -0.001)),
    ],
)
def test_userdefined_set_ultimate_strain_tuple(E, fy, eps_su):
    """Test UserDefined material set ultimate strain with a tuple."""
    x = [-3 * fy / E, 0, 3 * fy / E]
    y = [-fy, 0, fy]
    material = UserDefined(x=x, y=y)
    material.set_ultimate_strain(eps_su=eps_su)
    assert math.isclose(material.get_ultimate_strain()[0], min(eps_su))
    assert math.isclose(material.get_ultimate_strain()[1], max(eps_su))


@pytest.mark.parametrize(
    'E, fy, eps_su',
    [
        (210000, 350, (0.07, 0.02)),
        (210000, 400, (-0.03, -0.02)),
        (200000, 355, (-0.002, -0.001, 0.003)),
        (200000, 200, (-0.002,)),
        (210000, 300, [0.07, -0.07]),
    ],
)
def test_userdefined_set_ultimate_strain_tuple_error(E, fy, eps_su):
    """Test UserDefined material set ultimate strain raise error."""
    x = [-3 * fy / E, 0, 3 * fy / E]
    y = [-fy, 0, fy]
    material = UserDefined(x=x, y=y)
    with pytest.raises(ValueError):
        material.set_ultimate_strain(eps_su)


@pytest.mark.parametrize(
    'fc, eps_c1, eps_cu1, k',
    [
        (-12, -1.9e-3, -3.5e-3, 2.44),
        (-16, -2.0e-3, -3.5e-3, 2.36),
        (-20, -2.1e-3, -3.5e-3, 2.28),
        (-25, -2.2e-3, -3.5e-3, 2.15),
        (-30, -2.3e-3, -3.5e-3, 2.04),
        (
            -35,
            -2.3e-3,
            -3.5e-3,
            1.92,
        ),
        (-40, -2.4e-3, -3.5e-3, 1.82),
        (-45, -2.5e-3, -3.5e-3, 1.74),
        (-50, -2.6e-3, -3.4e-3, 1.66),
        (-55, -2.6e-3, -3.4e-3, 1.61),
        (-60, -2.7e-3, -3.3e-3, 1.55),
        (-70, -2.7e-3, -3.2e-3, 1.47),
        (-80, -2.8e-3, -3.1e-3, 1.41),
        (-90, -2.9e-3, -3.0e-3, 1.36),
        (-100, -3.0e-3, -3.0e-3, 1.32),
        (-110, -3.0e-3, -3.0e-3, 1.24),
        (-120, -3.0e-3, -3.0e-3, 1.18),
    ],
)
def test_sargin(fc, eps_c1, eps_cu1, k):
    """Test Sargin material."""
    law = Sargin(fc=fc, eps_c1=eps_c1, eps_cu1=eps_cu1, k=k)

    eps = np.linspace(0, eps_cu1, 20)

    # compute expected
    sig_expected = (
        fc
        * (k * eps / eps_c1 - (eps / eps_c1) ** 2)
        / (1 + (k - 2) * eps / eps_c1)
    )
    tan_expected = (
        fc
        / eps_c1
        * ((2 - k) * (eps / eps_c1) ** 2 - 2 * (eps / eps_c1) + k)
        / (1 + (k - 2) * eps / eps_c1) ** 2
    )

    # compute from Sargin
    sig_computed = law.get_stress(eps)
    tan_computed = law.get_tangent(eps)

    # Compare the two
    assert_allclose(sig_computed, sig_expected)
    assert_allclose(tan_computed, tan_expected)

    # Test getting ultimate strain
    eps_min, eps_max = law.get_ultimate_strain()
    assert math.isclose(eps_min, eps_cu1)
    assert math.isclose(eps_max, 100)

    eps_min, eps_max = law.get_ultimate_strain(yielding=True)
    assert math.isclose(eps_min, eps_c1)
    assert math.isclose(eps_max, 100)


@pytest.mark.parametrize(
    'fc, eps_c, eps_cu',
    [
        (12, 1.9e-3, 3.5e-3),
        (16, 2.0e-3, 3.5e-3),
        (20, 2.1e-3, 3.5e-3),
        (25, 2.2e-3, 3.5e-3),
        (30, 2.3e-3, 3.5e-3),
        (35, 2.3e-3, 3.5e-3),
        (40, 2.4e-3, 3.5e-3),
        (45, 2.5e-3, 3.5e-3),
        (50, 2.6e-3, 3.4e-3),
        (55, 2.6e-3, 3.4e-3),
        (60, 2.7e-3, 3.3e-3),
        (70, 2.7e-3, 3.2e-3),
        (80, 2.8e-3, 3.1e-3),
        (90, 2.9e-3, 3.0e-3),
        (100, 3.0e-3, 3.0e-3),
        (110, 3.0e-3, 3.0e-3),
        (120, 3.0e-3, 3.0e-3),
    ],
)
def test_popovics(fc, eps_c, eps_cu):
    """Test Popovics material."""
    Ec = 5000 * abs(fc) ** 0.5
    law = Popovics(fc=fc, eps_c=eps_c, eps_cu=eps_cu)

    eps = np.linspace(0, eps_cu, 20)

    # compute expected
    Esec = fc / eps_c
    n = Ec / (Ec - Esec)
    sig_expected = -fc * (eps / eps_c) * n / (n - 1 + (eps / eps_c) ** n)
    tan_expected = (
        fc
        * (1 - (eps / eps_c) ** n)
        * n
        * (n - 1)
        / (n - 1 + (eps / eps_c) ** n) ** 2
        / eps_c
    )

    # compute from Popovics
    sig_computed = law.get_stress(-eps)
    tan_computed = law.get_tangent(-eps)

    # Compare the two
    assert_allclose(sig_computed, sig_expected)
    assert_allclose(tan_computed, tan_expected)

    # Test getting ultimate strain
    eps_min, eps_max = law.get_ultimate_strain()
    assert math.isclose(eps_min, -eps_cu)
    assert math.isclose(eps_max, 100)

    eps_min, eps_max = law.get_ultimate_strain(yielding=True)
    assert math.isclose(eps_min, -eps_c)
    assert math.isclose(eps_max, 100)


@pytest.mark.parametrize(
    'fc, eps_c, eps_cu',
    [
        (12, 1.75e-3, 3.5e-3),
        (16, 1.75e-3, 3.5e-3),
        (20, 1.75e-3, 3.5e-3),
        (25, 1.75e-3, 3.5e-3),
        (30, 1.75e-3, 3.5e-3),
        (35, 1.75e-3, 3.5e-3),
        (40, 1.75e-3, 3.5e-3),
        (45, 1.75e-3, 3.5e-3),
        (50, 1.75e-3, 3.5e-3),
        (55, 1.8e-3, 3.1e-3),
        (60, 1.9e-3, 2.9e-3),
        (70, 2.0e-3, 2.7e-3),
        (80, 2.2e-3, 2.6e-3),
        (90, 2.3e-3, 2.6e-3),
    ],
)
def test_bilinearcompression(fc, eps_c, eps_cu):
    """Test BilinearCompression material."""
    law = BilinearCompression(fc=fc, eps_c=eps_c, eps_cu=eps_cu)

    eps = np.linspace(-2 * eps_cu, eps_c, 20)

    # compute expected
    E = fc / eps_c
    sig_expected = E * eps
    sig_expected[sig_expected < -fc] = -fc
    sig_expected[eps > 0] = 0
    sig_expected[eps < -eps_cu] = 0
    tan_expected = np.ones_like(eps) * E
    tan_expected[eps >= 0] = 0
    tan_expected[eps < -eps_c] = 0

    # compute from BilinearCompression
    sig_computed_array = law.get_stress(eps)
    sig_computed_scalar = np.array(
        [law.get_stress(eps_scalar) for eps_scalar in eps]
    )
    tan_computed_array = law.get_tangent(eps)
    tan_computed_scalar = np.array(
        [law.get_tangent(eps_scalar) for eps_scalar in eps]
    )

    # Compare the two
    assert_allclose(sig_computed_array, sig_expected)
    assert_allclose(sig_computed_scalar, sig_expected)
    assert_allclose(tan_computed_array, tan_expected)
    assert_allclose(tan_computed_scalar, tan_expected)

    # Test getting ultimate strain
    eps_min, eps_max = law.get_ultimate_strain()
    assert math.isclose(eps_min, -eps_cu)
    assert math.isclose(eps_max, 100)

    eps_min, eps_max = law.get_ultimate_strain(yielding=True)
    assert math.isclose(eps_min, -eps_c)
    assert math.isclose(eps_max, 100)


@pytest.mark.parametrize(
    'fy, eps_su',
    [
        (350, 0.07),
        (250, 0.03),
        (355, 0.002),
    ],
)
@pytest.mark.parametrize(
    'initial_strain', [0.0, 0.001, 0.005, 0.01, -0.001, -0.005]
)
def test_initstrain(fy, eps_su, initial_strain):
    """Test InitStrain constitutive law."""
    base_law = ElasticPlastic(fy=fy, E=200000, eps_su=eps_su)
    law = InitialStrain(
        constitutive_law=base_law, initial_strain=initial_strain
    )

    ult_strain_base = base_law.get_ultimate_strain()
    ult_strain = law.get_ultimate_strain()
    assert math.isclose(
        ult_strain_base[0] - initial_strain,
        ult_strain[0],
    )
    assert math.isclose(
        ult_strain_base[1] - initial_strain,
        ult_strain[1],
    )

    strain = np.linspace(-eps_su * 1.1, eps_su * 1.1, 100)

    # Test stresses
    sig_expected = base_law.get_stress(strain + initial_strain)
    sig_computed = law.get_stress(strain)

    assert_allclose(sig_computed, sig_expected)

    # Test tangents
    tan_expected = base_law.get_tangent(strain + initial_strain)
    tan_computed = law.get_tangent(strain)

    assert_allclose(tan_computed, tan_expected)


def test_initial_strain_not_constitutive_law():
    """Test if TypeError is raised if an InitialStrain law is initialized
    without a ConstitutiveLaw.
    """
    with pytest.raises(TypeError):
        InitialStrain(constitutive_law=None, initial_strain=1e-3)


def test_initial_strain_strain_compatibility_property():
    """Test the strain_compatibility property."""
    # Arrange
    constitutive_law = Elastic(E=30000)
    initial_strain_value = 2e-3

    initial_strain_no_strain_compatibility = InitialStrain(
        constitutive_law=constitutive_law,
        initial_strain=initial_strain_value,
        strain_compatibility=False,
    )

    initial_strain_strain_compatibility = InitialStrain(
        constitutive_law=constitutive_law,
        initial_strain=initial_strain_value,
        strain_compatibility=True,
    )

    # Act and assert
    assert not initial_strain_no_strain_compatibility.strain_compatibility
    assert initial_strain_strain_compatibility.strain_compatibility


def test_initial_strain_get_stress_tangent_no_strain_compatibility():
    """Test the get_stress and get_tangent methods without strain
    compatibility.
    """
    # Arrange
    constitutive_law = Elastic(E=30000)
    initial_strain_value = 2e-3

    initial_strain = InitialStrain(
        constitutive_law=constitutive_law,
        initial_strain=initial_strain_value,
        strain_compatibility=False,
    )

    initial_stress = constitutive_law.get_stress(initial_strain_value)
    initial_tangent = constitutive_law.get_tangent(0)

    strains = np.linspace(-2e-3, 2e-3, 10)

    # Act
    stresses = initial_strain.get_stress(strains)
    tangents = initial_strain.get_tangent(strains)

    # Assert
    assert np.allclose(stresses, initial_stress)
    assert np.allclose(tangents, initial_tangent * 1e-6)
