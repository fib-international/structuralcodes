"""Tests for the _concrete_material_properties module."""

import math

import numpy as np
import pytest
from numpy.testing import assert_allclose

from structuralcodes.materials.constitutive_laws import (
    Elastic,
    ElasticPlastic,
    ParabolaRectangle,
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
    assert math.isclose(Elastic(E).get_tangent(), E)
    assert math.isclose(Elastic(E).get_ultimate_strain()[0], 100)
    assert math.isclose(Elastic(E).get_ultimate_strain()[1], -100)


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
    assert math.isclose(material.get_ultimate_strain()[0], eps_su)
    assert math.isclose(material.get_ultimate_strain()[1], -eps_su)


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
    assert math.isclose(material.get_ultimate_strain()[0], max(eps_su))
    assert math.isclose(material.get_ultimate_strain()[1], min(eps_su))


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
    assert math.isclose(ElasticPlastic(E, fy).get_stress(strain)[0], expected)


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
    assert math.isclose(mat.get_stress(strain)[0], stress)
    assert math.isclose(mat.get_tangent(strain)[0], tangent)


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
        UserDefined(x, y, flag=flag).get_stress(strain)[0], expected
    )
    xmin = x[0] if x[0] < 0 else -x[-1]
    xmax = x[-1]
    assert math.isclose(UserDefined(x, y).get_ultimate_strain()[0], xmax)
    assert math.isclose(UserDefined(x, y).get_ultimate_strain()[1], xmin)

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
    assert math.isclose(material.get_ultimate_strain()[0], eps_su)
    assert math.isclose(material.get_ultimate_strain()[1], -eps_su)


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
    assert math.isclose(material.get_ultimate_strain()[0], max(eps_su))
    assert math.isclose(material.get_ultimate_strain()[1], min(eps_su))


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

    # compute stresses expected
    sig_expected = (
        fc
        * (k * eps / eps_c1 - (eps / eps_c1) ** 2)
        / (1 + (k - 2) * eps / eps_c1)
    )

    sig_computed = law.get_stress(eps)

    assert_allclose(sig_computed, sig_expected)
