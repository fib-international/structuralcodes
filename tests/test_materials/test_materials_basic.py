"""Tests for the constitutive_laws module."""

import math

import pytest

from structuralcodes.materials.basic import ElasticMaterial, GenericMaterial
from structuralcodes.materials.constitutive_laws import UserDefined


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
    mat = ElasticMaterial(E, density=7850)
    assert math.isclose(mat.constitutive_law.get_stress(strain), expected)
    assert math.isclose(mat.constitutive_law.get_tangent(strain), E)
    assert math.isclose(mat.constitutive_law.get_ultimate_strain()[0], -100)
    assert math.isclose(mat.constitutive_law.get_ultimate_strain()[1], 100)


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
    material = ElasticMaterial(E, density=7850, ultimate_strain=eps_su)
    assert math.isclose(
        material.constitutive_law.get_ultimate_strain()[0], -eps_su
    )
    assert math.isclose(
        material.constitutive_law.get_ultimate_strain()[1], eps_su
    )


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
    material = ElasticMaterial(E, density=7850, ultimate_strain=eps_su)
    assert math.isclose(
        material.constitutive_law.get_ultimate_strain()[0], min(eps_su)
    )
    assert math.isclose(
        material.constitutive_law.get_ultimate_strain()[1], max(eps_su)
    )


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
    mat = GenericMaterial(
        density=7850, constitutive_law=UserDefined(x, y, flag=flag)
    )
    assert math.isclose(mat.constitutive_law.get_stress(strain), expected)
    xmin = x[0] if x[0] < 0 else -x[-1]
    xmax = x[-1]
    assert math.isclose(mat.constitutive_law.get_ultimate_strain()[0], xmin)
    assert math.isclose(mat.constitutive_law.get_ultimate_strain()[1], xmax)


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
    material = GenericMaterial(
        density=7850, constitutive_law=UserDefined(x=x, y=y, eps_u=eps_su)
    )
    assert math.isclose(
        material.constitutive_law.get_ultimate_strain()[0], -eps_su
    )
    assert math.isclose(
        material.constitutive_law.get_ultimate_strain()[1], eps_su
    )


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
    material = GenericMaterial(
        density=7850, constitutive_law=UserDefined(x=x, y=y, eps_u=eps_su)
    )
    assert math.isclose(
        material.constitutive_law.get_ultimate_strain()[0], min(eps_su)
    )
    assert math.isclose(
        material.constitutive_law.get_ultimate_strain()[1], max(eps_su)
    )
