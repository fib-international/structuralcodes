"""Tests for the _concrete_material_properties module"""

import math
import numpy as np

import pytest

from structuralcodes.materials.constitutive_laws import Elastic
from structuralcodes.materials.constitutive_laws import ElasticPlastic
from structuralcodes.materials.constitutive_laws import ParabolaRectangle
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
    """Test the elastic material"""

    assert math.isclose(Elastic(E).get_stress(strain), expected)


def test_elastic_numpy():
    """Test the elastic material with numpy input"""
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
def test_elasticplastic_floats(E, fy, strain, expected):
    """Test the elasticPlastic material"""

    assert math.isclose(ElasticPlastic(E, fy).get_stress(strain)[0], expected)


@pytest.mark.parametrize(
    'fc, eps_0, eps_u, strain, expected',
    [
        (-30.0, -0.002, -0.0035, 0.001, 0.0),
        (-30.0, -0.002, -0.0035, -0.002, -30),
        (-30.0, -0.002, -0.0035, -0.003, -30),
        (-30.0, -0.002, -0.0035, -0.004, 0),
        (-30.0, -0.002, -0.0035, -0.001, -22.5),
        (-45.0, -0.002, -0.004, 0.001, 0.0),
        (-45.0, -0.002, -0.004, -0.002, -45),
        (-45.0, -0.002, -0.004, -0.0035, -45),
        (-45.0, -0.002, -0.004, -0.0045, 0.0),
        (-45.0, -0.002, -0.004, -0.001, -33.75),
    ],
)
def test_parabola_rectangle_floats(fc, eps_0, eps_u, strain, expected):
    """Test the parabola-rectangle material"""

    assert math.isclose(
        ParabolaRectangle(fc, eps_0, eps_u).get_stress(strain)[0], expected
    )


@pytest.mark.parametrize(
    'x, y, strain, expected',
    [
        ([0, 0.002, 0.005], [0, 20, 23], 0.001, 10),
        ([0, 0.002, 0.005], [0, 20, 23], 0.002, 20),
        ([0, 0.002, 0.005], [0, 20, 23], 0.003, 21),
        ([0, 0.002, 0.005], [0, 20, 23], 0.004, 22),
        ([0, 0.002, 0.005], [0, 20, 23], 0.005, 23),
        ([0, 0.002, 0.005], [0, 20, 23], 0.006, 0),
        ([0, 0.002, 0.005], [0, 20, 23], -0.001, -10),
        ([0, 0.002, 0.005], [0, 20, 23], -0.002, -20),
        ([0, 0.002, 0.005], [0, 20, 23], -0.003, -21),
        ([0, 0.002, 0.005], [0, 20, 23], -0.004, -22),
        ([0, 0.002, 0.005], [0, 20, 23], -0.005, -23),
        ([0, 0.002, 0.005], [0, 20, 23], -0.006, 0),
    ],
)
def test_user_defined_floats(x, y, strain, expected):
    """Test the parabola-rectangle material"""

    assert math.isclose(UserDefined(x, y).get_stress(strain)[0], expected)
