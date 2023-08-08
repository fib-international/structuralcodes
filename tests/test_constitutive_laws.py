"""Tests for the _concrete_material_properties module"""
import math

import pytest

from structuralcodes.materials.constitutive_laws import Elastic
from structuralcodes.materials.constitutive_laws import ElasticPlastic


@pytest.mark.parametrize(
    'E, strain, expected',
    [(210000, 0.003, 210000 * 0.003), (200000, 0.002, 200000 * 0.002)],
)
def test_elastic_floats(E, strain, expected):
    """Test the elastic material"""

    assert math.isclose(Elastic(E).get_stress(strain), expected)


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

    assert math.isclose(ElasticPlastic(E, fy).get_stress(strain), expected)
