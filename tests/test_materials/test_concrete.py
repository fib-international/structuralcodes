"""Tests for the concrete material class."""

import math

import pytest

from structuralcodes.materials.concrete import ConcreteMC2010
from structuralcodes.materials.constitutive_laws import (
    BilinearCompression,
    Elastic,
    ParabolaRectangle,
    Popovics,
    Sargin,
)


def test_constitutive_law_setter_valid():
    """Test the constitutive law setter, valid law."""
    # Arrange
    concrete = ConcreteMC2010(45)
    constitutive_law = Elastic(30000)

    # Act and assert
    concrete.constitutive_law = constitutive_law
    assert isinstance(concrete.constitutive_law, Elastic)


def test_constitutive_law_setter_factory():
    """Test the constitutive law setter, valid law."""
    # Arrange
    concrete = ConcreteMC2010(45)

    # Act and assert
    concrete.constitutive_law = 'elastic'
    assert isinstance(concrete.constitutive_law, Elastic)
    assert math.isclose(concrete.constitutive_law._E, concrete.Eci)

    # Act and assert
    concrete.constitutive_law = 'bilinearcompression'
    assert isinstance(concrete.constitutive_law, BilinearCompression)
    assert math.isclose(concrete.constitutive_law._fc, -concrete.fcd())
    assert math.isclose(
        concrete.constitutive_law._eps_cu, -abs(concrete.eps_cu3)
    )
    assert math.isclose(
        concrete.constitutive_law._eps_c, -abs(concrete.eps_c3)
    )

    # Act and assert
    concrete.constitutive_law = 'parabolarectangle'
    assert isinstance(concrete.constitutive_law, ParabolaRectangle)
    assert math.isclose(concrete.constitutive_law._fc, -concrete.fcd())
    assert math.isclose(
        concrete.constitutive_law._eps_0, -abs(concrete.eps_c2)
    )
    assert math.isclose(
        concrete.constitutive_law._eps_u, -abs(concrete.eps_cu2)
    )
    assert math.isclose(
        concrete.constitutive_law._n, concrete.n_parabolic_rectangular
    )

    # Act and assert
    concrete.constitutive_law = 'sargin'
    assert isinstance(concrete.constitutive_law, Sargin)
    assert math.isclose(concrete.constitutive_law._fc, -concrete.fcd())
    assert math.isclose(
        concrete.constitutive_law._eps_c1, -abs(concrete.eps_c1)
    )
    assert math.isclose(
        concrete.constitutive_law._eps_cu1, -abs(concrete.eps_cu1)
    )
    assert math.isclose(concrete.constitutive_law._k, concrete.k_sargin)

    # Act and assert
    concrete.constitutive_law = 'popovics'
    assert isinstance(concrete.constitutive_law, Popovics)
    assert math.isclose(concrete.constitutive_law._fc, -concrete.fcd())
    assert math.isclose(
        concrete.constitutive_law._eps_c, -abs(concrete.eps_c1)
    )
    assert math.isclose(
        concrete.constitutive_law._eps_cu, -abs(concrete.eps_cu1)
    )


def test_constitutive_law_setter_invalid():
    """Test the constitutive law setter, invalid law."""
    # Arrange
    concrete = ConcreteMC2010(45)

    # Act and assert
    with pytest.raises(ValueError):
        concrete.constitutive_law = 'elasticplastic'
