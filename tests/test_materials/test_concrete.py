"""Tests for the concrete material class."""

import math

import pytest

from structuralcodes.materials.concrete import (
    ConcreteEC2_2004,
    ConcreteEC2_2023,
    ConcreteMC2010,
)
from structuralcodes.materials.constitutive_laws import (
    BilinearCompression,
    Elastic,
    ElasticPlastic,
    ParabolaRectangle,
    Popovics,
    Sargin,
)


def test_constitutive_law_setter_valid():
    """Test the constitutive law setter, valid law."""
    # Arrange
    constitutive_law = Elastic(30000)
    concrete = ConcreteMC2010(45, constitutive_law=constitutive_law)

    # Act and assert
    assert isinstance(concrete.constitutive_law, Elastic)


def test_constitutive_law_setter_factory():
    """Test the constitutive law setter, valid law."""
    # Arrange
    concretes = [
        ConcreteMC2010(45),
        ConcreteEC2_2004(45),
        ConcreteEC2_2023(45),
    ]
    elastic_modulus_name = ['Eci', 'Ecm', 'Ecm']

    for concrete, E_name in zip(concretes, elastic_modulus_name):
        # Act and assert for elastic constitutive law
        new_concrete = type(concrete)(
            fck=concrete.fck, constitutive_law='elastic'
        )
        E = getattr(new_concrete, E_name)
        assert isinstance(new_concrete.constitutive_law, Elastic)
        assert math.isclose(new_concrete.constitutive_law._E, E)

        # Act and assert for parabolarectangle law
        new_concrete = type(concrete)(
            fck=concrete.fck, constitutive_law='parabolarectangle'
        )
        assert isinstance(new_concrete.constitutive_law, ParabolaRectangle)
        assert math.isclose(
            new_concrete.constitutive_law._fc, -new_concrete.fcd()
        )
        assert math.isclose(
            new_concrete.constitutive_law._eps_0, -abs(new_concrete.eps_c2)
        )
        assert math.isclose(
            new_concrete.constitutive_law._eps_u, -abs(new_concrete.eps_cu2)
        )
        assert math.isclose(
            new_concrete.constitutive_law._n,
            new_concrete.n_parabolic_rectangular,
        )

        # Act and assert for Sargin law
        new_concrete = type(concrete)(
            fck=concrete.fck, constitutive_law='sargin'
        )
        assert isinstance(new_concrete.constitutive_law, Sargin)
        assert math.isclose(
            new_concrete.constitutive_law._fc, -new_concrete.fcd()
        )
        assert math.isclose(
            new_concrete.constitutive_law._eps_c1, -abs(new_concrete.eps_c1)
        )
        assert math.isclose(
            new_concrete.constitutive_law._eps_cu1, -abs(new_concrete.eps_cu1)
        )
        assert math.isclose(
            new_concrete.constitutive_law._k, new_concrete.k_sargin
        )

        # Act and assert for Popovics law
        new_concrete = type(concrete)(
            fck=concrete.fck, constitutive_law='popovics'
        )
        assert isinstance(new_concrete.constitutive_law, Popovics)
        assert math.isclose(
            new_concrete.constitutive_law._fc, -new_concrete.fcd()
        )
        assert math.isclose(
            new_concrete.constitutive_law._eps_c, -abs(new_concrete.eps_c1)
        )
        assert math.isclose(
            new_concrete.constitutive_law._eps_cu, -abs(new_concrete.eps_cu1)
        )

    # Test bilinear law only for MC2010 and EC2_2004
    concretes = concretes[:-1]
    elastic_modulus_name = elastic_modulus_name[:-1]

    # Act and assert
    for concrete, E_name in zip(concretes, elastic_modulus_name):
        new_concrete = type(concrete)(
            fck=concrete.fck, constitutive_law='bilinearcompression'
        )
        assert isinstance(new_concrete.constitutive_law, BilinearCompression)
        assert math.isclose(
            new_concrete.constitutive_law._fc, -new_concrete.fcd()
        )
        assert math.isclose(
            new_concrete.constitutive_law._eps_cu, -abs(new_concrete.eps_cu3)
        )
        assert math.isclose(
            new_concrete.constitutive_law._eps_c, -abs(new_concrete.eps_c3)
        )


def test_constitutive_law_setter_invalid():
    """Test the constitutive law setter, invalid law."""
    # Act and assert
    with pytest.raises(ValueError):
        ConcreteMC2010(45, constitutive_law='elasticplastic')


@pytest.mark.parametrize(
    'concrete_type', (ConcreteMC2010, ConcreteEC2_2004, ConcreteEC2_2023)
)
def test_invalid_constitutive_law(concrete_type):
    """Test initializing a concrete object with an invalid constitutive law."""
    invalid_constitutive_law = ElasticPlastic(
        E=30000, fy=45, Eh=0, eps_su=3.5e-3
    )
    with pytest.raises(ValueError):
        concrete_type(fck=45, constitutive_law=invalid_constitutive_law)
