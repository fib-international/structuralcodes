"""Tests for the concrete material class."""

import math

import pytest

from structuralcodes.materials.concrete import (
    Concrete,
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


@pytest.mark.parametrize(
    'concrete_type, E_name',
    [
        (ConcreteMC2010, 'Eci'),
        (ConcreteEC2_2004, 'Ecm'),
        (ConcreteEC2_2023, 'Ecm'),
    ],
)
@pytest.mark.parametrize('fck', [20, 30, 40, 50, 60])
def test_constitutive_law_setter_factory(
    concrete_type: Concrete, fck: float, E_name: str
):
    """Test the constitutive law setter, valid law."""
    # Act and assert for elastic law
    concrete = concrete_type(fck=fck, constitutive_law='elastic')

    E = getattr(concrete, E_name)
    assert isinstance(concrete.constitutive_law, Elastic)
    assert math.isclose(concrete.constitutive_law._E, E)

    # Act and assert for parabolarectangle law
    concrete = concrete_type(fck=fck, constitutive_law='parabolarectangle')
    assert isinstance(concrete.constitutive_law, ParabolaRectangle)
    assert math.isclose(concrete.constitutive_law._fc, -concrete.fcd())
    assert math.isclose(
        concrete.constitutive_law._eps_0, -abs(concrete.eps_c2)
    )
    assert math.isclose(
        concrete.constitutive_law._eps_u, -abs(concrete.eps_cu2)
    )
    assert math.isclose(
        concrete.constitutive_law._n,
        concrete.n_parabolic_rectangular,
    )

    # Act and assert for Sargin law
    concrete = concrete_type(fck=concrete.fck, constitutive_law='sargin')
    assert isinstance(concrete.constitutive_law, Sargin)
    assert math.isclose(concrete.constitutive_law._fc, -concrete.fcm)
    assert math.isclose(
        concrete.constitutive_law._eps_c1, -abs(concrete.eps_c1)
    )
    assert math.isclose(
        concrete.constitutive_law._eps_cu1, -abs(concrete.eps_cu1)
    )
    assert math.isclose(concrete.constitutive_law._k, concrete.k_sargin)

    # Act and assert for Popovics law
    concrete = concrete_type(fck=concrete.fck, constitutive_law='popovics')
    assert isinstance(concrete.constitutive_law, Popovics)
    assert math.isclose(concrete.constitutive_law._fc, -concrete.fcd())
    assert math.isclose(
        concrete.constitutive_law._eps_c, -abs(concrete.eps_c1)
    )
    assert math.isclose(
        concrete.constitutive_law._eps_cu, -abs(concrete.eps_cu1)
    )

    # Test bilinear law only for MC2010 and EC2_2004
    if type(concrete) in (ConcreteMC2010, ConcreteEC2_2004):
        concrete = concrete_type(
            fck=concrete.fck, constitutive_law='bilinearcompression'
        )
        assert isinstance(concrete.constitutive_law, BilinearCompression)
        assert math.isclose(concrete.constitutive_law._fc, -concrete.fcd())
        assert math.isclose(
            concrete.constitutive_law._eps_cu, -abs(concrete.eps_cu3)
        )
        assert math.isclose(
            concrete.constitutive_law._eps_c, -abs(concrete.eps_c3)
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
