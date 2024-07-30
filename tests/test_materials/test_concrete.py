"""Tests for the concrete material class."""

import pytest

from structuralcodes.materials.concrete import ConcreteMC2010
from structuralcodes.materials.constitutive_laws import Elastic, ElasticPlastic


def test_constitutive_law_setter_valid():
    """Test the constitutive law setter, valid law."""
    # Arrange
    concrete = ConcreteMC2010(45)
    constitutive_law = Elastic(30000)

    # Act and assert
    concrete.constitutive_law = constitutive_law
    assert isinstance(concrete.constitutive_law, Elastic)


def test_constitutive_law_setter_invalid():
    """Test the constitutive law setter, invalid law."""
    # Arrange
    concrete = ConcreteMC2010(45)
    constitutive_law = ElasticPlastic(30000, 45)

    # Act and assert
    with pytest.raises(ValueError):
        concrete.constitutive_law = constitutive_law
