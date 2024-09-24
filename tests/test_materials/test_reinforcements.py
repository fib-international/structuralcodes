"""Tests for the code specific reinforcement classes."""

import math

import pytest

from structuralcodes.materials.constitutive_laws import (
    Elastic,
    ElasticPlastic,
    ParabolaRectangle,
)
from structuralcodes.materials.reinforcement import (
    ReinforcementEC2_2004,
    ReinforcementEC2_2023,
    ReinforcementMC2010,
)


@pytest.mark.parametrize(
    'reinforcement_material',
    (
        ReinforcementMC2010,
        ReinforcementEC2_2004,
        ReinforcementEC2_2023,
    ),
)
def test_reinforcements(reinforcement_material):
    """Basic test for the code specific reinforcement materials."""
    # Arrange
    fyk = 500
    Es = 200000
    ftk = 1.15 * fyk
    epsuk = 7.5e-2
    reinf = reinforcement_material(fyk=fyk, Es=Es, ftk=ftk, epsuk=epsuk)

    # Assert
    assert reinf.fyk == fyk
    assert reinf.Es == Es
    assert math.isclose(reinf.fyd(), fyk / 1.15)

    # Set new values
    reinf.fyk = 1.05 * reinf.fyk
    reinf.Es = 1.05 * reinf.Es
    reinf.ftk = 1.05 * reinf.ftk
    reinf.epsuk = 0.5 * reinf.epsuk

    # Assert
    assert reinf.fyk == 1.05 * fyk
    assert reinf.Es == 1.05 * Es
    assert reinf.ftk == 1.05 * ftk
    assert reinf.epsuk == 0.5 * epsuk


def test_constitutive_law_setter_valid():
    """Test the constitutive law setter, valid law."""
    # Arrange
    steel = ReinforcementEC2_2004(
        fyk=500, Es=200000, density=7850, ftk=550, epsuk=0.07
    )
    constitutive_law = Elastic(200000)

    # Act and assert
    steel.constitutive_law = constitutive_law
    assert isinstance(steel.constitutive_law, Elastic)


def test_constitutive_law_setter_factory():
    """Test the constitutive law setter, valid law."""
    # Arrange
    steel = ReinforcementEC2_2004(
        fyk=500, Es=200000, density=7850, ftk=550, epsuk=0.07
    )

    # Act and assert
    steel.constitutive_law = 'elastic'
    assert isinstance(steel.constitutive_law, Elastic)
    assert math.isclose(steel.constitutive_law._E, steel.Es)

    # Act and assert
    steel.constitutive_law = 'elasticplastic'
    assert isinstance(steel.constitutive_law, ElasticPlastic)
    assert math.isclose(steel.constitutive_law._E, steel.Es)
    assert math.isclose(steel.constitutive_law._fy, steel.fyd())
    assert math.isclose(steel.constitutive_law._eps_su, steel.epsud())
    assert math.isclose(steel.constitutive_law._eps_sy, steel.epsyd)

    # Act and assert
    steel.constitutive_law = 'elasticperfectlyplastic'
    assert isinstance(steel.constitutive_law, ElasticPlastic)
    assert math.isclose(steel.constitutive_law._E, steel.Es)
    assert math.isclose(steel.constitutive_law._fy, steel.fyd())
    assert math.isclose(steel.constitutive_law._Eh, 0.0)
    assert math.isclose(steel.constitutive_law._eps_su, steel.epsud())
    assert math.isclose(steel.constitutive_law._eps_sy, steel.epsyd)


def test_constitutive_law_setter_invalid():
    """Test the constitutive law setter, invalid law."""
    # Arrange
    steel = ReinforcementEC2_2004(
        fyk=500, Es=200000, density=7850, ftk=550, epsuk=0.07
    )
    constitutive_law = ParabolaRectangle(500)

    # Act and assert
    with pytest.raises(ValueError):
        steel.constitutive_law = constitutive_law
        steel.constitutive_law = 'parabolarectangle'
