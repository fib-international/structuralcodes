"""Tests for the code specific reinforcement classes."""

import math

import pytest

from structuralcodes.materials.constitutive_laws import (
    BilinearCompression,
    Elastic,
    ElasticPlastic,
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


def test_constitutive_law_setter_valid():
    """Test the constitutive law setter, valid law."""
    # Arrange
    steel = ReinforcementEC2_2004(
        fyk=500,
        Es=200000,
        density=7850,
        ftk=550,
        epsuk=0.07,
        constitutive_law=Elastic(200000),
    )

    # Assert
    assert isinstance(steel.constitutive_law, Elastic)


@pytest.mark.parametrize(
    'reinforcement_material',
    (ReinforcementEC2_2004, ReinforcementEC2_2023, ReinforcementMC2010),
)
@pytest.mark.parametrize(
    'constitutive_law',
    ('elastic', 'elasticplastic', 'elasticperfectlyplastic'),
)
def test_constitutive_law_setter_factory(
    constitutive_law, reinforcement_material
):
    """Test the constitutive law setter, valid law."""
    # Arrange
    steel = reinforcement_material(
        fyk=500,
        Es=200000,
        density=7850,
        ftk=550,
        epsuk=0.07,
        constitutive_law=constitutive_law,
    )

    # Assert
    if constitutive_law == 'elastic':
        assert isinstance(steel.constitutive_law, Elastic)
        assert math.isclose(steel.constitutive_law._E, steel.Es)
    elif constitutive_law == 'elasticplastic':
        assert isinstance(steel.constitutive_law, ElasticPlastic)
        assert math.isclose(steel.constitutive_law._E, steel.Es)
        assert math.isclose(steel.constitutive_law._fy, steel.fyd())
        assert math.isclose(steel.constitutive_law._eps_su, steel.epsud())
        assert math.isclose(steel.constitutive_law._eps_sy, steel.epsyd)
    elif constitutive_law == 'elasticperfectlyplastic':
        assert isinstance(steel.constitutive_law, ElasticPlastic)
        assert math.isclose(steel.constitutive_law._E, steel.Es)
        assert math.isclose(steel.constitutive_law._fy, steel.fyd())
        assert math.isclose(steel.constitutive_law._Eh, 0.0)
        assert math.isclose(steel.constitutive_law._eps_su, steel.epsud())
        assert math.isclose(steel.constitutive_law._eps_sy, steel.epsyd)


def test_constitutive_law_setter_invalid():
    """Test the constitutive law setter, invalid law."""
    # Act and assert
    with pytest.raises(ValueError):
        ReinforcementEC2_2004(
            fyk=500,
            Es=200000,
            density=7850,
            ftk=550,
            epsuk=0.07,
            constitutive_law='parabolarectangle',
        )


@pytest.mark.parametrize(
    'reinforcement_type',
    (ReinforcementMC2010, ReinforcementEC2_2004, ReinforcementEC2_2023),
)
def test_invalid_constitutive_law(reinforcement_type):
    """Test initializing a reinforcement object with an invalid constitutive
    law.
    """
    invalid_constitutive_law = BilinearCompression(500, 500 / 200000, 6e-2)
    with pytest.raises(ValueError):
        reinforcement_type(
            fyk=500,
            ftk=500,
            Es=200000,
            epsuk=6e-2,
            constitutive_law=invalid_constitutive_law,
        )


@pytest.mark.parametrize(
    'initial_strain, initial_stress, strain_compatibility',
    [
        (None, 300, True),
        (0.0015, None, None),
        (None, 450, None),
        (0.0025, None, True),
        (0.0025, None, False),
    ],
)
def test_initial_strain_and_stress(
    initial_strain, initial_stress, strain_compatibility
):
    """Test initializing reinforcement with initial strain and stress."""
    # Arrange
    fyk = 500
    Es = 200000
    ftk = 1.15 * fyk
    epsuk = 7.5e-2

    # Act
    reinf = ReinforcementEC2_2004(
        fyk=fyk,
        Es=Es,
        ftk=ftk,
        epsuk=epsuk,
        initial_strain=initial_strain,
        initial_stress=initial_stress,
        strain_compatibility=strain_compatibility,
    )

    # Assert
    assert reinf.strain_compatibility == strain_compatibility
    assert reinf.fyk == fyk
    assert reinf.Es == Es
    assert math.isclose(reinf.fyd(), fyk / 1.15)
    if initial_strain is not None:
        expected_stress = (
            initial_strain * Es
            if initial_strain < reinf.epsyd
            else reinf.constitutive_law.wrapped_law.get_stress(initial_strain)
        )
        # Check that stra_compatibility is not None and is false
        if (
            reinf.strain_compatibility is not None
            and not reinf._strain_compatibility
        ):
            expected_stress *= 0
            expected_stress += reinf.initial_stress
        assert math.isclose(reinf.initial_strain, initial_strain)
        assert math.isclose(reinf.initial_stress, expected_stress)
    if initial_stress is not None:
        expected_strain = (
            initial_stress / Es
            if initial_stress < reinf.fyd()
            else (initial_stress - reinf.fyd())
            / reinf.constitutive_law.wrapped_law._Eh
            + reinf.epsyd
        )
        if (
            reinf.strain_compatibility is not None
            and not reinf._strain_compatibility
        ):
            expected_stress *= 0
        assert math.isclose(reinf.initial_strain, expected_strain)
        assert math.isclose(reinf.initial_stress, initial_stress)


def test_initial_strain_and_stress_invalid():
    """Test initializing reinforcement with initial strain and stress."""
    # Arrange
    fyk = 500
    Es = 200000
    ftk = 1.15 * fyk
    epsuk = 7.5e-2

    # Act
    with pytest.raises(ValueError):
        ReinforcementEC2_2004(
            fyk=fyk,
            Es=Es,
            ftk=ftk,
            epsuk=epsuk,
            initial_strain=0.002,
            initial_stress=400,
        )
