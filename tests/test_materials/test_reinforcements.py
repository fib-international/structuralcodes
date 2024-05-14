"""Tests for the code specific reinforcement classes."""

import math

import pytest

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
    assert math.isclose(reinf.fyd, fyk / 1.15)

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
