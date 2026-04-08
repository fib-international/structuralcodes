"""Tests for the ReinforcementACI318 material class."""

import math

import pytest

import structuralcodes
from structuralcodes.materials.reinforcement import (
    ReinforcementACI318,
    create_reinforcement,
)


@pytest.fixture(autouse=True)
def _reset_design_code():
    """Reset the global design code after each test."""
    yield
    structuralcodes.set_design_code(None)


def test_create_via_factory():
    """Test creating reinforcement via the factory function."""
    r = create_reinforcement(
        fyk=420,
        Es=200000,
        ftk=550,
        epsuk=0.05,
        design_code='aci318',
    )
    assert isinstance(r, ReinforcementACI318)


def test_default_name():
    """Test default name generation."""
    r = ReinforcementACI318(
        fyk=420,
        Es=200000,
        ftk=550,
        epsuk=0.05,
    )
    assert r.name == 'Reinforcement420'


def test_gamma_s_default():
    """Test default gamma_s is 1.0 for ACI."""
    r = ReinforcementACI318(
        fyk=420,
        Es=200000,
        ftk=550,
        epsuk=0.05,
    )
    assert r.gamma_s == 1.0


def test_fyd():
    """Test fyd returns unreduced fy (gamma_s=1.0)."""
    r = ReinforcementACI318(
        fyk=420,
        Es=200000,
        ftk=550,
        epsuk=0.05,
    )
    assert math.isclose(r.fyd(), 420)


def test_ftd():
    """Test ftd returns unreduced ftk (gamma_s=1.0)."""
    r = ReinforcementACI318(
        fyk=420,
        Es=200000,
        ftk=550,
        epsuk=0.05,
    )
    assert math.isclose(r.ftd(), 550)


def test_epsud():
    """Test epsud returns epsuk (no reduction for ACI)."""
    r = ReinforcementACI318(
        fyk=420,
        Es=200000,
        ftk=550,
        epsuk=0.05,
    )
    assert r.epsud() == 0.05


def test_constitutive_law_elastic():
    """Test elastic constitutive law creation."""
    r = ReinforcementACI318(
        fyk=420,
        Es=200000,
        ftk=550,
        epsuk=0.05,
        constitutive_law='elastic',
    )
    assert r.constitutive_law is not None


def test_constitutive_law_elasticplastic():
    """Test elastic-plastic constitutive law creation."""
    r = ReinforcementACI318(
        fyk=420,
        Es=200000,
        ftk=550,
        epsuk=0.05,
        constitutive_law='elasticplastic',
    )
    assert r.constitutive_law is not None


def test_constitutive_law_perfectly_plastic():
    """Test elastic perfectly plastic constitutive law."""
    r = ReinforcementACI318(
        fyk=420,
        Es=200000,
        ftk=550,
        epsuk=0.05,
        constitutive_law='elasticperfectlyplastic',
    )
    assert r.constitutive_law is not None
