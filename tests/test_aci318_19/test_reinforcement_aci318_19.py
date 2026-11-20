"""Tests for the ReinforcementACI318_19 material class."""

import math

import pytest

import structuralcodes
from structuralcodes.codes import aci318_19
from structuralcodes.materials.reinforcement import (
    ReinforcementACI318_19,
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
        design_code='aci318_19',
    )
    assert isinstance(r, ReinforcementACI318_19)


def test_default_name():
    """Test default name generation."""
    r = ReinforcementACI318_19(
        fyk=420,
        Es=200000,
        ftk=550,
        epsuk=0.05,
    )
    assert r.name == 'Reinforcement420'


def test_gamma_s_default():
    """Test default gamma_s is 1.0 for ACI."""
    r = ReinforcementACI318_19(
        fyk=420,
        Es=200000,
        ftk=550,
        epsuk=0.05,
    )
    assert r.gamma_s == 1.0


def test_from_grade():
    """Test creating reinforcement from an ASTM grade designation."""
    r = ReinforcementACI318_19.from_grade('60', epsuk=0.05)
    assert r.name == 'Grade 60'
    assert math.isclose(r.fyk, 420)
    assert math.isclose(r.Es, 200000)
    assert math.isclose(r.ftk, 550)
    assert math.isclose(r.epsuk, 0.05)
    assert math.isclose(r.density, aci318_19.pcf_to_kg_per_m3(490))


def test_from_ksi():
    """Test creating reinforcement from US customary stress inputs."""
    r = ReinforcementACI318_19.from_ksi(
        fy_ksi=60,
        fu_ksi=80,
        epsuk=0.05,
    )
    assert r.name == 'Reinforcement60ksi'
    assert math.isclose(r.fyk, aci318_19.ksi_to_mpa(60))
    assert math.isclose(r.Es, aci318_19.ksi_to_mpa(29000))
    assert math.isclose(r.ftk, aci318_19.ksi_to_mpa(80))
    assert math.isclose(r.epsuk, 0.05)


def test_invalid_gamma_s():
    """Test gamma_s values other than 1.0 are rejected."""
    with pytest.raises(ValueError):
        ReinforcementACI318_19(
            fyk=420,
            Es=200000,
            ftk=550,
            epsuk=0.05,
            gamma_s=1.15,
        )


def test_fyd():
    """Test fyd returns unreduced fy (gamma_s=1.0)."""
    r = ReinforcementACI318_19(
        fyk=420,
        Es=200000,
        ftk=550,
        epsuk=0.05,
    )
    assert math.isclose(r.fyd(), 420)


def test_ftd():
    """Test ftd returns unreduced ftk (gamma_s=1.0)."""
    r = ReinforcementACI318_19(
        fyk=420,
        Es=200000,
        ftk=550,
        epsuk=0.05,
    )
    assert math.isclose(r.ftd(), 550)


def test_epsud():
    """Test epsud returns epsuk (no reduction for ACI)."""
    r = ReinforcementACI318_19(
        fyk=420,
        Es=200000,
        ftk=550,
        epsuk=0.05,
    )
    assert r.epsud() == 0.05


def test_constitutive_law_elastic():
    """Test elastic constitutive law creation."""
    r = ReinforcementACI318_19(
        fyk=420,
        Es=200000,
        ftk=550,
        epsuk=0.05,
        constitutive_law='elastic',
    )
    assert r.constitutive_law is not None


def test_constitutive_law_elasticplastic():
    """Test elastic-plastic constitutive law creation."""
    r = ReinforcementACI318_19(
        fyk=420,
        Es=200000,
        ftk=550,
        epsuk=0.05,
        constitutive_law='elasticplastic',
    )
    assert r.constitutive_law is not None


def test_constitutive_law_perfectly_plastic():
    """Test elastic perfectly plastic constitutive law."""
    r = ReinforcementACI318_19(
        fyk=420,
        Es=200000,
        ftk=550,
        epsuk=0.05,
        constitutive_law='elasticperfectlyplastic',
    )
    assert r.constitutive_law is not None
