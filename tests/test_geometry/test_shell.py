"""Tests for the Geometry."""

import numpy as np
import pytest

from structuralcodes.geometry import (
    PointGeometry,
    ShellGeometry,
    ShellReinforcement,
)
from structuralcodes.materials.concrete import ConcreteEC2_2004
from structuralcodes.materials.constitutive_laws import (
    Elastic2D,
    ElasticPlastic,
)


def test_shell_geometry():
    """Test the ShellGeometry class."""
    # Create a material to use
    concrete = ConcreteEC2_2004(fck=35)
    # Choose a constitutive law to use
    const = Elastic2D(E=200000, nu=0.2)

    for mat in (concrete, const):
        shell = ShellGeometry(
            thickness=200, material=mat, name='Shell', group_label='Group1'
        )
        assert shell.thickness == 200
        assert isinstance(shell.material, (ConcreteEC2_2004, Elastic2D))
        assert shell.name == 'Shell'
        assert shell.group_label == 'Group1'


def test_negative_thickness_raises():
    """Test that a negative thickness raises a ValueError."""
    with pytest.raises(ValueError):
        ShellGeometry(thickness=-200, material=Elastic2D(E=200000, nu=0.2))


def test_shell_reinforcement():
    """Test the ShellReinforcement class."""
    z = -60
    n_bars = 4
    cc_bars = 500
    d = 16
    material = Elastic2D(E=200000, nu=0.3)
    phi = np.pi / 4
    shell_reinforcement = ShellReinforcement(
        z=z,
        n_bars=n_bars,
        cc_bars=cc_bars,
        diameter_bar=d,
        material=material,
        phi=phi,
        name='rebars_1',
        group_label='group_A',
    )
    assert shell_reinforcement.z == z
    assert shell_reinforcement.n_bars == n_bars
    assert shell_reinforcement.cc_bars == cc_bars
    assert shell_reinforcement.diameter_bar == d
    assert shell_reinforcement.material == material
    assert shell_reinforcement.phi == phi
    assert shell_reinforcement.name == 'rebars_1'
    assert shell_reinforcement.group_label == 'group_A'


def test_add_reinforcement():
    """Test the add_reinforcement function."""
    # Create a shell geometry
    shell = ShellGeometry(thickness=200, material=Elastic2D(E=200000, nu=0.2))

    # Create a reinforcement
    reinf_1 = ShellReinforcement(
        z=-60,
        n_bars=4,
        cc_bars=500,
        diameter_bar=16,
        material=Elastic2D(E=200000, nu=0.3),
        phi=0,
    )

    reinf_2 = ShellReinforcement(
        z=60,
        n_bars=6,
        cc_bars=600,
        diameter_bar=20,
        material=Elastic2D(E=200000, nu=0.3),
        phi=np.pi / 4,
    )

    # Add a single reinforcement
    shell.add_reinforcement(reinf_1)
    assert len(shell.reinforcement) == 1
    assert shell.reinforcement[0] is reinf_1

    # Add reinforcement as a list
    shell.add_reinforcement([reinf_1, reinf_2])
    assert len(shell.reinforcement) == 3


def test_add_reinforcement_invalid_type():
    """Test that adding a reinforcement of invalid type raises an error."""
    # Create a shell geometry
    shell = ShellGeometry(thickness=200, material=Elastic2D(E=200000, nu=0.2))

    # Create a reinforcement
    reinf_1 = PointGeometry(
        np.array([2, 3]), 12, ElasticPlastic(E=200000, fy=500), name='Rebar'
    )
    with pytest.raises(TypeError):
        shell.add_reinforcement(reinf_1)


@pytest.mark.parametrize(
    'invalid_z',
    [
        95,  # barely outside of shell (half_thickness=100, d/2=6-> max=94)
        250,  # clearly outside the shell
    ],
)
def test_add_reinforcement_invalid_z_value(invalid_z):
    """Test that adding a reinforcement outside the shell raises an error."""
    shell = ShellGeometry(thickness=200, material=Elastic2D(E=200000, nu=0.2))
    reinf = ShellReinforcement(
        z=invalid_z,
        n_bars=1,
        cc_bars=100,
        diameter_bar=12,
        material=Elastic2D(E=200000, nu=0.3),
        phi=0,
    )
    with pytest.raises(ValueError):
        shell.add_reinforcement([reinf])


def test_repr_svg():
    """Test the SVG representation of the shell geometry."""
    # Create a shell geometry
    shell = ShellGeometry(thickness=200, material=Elastic2D(E=200000, nu=0.2))

    # Add a reinforcement
    reinf = ShellReinforcement(
        z=-60,
        n_bars=4,
        cc_bars=500,
        diameter_bar=16,
        material=Elastic2D(E=200000, nu=0.3),
        phi=0,
    )
    shell.add_reinforcement(reinf)
    # Get the SVG representation
    svg = shell._repr_svg_()
    assert '<svg' in svg
