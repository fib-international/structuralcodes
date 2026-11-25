"""Tests for the Annular Geometry."""

import math

import numpy as np
import pytest

from structuralcodes.geometry import AnnularGeometry, add_reinforcement_circle
from structuralcodes.materials.basic import ElasticMaterial
from structuralcodes.materials.concrete import ConcreteMC2010
from structuralcodes.materials.reinforcement import ReinforcementMC2010


# Test create a annular geometry
@pytest.mark.parametrize(
    'outer_diameter, inner_diameter',
    [(100, 80), (200, 160), (300, 250), (400, 300), (500, 320)],
)
def test_create_annular_geometry(outer_diameter, inner_diameter):
    """Test creating a AnnularGeometry."""
    mat = ElasticMaterial(E=300000, density=2450)
    annular = AnnularGeometry(outer_diameter, inner_diameter, mat)

    assert math.isclose(annular.outer_diameter, outer_diameter)
    assert math.isclose(annular.outer_radius, outer_diameter / 2.0)
    assert math.isclose(annular.inner_diameter, inner_diameter)
    assert math.isclose(annular.inner_radius, inner_diameter / 2.0)
    assert math.isclose(
        annular.thickness, (outer_diameter - inner_diameter) / 2.0
    )

    area = annular.area
    expected_area = np.pi * (outer_diameter**2 - inner_diameter**2) / 4.0
    assert math.isclose(area, expected_area, rel_tol=0.05)


# Test wrong input
@pytest.mark.parametrize(
    'outer_diameter, inner_diameter',
    [
        (-100, -100),  # Wrong: both negative
        (-100, 100),  # Wrong: at least one negative
        (100, -50),  # Wrong: at least one negative
        (200, 300),  # Wrong: inner > outer
        (0, 100),  # Wrong: at least one zero
        (100, 0),  # Wrong: at least one zero
        (100, 100),  # Wrong: outer == inner
    ],
)
def test_create_annular_geometry_exception(outer_diameter, inner_diameter):
    """Test raising exception when inputing wrong value."""
    mat = ElasticMaterial(E=300000, density=2450)
    with pytest.raises(ValueError):
        AnnularGeometry(outer_diameter, inner_diameter, mat)


# Test providing origin
@pytest.mark.parametrize(
    'origin',
    [(0, 0), (30, 600), (1.0,), None],
)
def test_annulus_with_origin(origin):
    """Test creating an annulus with an origin."""
    # Arrange
    mat = ElasticMaterial(E=300000, density=2450)
    outer_diameter = 100
    inner_diameter = 80
    thickness = (outer_diameter - inner_diameter) / 2.0

    # Act and assert
    if origin is not None and len(origin) == 2:
        geom = AnnularGeometry(
            outer_diameter=outer_diameter,
            inner_diameter=inner_diameter,
            material=mat,
            origin=origin,
        )
        assert np.allclose(
            (geom.polygon.centroid.xy[0][0], geom.polygon.centroid.xy[1][0]),
            origin,
        )
        assert math.isclose(thickness, geom.thickness)
    elif origin is None:
        geom = AnnularGeometry(
            outer_diameter=outer_diameter,
            inner_diameter=inner_diameter,
            material=mat,
            origin=origin,
        )
        assert np.allclose(
            (geom.polygon.centroid.xy[0][0], geom.polygon.centroid.xy[1][0]),
            (0, 0),
        )
        assert math.isclose(thickness, geom.thickness)
    else:
        with pytest.raises(ValueError):
            geom = AnnularGeometry(
                outer_diameter=outer_diameter,
                inner_diameter=inner_diameter,
                material=mat,
                origin=origin,
            )


# Test adding reinforcement in circular pattern
@pytest.mark.parametrize(
    'outer_diameter, inner_diameter, cover, n',
    [(500, 300, 30, 16), (600, 300, 40, 32)],
)
def test_reinforced_annular_geometry(outer_diameter, inner_diameter, cover, n):
    """Test adding reinforcemnet to a AnnularGeometry."""
    # Create materials
    concrete = ConcreteMC2010(40)
    steel = ReinforcementMC2010(450, 200000, 450, 0.065)

    # create the annular geometry
    annular = AnnularGeometry(outer_diameter, inner_diameter, concrete)

    assert annular.concrete

    # add reinforcement in the circular pattern
    rc = add_reinforcement_circle(
        geo=annular,
        center=(0, 0),
        radius=annular.outer_radius - cover,
        diameter=8,
        material=steel,
        n=n,
    )
    rc = add_reinforcement_circle(
        geo=rc,
        center=(0, 0),
        radius=annular.inner_radius + cover,
        diameter=10,
        material=steel,
        n=n,
    )

    assert len(rc.point_geometries) == 2 * n


def test_annular_geometry_name_group_label():
    """Test the name and group label attribute of a AnnularGeometry."""
    # Arrange
    outer_diameter = 500
    inner_diameter = 320
    concrete = ConcreteMC2010(fck=35)
    name = 'concrete_geometry'
    group_label = 'concrete'

    # Act
    geometry = AnnularGeometry(
        outer_diameter=outer_diameter,
        inner_diameter=inner_diameter,
        material=concrete,
        name=name,
        group_label=group_label,
    )

    # Assert
    assert geometry.name == name
    assert geometry.group_label == group_label
