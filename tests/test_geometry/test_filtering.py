"""Tests for geometry filtering methods."""

import pytest

from structuralcodes.geometry import (
    RectangularGeometry,
    add_reinforcement_line,
)
from structuralcodes.materials.basic import ElasticMaterial


@pytest.fixture
def reinforced_t_section_geometry():
    """Create a T section with two named surfaces and grouped reinforcement."""
    concrete = ElasticMaterial(E=30000, density=2450)
    reinforcement = ElasticMaterial(E=200000, density=7850)

    geo = RectangularGeometry(
        width=200,
        height=400,
        material=concrete,
        name='web',
    )
    geo += RectangularGeometry(
        width=800,
        height=150,
        material=concrete,
        origin=(0, 275),
        name='flange',
    )

    geo = add_reinforcement_line(
        geo,
        coords_i=(-60, -160),
        coords_j=(60, -160),
        diameter=16,
        n=4,
        material=reinforcement,
        group_label='bottom_reinforcement',
    )
    return add_reinforcement_line(
        geo,
        coords_i=(-360, 310),
        coords_j=(360, 310),
        diameter=12,
        n=12,
        material=reinforcement,
        group_label='top_reinforcement',
    )


def test_name_filter_none_returns_all_geometries_flat(
    reinforced_t_section_geometry,
):
    """Test filtering by name with pattern=None and flat return mode."""
    geo = reinforced_t_section_geometry

    filtered = geo.name_filter(None)

    assert len(filtered) == 18
    assert len(geo.geometries) == 2
    assert len(geo.point_geometries) == 16


def test_name_filter_none_returns_split_dict(reinforced_t_section_geometry):
    """Test filtering by name with pattern=None and split return mode."""
    geo = reinforced_t_section_geometry

    filtered = geo.name_filter(None, return_mode='split')

    assert set(filtered.keys()) == {'surfaces', 'points'}
    assert len(filtered['surfaces']) == 2
    assert len(filtered['points']) == 16


def test_name_filter_exact_surface_name(reinforced_t_section_geometry):
    """Test filtering by exact surface name."""
    geo = reinforced_t_section_geometry

    filtered = geo.name_filter('web', return_mode='split')

    assert len(filtered['surfaces']) == 1
    assert filtered['surfaces'][0].name == 'web'
    assert len(filtered['points']) == 0


def test_name_filter_wildcard_case_sensitive(reinforced_t_section_geometry):
    """Test wildcard filtering by name in case-sensitive mode."""
    geo = reinforced_t_section_geometry

    filtered = geo.name_filter('we*')

    assert len(filtered) == 1
    assert filtered[0].name == 'web'

    filtered = geo.name_filter('we?')

    assert len(filtered) == 1
    assert filtered[0].name == 'web'


def test_name_filter_wildcard_case_insensitive(reinforced_t_section_geometry):
    """Test wildcard filtering by name in case-insensitive mode."""
    geo = reinforced_t_section_geometry

    filtered = geo.name_filter('*geo*', case_sensitive=False)

    assert len(filtered) == 16
    assert all(g.name.startswith('Geometry_') for g in filtered)


def test_group_filter_wildcards(reinforced_t_section_geometry):
    """Test wildcard filtering by reinforcement group labels."""
    geo = reinforced_t_section_geometry

    all_reinforcement = geo.group_filter('*reinf*')
    bottom_reinforcement = geo.group_filter('bottom*')
    top_reinforcement = geo.group_filter('top*')

    assert len(all_reinforcement) == 16
    assert len(bottom_reinforcement) == 4
    assert len(top_reinforcement) == 12


def test_group_filter_none_returns_split_dict(reinforced_t_section_geometry):
    """Test filtering by group_label with pattern=None and split return."""
    geo = reinforced_t_section_geometry

    filtered = geo.group_filter(None, return_mode='split')

    assert set(filtered.keys()) == {'surfaces', 'points'}
    assert len(filtered['surfaces']) == 2
    assert len(filtered['points']) == 16


def test_group_filter_wildcards_returns_split_dict(
    reinforced_t_section_geometry,
):
    """Test filtering by group_labele and split return."""
    geo = reinforced_t_section_geometry

    filtered = geo.group_filter('top*', return_mode='split')

    assert set(filtered.keys()) == {'surfaces', 'points'}
    assert len(filtered['surfaces']) == 0
    assert len(filtered['points']) == 12


def test_group_filter_none_returns_all_geometries(
    reinforced_t_section_geometry,
):
    """Test filtering by group_label with pattern=None."""
    geo = reinforced_t_section_geometry

    filtered = geo.group_filter(None)

    assert len(filtered) == 18


def test_group_filter_case_sensitive_option(reinforced_t_section_geometry):
    """Test case sensitivity option when filtering by group_label."""
    geo = reinforced_t_section_geometry

    filtered_sensitive = geo.group_filter('*REINF*')
    filtered_insensitive = geo.group_filter('*REINF*', case_sensitive=False)

    assert len(filtered_sensitive) == 0
    assert len(filtered_insensitive) == 16


def test_combined_name_and_group_filters_intersection(
    reinforced_t_section_geometry,
):
    """Test combining name and group filters through intersection."""
    geo = reinforced_t_section_geometry

    geom_name = geo.name_filter('we*')
    geom_group = geo.group_filter(None)

    filtered = list(set(geom_name) & set(geom_group))

    assert len(filtered) == 1
    assert filtered[0].name == 'web'
