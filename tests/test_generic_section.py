"""Tests for the Geometry."""

import math

import numpy as np

import pytest

from structuralcodes.geometry import (
    PointGeometry,
    SurfaceGeometry,
)
from structuralcodes.sections._reinforcement import (
    add_reinforcement_line,
)
from structuralcodes.materials.concrete import ConcreteMC2010
from structuralcodes.materials.constitutive_laws import (
    ElasticPlastic,
    ParabolaRectangle,
)
from structuralcodes.sections._generic import GenericSection
from shapely import Polygon


# Test rectangular section
def test_rectangular_section():
    """Test rectangular section."""

    # Create materials to use
    concrete = ConcreteMC2010(25)
    steel = ElasticPlastic(210000, 450, eps_su=0.0675)

    # The section
    poly = Polygon(((0, 0), (200, 0), (200, 400), (0, 400)))
    geo = SurfaceGeometry(poly, concrete)
    geo = add_reinforcement_line(geo, (40, 40), (160, 40), 16, steel, n=4)
    geo = add_reinforcement_line(geo, (40, 360), (160, 360), 16, steel, n=4)
    geo = geo.translate(-100, -200)
    assert geo.geometries[0].centroid[0] == 0
    assert geo.geometries[0].centroid[1] == 0

    # Create the section (default Marin integrator)
    sec = GenericSection(geo)
    assert sec.name == 'GenericSection'

    assert math.isclose(sec.gross_properties.area, 200 * 400)

    # Compute bending strength
    res_marin = sec.section_analyzer.calculate_bending_strength(theta=0, n=0)

    # Compute moment curvature
    res_mc_marin = sec.section_analyzer.calculate_moment_curvature(
        theta=0, n=0
    )

    # Use fiber integration
    sec = GenericSection(geo, integrator='Fiber', mesh_size=0.0001)
    assert math.isclose(sec.gross_properties.area, 200 * 400)

    # Compute bending strength
    res_fiber = sec.section_analyzer.calculate_bending_strength(theta=0, n=0)

    assert math.isclose(res_marin.m_x, res_fiber.m_x, rel_tol=1e-2)

    # Compute moment curvature
    res_mc_fiber = sec.section_analyzer.calculate_moment_curvature(
        theta=0, n=0
    )

    assert math.isclose(res_mc_marin.my[-1], res_mc_fiber.my[-1], rel_tol=1e-2)
