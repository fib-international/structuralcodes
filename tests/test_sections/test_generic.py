"""Test initializing the deprecated GenericSection class."""

import pytest
from shapely import Polygon

from structuralcodes.geometry import SurfaceGeometry
from structuralcodes.materials.concrete import ConcreteMC2010
from structuralcodes.sections import BeamSection, GenericSection


# Test initializing a GenericSection
def test_deprecated_generic_section():
    """Test initializing the deprecated GenericSection."""
    # Create materials to use
    concrete = ConcreteMC2010(25)

    # The section
    poly = Polygon(((0, 0), (200, 0), (200, 400), (0, 400)))
    geo = SurfaceGeometry(poly, concrete)

    # Create the section
    with pytest.warns(DeprecationWarning):
        sec = GenericSection(geo)

    assert isinstance(sec, BeamSection)
