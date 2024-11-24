"""Tests for the Specific Section."""

import math

from sections.specific import RectangularSection

from structuralcodes.materials.concrete import ConcreteEC2_2004


# Test rectangular section
def test_specific_rectangular_section():
    """Test rectangular section."""
    # Create materials to use
    concrete = ConcreteEC2_2004(35)

    width, height = 250, 300
    sec = RectangularSection(width, height, concrete)

    assert sec.name == 'RectangularSection'

    # Testing and comparing area
    area: float = width * height
    assert math.isclose(sec.gross_properties.area, area)

    # Compute max / min axial load
    n_min_marin = sec.section_calculator.n_min

    # Testing max axial tension based off f_cd
    f_cd: float = -concrete.fcd()
    compression_capacity = area * f_cd
    assert math.isclose(compression_capacity, n_min_marin)
