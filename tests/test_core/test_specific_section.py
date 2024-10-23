"""Tests for the Specific Section."""
from structuralcodes.materials.concrete import ConcreteEC2_2004
from structuralcodes.materials.reinforcement import ReinforcementEC2_2004


# Test rectangular section
def test_specific_rectangular_section():
    """Test rectangular section."""
    # Create materials to use
    concrete = ConcreteEC2_2004(35)
    steel = ReinforcementEC2_2004(fyk=500, Es=200000, ftk=450, epsuk=0.0675)