"""Tests for the reinforcement factory."""

import pytest

from structuralcodes.codes import set_design_code
from structuralcodes.materials.reinforcement import (
    ReinforcementEC2_2004,
    ReinforcementEC2_2023,
    ReinforcementMC2010,
    create_reinforcement,
)


@pytest.mark.parametrize(
    'design_code, expected_reinforcement',
    [
        (None, None),
        ('mc2010', ReinforcementMC2010),
        ('ec2_2004', ReinforcementEC2_2004),
        ('ec2_2023', ReinforcementEC2_2023),
    ],
)
def test_reinforcement_factory(design_code, expected_reinforcement):
    """Test the reinforcement factory."""
    # Arrange
    fyk = 500
    Es = 200000

    # Act and assert
    if design_code is not None:
        set_design_code(design_code=design_code)
        reinforcement = create_reinforcement(fyk=fyk, Es=Es)
        assert isinstance(reinforcement, expected_reinforcement)
    else:
        set_design_code(design_code=design_code)
        with pytest.raises(ValueError):
            create_reinforcement(fyk=fyk, Es=Es)
