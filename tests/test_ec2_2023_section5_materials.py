"""Tests for _section5_materials module"""
import math

import pytest

from structuralcodes.codes.ec2_2023 import _section5_materials


@pytest.mark.parametrize(
    'test_input1, expected',
    [
        (12, 1.572),
        (16, 1.905),
        (20, 2.21),
        (25, 2.565),
        (30, 2.896),
        (35, 3.21),
        (40, 3.509),
        (45, 3.795),
        (50, 4.072),
        (55, 4.183),
        (60, 4.306),
        (70, 4.533),
        (80, 4.74),
        (90, 4.93),
        (100, 5.106),
    ],
)
def test_fctm(test_input1, expected):
    """Test the fctm function."""
    assert math.isclose(
        _section5_materials.fctm(test_input1), expected, rel_tol=0.05
    )
