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
        _section5_materials.fctm(test_input1), expected, rel_tol=0.01
    )


@pytest.mark.parametrize(
    'test_input1, test_input2, expected',
    [
        (20, 9500, 25787),
        (24, 9500, 27403),
        (24, 14000, 40383),
        (33, 9500, 30472),
        (38, 9500, 31939),
        (43, 9500, 33282),
        (48, 9500, 34525),
        (53, 9500, 35685),
        (58, 9500, 36773),
        (63, 9500, 37801),
        (68, 9500, 38776),
        (78, 9500, 40590),
        (88, 9500, 42256),
        (98, 9500, 43799),
        (108, 9500, 45241),
    ],
)
def test_Ecm(test_input1, test_input2, expected):
    """Test the fctm function."""
    assert math.isclose(
        _section5_materials.Ecm(test_input1, test_input2),
        expected,
        rel_tol=0.01,
    )
