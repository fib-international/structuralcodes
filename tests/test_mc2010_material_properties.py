"""Tests for the _concrete_material_properties module"""
import math

import pytest

from structuralcodes.codes.mc2010 import _concrete_material_properties


@pytest.mark.parametrize(
    'test_input, expect',
    [(12, 20), (35, 43), (55, 63), (90, 98), (120, 128)],
)
def test_fcm(test_input, expect):
    """Test the fcm function."""
    assert math.isclose(_concrete_material_properties.fcm(test_input), expect)


@pytest.mark.parametrize(
    'test_input, expected',
    [
        (12, 1.6),
        (16, 1.9),
        (20, 2.2),
        (25, 2.6),
        (30, 2.9),
        (35, 3.2),
        (40, 3.5),
        (45, 3.8),
        (50, 4.1),
        (55, 4.2),
        (60, 4.4),
        (70, 4.6),
        (80, 4.8),
        (90, 5.0),
        (100, 5.2),
        (110, 5.4),
        (120, 5.6),
    ],
)
def test_fctm(test_input, expected):
    """Test the fctm function."""
    assert math.isclose(
        _concrete_material_properties.fctm(test_input), expected, rel_tol=0.02
    )


@pytest.mark.parametrize(
    'test_input, expected',
    [
        (12, 1.1),
        (16, 1.3),
        (20, 1.5),
        (25, 1.8),
        (30, 2.0),
        (35, 2.2),
        (40, 2.5),
        (45, 2.7),
        (50, 2.9),
        (55, 3.0),
        (60, 3.1),
        (70, 3.2),
        (80, 3.4),
        (90, 3.5),
        (100, 3.7),
        (110, 3.8),
        (120, 3.9),
    ],
)
def test_fctkmin(test_input, expected):
    """Test the fctkmin function."""
    assert math.isclose(
        _concrete_material_properties.fctkmin(
            _concrete_material_properties.fctm(test_input)
        ),
        expected,
        rel_tol=0.031,
    )


@pytest.mark.parametrize(
    'test_input, expected',
    [
        (12, 2.0),
        (16, 2.5),
        (20, 2.9),
        (25, 3.3),
        (30, 3.8),
        (35, 4.2),
        (40, 4.6),
        (45, 4.8),
        (50, 5.3),
        (55, 5.5),
        (60, 5.7),
        (70, 6.0),
        (80, 6.3),
        (90, 6.6),
        (100, 6.8),
        (110, 7.0),
        (120, 7.2),
    ],
)
def test_fctkmax(test_input, expected):
    """Test the fctkmax function."""
    assert math.isclose(
        _concrete_material_properties.fctkmax(
            _concrete_material_properties.fctm(test_input)
        ),
        expected,
        rel_tol=0.028,
    )


@pytest.mark.parametrize(
    'test_input, expected',
    [
        (12, 125.172),
        (35, 143.664),
        (55, 153.888),
        (90, 166.626),
        (120, 174.832),
    ],
)
def test_Gf(test_input, expected):
    """Test the Gf function."""
    assert math.isclose(
        _concrete_material_properties.Gf(test_input),
        expected,
        rel_tol=1e-5,
    )
