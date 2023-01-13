"""Tests for the _section9_sls module"""
import math

import pytest

from structuralcodes.codes.ec2_2023 import _section9_sls


@pytest.mark.parametrize(
    'test_input1, test_input2, expected',
    [(33, 1.5, 12188.63), (33, 2, 10157.19), (33, 2.5, 8706.16)],
)
def test_Ec_eff(test_input1, test_input2, expected):
    """Test the Ec_eff function."""
    assert math.isclose(
        _section9_sls.Ec_eff(test_input1, test_input2), expected, rel_tol=0.05
    )


@pytest.mark.parametrize(
    'test_input1, test_input2, test_input3, test_input4, test_input5,'
    + 'expected1, expected2',
    [
        (-1000, 1, 0.15, 2.565, 500, 0, 0),
        (-500, 1, 0.15, 2.565, 500, 0, 0),
        (0, 1, 0.15, 2.565, 500, 1.231, 0),
        (500, 1, 0.15, 2.565, 500, 3.078, 3.078),
        (3000, 1, 0.15, 2.565, 500, 3.078, 3.078),
        (-1000, 1, 0.5, 2.565, 400, 0, 0),
        (-500, 1, 0.5, 2.565, 400, 0.611, 0),
        (0, 1, 0.5, 2.565, 400, 4.361, 0),
        (500, 1, 0.5, 2.565, 400, 8.111, 4.389),
        (1000, 1, 0.5, 2.565, 400, 10.901, 10.901),
        (-1000, 1, 0.5, 3.795, 500, 0, 0),
        (-500, 1, 0.8, 3.795, 500, 3.072, 0),
        (0, 1, 0.8, 3.795, 500, 6.072, 0),
        (500, 1, 0.8, 3.795, 500, 9.072, 0.928),
        (1000, 1, 0.8, 3.795, 500, 12.072, 7.928),
    ],
)
def test_As_min_y(
    test_input1,
    test_input2,
    test_input3,
    test_input4,
    test_input5,
    expected1,
    expected2,
):
    """Test the As_min_y function."""
    (As_min_y1, As_min_y2) = _section9_sls.As_min_y(
        test_input1, test_input2, test_input3, test_input4, test_input5
    )

    assert math.isclose(As_min_y1, expected1, rel_tol=0.05)
    assert math.isclose(As_min_y2, expected2, rel_tol=0.05)


@pytest.mark.parametrize(
    'test_input1, test_input2, expected',
    [
        (0.15, 0.15, 0.8),
        (0.15, 0.25, 0.8),
        (0.25, 0.15, 0.8),
        (0.25, 0.25, 0.8),
        (0.3, 0.3, 0.8),
        (0.4, 0.4, 0.74),
        (0.5, 0.5, 0.68),
        (0.6, 0.6, 0.62),
        (0.7, 0.7, 0.56),
        (0.8, 0.8, 0.5),
        (0.9, 0.9, 0.5),
        (1, 1, 0.5),
        (1.1, 1.1, 0.5),
    ],
)
def test_kh(test_input1, test_input2, expected):
    """Test the kh function."""
    assert math.isclose(
        _section9_sls.kh(test_input1, test_input2), expected, rel_tol=0.05
    )
