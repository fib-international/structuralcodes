"""Tests for _AnnexB_TimeDependent module"""
import math

import pytest

from structuralcodes.codes.ec2_2023 import _annexB_time_dependent


@pytest.mark.parametrize(
    'test_input1, expected',
    [
        (20, 1.183),
        (24, 1.17),
        (28, 1.158),
        (33, 1.143),
        (38, 1.128),
        (43, 1.114),
        (48, 1.1),
        (53, 1.086),
        (58, 1.073),
        (63, 1.06),
        (68, 1.048),
        (78, 1.023),
        (88, 1),
        (98, 1),
        (108, 1),
    ],
)
def test_alphac(test_input1, expected):
    """Test the fctm function."""
    assert math.isclose(
        _annexB_time_dependent.alpha_c(test_input1), expected, rel_tol=0.01
    )
