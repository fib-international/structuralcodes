import math

import pytest

from structuralcodes.code.mc2010 import _concrete_shear


@pytest.mark.parametrize(
    'test_input1, testinput2, test, expected',
    [
        (12, 150, 50, 2625.424),
        (14, 150, 50, 2835.782),
        (16, 150, 50, 3031.579),
        (20, 150, 50, 3389.408),
        (25, 150, 50, 3789.474),
        (30, 150, 50, 4151.160)
    ]
)
def test_fcm(test_input1, test_input2, test_input3, expected):
    """Test the fcm function."""
    assert math.isclose(
        _concrete_shear.vrdc(test_input1, test_input2, test_input3, expected)
    )

