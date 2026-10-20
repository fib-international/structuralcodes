"""Test for Functions from Section 8.6 of EN 1992-1-1:2023."""

import pytest

from structuralcodes.codes.ec2_2023 import _section_8_6_partially_loaded_areas


@pytest.mark.parametrize(
    'f_cd, a0, b0, a1, b, ea, eb, nu_part, expected',
    [
        (25, 300, 200, 350, 500, 0, 0, 3.0, 30.19),
        (30, 400, 300, 450, 600, 20, 10, 3.0, 37.5),
        (40, 500, 400, 550, 700, 30, 20, 3.0, 50.0),
    ],
)
def test_calculate_design_resistance(
    f_cd, a0, b0, a1, b, ea, eb, nu_part, expected
):
    """Test calculate_design_resistance with various
    inputs to verify correctness.
    """
    assert (
        pytest.approx(
            _section_8_6_partially_loaded_areas.sigma_Rd_u(
                f_cd, a0, b0, a1, b, ea, eb, nu_part
            ),
            rel=1e-2,
        )
        == expected
    )
