"""Tests for EUROCODE 1992-3:2006."""

import math

import pytest

from structuralcodes.codes.ec2_2004 import (
    _cracking_restraint_imposed_deformations,
)


@pytest.mark.parametrize(
    'alpha_e, rho_p_eff, kc, k, fct_eff, Es, expected',
    [
        (5.869, 0.01478, 1, 1, 3.2, 200000, 0.000590),
        (5.869, 0.009690, 1, 1, 3.2, 200000, 0.000875),
        (5.512, 0.01478, 1, 1, 3.8, 200000, 0.000694),
        (5.512, 0.009690, 1, 1, 3.8, 200000, 0.001032),
    ],
)
def test_eps_sm_eps_cm_m1_returns_expected_values(
    alpha_e, kc, k, fct_eff, rho_p_eff, Es, expected
):
    """Test that eps_sm_cm_m1 returns expected values."""
    # Assert
    assert math.isclose(
        _cracking_restraint_imposed_deformations.eps_sm_eps_cm_restraint_end(
            alpha_e, rho_p_eff, kc, k, fct_eff, Es
        ),
        expected,
        abs_tol=10e-5,
    )
