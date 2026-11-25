"""Tests for ACI 318-19 one-way shear functions (Sec. 22.5)."""

import pytest

from structuralcodes.codes.aci318_19 import _shear as sh

IN = 25.4  # mm per in


class TestLambdaS:
    """Eq. (22.5.5.1.3): sqrt(2 / (1 + d/10)) <= 1, d in inches."""

    def test_d_20_in(self):
        """Hand calculation: sqrt(2 / 3) = 0.8165."""
        assert sh.lambda_s(20.0 * IN) == pytest.approx(0.8165, rel=1e-4)

    def test_d_40_in(self):
        """Hand calculation: sqrt(2 / 5) = 0.6325."""
        assert sh.lambda_s(40.0 * IN) == pytest.approx(0.6325, rel=1e-4)

    def test_d_10_in_is_unity(self):
        """Hand calculation: sqrt(2 / 2) = 1."""
        assert sh.lambda_s(10.0 * IN) == pytest.approx(1.0, rel=1e-9)

    def test_shallow_member_capped_at_one(self):
        """With d = 5 in the expression gives 1.155, limited to 1."""
        assert sh.lambda_s(5.0 * IN) == 1.0

    def test_non_positive_depth_raises(self):
        """A depth of zero is not a section."""
        with pytest.raises(ValueError):
            sh.lambda_s(0.0)
