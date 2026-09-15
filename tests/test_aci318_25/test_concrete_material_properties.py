"""Tests for ACI 318-25 concrete material property functions (Ch. 19)."""

import math

import pytest

from structuralcodes.codes.aci318_25 import (
    _concrete_material_properties as cmp,
)


class TestEc:
    """Tests for the modulus of elasticity Ec (Table 19.2.2.1)."""

    def test_normalweight_4000psi(self):
        """4000 psi = 27.58 MPa, wc = 2320 kg/m3 (normalweight)."""
        fc = 27.58
        wc = 2320.0
        expected = (wc**1.5) * 0.043 * math.sqrt(fc)
        assert math.isclose(cmp.Ec(fc, wc), expected, rel_tol=1e-6)

    def test_normalweight_28mpa(self):
        """Fc = 28 MPa, default wc."""
        expected = (2320.0**1.5) * 0.043 * math.sqrt(28.0)
        assert math.isclose(cmp.Ec(28.0), expected, rel_tol=1e-6)

    def test_custom_wc(self):
        """Custom unit weight wc = 1800 kg/m3."""
        fc = 30.0
        wc = 1800.0
        expected = (wc**1.5) * 0.043 * math.sqrt(fc)
        assert math.isclose(cmp.Ec(fc, wc), expected, rel_tol=1e-6)

    def test_invalid_fc_zero(self):
        """Fc = 0 should raise ValueError."""
        with pytest.raises(ValueError):
            cmp.Ec(0.0)

    def test_invalid_fc_negative(self):
        """Negative fc should raise ValueError."""
        with pytest.raises(ValueError):
            cmp.Ec(-10.0)

    def test_invalid_wc_too_low(self):
        """Wc below 1440 should raise ValueError."""
        with pytest.raises(ValueError):
            cmp.Ec(28.0, wc=1400.0)

    def test_invalid_wc_too_high(self):
        """Wc above 2560 should raise ValueError."""
        with pytest.raises(ValueError):
            cmp.Ec(28.0, wc=2600.0)


class TestFr:
    """Tests for the modulus of rupture fr (Eq. 19.2.3.1)."""

    def test_normalweight(self):
        """Normalweight concrete, lambda_s = 1.0."""
        fc = 28.0
        expected = 0.62 * 1.0 * math.sqrt(fc)
        assert math.isclose(cmp.fr(fc), expected, rel_tol=1e-6)

    def test_lightweight(self):
        """Lightweight concrete, lambda_s = 0.75."""
        fc = 28.0
        lambda_s = 0.75
        expected = 0.62 * lambda_s * math.sqrt(fc)
        assert math.isclose(
            cmp.fr(fc, lambda_s=lambda_s), expected, rel_tol=1e-6
        )

    def test_invalid_fc(self):
        """Fc <= 0 should raise ValueError."""
        with pytest.raises(ValueError):
            cmp.fr(0.0)

    def test_invalid_lambda_zero(self):
        """lambda_s = 0 should raise ValueError (must be > 0)."""
        with pytest.raises(ValueError):
            cmp.fr(28.0, lambda_s=0.0)

    def test_invalid_lambda_above_one(self):
        """lambda_s > 1 should raise ValueError."""
        with pytest.raises(ValueError):
            cmp.fr(28.0, lambda_s=1.1)


class TestBeta1:
    """Tests for the Whitney stress block factor beta1 (Table 22.2.2.4.3)."""

    @pytest.mark.parametrize(
        'fc, expected',
        [
            (21.0, 0.85),
            (27.58, 0.85),
            (34.47, 0.85 - 0.05 * (34.47 - 28) / 7),
            (55.16, 0.65),
            (68.95, 0.65),
        ],
    )
    def test_beta1_parametric(self, fc, expected):
        """Test beta1 for various fc values."""
        assert math.isclose(cmp.beta1(fc), expected, rel_tol=1e-6)

    def test_invalid_fc(self):
        """Fc <= 0 should raise ValueError."""
        with pytest.raises(ValueError):
            cmp.beta1(0.0)


class TestEpsCu:
    """Tests for the ultimate concrete strain eps_cu (Sec. 22.2.2.1)."""

    def test_value(self):
        """Ultimate strain must be 0.003."""
        assert cmp.eps_cu() == 0.003


class TestAlpha1:
    """Tests for the stress block intensity factor alpha1 (Sec. 22.2.2.4.1)."""

    def test_value(self):
        """Stress block intensity must be 0.85."""
        assert cmp.alpha1() == 0.85


class TestFct:
    """Tests for the splitting tensile strength fct (Sec. 19.2.4.3)."""

    def test_normalweight(self):
        """Normalweight concrete, lambda_s = 1.0."""
        fc = 28.0
        expected = 0.56 * 1.0 * math.sqrt(fc)
        assert math.isclose(cmp.fct(fc), expected, rel_tol=1e-6)

    def test_invalid_fc(self):
        """Fc <= 0 should raise ValueError."""
        with pytest.raises(ValueError):
            cmp.fct(0.0)


class TestLambdaFactor:
    """Tests for the lightweight concrete factor lambda_factor.

    Reference: Table 19.2.4.2.
    """

    @pytest.mark.parametrize(
        'concrete_type, expected',
        [
            ('normalweight', 1.0),
            ('sand-lightweight', 0.85),
            ('all-lightweight', 0.75),
        ],
    )
    def test_lambda_factor_parametric(self, concrete_type, expected):
        """Test lambda_factor for all defined concrete types."""
        assert math.isclose(
            cmp.lambda_factor(concrete_type), expected, rel_tol=1e-9
        )

    def test_invalid_type(self):
        """Unknown concrete type should raise ValueError."""
        with pytest.raises(ValueError):
            cmp.lambda_factor('medium-lightweight')
