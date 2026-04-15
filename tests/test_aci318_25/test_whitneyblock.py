"""Tests for WhitneyBlock constitutive law."""

import math

import numpy as np
import pytest

from structuralcodes.materials.constitutive_laws._whitneyblock import (
    WhitneyBlock,
)


@pytest.fixture
def wb():
    """4000 psi concrete: fc = 0.85 * 27.58 = 23.44, beta1 = 0.85."""
    return WhitneyBlock(fc=23.44, beta1=0.85, eps_cu=0.003)


class TestGetStress:
    """Tests for get_stress method."""

    def test_in_active_zone(self, wb):
        """Strain in the active zone returns compressive stress."""
        sig = wb.get_stress(-0.002)
        assert math.isclose(sig, -23.44, rel_tol=1e-10)

    def test_in_zero_zone(self, wb):
        """Strain between transition and zero returns 0.

        eps_transition = -0.003 * (1 - 0.85) = -0.00045
        eps = -0.0002 is between -0.00045 and 0, so outside active zone.
        """
        sig = wb.get_stress(-0.0002)
        assert sig == 0.0

    def test_positive_strain(self, wb):
        """Positive strain (tension) returns 0."""
        sig = wb.get_stress(0.001)
        assert sig == 0.0

    def test_beyond_ultimate(self, wb):
        """Strain beyond ultimate (more negative than eps_cu) returns 0."""
        sig = wb.get_stress(-0.004)
        assert sig == 0.0

    def test_array_input(self, wb):
        """Array of strains returns correct stress array."""
        eps = np.array([-0.004, -0.002, -0.0002, 0.001])
        sig = wb.get_stress(eps)
        expected = np.array([0.0, -23.44, 0.0, 0.0])
        np.testing.assert_allclose(sig, expected, rtol=1e-10)

    def test_at_eps_cu_boundary(self, wb):
        """Strain exactly at eps_cu should be in active zone."""
        sig = wb.get_stress(-0.003)
        assert math.isclose(sig, -23.44, rel_tol=1e-10)

    def test_at_transition_boundary(self, wb):
        """Strain exactly at eps_transition should be in active zone."""
        # Use the actual computed transition value to avoid fp issues
        sig = wb.get_stress(wb._eps_transition)
        assert math.isclose(sig, -23.44, rel_tol=1e-10)


class TestGetUltimateStrain:
    """Tests for get_ultimate_strain method."""

    def test_ultimate_strain(self, wb):
        """Ultimate strain returns (eps_cu, 0.0)."""
        eps_min, eps_max = wb.get_ultimate_strain()
        assert math.isclose(eps_min, -0.003, rel_tol=1e-10)
        assert eps_max == 0.0

    def test_ultimate_strain_yielding(self, wb):
        """Yielding flag does not change result."""
        eps_min, eps_max = wb.get_ultimate_strain(yielding=True)
        assert math.isclose(eps_min, -0.003, rel_tol=1e-10)
        assert eps_max == 0.0


class TestGetTangent:
    """Tests for get_tangent method."""

    def test_in_active_zone(self, wb):
        """Tangent in active zone is 0."""
        tangent = wb.get_tangent(-0.002)
        assert tangent == 0.0

    def test_at_zero(self, wb):
        """Tangent at zero strain is 0."""
        tangent = wb.get_tangent(0.0)
        assert tangent == 0.0

    def test_array_input(self, wb):
        """Tangent for array input is all zeros."""
        eps = np.array([-0.003, -0.001, 0.0, 0.001])
        tangent = wb.get_tangent(eps)
        np.testing.assert_array_equal(tangent, np.zeros(4))


class TestMarin:
    """Tests for __marin__ method."""

    def test_uniform_strain_in_active_zone(self, wb):
        """Uniform strain in active zone returns fc coefficient."""
        strains, coeff = wb.__marin__([-0.002, 0])
        assert strains is None
        assert len(coeff) == 1
        assert math.isclose(coeff[0][0], -23.44, rel_tol=1e-10)

    def test_uniform_strain_outside_active_zone(self, wb):
        """Uniform strain outside active zone returns zero coefficient."""
        strains, coeff = wb.__marin__([0.001, 0])
        assert strains is None
        assert len(coeff) == 1
        assert coeff[0][0] == 0.0

    def test_varying_strain(self, wb):
        """Varying strain returns two zones with correct limits."""
        strains, coeff = wb.__marin__([-0.003, 0.001])
        assert strains is not None
        assert len(strains) == 2
        assert len(coeff) == 2

        # First zone: active stress zone
        assert math.isclose(strains[0][0], -0.003, rel_tol=1e-10)
        assert math.isclose(strains[0][1], -0.00045, rel_tol=1e-10)
        assert math.isclose(coeff[0][0], -23.44, rel_tol=1e-10)

        # Second zone: zero stress zone
        assert math.isclose(strains[1][0], -0.00045, rel_tol=1e-10)
        assert strains[1][1] == 0
        assert coeff[1][0] == 0.0


class TestMarinTangent:
    """Tests for __marin_tangent__ method."""

    def test_uniform_strain(self, wb):
        """Uniform strain returns zero tangent."""
        strains, coeff = wb.__marin_tangent__([-0.002, 0])
        assert strains is None
        assert len(coeff) == 1
        assert coeff[0][0] == 0.0

    def test_varying_strain(self, wb):
        """Varying strain returns single zone with zero tangent."""
        strains, coeff = wb.__marin_tangent__([-0.003, 0.001])
        assert strains is not None
        assert len(strains) == 1
        assert len(coeff) == 1
        assert math.isclose(strains[0][0], -0.003, rel_tol=1e-10)
        assert strains[0][1] == 0
        assert coeff[0][0] == 0.0


class TestConstructor:
    """Tests for constructor behavior."""

    def test_default_name(self):
        """Default name is WhitneyBlock."""
        wb = WhitneyBlock(fc=23.44, beta1=0.85)
        assert wb.name == 'WhitneyBlock'

    def test_custom_name(self):
        """Custom name is used."""
        wb = WhitneyBlock(fc=23.44, beta1=0.85, name='MyBlock')
        assert wb.name == 'MyBlock'

    def test_fc_stored_negative(self):
        """Fc is stored as negative (compression)."""
        wb = WhitneyBlock(fc=23.44, beta1=0.85)
        assert wb._fc < 0
        assert math.isclose(wb._fc, -23.44, rel_tol=1e-10)

    def test_eps_cu_stored_negative(self):
        """eps_cu is stored as negative (compression)."""
        wb = WhitneyBlock(fc=23.44, beta1=0.85, eps_cu=0.003)
        assert wb._eps_cu < 0
        assert math.isclose(wb._eps_cu, -0.003, rel_tol=1e-10)

    def test_eps_transition(self):
        """eps_transition is computed correctly."""
        wb = WhitneyBlock(fc=23.44, beta1=0.85, eps_cu=0.003)
        expected = -0.003 * (1.0 - 0.85)
        assert math.isclose(wb._eps_transition, expected, rel_tol=1e-10)

    def test_default_eps_cu(self):
        """Default eps_cu is 0.003."""
        wb = WhitneyBlock(fc=23.44, beta1=0.85)
        assert math.isclose(wb._eps_cu, -0.003, rel_tol=1e-10)

    def test_materials_attribute(self):
        """__materials__ includes concrete."""
        assert 'concrete' in WhitneyBlock.__materials__
