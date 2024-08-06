"""Test for functions from Section 11 of EN 1992-1-1:2023."""

import pytest

from structuralcodes.codes.ec2_2023 import _section11_detailing


@pytest.mark.parametrize(
    'phi, D_upper, expected',
    [
        (12, 15, 20),
        (10, 14, 20),
        (25, 10, 25),
    ],
)
def test_clear_distance_between_bars(phi, D_upper, expected):
    """Test clear distance between parallel bars."""
    assert (
        _section11_detailing.min_clear_distance_between_bars(phi, D_upper)
        == expected
    )


@pytest.mark.parametrize(
    'surface_roughness, min_bond_distance, expected',
    [
        (True, 8, 5),
        (False, 8, 8),
    ],
)
def test_clear_distance_to_poured_concrete(
    surface_roughness, min_bond_distance, expected
):
    """Test clear distance between poured concrete and parallel bar."""
    assert (
        _section11_detailing.min_clear_distance_to_poured_concrete(
            surface_roughness, min_bond_distance
        )
        == expected
    )


@pytest.mark.parametrize(
    'phi, expected',
    [
        (10, 40),
        (20, 140),
    ],
)
def test_mandrel_diameter(phi, expected):
    """Test minimum mandrel diameter for bending bars."""
    assert _section11_detailing.min_mandrel_phi(phi) == expected


@pytest.mark.parametrize(
    'fck, gamma_c, d_g, phi, phi_mand, c_d, k_bend, expected',
    [
        (30, 1.5, 16, 10, 40, 20, 32, 423.559),
        (40, 1.5, 20, 15, 40, 20, 32, 334.359),
    ],
)
def test_steel_stress_limit(
    fck, gamma_c, d_g, phi, phi_mand, c_d, k_bend, expected
):
    """Test steel stress limit to avoid concrete failure inside the bend."""
    result = _section11_detailing.sigma_sd_bend_concrete(
        fck, gamma_c, d_g, phi, phi_mand, c_d, k_bend
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'n_trans, phi, phi_mand, phi_trans, alpha_bend, expected',
    [
        (2, 10, 40, 10, 45, 3),
        (5, 20, 50, 15, 45, 5.5),
    ],
)
def test_k_trans_factor(
    n_trans, phi, phi_mand, phi_trans, alpha_bend, expected
):
    """Test increase factor for steel stress limit with transverse bars."""
    result = _section11_detailing.k_trans(
        n_trans, phi, phi_mand, phi_trans, alpha_bend
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'alpha_bend, expected',
    [
        (45, 32),
        (30, 48),
    ],
)
def test_k_bend(alpha_bend, expected):
    """Test test_k_bend."""
    result = _section11_detailing.k_bend(alpha_bend)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'fck, phi, bond_condition, expected',
    [
        (20, 12, 'good', 564),
        (30, 16, 'good', 672),
        (40, 25, 'good', 1075),
        (25, 28, 'poor', 1881.60),
        (28, 14, 'good', 590.800),
        (50, 32, 'good', 1312),
        (35, 20, 'poor', 1008.0),
    ],
)
def test_lbd_simple(fck, phi, bond_condition, expected):
    """Test anchorage length calculation with
    various conditions including bilinear interpolation.
    """
    result = _section11_detailing.lbd_simple(fck, phi, bond_condition)
    assert pytest.approx(result, rel=1e-2) == expected


@pytest.mark.parametrize(
    'phi, fck, sigma_sd, cd, c, bond_condition, design_situation, expected',
    [
        (20, 30, 400, 30, 20, 'good', 'persistent', 1138.362),
        (25, 40, 435, 40, 20, 'poor', 'persistent', 1749.1829),
        (32, 50, 400, 50, 20, 'bentonite', 'persistent', 2530.665),
        (28, 25, 350, 28, 20, 'good', 'accidental', 1370.4532),
    ],
)
def test_lbd(
    phi, fck, sigma_sd, cd, c, bond_condition, design_situation, expected
):
    """Test design anchorage length calculation with various conditions."""
    result = _section11_detailing.lbd(
        phi, fck, sigma_sd, cd, c, bond_condition, design_situation
    )
    assert pytest.approx(result, rel=1e-2) == expected
