"""Test for the function of _concrete_punching."""

import math

import pytest

from structuralcodes.codes.mc2010 import _concrete_punching


def test_b0():
    """Test the b_0 function."""
    assert math.isclose(
        _concrete_punching.b_0(v_ed=100, v_prep_d_max=2.3),
        43.47,  # math check
        rel_tol=0.001,
    )
    assert _concrete_punching.b_0(v_ed=1000, v_prep_d_max=0.001) == 1000000
    assert _concrete_punching.b_0(v_ed=0, v_prep_d_max=1) == 0


@pytest.mark.parametrize(
    'v_ed, e_u, l_x, l_y, inner, edge_par, edge_per, corner, expected',
    [
        (10e3, 20, 2e3, 2e3, True, False, False, False, 1401),  # Inner
        (10e3, 20, 2e3, 3e3, False, True, False, False, 2500),  # Edge par
        (10e3, 20, 2e3, 3e3, False, False, True, False, 1497),  # Edge perp
        (10e3, 20, 2e3, 3e3, False, False, False, True, 5000),  # Corner
        # Edge cases
        (10e3, 0, 2e3, 2e3, True, False, False, False, 1250),  # No ecc
        (10e3, 100, 2e3, 2e3, True, False, False, False, 2007.58),  # Large ecc
    ],
)
def test_m_ed(
    v_ed, e_u, l_x, l_y, inner, edge_par, edge_per, corner, expected
):
    """Test the m_ed function."""
    assert math.isclose(
        _concrete_punching.m_ed(
            v_ed, e_u, l_x, l_y, inner, edge_par, edge_per, corner
        ),
        expected,
        rel_tol=0.001,
    )


def test_m_ed_invalid():
    """Test m_ed with invalid inputs."""
    with pytest.raises(ValueError, match='Placement is not defined'):
        _concrete_punching.m_ed(
            10e3, 20, 2e3, 2e3, False, False, False, False
        )


@pytest.mark.parametrize(
    (
        'l_x, l_y, f_yd, d_eff, e_s, approx_lvl_p, m_rd, m_ed, '
        'x_direction, expected'
    ),
    [
        # Level 1 approximation
        (2e3, 3e3, 434, 160, 200e3, 1, 140, 1400, True, 0.013427),
        # Level 2 approximation
        (2e3, 3e3, 434, 160, 200e3, 2, 140, 1500, False, 0.47089),
        # Edge cases
        (1e3, 1e3, 434, 160, 200e3, 1, 140, 1400, True, 0.004476),
        (2e3, 3e3, 434, 50, 200e3, 1, 140, 1400, True, 0.042966),
    ],
)
def test_psi_punching(
    l_x, l_y, f_yd, d_eff, e_s, approx_lvl_p, m_rd, m_ed, x_direction, expected
):
    """Test the psi_punching function."""
    assert math.isclose(
        _concrete_punching.psi_punching(
            l_x, l_y, f_yd, d_eff, e_s, approx_lvl_p, m_rd, m_ed, x_direction
        ),
        expected,
        rel_tol=0.001,
    )


@pytest.mark.parametrize(
    'psi_punching_val, d_eff, d_g, expected',
    [
        (0.015, 160, 32, 0.32051),  # Normal case
        (0.05, 160, 16, 0.11494),  # Smaller aggregate
        (0.001, 160, 32, 0.6),  # Small rotation
        (0.1, 160, 8, 0.04831),  # Very small aggregate
    ],
)
def test_k_psi(psi_punching_val, d_eff, d_g, expected):
    """Test the k_psi function."""
    assert math.isclose(
        _concrete_punching.k_psi(psi_punching_val, d_eff, d_g),
        expected,
        rel_tol=0.001,
    )


@pytest.mark.parametrize(
    'k_psi_val, b_0, d_v, f_ck, gamma_c, expected',
    [
        (0.6, 1000, 160, 30, 1.5, 350542.44),  # Normal case
        (0.3, 1000, 160, 30, 1.5, 175271.22),  # Lower k_psi
        (0.6, 1000, 160, 90, 1.5, 607157.31),  # High strength concrete
        (0.6, 1000, 160, 30, 2.0, 262906.83),  # Different safety factor
    ],
)
def test_v_rdc_punching(k_psi_val, b_0, d_v, f_ck, gamma_c, expected):
    """Test the v_rdc_punching function."""
    assert math.isclose(
        _concrete_punching.v_rdc_punching(k_psi_val, b_0, d_v, f_ck, gamma_c),
        expected,
        rel_tol=0.001,
    )


def test_v_rds_punching():
    """Test the v_rds_punching function."""
    # Normal case
    result = _concrete_punching.v_rds_punching(
        f_ywk=500,
        gamma_s=1.15,
        e_u=20,
        b_u=1000,
        e_s=200e3,
        psi_punching=0.015,
        alpha=90,
        f_bd=3,
        d_eff=160,
        phi_w=8,
        a_sw=100,
        v_ed=1000,
    )
    assert result > 0

    # Test warning for insufficient reinforcement
    warn_msg = 'Consider increasing punching shear reinforcement'
    with pytest.warns(UserWarning, match=warn_msg):
        _concrete_punching.v_rds_punching(
            f_ywk=500,
            gamma_s=1.15,
            e_u=20,
            b_u=1000,
            e_s=200e3,
            psi_punching=0.015,
            alpha=90,
            f_bd=3,
            d_eff=160,
            phi_w=8,
            a_sw=1,  # Very small area
            v_ed=1000,
        )


@pytest.mark.parametrize(
    ('d_v, f_ck, d_head, stirrups_compression, '
     'b0_val, k_psi_val, gamma_c, expected'),
    [
        # Normal cases
        (160, 30, True, True, 500, 0.6, 1.5, 292118.70),
        (160, 30, True, False, 500, 0.6, 1.5, 292118.70),
        (160, 30, False, True, 500, 0.6, 1.5, 292118.70),
        (160, 30, False, False, 500, 0.6, 1.5, 292118.70),
        # tests to trigger "k_sys * k_psi_val * base_resistance"
        (160, 30, False, False, 500, 0.1, 1.5, 58423),  # k = 2
        (160, 30, True, True, 500, 0.1, 1.5, 81793),  # k = 2.8
        (160, 30, False, True, 500, 0.1, 1.5, 70108),  # k = 2.4
        (160, 30, True, False, 500, 0.1, 1.5, 81793),  # k = 2.8
    ],
)
def test_v_rd_max_punching(
    d_v, f_ck, d_head, stirrups_compression,
    b0_val, k_psi_val, gamma_c, expected
):
    """Test the v_rd_max_punching function."""
    assert math.isclose(
        _concrete_punching.v_rd_max_punching(
            d_v, f_ck, d_head, stirrups_compression,
            b0_val, k_psi_val, gamma_c
        ),
        expected,
        rel_tol=0.001,
    )


@pytest.mark.parametrize(
    'v_rd_c, v_rd_s, v_rd_max, expected',
    [
        (100, 50, 200, 150),
        (150, 100, 200, 200),
        (100, 0, 150, 100),
        (0, 100, 150, 100),
        (120, 80, 200, 200),
    ],
)
def test_v_rd_punching(v_rd_c, v_rd_s, v_rd_max, expected):
    """Test the v_rd_punching function."""
    assert math.isclose(
        _concrete_punching.v_rd_punching(v_rd_c, v_rd_s, v_rd_max),
        expected,
        rel_tol=0.001,
    )
