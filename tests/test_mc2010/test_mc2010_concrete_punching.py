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
    'v_ed, e_u, b_s, inner, edge_par, edge_per, corner, expected',
    [
        (10e3, 20, 660, True, False, False, False, 1401),  # Inner
        (10e3, 20, 660, False, True, False, False, 2500),  # Edge par
        (10e3, 20, 660, False, False, True, False, 1553.03),  # Edge perp
        (10e3, 20, 660, False, False, False, True, 5000),  # Corner
        # Edge cases
        (10e3, 0, 660, True, False, False, False, 1250),  # No ecc
        (10e3, 100, 660, True, False, False, False, 2007.58),  # Large ecc
    ],
)
def test_m_ed(v_ed, e_u, b_s, inner, edge_par, edge_per, corner, expected):
    """Test the m_ed function."""
    assert math.isclose(
        _concrete_punching.m_ed(
            v_ed, e_u, b_s, inner, edge_par, edge_per, corner
        ),
        expected,
        rel_tol=0.001,
    )


def test_m_ed_invalid():
    """Test m_ed with invalid inputs."""
    with pytest.raises(ValueError, match='Placement is not defined'):
        _concrete_punching.m_ed(10e3, 20, 660, False, False, False, False)


@pytest.mark.parametrize(
    (
        'psi_punching_level_one, psi_punching_level_two, '
        'psi_punching_level_three, approx_lvl_p, expected'
    ),
    [
        (0.013427, 0.47089, 0.47089, 1, 0.013427),  # Level 1
        (0.013427, 0.47089, 0.47089, 2, 0.47089),  # Level 2
        (0.013427, 0.47089, 0.37671, 3, 0.37671),  # Level 3
    ],
)
def test_psi_punching(
    psi_punching_level_one,
    psi_punching_level_two,
    psi_punching_level_three,
    approx_lvl_p,
    expected,
):
    """Test the psi_punching function."""
    assert math.isclose(
        _concrete_punching.psi_punching(
            psi_punching_level_one,
            psi_punching_level_two,
            psi_punching_level_three,
            approx_lvl_p,
        ),
        expected,
        rel_tol=0.001,
    )


def test_psi_punching_invalid():
    """Test psi_punching with invalid inputs."""
    with pytest.raises(ValueError, match='Approximation level is not defined'):
        _concrete_punching.psi_punching(0.01, 0.02, 0.03, 4)


@pytest.mark.parametrize(
    'k_dg_val, d_eff, psi_punching_val, expected',
    [
        (0.75, 160, 0.015, 0.3205),  # Normal case
        (0.75, 160, 0.05, 0.145),  # Larger rotation
        (0.75, 160, 0.001, 0.6),  # Small rotation
        (1.0, 160, 0.015, 0.273),  # Different k_dg
    ],
)
def test_k_psi(k_dg_val, d_eff, psi_punching_val, expected):
    """Test the k_psi function."""
    assert math.isclose(
        _concrete_punching.k_psi(k_dg_val, d_eff, psi_punching_val),
        expected,
        rel_tol=0.001,
    )


@pytest.mark.parametrize(
    'd_g, expected',
    [
        (32, 0.75),  # Large aggregate, limited by 0.75
        (16, 1.0),  # Medium aggregate
        (8, 1.333),  # Small aggregate
        (0, 2.0),  # No aggregate
    ],
)
def test_k_dg(d_g, expected):
    """Test the k_dg function."""
    assert math.isclose(
        _concrete_punching.k_dg(d_g),
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
        f_ywd=434.78,  # 500/1.15
        e_u=20,
        b_u=1000,
        alpha=90,
        sigma_swd=400,
        a_sw=100,
        v_ed=1000,
    )
    assert math.isclose(result, 39215.686, rel_tol=0.001)

    # Test warning for insufficient reinforcement
    warn_msg = 'Consider increasing punching shear reinforcement'
    with pytest.warns(UserWarning, match=warn_msg):
        _concrete_punching.v_rds_punching(
            f_ywd=434.78,  # 500/1.15
            e_u=20,
            b_u=1000,
            alpha=90,
            sigma_swd=400,
            a_sw=1,  # Very small area
            v_ed=1000,
        )


@pytest.mark.parametrize(
    (
        'd_v, f_ck, d_head, stirrups_compression, '
        'b0_val, k_psi_val, gamma_c, expected'
    ),
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
    d_v,
    f_ck,
    d_head,
    stirrups_compression,
    b0_val,
    k_psi_val,
    gamma_c,
    expected,
):
    """Test the v_rd_max_punching function."""
    assert math.isclose(
        _concrete_punching.v_rd_max_punching(
            d_v, f_ck, d_head, stirrups_compression, b0_val, k_psi_val, gamma_c
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


@pytest.mark.parametrize(
    'l_x, l_y, expected',
    [
        (2000, 2000, 660),
        (2000, 3000, 808.332),
        (1000, 1000, 330.0),
        (5000, 500, 500),
    ],
)
def test_b_s(l_x, l_y, expected):
    """Test the b_s function."""
    assert math.isclose(
        _concrete_punching.b_s(l_x=l_x, l_y=l_y),
        expected,
        rel_tol=0.001,
    )


@pytest.mark.parametrize(
    'l_x, l_y, expected',
    [
        (2000, 3000, 2000),
        (5000, 500, 500),
    ],
)
def test_b_sr(l_x, l_y, expected):
    """Test the b_sr function."""
    assert math.isclose(
        _concrete_punching.b_sr(l_x=l_x, l_y=l_y),
        expected,
        rel_tol=0.001,
    )


@pytest.mark.parametrize(
    'l_x, l_y, f_yd, d_eff, e_s, expected',
    [
        (2000, 3000, 434, 160, 200e3, 0.013427),  # Normal case
        (1000, 1000, 434, 160, 200e3, 0.004476),  # Equal spans
        (2000, 3000, 434, 50, 200e3, 0.042966),  # Small effective depth
        (2000, 3000, 500, 160, 200e3, 0.015476),  # Higher steel strength
    ],
)
def test_psi_punching_level_one(l_x, l_y, f_yd, d_eff, e_s, expected):
    """Test the psi_punching_level_one function."""
    assert math.isclose(
        _concrete_punching.psi_punching_level_one(l_x, l_y, f_yd, d_eff, e_s),
        expected,
        rel_tol=0.001,
    )


@pytest.mark.parametrize(
    'r_s, f_yd, d_eff, e_s, m_ed, m_rd, m_Pd, expected',
    [
        (440, 434, 160, 200e3, 1500, 140, 0, 0.3139),  # Normal case
        (440, 434, 160, 200e3, 1400, 140, 0, 0.2831),  # Lower m_ed
        (440, 434, 160, 200e3, 1500, 200, 0, 0.1839),  # Higher m_rd
        (440, 434, 160, 200e3, 1500, 140, 50, 0.5789),  # With prestressing
    ],
)
def test_psi_punching_level_two(
    r_s, f_yd, d_eff, e_s, m_ed, m_rd, m_Pd, expected
):
    """Test the psi_punching_level_two function."""
    assert math.isclose(
        _concrete_punching.psi_punching_level_two(
            r_s, f_yd, d_eff, e_s, m_ed, m_rd, m_Pd
        ),
        expected,
        rel_tol=0.001,
    )


@pytest.mark.parametrize(
    (
        'psi_punching_level_two, is_uncracked_model, '
        'is_moment_from_uncracked_model, expected'
    ),
    [
        (0.47089, False, False, 0.47089),  # No change
        (0.47089, True, True, 0.37671),  # Both conditions true
        (0.47089, True, False, 0.47089),  # Only one condition true
        (0.47089, False, True, 0.47089),  # Only one condition true
    ],
)
def test_psi_punching_level_three(
    psi_punching_level_two,
    is_uncracked_model,
    is_moment_from_uncracked_model,
    expected,
):
    """Test the psi_punching_level_three function."""
    assert math.isclose(
        _concrete_punching.psi_punching_level_three(
            psi_punching_level_two,
            is_uncracked_model,
            is_moment_from_uncracked_model,
        ),
        expected,
        rel_tol=0.001,
    )


@pytest.mark.parametrize(
    (
        'l_x, l_y, x_direction, is_level_three_approximation, '
        'column_edge_or_corner, b_sr_val, expected'
    ),
    [
        (2000, 3000, True, False, False, -1, 440),  # x direction
        (2000, 3000, False, False, False, -1, 660),  # y direction
        (2000, 3000, True, True, True, 2000, 1340),  # Level 3, edge/corner,
        # x direction
        (2000, 3000, False, True, True, 2000, 1340),  # Level 3, edge/corner,
        # y direction
    ],
)
def test_r_s(
    l_x,
    l_y,
    x_direction,
    is_level_three_approximation,
    column_edge_or_corner,
    b_sr_val,
    expected,
):
    """Test the r_s function."""
    assert math.isclose(
        _concrete_punching.r_s(
            l_x,
            l_y,
            x_direction,
            is_level_three_approximation,
            column_edge_or_corner,
            b_sr_val,
        ),
        expected,
        rel_tol=0.001,
    )


def test_r_s_error():
    """Test r_s with invalid inputs."""
    with pytest.raises(
        ValueError, match='b_sr is not defined for Level 3 of Approximation'
    ):
        _concrete_punching.r_s(2000, 3000, True, True, True, None)


@pytest.mark.parametrize(
    'f_ywk, gamma_s, expected',
    [
        (500, 1.15, 434.78),  # Normal case
        (600, 1.15, 521.74),  # Higher steel strength
        (500, 1.0, 500),  # No safety factor
    ],
)
def test_f_ywd(f_ywk, gamma_s, expected):
    """Test the f_ywd function."""
    assert math.isclose(
        _concrete_punching.f_ywd(f_ywk, gamma_s),
        expected,
        rel_tol=0.001,
    )


@pytest.mark.parametrize(
    'e_s, psi_punching, alpha, f_bd, d_eff, f_ywd, phi_w, expected',
    [
        (200e3, 0.015, 90, 3, 160, 434.78, 8, 434.78),  # Limited by f_ywd
        (200e3, 0.001, 90, 3, 160, 434.78, 8, 37.9334),  # Small rotation
        (200e3, 0.015, 45, 3, 160, 434.78, 8, 434.78),  # Different angle
        (200e3, 0.015, 90, 3, 160, 300, 8, 300),  # Lower f_ywd
    ],
)
def test_sigma_swd(
    e_s, psi_punching, alpha, f_bd, d_eff, f_ywd, phi_w, expected
):
    """Test the sigma_swd function."""
    assert math.isclose(
        _concrete_punching.sigma_swd(
            e_s, psi_punching, alpha, f_bd, d_eff, f_ywd, phi_w
        ),
        expected,
        rel_tol=0.001,
    )
