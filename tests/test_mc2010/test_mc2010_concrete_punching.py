"""Test for the function of _concrete_punching."""

import math

import pytest

from structuralcodes.codes.mc2010._concrete_punching import (
    ConcretePunching,
)


@pytest.mark.parametrize(
    'Ved, e_u, l_x, l_y, inner, edge_par, edge_per, corner, expected',
    [
        (10e3, 20, 2e3, 2e3, True, False, False, False, 1401),
        (10e3, 20, 2e3, 3e3, False, True, False, False, 2500),
        (10e3, 20, 2e3, 3e3, False, False, True, False, 1497),
        (10e3, 20, 2e3, 3e3, False, False, False, True, 5000),
    ],
)
def test_m_ed(Ved, e_u, l_x, l_y, inner, edge_par, edge_per, corner, expected):
    """Test the m_ed function."""
    # Create an instance with some default
    # values that don't affect m_ed calculation
    punching = ConcretePunching(
        l_x=l_x,
        l_y=l_y,
        f_yd=-1,  # does not effect m_ed calculation
        d=-1,  # does not effect m_ed calculation
        e_s=-1,  # does not effect m_ed calculation
        v_ed=Ved,
        e_u=e_u,
        inner=inner,
        edge_par=edge_par,
        edge_per=edge_per,
        corner=corner,
        m_rd=-1,  # does not effect m_ed calculation
        m_pd=-1,  # does not effect m_ed calculation
        x_direction=-1,  # does not effect m_ed calculation
    )
    assert math.isclose(
        punching.m_ed(),
        expected,
        rel_tol=0.001,
    )


@pytest.mark.parametrize(
    (
        'l_x, l_y, f_yd, d, e_s, approx_lvl_p, Ved, e_u, inner, edge_par, '
        'edge_per, corner, m_rd, x_direction, expected'
    ),
    [
        (
            2e3,
            3e3,
            434,
            160,
            200e3,
            1,
            50e3,
            20,
            True,
            False,
            False,
            False,
            140,
            True,
            0.013426875,
        ),
        (
            2e3,
            3e3,
            434,
            160,
            200e3,
            2,
            10e3,
            20,
            True,
            False,
            False,
            False,
            140,
            False,
            0.41269218238,
        ),
    ],
)
def test_psi_punching(
    l_x,
    l_y,
    f_yd,
    d,
    e_s,
    approx_lvl_p,
    Ved,
    e_u,
    inner,
    edge_par,
    edge_per,
    corner,
    m_rd,
    x_direction,
    expected,
):
    """Test the psi_punching function."""
    punching = ConcretePunching(
        l_x=l_x,
        l_y=l_y,
        f_yd=f_yd,
        d=d,
        e_s=e_s,
        v_ed=Ved,
        e_u=e_u,
        inner=inner,
        edge_par=edge_par,
        edge_per=edge_per,
        corner=corner,
        m_rd=m_rd,
        m_pd=-1,  # does not effect psi_punching calculation
        x_direction=x_direction,
    )
    assert math.isclose(
        punching.psi_punching(approx_lvl_p),
        expected,
        rel_tol=0.001,
    )
