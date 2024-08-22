"""Test for functions from Section 13 of EN 1992-1-1:2023."""

import pytest

from structuralcodes.codes.ec2_2023 import _section13_precast


@pytest.mark.parametrize(
    'f_cmp, f_cm_tref, t, tref, tp, expected',
    [
        (40, 50, 7, 28, 1, 45.83),
        (30, 45, 14, 28, 2, 41.67),
    ],
)
def test_fcm_t_precast(f_cmp, f_cm_tref, t, tref, tp, expected):
    """Test fcm_t_precast."""
    assert _section13_precast.fcm_t_precast(
        f_cmp, f_cm_tref, t, tref, tp
    ) == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'temperatures, time_intervals, tmax, expected',
    [
        ([50, 55], [2, 3], 30, 61.16),
        ([60, 65, 70], [1, 1, 2], 50, 314.19),
    ],
)
def test_t_eq_precast(temperatures, time_intervals, tmax, expected):
    """Test t_eq_precast."""
    assert _section13_precast.t_eq_precast(
        temperatures, time_intervals, tmax
    ) == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'Ep, alpha_c_th, Tmax, T0, expected',
    [
        (200000, 1e-5, 60, 20, 40.0),
        (210000, 1.2e-5, 70, 30, 50.4),
    ],
)
def delta_signma_th_precast(Ep, alpha_c_th, Tmax, T0, expected):
    """Test specific thermal loss calculation."""
    assert _section13_precast.delta_signma_th_precast(
        Ep, alpha_c_th, Tmax, T0
    ) == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'csx, phi_p, Dupper, expected',
    [
        (30, 12, 10, 24),
        (40, 10, 20, 25),
    ],
)
def test_min_csx_precast(csx, phi_p, Dupper, expected):
    """Test min_csx_precast."""
    result = _section13_precast.min_csx_precast(csx, phi_p, Dupper)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'csy, phi_p, Dupper, expected',
    [
        (30, 12, 10, 24),
        (40, 10, 20, 25),
    ],
)
def test_min_csy_precast(csy, phi_p, Dupper, expected):
    """Test min_csy_precast."""
    result = _section13_precast.min_csy_precast(csy, phi_p, Dupper)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    's, phi_p, tendon_type, expected',
    [
        (24, 12, 'strand', 36.0),
        (30, 10, 'indented_wire', 40.0),
        (30, 12, 'strand', 30.0),
        (30, 10, 'indented_wire', 40.0),
    ],
)
def test_cmin_b_precast(s, phi_p, tendon_type, expected):
    """Test cmin_b_precast."""
    assert _section13_precast.cmin_b_precast(
        s, phi_p, tendon_type
    ) == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'gamma_c, release, tendon, sigma_pm0, position, fck_t, phi_p, expected',
    [
        (1.5, 'gradual', 'indented_wires', 120, 'favorable', 25, 12, 115.19),
        (1.5, 'sudden', 'stands', 130, 'unfavorable', 25, 15, 181.07),
    ],
)
def test_lpt_precast(
    gamma_c, release, tendon, sigma_pm0, position, fck_t, phi_p, expected
):
    """Test lpt_precast."""
    assert _section13_precast.lpt_precast(
        gamma_c, release, tendon, sigma_pm0, position, fck_t, phi_p
    ) == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'lpt, expected',
    [
        (15, 12.0),
        (20, 16.0),
    ],
)
def test_lpt1_precast(lpt, expected):
    """Test lpt1_precast."""
    assert _section13_precast.lpt1_precast(lpt) == pytest.approx(
        expected, rel=1e-2
    )


@pytest.mark.parametrize(
    'lpt, expected',
    [
        (15, 18.0),
        (20, 24.0),
    ],
)
def test_lpt2_precast(lpt, expected):
    """Test lpt2_precast."""
    assert _section13_precast.lpt2_precast(lpt) == pytest.approx(
        expected, rel=1e-2
    )


@pytest.mark.parametrize(
    'lpt, d, expected',
    [
        (15, 10, 18.03),
        (20, 15, 25.00),
    ],
)
def test_ldisp_precast(lpt, d, expected):
    """Test ldisp_precast."""
    assert _section13_precast.ldisp_precast(lpt, d) == pytest.approx(
        expected, rel=1e-2
    )


@pytest.mark.parametrize(
    'lpt2, gamma_c, tendon, fatigue, position, sigma_pd, '
    + 'sigma_pm_inf, fck, phi_p, expected',
    [
        (
            15,
            1.5,
            'indented_wires',
            True,
            'favorable',
            100,
            90,
            25,
            12,
            253.8,
        ),
        (
            20,
            1.5,
            'stands',
            False,
            'unfavorable',
            100,
            90,
            25,
            15,
            422.28,
        ),
    ],
)
def test_anchorage_length_lbpd(
    lpt2,
    gamma_c,
    tendon,
    fatigue,
    position,
    sigma_pd,
    sigma_pm_inf,
    fck,
    phi_p,
    expected,
):
    """Test total anchorage length lbpd calculation."""
    assert _section13_precast.lbpd_precast(
        lpt2,
        gamma_c,
        tendon,
        fatigue,
        position,
        sigma_pd,
        sigma_pm_inf,
        fck,
        phi_p,
    ) == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'sigma_xEd_y, tau_Ed_y, expected',
    [
        (10.0, 5.0, 12.071),
        (15.0, 7.0, 17.759),
    ],
)
def test_sigma_1_Ed_precast(sigma_xEd_y, tau_Ed_y, expected):
    """Test sigma_1_Ed_precast."""
    assert _section13_precast.sigma_1_Ed_precast(
        sigma_xEd_y, tau_Ed_y
    ) == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'V_Ed, S_y, In, b_y, expected',
    [
        (100, 5000, 3000000, 300, 0.5555),
        (200, 8000, 4000000, 400, 1.0),
    ],
)
def test_tau_Ed_y_precast(V_Ed, S_y, In, b_y, expected):
    """Test tau_Ed_y_precast."""
    assert _section13_precast.tau_Ed_y_precast(
        V_Ed, S_y, In, b_y
    ) == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'f_ctk_005, gamma_c, expected',
    [
        (3.0, 1.5, 2.0),
        (2.5, 1.4, 1.785),
    ],
)
def test_sigma_1_Rd_precast(f_ctk_005, gamma_c, expected):
    """Test sigma_1_Rd_precast."""
    assert _section13_precast.sigma_1_Rd_precast(
        f_ctk_005, gamma_c
    ) == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'qEd, be, expected',
    [
        (10, 2, 6.67),
        (15, 3, 15.0),
    ],
)
def test_vEd_joints_precast(qEd, be, expected):
    """Test vEd_joints_precast."""
    assert _section13_precast.vEd_joints_precast(qEd, be) == pytest.approx(
        expected, rel=1e-2
    )


@pytest.mark.parametrize(
    'sL, lL, h, load_type, expected',
    [
        (500, 4000, 200, 'residential_snow', float('inf')),
        (700, 5600, 250, 'other', 2500.0),
        (400, 3200, 150, 'other', 1500.0),
    ],
)
def test_sT_transverse_ribs(sL, lL, h, load_type, expected):
    """Test sT_transverse_ribs."""
    assert _section13_precast.sT_transverse_ribs(
        sL, lL, h, load_type
    ) == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'FEd, fyd, t, h, expected',
    [
        (100, 500, 20, 200, 5.0),
        (500, 500, 25, 300, 20.83),
        (0, 500, 25, 300, 0),
    ],
)
def test_calculate_reinforcement_area(FEd, fyd, t, h, expected):
    """Test calculate_reinforcement_area with a range of values."""
    assert _section13_precast.As_conn_precast(FEd, fyd, t, h) == pytest.approx(
        expected, rel=1e-2
    )


# Pytest for check_wall_vertical_load
@pytest.mark.parametrize(
    'FEd, h, fcd, reinforced, expected',
    [
        (100, 200, 30, False, 3000.0),
        (100, 200, 30, True, 3600.0),
        (70, 200, 30, False, 3000.0),
        (100, 100, 30, False, 1500.0),
    ],
)
def test_check_wall_vertical_load(FEd, h, fcd, reinforced, expected):
    """Test check_wall_vertical_load with a range of values."""
    assert (
        _section13_precast.FRd_conn_precast(FEd, h, fcd, reinforced)
        == expected
    )


@pytest.mark.parametrize(
    'hcol, MEd, NEd, expected',
    [
        (300, 50, 500, 383.78),
        (300, 2000, 500, 600),
        (300, 40, 500, 375.13),
    ],
)
def test_l_emb_min_precast(hcol, MEd, NEd, expected):
    """Test calculate_minimum_embedded_length with a range of values."""
    assert _section13_precast.l_emb_min_precast(
        hcol, MEd, NEd
    ) == pytest.approx(expected, rel=1e-2)
