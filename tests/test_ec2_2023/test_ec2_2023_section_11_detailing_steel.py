"""Test for functions from Section 11 of EN 1992-1-1:2023."""

import pytest

from structuralcodes.codes.ec2_2023 import _section11_detailing_steel


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
    result = _section11_detailing_steel.min_clear_distance_between_bars(
        phi, D_upper
    )
    assert result == expected


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
    result = _section11_detailing_steel.min_clear_distance_to_poured_concrete(
        surface_roughness, min_bond_distance
    )
    assert result == expected


@pytest.mark.parametrize(
    'phi, expected',
    [
        (10, 40),
        (20, 140),
    ],
)
def test_mandrel_diameter(phi, expected):
    """Test minimum mandrel diameter for bending bars."""
    assert _section11_detailing_steel.min_mandrel_phi(phi) == expected


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
    result = _section11_detailing_steel.sigma_sd_bend_concrete(
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
    result = _section11_detailing_steel.k_trans(
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
    result = _section11_detailing_steel.k_bend(alpha_bend)
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
    result = _section11_detailing_steel.lbd_simple(fck, phi, bond_condition)
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
    result = _section11_detailing_steel.lbd(
        phi, fck, sigma_sd, cd, c, bond_condition, design_situation
    )
    assert pytest.approx(result, rel=1e-2) == expected


@pytest.mark.parametrize(
    'cx, cy, phi_t, st, cs, phi, rho_conf, sigma_ccd, f_ck, expected',
    [
        (40, 40, 15, 200, 50, 16, 0.01, 5, 30, 96),
        (50, 45, 15, 200, 55, 20, 0.02, 7, 35, 120),
    ],
)
def test_cd_conf(
    cx, cy, phi_t, st, cs, phi, rho_conf, sigma_ccd, f_ck, expected
):
    """Test calculate_cd_conf function with example values."""
    result = _section11_detailing_steel.cd_conf(
        cx, cy, phi_t, st, cs, phi, rho_conf, sigma_ccd, f_ck
    )
    assert pytest.approx(result, rel=1e-2) == expected


@pytest.mark.parametrize(
    'nc, phi_c, nb, phi, sc, expected',
    [
        (4, 10, 2, 16, 150, 0.0654),
        (3, 12, 1, 20, 100, 0.1696),
    ],
)
def test_rho_conf(nc, phi_c, nb, phi, sc, expected):
    """Test calculate_rho_conf function with example values."""
    result = _section11_detailing_steel.rho_conf(nc, phi_c, nb, phi, sc)
    assert pytest.approx(result, rel=1e-3) == expected


@pytest.mark.parametrize(
    'phi, n1, n2, expected_Ast, expected_Asc',
    [
        (16, 2, 3, 102.4, 153.6),
        (20, 1, 2, 80.0, 160.0),
    ],
)
def test_calculate_additional_reinforcement(
    phi, n1, n2, expected_Ast, expected_Asc
):
    """Test calculate_additional_reinforcement function with example values."""
    Ast, Asc = (
        _section11_detailing_steel.additional_transverse_confinement_reinf(
            phi, n1, n2
        )
    )
    assert pytest.approx(Ast, rel=1e-3) == expected_Ast
    assert pytest.approx(Asc, rel=1e-3) == expected_Asc


@pytest.mark.parametrize(
    'area_bars, expected',
    [
        (804, 32.0),
        (1256.64, 40.0),
    ],
)
def test_phi_b_anchor(area_bars, expected):
    """Test calculate_equivalent_diameter function with example values."""
    result = _section11_detailing_steel.phi_b_anchor(area_bars)
    assert pytest.approx(result, rel=1e-3) == expected


@pytest.mark.parametrize(
    'cs, cx, cy, cyb, expected',
    [
        (20, 25, 30, 15, 10),
        (30, 40, 20, 25, 15),
        (50, 35, 45, 20, 20),
    ],
)
def test_cd(cs, cx, cy, cyb, expected):
    """Test nominal cover cd calculation."""
    assert _section11_detailing_steel.cd(cs, cx, cy, cyb) == expected


@pytest.mark.parametrize(
    'lbd, phi, expected',
    [
        (200, 10, 100),
        (180, 12, 120),
        (220, 8, 100),
    ],
)
def test_lb_anchor_tension(lbd, phi, expected):
    """Test design anchorage length for tension."""
    assert _section11_detailing_steel.lb_anchor_tension(lbd, phi) == expected


@pytest.mark.parametrize(
    'lbd, phi, perpendicular_distance, expected',
    [
        (200, 10, 40, 100),
        (180, 12, 36, 180),
        (220, 8, 28, 100),
    ],
)
def test_lb_anchor_compression(lbd, phi, perpendicular_distance, expected):
    """Test design anchorage length for compression."""
    assert (
        _section11_detailing_steel.lb_anchor_compression(
            lbd, phi, perpendicular_distance
        )
        == expected
    )


@pytest.mark.parametrize(
    'lbd, phi, phi_t, s, number_of_transverse_bars, expected',
    [
        (200, 10, 6, 75, 2, 50),
        (180, 12, 8, 80, 1, 60),
        (250, 16, 10, 60, 2, 80),
    ],
)
def test_anchorage_length_welded_transverse_reinforcement(
    lbd, phi, phi_t, s, number_of_transverse_bars, expected
):
    """Test design anchorage length calculation
    for welded transverse reinforcement.
    """
    assert (
        _section11_detailing_steel.lb_anchor_weld_trans(
            lbd, phi, phi_t, s, number_of_transverse_bars
        )
        == expected
    )


@pytest.mark.parametrize(
    'lbd, phi, expected',
    [
        (250, 12, 120),
        (300, 16, 160),
        (200, 10, 100),
    ],
)
def test_anchorage_length_u_bar_loops(lbd, phi, expected):
    """Test design anchorage length calculation for U-bar loops."""
    assert (
        _section11_detailing_steel.lb_anchor_ubar_loops_tension(lbd, phi)
        == expected
    )


@pytest.mark.parametrize(
    'Ah, expected',
    [
        (2800, 59.708),
        (4400, 75),
        (6300, 90),
        (0, 0),
    ],
)
def test_phi_h(Ah, expected):
    """Test calculation of equivalent diameter for circular head."""
    result = _section11_detailing_steel.phi_h(Ah)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'fck, phi, phi_h, th, ay, ax, sx, Ah, concrete_cracked, expected',
    [
        (30, 20, 60, 15, 60, 130, 240, 2827.43, False, False),
        (28, 25, 75, 18, 75, 160, 320, 4417.86, True, False),
        (25, 22, 55, 20, 80, 100, 180, 2375.08, True, False),
        (40, 20, 60, 15, 50, 120, 200, 2827.43, False, False),
        (35, 30, 90, 20, 90, 190, 360, 6361.73, True, False),
    ],
)
def test_check_anchorage_headed_bars_in_tension(
    fck, phi, phi_h, th, ay, ax, sx, Ah, concrete_cracked, expected
):
    """Test if the headed bar configuration
    satisfies anchorage requirements.
    """
    assert (
        _section11_detailing_steel.check_anchorage_headed_bars_in_tension(
            fck, phi, phi_h, th, ay, ax, sx, Ah, concrete_cracked
        )
        == expected
    )


@pytest.mark.parametrize(
    'phi_h, phi, expected',
    [
        (60, 20, 8.0),
        (75, 25, 8.0),
        (80, 20, 15.0),
    ],
)
def test_kh_A(phi_h, phi, expected):
    """Test calculation of kh_A."""
    result = _section11_detailing_steel.kh_A(phi_h, phi)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'ay, ax, sx, phi, phi_h, expected',
    [
        (60, 130, 240, 20, 60, 44.5),
        (75, 100, 200, 25, 75, 40),
        (80, 100, 180, 22, 55, 48.5),
        (50, 70, 150, 18, 50, 27.5),
    ],
)
def test_a_d(ay, ax, sx, phi, phi_h, expected):
    """Test calculation of ad."""
    result = _section11_detailing_steel.a_d(ay, ax, sx, phi, phi_h)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'fck, gamma_c, phi, phi_h, a_d, kh_A, nu_part, ddg, expected',
    [
        (30, 1.5, 20, 60, 60, 8.0, 11.0, 32, 512.066),
        (28, 1.5, 25, 75, 75, 8.25, 8.0, 32, 383.6342),
        (25, 1.5, 22, 55, 58.75, 15.0, 8.0, 32, 423.146),
    ],
)
def test_sigma_sd_prime(
    fck, gamma_c, phi, phi_h, a_d, kh_A, nu_part, ddg, expected
):
    """Test calculation of maximum tensile stress sigma_sd_prime."""
    result = _section11_detailing_steel.sigma_sd_prime(
        fck, gamma_c, phi, phi_h, a_d, kh_A, nu_part, ddg
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'lbd_sigma_sd, lbd_sigma_sd_prime, expected',
    [
        (200, 150, 55),
        (300, 100, 220),
        (180, 180, 0),
    ],
)
def test_calculate_lbd(lbd_sigma_sd, lbd_sigma_sd_prime, expected):
    """Test calculation of design length lbd."""
    result = _section11_detailing_steel.lbd_bar_tension(
        lbd_sigma_sd, lbd_sigma_sd_prime
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'phi, lbd_pi, drilling_method, use_drilling_aid, expected',
    [
        (20, 200, 'rotary_percussion', False, 42),
        (26, 300, 'diamond', True, 46),
        (24, 150, 'compressed_air', False, 62),
        (26, 250, 'compressed_air', True, 65),
        (20, 200, 'rotary_percussion', True, 34),
    ],
)
def test_cmin_b(phi, lbd_pi, drilling_method, use_drilling_aid, expected):
    """Test calculation of minimum concrete cover cmin,b."""
    result = _section11_detailing_steel.cmin_b(
        phi, lbd_pi, drilling_method, use_drilling_aid
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'lbd, kb_pi, phi, alpha_lb, expected',
    [
        (200, 1.0, 20, 1.5, 300),
        (300, 1.5, 25, 1.5, 375),
        (250, 1.2, 22, 1.5, 330),
    ],
)
def test_lbd_pi(lbd, kb_pi, phi, alpha_lb, expected):
    """Test calculation of design anchorage length lbd,pi."""
    assert (
        _section11_detailing_steel.lbd_pi(lbd, kb_pi, phi, alpha_lb)
        == expected
    )


@pytest.mark.parametrize(
    'type_of_lap, state, lbd, phi, kls, lbd_pi, alpha_lb, phi_mand, expected',
    [
        ('straight_bars', 'tension', 200, 10, 1.2, 0, 1.0, 0, 240),
        ('bends_and_hooks', 'tension', 250, 15, 1.2, 0, 1.0, 0, 300),
        ('loops', 'tension', 0, 12, 1.0, 0, 1.0, 50, 98),
        ('headed_bars', 'compression', 300, 0, 1.0, 0, 1.0, 0, 300),
        ('intermeshed_fabric', 'tension', 220, 20, 1.2, 0, 1.0, 0, 300),
        ('layered_fabric', 'compression', 220, 20, 1.2, 0, 1.0, 0, 300),
        ('bonded_postinstalled', 'tension', 0, 15, 1.0, 280, 1.5, 0, 337.5),
    ],
)
def test_lsd(
    type_of_lap, state, lbd, phi, kls, lbd_pi, alpha_lb, phi_mand, expected
):
    """Test calculation of design lap length
    lsd for various types of lap splices.
    """
    result = _section11_detailing_steel.lsd(
        type_of_lap, state, lbd, phi, kls, lbd_pi, alpha_lb, phi_mand
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'phi, expected',
    [
        (10, 20),
        (15, 30),
        (8, 20),
    ],
)
def test_cs_laps(phi, expected):
    """Test calculation of minimum clear distance between adjacent laps."""
    assert _section11_detailing_steel.cs_laps(phi) == expected


@pytest.mark.parametrize(
    'phi, expected',
    [
        (10, 40),
        (15, 50),
        (8, 32),
    ],
)
def test_cs_bars(phi, expected):
    """Test calculation of clear distance between lapping bars."""
    assert _section11_detailing_steel.cs_bars(phi) == expected


@pytest.mark.parametrize(
    'phi, expected',
    [
        (10, 40),
        (15, 50),
        (8, 32),
    ],
)
def test_calculate_min_clear_distance(phi, expected):
    """Test calculation of minimum clear distance between lapped bars."""
    assert _section11_detailing_steel.min_clear_distante_laps(phi) == expected


@pytest.mark.parametrize(
    'phi_mand, phi, lsd, expected',
    [
        (50, 20, 200, 12971.0),
        (40, 15, 150, 7614.75),
    ],
)
def test_Ac_u_bars(phi_mand, phi, lsd, expected):
    """Test calculation of total effective concrete area Ac."""
    result = _section11_detailing_steel.Ac_u_bars(phi_mand, phi, lsd)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'Ast, fyd, Ac, lsd, fcd, ddg, expected',
    [
        (100, 500, 5000, 200, 30, 16, 1),
        (120, 400, 15000, 150, 25, 16, 0.867),
    ],
)
def test_kst(Ast, fyd, Ac, lsd, fcd, ddg, expected):
    """Test calculation of resistance factor kst."""
    result = _section11_detailing_steel.kst_u_bar_loops(
        Ast, fyd, Ac, fcd, lsd, ddg
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'fcd, lsd, ddg, kst, cs, Ac, expected',
    [
        (30, 200, 16, 0.8, 50, 10000, 17.546),
        (25, 150, 16, 0.9, 40, 12500, 21.305),
    ],
)
def test_TRd_c(fcd, lsd, ddg, kst, cs, Ac, expected):
    """Test calculation of resistance of a single lap splice TRd_c."""
    result = _section11_detailing_steel.TRd_c_u_bar_loops(
        fcd, lsd, ddg, kst, cs, Ac
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'fck, Ac, fyk, expected',
    [
        (30, 7640, 500, 41.846),
        (25, 4359, 400, 27.243),
    ],
)
def test_Ast_min_u_bar_loops(fck, Ac, fyk, expected):
    """Test calculation of minimum amount of confinement reinforcement."""
    result = _section11_detailing_steel.Ast_min_u_bar_loops(fck, Ac, fyk)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'lsd, phi, bh1, expected',
    [
        (200, 20, 20, 3200),
        (150, 15, 15, 1800),
    ],
)
def test_Ac_headed_laps(lsd, phi, bh1, expected):
    """Test calculation of effective concrete area Ac for headed bars."""
    result = _section11_detailing_steel.Ac_headed_laps(lsd, phi, bh1)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'phi_h, expected',
    [
        (50, 44.4),
        (40, 35.52),
    ],
)
def test_bh1_headed_laps(phi_h, expected):
    """Test calculation of effective width bh1 for circular heads."""
    result = _section11_detailing_steel.bh1_headed_laps(phi_h)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'Ast, fyd, lsd, fcd, phi, ddg, Ac, expected',
    [
        (100, 500, 200, 30, 20, 16, 10000, 0.8),
        (120, 400, 150, 25, 15, 16, 12500, 0.711),
    ],
)
def test_kst_headed_bars(Ast, fyd, lsd, fcd, phi, ddg, Ac, expected):
    """Test calculation of resistance factor kst
    for transverse reinforcement.
    """
    result = _section11_detailing_steel.kst_headed_bars(
        Ast, fyd, fcd, Ac, lsd, phi, ddg
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'fcd, lsd, ddg, kst, cs, phi, bh1, Ac, expected',
    [
        (30, 200, 16, 0.8, 50, 20, 50, 10000, 53.048),
        (25, 150, 16, 0.9, 40, 15, 40, 12500, 64.389),
    ],
)
def test_TRd_c_headed(fcd, lsd, ddg, kst, cs, phi, bh1, Ac, expected):
    """Test calculation of resistance of a single headed bar lap TRd_c."""
    result = _section11_detailing_steel.TRd_c_headed(
        fcd, lsd, ddg, kst, cs, phi, bh1, Ac
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'fck, Ac, fyk, phi, expected',
    [
        (30, 3600, 500, 20, 157.07),
        (25, 1575, 400, 15, 88.357),
    ],
)
def test_Ast_min_headed_bars(fck, Ac, fyk, phi, expected):
    """Test calculation of minimum amount of transverse reinforcement."""
    result = _section11_detailing_steel.Ast_min_headed_bars(fck, Ac, fyk, phi)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'phi, expected',
    [
        (20, 48.0),
        (15, 27.0),
    ],
)
def test_Ast_min_headed_single_bars(phi, expected):
    """Test calculation of minimum area of tie down reinforcement."""
    result = _section11_detailing_steel.Ast_min_headed_single_bars(phi)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'Dupper, max_bar_diameter, expected',
    [
        (30, 25, 35),
        (40, 45, 45),
        (50, 55, 55),
    ],
)
def test_min_clear_distante_laps_couplers(Dupper, max_bar_diameter, expected):
    """Test calculation of minimum clear distance between
    couplers and between couplers and adjacent bars.
    """
    assert (
        _section11_detailing_steel.min_clear_distante_laps_couplers(
            Dupper, max_bar_diameter
        )
        == expected
    )


@pytest.mark.parametrize(
    'Dupper, phi_duct, expected',
    [
        (30, 25, 50),
        (40, 45, 50),
        (60, 55, 65),
    ],
)
def test_csx_duct(Dupper, phi_duct, expected):
    """Test calculation of minimum horizontal
    spacing csx for post-tensioning tendons.
    """
    assert _section11_detailing_steel.csx_duct(Dupper, phi_duct) == expected


@pytest.mark.parametrize(
    'Dupper, phi_duct, expected',
    [
        (30, 25, 40),
        (40, 45, 45),
        (60, 55, 60),
    ],
)
def test_csy_duct(Dupper, phi_duct, expected):
    """Test calculation of minimum vertical spacing
    csy for post-tensioning tendons.
    """
    assert _section11_detailing_steel.csy_duct(Dupper, phi_duct) == expected


def test_min_spacing_bundled_tendons():
    """Test calculation of minimum spacing between bundled tendons."""
    assert _section11_detailing_steel.min_spacing_bundled_tendons() == 100.0


@pytest.mark.parametrize(
    'sigma_pd, pRd, Ap, expected',
    [
        (1000, 15, 200, 942.81),
        (1200, 70, 150, 209.96),
        (1100, 30, 180, 491.93),
    ],
)
def test_Rmin_unconf_concrete(sigma_pd, pRd, Ap, expected):
    """Test calculation of minimum radius of curvature Rmin for tendons."""
    result = _section11_detailing_steel.Rmin_unconf_concrete(sigma_pd, pRd, Ap)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'case, fcd, expected',
    [
        ('internal_corrugated', 20, 15),
        ('internal_U_shape', 30, 70),
        ('external_smooth', 25, 30),
    ],
)
def test_max_pRd(case, fcd, expected):
    """Test getting the maximum transverse bearing
    stress pRd based on the case.
    """
    assert _section11_detailing_steel.max_pRd(case, fcd) == expected


@pytest.mark.parametrize(
    'cs, cy, phi, expected',
    [
        (100, 10, 12, 55.43),
        (30, 15, 20, 30),
    ],
)
def test_cu(cs, cy, phi, expected):
    """Test calculation of parameter cu."""
    result = _section11_detailing_steel.cu(cs, cy, phi)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'Ftd, r, cu, fck, gamma_c, expected',
    [
        (100, 2e6, 50, 30, 1.5, True),
        (200, 1e4, 30, 25, 1.5, False),
    ],
)
def test_verify_deviation_forces(Ftd, r, cu, fck, gamma_c, expected):
    """Test verification if the deviation forces can
    be resisted by the concrete in tension.
    """
    assert (
        _section11_detailing_steel.verify_deviation_forces(
            Ftd, r, cu, fck, gamma_c
        )
        == expected
    )


@pytest.mark.parametrize(
    'Ftd, r, cu, fck, gamma_c, lsd, ls, expected',
    [
        (100, 2e6, 50, 30, 1.5, 100, 120, True),
        (200, 1e4, 30, 25, 1.5, 80, 100, False),
    ],
)
def test_verify_interaction_bond(Ftd, r, cu, fck, gamma_c, lsd, ls, expected):
    """Test verification of the interaction between bond
    and transverse tensile stresses due to deviation forces.
    """
    assert (
        _section11_detailing_steel.verify_interaction_bond(
            Ftd, r, cu, fck, gamma_c, lsd, ls
        )
        == expected
    )
