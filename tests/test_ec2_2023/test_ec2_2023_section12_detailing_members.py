"""Test for functions from Section 12 of EN 1992-1-1:2023."""

import pytest

from structuralcodes.codes.ec2_2023 import _section12_detailing_members


@pytest.mark.parametrize(
    'A_c, f_ctm, f_yk, expected',
    [
        (1000, 2.5, 500, 5.0),
        (1500, 3.0, 600, 7.5),
        (1200, 2.0, 400, 6.0),
        (-1000, 2.5, 500, None),
    ],
)
def test_As_min(A_c, f_ctm, f_yk, expected):
    """Test As_min."""
    if A_c < 0 or f_ctm < 0 or f_yk < 0:
        with pytest.raises(ValueError):
            _section12_detailing_members.As_min(A_c, f_ctm, f_yk)
    else:
        assert (
            _section12_detailing_members.As_min(A_c, f_ctm, f_yk) == expected
        )


@pytest.mark.parametrize(
    'f_ck, f_yk, alpha, b_w, s, ductility_class, expected',
    [
        (30, 500, 45, 200, 100, 'C', 9.91),
        (25, 400, 30, 150, 80, 'B', 5.4),
        (40, 600, 60, 250, 120, 'A', 21.91),
    ],
)
def test_As_w_min(f_ck, f_yk, alpha, b_w, s, ductility_class, expected):
    """Test As_w_min."""
    assert (
        pytest.approx(
            _section12_detailing_members.As_w_min(
                f_ck, f_yk, alpha, b_w, s, ductility_class
            ),
            rel=1e-2,
        )
        == expected
    )


@pytest.mark.parametrize(
    'M_ed, ductility_class, expected',
    [
        (100, 'C', 100),
        (80, 'B', 88),
        (60, 'A', 78),
    ],
)
def test_MRd_min(M_ed, ductility_class, expected):
    """Test the minimum moment resistance."""
    assert (
        _section12_detailing_members.MRd_min(M_ed, ductility_class) == expected
    )


@pytest.mark.parametrize(
    'd, alpha, expected',
    [
        (500, 45, 750),
        (400, 30, 819.62),
    ],
)
def test_sl_max(d, alpha, expected):
    """Test maximum longitudinal spacing of shear assemblies."""
    assert _section12_detailing_members.sl_max(d, alpha) == pytest.approx(
        expected, rel=1e-2
    )


@pytest.mark.parametrize(
    'd, alpha, expected',
    [
        (500, 45, 600.0),
        (400, 30, 655.69),
    ],
)
def test_sbu_max(d, alpha, expected):
    """Test maximum longitudinal spacing of bent-up bars."""
    assert _section12_detailing_members.sbu_max(d, alpha) == pytest.approx(
        expected, rel=1e-2
    )


@pytest.mark.parametrize(
    'd, expected',
    [
        (500, 375.0),
        (400, 300.0),
    ],
)
def test_str_max(d, expected):
    """Test maximum transverse spacing of shear legs."""
    assert _section12_detailing_members.str_max(d) == pytest.approx(
        expected, rel=1e-2
    )


@pytest.mark.parametrize(
    'rho_w_req, expected',
    [
        (0.02, 0.01),
        (0.04, 0.02),
    ],
)
def test_rho_w_stir_min(rho_w_req, expected):
    """Test minimum ratio of shear reinforcement in the form of stirrups."""
    assert _section12_detailing_members.rho_w_stir_min(
        rho_w_req
    ) == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'rho_w_req, expected',
    [
        (0.02, 0.004),
        (0.04, 0.008),
    ],
)
def test_rho_w_tors_min(rho_w_req, expected):
    """Test minimum ratio of torsion reinforcement in
    the form of closed stirrups.
    """
    assert _section12_detailing_members.rho_w_tors_min(
        rho_w_req
    ) == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'u, b, h, expected',
    [
        (1600, 300, 500, 200.0),
        (2400, 250, 400, 250.0),
    ],
)
def test_sstir_max(u, b, h, expected):
    """Test maximum spacing for torsion assemblies/stirrups."""
    assert _section12_detailing_members.s_stir_max(u, b, h) == pytest.approx(
        expected, rel=1e-2
    )


def test_as_web_min():
    """Test minimum area of longitudinal surface
    reinforcement in beams with a downstand ≥ 600 mm.
    """
    assert _section12_detailing_members.sl_surf_max() == 300.0


@pytest.mark.parametrize(
    'z, theta, alpha, expected',
    [
        (500, 30, 45, 183.01),
        (400, 45, 60, 84.53),
    ],
)
def test_al_with_shear_reinforcement(z, theta, alpha, expected):
    """Test shift distance calculation for members with shear reinforcement."""
    assert _section12_detailing_members.al_with_shear_reinforcement(
        z, theta, alpha
    ) == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'd, expected',
    [
        (500, 500),
        (400, 400),
    ],
)
def test_al_without_shear_reinforcement(d, expected):
    """Test shift distance calculation for
    members without shear reinforcement.
    """
    assert _section12_detailing_members.al_without_shear_reinforcement(
        d
    ) == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'lbd, zone, expected',
    [
        (200, 'tension', 260.0),
        (150, 'compression', 105.0),
    ],
)
def test_anchorage_length_bent_up_bar(lbd, zone, expected):
    """Test anchorage length calculation for bent-up bars."""
    assert _section12_detailing_members.min_anchorage_length_bent_up_bar(
        lbd, zone
    ) == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'phi, expected',
    [
        (16, 160),
        (12, 120),
    ],
)
def test_bottom_reinforcement_extension(phi, expected):
    """Test minimum extension of bottom
    reinforcement at intermediate supports.
    """
    assert _section12_detailing_members.bottom_reinforcement_extension(
        phi
    ) == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'as_req_span, expected',
    [
        (1000, 200),
        (1200, 240),
    ],
)
def test_minimum_secondary_reinforcement(as_req_span, expected):
    """Test minimum secondary reinforcement calculation in slabs."""
    assert (
        _section12_detailing_members.As_slab_min_secondary_reinforcement(
            as_req_span
        )
        == expected
    )


@pytest.mark.parametrize(
    'as_req_span, expected',
    [
        (1000, 250),
        (1200, 300),
    ],
)
def test_minimum_bottom_reinforcement_inner_supports(as_req_span, expected):
    """Test minimum bottom reinforcement at inner supports in slabs."""
    assert (
        _section12_detailing_members.As_slab_min_bottom_reinforcement_inner_supports(
            as_req_span
        )
        == expected
    )


@pytest.mark.parametrize(
    'as_req_span, expected',
    [
        (1000, 250),
        (1200, 300),
    ],
)
def test_minimum_bottom_reinforcement_end_supports(as_req_span, expected):
    """Test minimum bottom reinforcement at end supports in slabs."""
    assert (
        _section12_detailing_members.As_slab_min_bottom_reinforcement_end_supports(
            as_req_span
        )
        == expected
    )


@pytest.mark.parametrize(
    'as_req_span, as_min, expected',
    [
        (1000, 150, 250),
        (800, 250, 250),
    ],
)
def test_minimum_top_reinforcement_end_supports(as_req_span, as_min, expected):
    """Test minimum top reinforcement at end supports in slabs."""
    assert (
        _section12_detailing_members.As_slab_min_top_reinforcement_end_supports(
            as_req_span, as_min
        )
        == expected
    )


@pytest.mark.parametrize(
    'h, expected',
    [
        (150, 400),
        (100, 300),
    ],
)
def test_sslab_max(h, expected):
    """Test maximum spacing of bars in slabs."""
    assert _section12_detailing_members.s_slab_max(h) == expected


@pytest.mark.parametrize(
    'd, alpha, expected',
    [
        (200, 45, 300),
        (180, 30, 368.82),
    ],
)
def test_sl_max_slab(d, alpha, expected):
    """Test maximum longitudinal spacing of shear assemblies in slabs."""
    assert _section12_detailing_members.sl_max_slab(d, alpha) == pytest.approx(
        expected, rel=1e-2
    )


@pytest.mark.parametrize(
    'd, expected',
    [
        (200, 200),
        (180, 180),
    ],
)
def test_sbu_max_slab(d, expected):
    """Test maximum longitudinal spacing of bent-up bars in slabs."""
    assert _section12_detailing_members.sbu_max_slab(d) == expected


@pytest.mark.parametrize(
    'd, expected',
    [
        (200, 300),
        (180, 270),
    ],
)
def test_str_max_slab(d, expected):
    """Test maximum transverse spacing of shear legs in slabs."""
    assert _section12_detailing_members.str_max_slab(d) == expected


@pytest.mark.parametrize(
    'tau_ed, fcd, alignment, expected',
    [
        (2.4, 30, 0, True),
        (4.8, 30, 45, True),
        (2.5, 30, 0, False),
        (5.0, 30, 30, False),
        (5.5, 30, 30, False),
    ],
)
def test_check_shear_stress_conditions(tau_ed, fcd, alignment, expected):
    """Test shear stress condition check with
    interpolation based on alignment angle.
    """
    assert (
        _section12_detailing_members.check_slab_shear_stress_conditions(
            tau_ed, fcd, alignment
        )
        == expected
    )


@pytest.mark.parametrize(
    'h, lbd, expected',
    [
        (40, 30, 80),
        (50, 30, 100),
        (60, 30, 120),
        (10, 30, 30),
        (25, 30, 50),
    ],
)
def test_min_slab_reinforcement_along_free_edge(h, lbd, expected):
    """Test min_slab_reinforcement_along_free_edge."""
    assert _section12_detailing_members.min_slab_reinforcement_along_free_edge(
        h, lbd
    ) == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'h, reinforcement_type, expected',
    [
        (220, 'stirrups', True),
        (180, 'bent_up_bars', True),
        (190, 'stirrups', False),
        (150, 'bent_up_bars', False),
    ],
)
def test_check_slab_depth_for_shear_reinforcement(
    h, reinforcement_type, expected
):
    """Test slab depth check for various shear reinforcement types."""
    assert (
        _section12_detailing_members.check_slab_depth_for_shear_reinforcement(
            h, reinforcement_type
        )
        == expected
    )


@pytest.mark.parametrize(
    'd, reinforcement_type, expected',
    [
        (200, 'single_leg', 10),
        (300, 'closed_stirrups', 13.47),
        (250, 'bent_up_bars', 17.88),
    ],
)
def test_max_shear_reinforcement_diameter(d, reinforcement_type, expected):
    """Test maximum effective diameter of shear reinforcement."""
    assert (
        _section12_detailing_members.slab_column_max_leg_shear_reinforcement_diameter(
            d, reinforcement_type
        )
        == pytest.approx(expected, rel=1e-2)
    )


@pytest.mark.parametrize(
    'dv, distance_to_column_edge, slab_type, expected',
    [
        (200, 300, 'column_edge', 300.0),
        (180, 360, 'flat_slab', 135.0),
        (200, 100, 'column_base', 100.0),
    ],
)
def test_tangential_spacing_limit(
    dv, distance_to_column_edge, slab_type, expected
):
    """Test tangential spacing limit based on
    the distance to the column edge.
    """
    assert _section12_detailing_members.slab_column_tangential_spacing_limit(
        dv, distance_to_column_edge, slab_type
    ) == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'As_int, fyd, kint, VEd, expected',
    [
        (4000, 500, 'B', 730, 740.0),
        (3000, 600, 'C', 900, 900),
    ],
)
def test_integrity_reinforcement_vrd_int(As_int, fyd, kint, VEd, expected):
    """Test integrity reinforcement resistance VRd,int."""
    assert _section12_detailing_members.Vrd_int_flat_slab(
        As_int, fyd, kint, VEd
    ) == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'rho_w, fywd, b0_5, dv, VEd, expected',
    [
        (0.002, 500, 3000, 200, 800, 600.0),
        (0.0015, 600, 3500, 180, 400, 400.0),
    ],
)
def test_integrity_reinforcement_vrd_w_int(
    rho_w, fywd, b0_5, dv, VEd, expected
):
    """Test integrity reinforcement resistance VRd,w,int."""
    assert _section12_detailing_members.Vrd_w_int_flat_slab(
        rho_w, fywd, b0_5, dv, VEd
    ) == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'nhog, fck, gamma_c, phi, s, c, expected',
    [
        (4, 30, 1.5, 16, 150, 30, 22.43),
        (6, 40, 1.5, 20, 180, 35, 60.72),
    ],
)
def test_hogging_reinforcement_vrd_hog(
    nhog, fck, gamma_c, phi, s, c, expected
):
    """Test hogging reinforcement resistance VRd,hog."""
    assert _section12_detailing_members.Vrd_hog_flat_slab(
        nhog, fck, gamma_c, phi, s, c
    ) == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'NEd, fyd, Ac, expected',
    [
        (1000, 500, 400000, 800.0),
        (2000, 600, 500000, 1000.0),
        (1500, 400, 300000, 600.0),
    ],
)
def test_As_column_min(NEd, fyd, Ac, expected):
    """Test minimum longitudinal reinforcement calculation."""
    assert _section12_detailing_members.As_column_min(
        NEd, fyd, Ac
    ) == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'h, b, expected',
    [
        (600, 300, 300.0),
        (500, 500, 400.0),
        (800, 450, 400.0),
    ],
)
def test_s_max_poly_col(h, b, expected):
    """Test maximum longitudinal spacing for polygonal cross-sections."""
    assert _section12_detailing_members.s_max_poly_col(h, b) == pytest.approx(
        expected, rel=1e-2
    )


@pytest.mark.parametrize(
    'n_bars, diameter, expected',
    [
        (8, 300, 117.8),
        (6, 400, 209.4),
        (10, 500, 157.1),
    ],
)
def test_s_max_circular_col(n_bars, diameter, expected):
    """Test s_max_circular_col."""
    assert _section12_detailing_members.s_max_circular_col(
        n_bars, diameter
    ) == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'phi_l_max, h, b, long_bars_res, expected',
    [
        (16, 600, 300, True, 300.0),
        (16, 600, 300, False, 300.0),
        (20, 500, 500, True, 400.0),
    ],
)
def test_s_max_col_int(phi_l_max, h, b, long_bars_res, expected):
    """Test s_max_col_int."""
    assert _section12_detailing_members.s_max_col_int(
        phi_l_max, h, b, long_bars_res
    ) == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    's_max_col, expected',
    [
        (300, 180.0),
        (250, 150.0),
    ],
)
def test_s_max_col_end(s_max_col, expected):
    """Test s_max_col_end."""
    assert _section12_detailing_members.s_max_col_end(
        s_max_col
    ) == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'Ac, fctm, fyk, design_case, expected',
    [
        (200000, 2.9, 500, 'in_plane_stress', 290.0),
        (200000, 2.9, 500, 'compression_bending', 200.0),
    ],
)
def test_As_wall_min_v(Ac, fctm, fyk, design_case, expected):
    """Test As_wall_min_v."""
    assert _section12_detailing_members.As_wall_min_v(
        Ac, fctm, fyk, design_case
    ) == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'Ac, As_v, fctm, fyk, design_case, expected',
    [
        (200000, 300, 2.9, 500, 'in_plane_stress', 290.0),
        (200000, 300, 2.9, 500, 'compression_bending', 75.0),
    ],
)
def test_As_wall_min_h(Ac, As_v, fctm, fyk, design_case, expected):
    """Test As_wall_min_h."""
    assert _section12_detailing_members.As_wall_min_h(
        Ac, As_v, fctm, fyk, design_case
    ) == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'h, expected',
    [
        (100, 300.0),
        (150, 400.0),
    ],
)
def test_s_max_wall_v(h, expected):
    """Test s_max_wall_v."""
    assert _section12_detailing_members.s_max_wall_v(h) == pytest.approx(
        expected, rel=1e-2
    )


def test_s_max_wall_h():
    """Test s_max_wall_h."""
    assert _section12_detailing_members.s_max_wall_h() == pytest.approx(
        400.0, rel=1e-2
    )


@pytest.mark.parametrize(
    'ch_i, delta_a, r_i, loop_type, expected',
    [
        (40, 10, 'horizontal_loops', 0, 50.0),
        (40, 10, 'vertical_bent', 20, 70.0),
    ],
)
def test_min_di_support_and_joint(ch_i, delta_a, r_i, loop_type, expected):
    """Test min_di_support_and_joint."""
    assert _section12_detailing_members.min_di_support_and_joint(
        ch_i, delta_a, r_i, loop_type
    ) == pytest.approx(expected, rel=1e-2)
