"""Test for functions from Section 14 of EN 1992-1-1:2023."""

import pytest

from structuralcodes.codes.ec2_2023 import _section14_plain_and_light


@pytest.mark.parametrize(
    'kc_pl, fcd, expected', [(0.8, 30, 24.0), (0.7, 25, 17.5), (1.0, 40, 40.0)]
)
def test_fcd_pl(kc_pl, fcd, expected):
    """Test fcd_pl function with various input values."""
    assert _section14_plain_and_light.fcd_pl(kc_pl, fcd) == pytest.approx(
        expected, rel=1e-2
    )


@pytest.mark.parametrize('kc_pl, fcd', [(-0.8, 30), (0.8, -30), (-0.8, -30)])
def test_fcd_pl_value_error(kc_pl, fcd):
    """Test fcd_pl function raises ValueError for negative inputs."""
    with pytest.raises(ValueError):
        _section14_plain_and_light.fcd_pl(kc_pl, fcd)


@pytest.mark.parametrize(
    'kt_pl, fctd, expected',
    [(0.8, 2.5, 2.0), (0.7, 3.0, 2.1), (1.0, 3.5, 3.5)],
)
def test_fctd_pl(kt_pl, fctd, expected):
    """Test fctd_pl function with various input values."""
    assert _section14_plain_and_light.fctd_pl(kt_pl, fctd) == pytest.approx(
        expected, rel=1e-2
    )


@pytest.mark.parametrize(
    'kt_pl, fctd', [(-0.8, 2.5), (0.8, -2.5), (-0.8, -2.5)]
)
def test_fctd_pl_value_error(kt_pl, fctd):
    """Test fctd_pl function raises ValueError for negative inputs."""
    with pytest.raises(ValueError):
        _section14_plain_and_light.fctd_pl(kt_pl, fctd)


@pytest.mark.parametrize(
    'fcd_pl, b, h, e, expected',
    [
        (24, 300, 600, 50, 3600),
        (20, 250, 500, 30, 2200),
        (30, 400, 800, 100, 7200),
    ],
)
def test_NRd_pl(fcd_pl, b, h, e, expected):
    """Test n_rd function with various input values."""
    assert _section14_plain_and_light.Nrd_pl(fcd_pl, b, h, e) == expected


@pytest.mark.parametrize(
    'fcd_pl, b, h, e',
    [
        (-24, 300, 600, 50),
        (24, -300, 600, 50),
        (24, 300, -600, 50),
        (24, 300, 600, -50),
        (24, 300, 600, 400),
    ],
)
def test_n_rd_value_error(fcd_pl, b, h, e):
    """Test n_rd function raises ValueError for invalid inputs."""
    with pytest.raises(ValueError):
        _section14_plain_and_light.Nrd_pl(fcd_pl, b, h, e)


@pytest.mark.parametrize(
    'NEd, Acc, expected',
    [
        (100, 10000, 10),
        (200, 20000, 10),
        (150, 15000, 10),
    ],
)
def test_sigma_cp(NEd, Acc, expected):
    """Test sigma_cp function with various input values."""
    result = _section14_plain_and_light.sigma_cp_pl(NEd, Acc)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'VEd, Acc, expected',
    [
        (100, 10000, 15),
        (200, 20000, 15),
        (150, 15000, 15),
    ],
)
def test_tau_cp(VEd, Acc, expected):
    """Test tau_cp function with various input values."""
    result = _section14_plain_and_light.tau_cp_pl(VEd, Acc)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'fcd_pl, fctd_pl, expected',
    [
        (30, 2, 14),
        (25, 2.5, 8.416),
        (35, 1.5, 20.201),
    ],
)
def test_sigma_c_lim(fcd_pl, fctd_pl, expected):
    """Test sigma_c_lim function with various input values."""
    result = _section14_plain_and_light.sigma_c_lim_pl(fcd_pl, fctd_pl)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'fctd_pl, sigma_cp, sigma_c_lim, expected',
    [
        (2, 0.01, 26.054, 2.0),
        (2.5, 0.02, 20.920, 2.51),
        (1.5, 0.015, 31.054, 1.51),
    ],
)
def test_tau_Rd_pl(fctd_pl, sigma_cp, sigma_c_lim, expected):
    """Test tau_Rd_pl function with various input values."""
    result = _section14_plain_and_light.tau_Rd_pl(
        fctd_pl, sigma_cp, sigma_c_lim
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'b, lw, num_sides, expected',
    [
        (1000, 3000, 3, 0.5),
        (2000, 3000, 4, 0.333),
        (4000, 4000, 4, 0.5),
    ],
)
def test_beta_Eul(b, lw, num_sides, expected):
    """Test beta_Eul function with various input values."""
    result = _section14_plain_and_light.beta_Eul(b, lw, num_sides)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'beta_Eul, lw, factor, expected',
    [
        (0.979, 3000, 1.0, 2937.0),
        (0.667, 3000, 1.0, 2001.0),
        (0.5, 4000, 0.85, 1700.0),
    ],
)
def test_l0(beta_Eul, lw, factor, expected):
    """Test l0 function with various input values."""
    result = _section14_plain_and_light.l0(beta_Eul, lw, factor)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'e_tot, l0, h, fcd_pl, phi_eff, expected',
    [
        (10, 3000, 500, 20, 100, 0.251),
        (15, 4000, 600, 25, 150, 0.115),
        (20, 5000, 700, 30, 200, 0.066),
    ],
)
def test_phi(e_tot, l0, h, fcd_pl, phi_eff, expected):
    """Test phi function with various input values."""
    result = _section14_plain_and_light.phi(e_tot, l0, h, fcd_pl, phi_eff)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'b, h, fcd_pl, phi, expected',
    [
        (300, 500, 20, 0.619, 1857.0),
        (400, 600, 25, 0.551, 3306.0),
        (500, 700, 30, 0.489, 5143.5),
    ],
)
def test_NRd_pl_simp(b, h, fcd_pl, phi, expected):
    """Test n_rd function with various input values."""
    result = _section14_plain_and_light.NRd_pl_simp(b, h, fcd_pl, phi)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'e0, ei, expected',
    [
        (10, 5, 15),
        (20, 10, 30),
        (15, 25, 40),
    ],
)
def test_e_tot(e0, ei, expected):
    """Test e_tot function with various input values."""
    result = _section14_plain_and_light.e_tot(e0, ei)
    assert result == pytest.approx(expected, rel=1e-2)


def test_min_tw():
    """Test min_tw."""
    assert _section14_plain_and_light.min_tw() == 120


@pytest.mark.parametrize(
    'sigma_gd, fctd_pl, expected',
    [
        (0.3, 2.5, 0.7058),
        (0.5, 3.0, 0.8318),
        (0.4, 2.0, 0.9112),
        (0.2, 2.5, 0.5763),
    ],
)
def test_min_footing_depth_ratio(sigma_gd, fctd_pl, expected):
    """Test min_footing_depth_ratio."""
    result = _section14_plain_and_light.min_footing_depth_ratio(
        sigma_gd, fctd_pl
    )
    assert result == pytest.approx(expected, rel=1e-2)


def test_min_footing_depth_ratio_simp():
    """Test min_footing_depth_ratio_simp."""
    assert _section14_plain_and_light.min_footing_depth_ratio_simp() == 2.0
