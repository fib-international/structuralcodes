"""EC2-2023 Section 8.1 Bending Test."""

from typing import List

import pytest

from structuralcodes.codes.ec2_2023 import _section_8_1_bending


@pytest.mark.parametrize(
    'h, NEd, expected',
    [
        (3000, 500, 50.0),
        (600, 1000, 20.0),
        (450, 250, 5.0),
        (1500, 800, 40.0),
        (2000, 300, 20.0),
    ],
)
def test_MEd_min(h, NEd, expected):
    """Test minimum beding due to imperfections."""
    assert _section_8_1_bending.MEd_min(h, NEd) == pytest.approx(
        expected, rel=1e-2
    )


@pytest.mark.parametrize(
    'h, expected',
    [
        (600, 0.02),
        (900, 0.03),
        (3000, 0.1),
        (30, 0.02),
        (1500, 0.05),
        (20, 0.02),
    ],
)
def test_e_min(h, expected):
    """Test minimum eccentricity for geometric imperfections."""
    assert _section_8_1_bending.e_min(h) == pytest.approx(expected, rel=1e-5)


@pytest.mark.parametrize(
    'Ac, fcd, As, fyd, expected',
    [
        (10000, 30, 1000, 500, 800),
        (20000, 25, 1500, 400, 1100),
        (15000, 20, 2000, 450, 1200),
        (25000, 40, 2500, 600, 2500),
        (5000, 50, 500, 250, 375),
    ],
)
def test_NRd0(Ac, fcd, As, fyd, expected):
    """Test resistance compressive strength."""
    assert _section_8_1_bending.NRd0(Ac, fcd, As, fyd) == pytest.approx(
        expected, rel=1e-2
    )


@pytest.mark.parametrize(
    'Ac, fcd, As, fyd',
    [
        (-10000, 30, 1000, 500),  # Negative concrete area
        (10000, -30, 1000, 500),  # Negative concrete resistance
        (10000, 30, -1000, 500),  # Negative steel area
        (10000, 30, 1000, -500),  # Negative steel resistance
    ],
)
def test_NRd0_errors(Ac, fcd, As, fyd):
    """Test resistance compressive strength with negative values."""
    with pytest.raises(ValueError):
        _section_8_1_bending.NRd0(Ac, fcd, As, fyd)


@pytest.mark.parametrize(
    'MEdz_MRdz, MEdy_MRdy, Ned_NRd, section_type, expected',
    [
        # Cases for elliptical and circular sections
        (0.5, 0.5, 0.5, 'circular', 0.5**2 + 0.5**2),  # an = 2.0
        (0.3, 0.4, 0.6, 'elliptical', 0.3**2 + 0.4**2),  # an = 2.0
        # Cases for rectangular sections with different Ned_NRd
        (
            0.2,
            0.3,
            0.1,
            'rectangular',
            0.2**1.0 + 0.3**1.0,
        ),  # an = 1.0 (Ned_NRd <= 0.1)
        (
            0.3,
            0.4,
            1.0,
            'rectangular',
            0.3706,
        ),  # an = 2.0 (Ned_NRd >= 1.0)
        (
            0.1,
            0.2,
            0.5,
            'rectangular',
            0.1633,
        ),  # an interpolated between 1.0 and 2.0
    ],
)
def test_biaxial_resistant_ratio(
    MEdz_MRdz, MEdy_MRdy, Ned_NRd, section_type, expected
):
    """Computes the biaxial resistant ratio."""
    result = _section_8_1_bending.biaxial_resistant_ratio(
        MEdz_MRdz, MEdy_MRdy, Ned_NRd, section_type
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'fcd, eps_c, expected',
    [
        (30.0, 0.0, 0.0),  # eps_c <= 0
        (30.0, 0.001, 22.5),  # eps_c <= 0.002
        (30.0, 0.002, 30.0),  # eps_c == 0.002
        (30.0, 0.003, 30.0),  # 0.002 < eps_c <= 0.0035
        (30.0, 0.0035, 30.0),  # eps_c == 0.0035
    ],
)
def test_sigma_cd(fcd, eps_c, expected):
    """Tests the concrete stress calculation."""
    assert _section_8_1_bending.sigma_cd(fcd, eps_c) == pytest.approx(
        expected, rel=1e-2
    )


@pytest.mark.parametrize(
    'fcd, eps_c',
    [
        (30.0, 0.004),  # eps_c > 0.0035
    ],
)
def test_sigma_cd_errors(fcd, eps_c):
    """Tests the compressive stress concrete calculation raises Errors."""
    with pytest.raises(ValueError):
        _section_8_1_bending.sigma_cd(fcd, eps_c)


@pytest.mark.parametrize(
    'sigma_c2d, f_cd, ddg, expected',
    [
        (10.0, 30.0, 32.0, 40.0),  # Test case with sigma_c2d <= 0.6 * f_cd
        (20.0, 30.0, 32.0, 77.467),  # Test case with sigma_c2d > 0.6 * f_cd
        (10.0, 30.0, 16.0, 20.0),  # Test case with ddg < 32
    ],
)
def test_delta_fcd_confined(
    sigma_c2d: float, f_cd: float, ddg: float, expected: float
):
    """Test delta_fcd_confined."""
    result = _section_8_1_bending.delta_fcd_confined(sigma_c2d, f_cd, ddg)
    assert pytest.approx(result, 0.01) == expected


@pytest.mark.parametrize(
    'A_s_conf, f_yd, b_cs, s, expected',
    [
        (500.0, 500.0, 200.0, 150.0, 16.67),
    ],
)
def test_confinement_sigma_c2d_circular_square(
    A_s_conf: float, f_yd: float, b_cs: float, s: float, expected: float
):
    """Test confinement_sigma_c2d_circular_square."""
    result = _section_8_1_bending.confinement_sigma_c2d_circular_square(
        A_s_conf, f_yd, b_cs, s
    )
    assert pytest.approx(result, 0.01) == expected


@pytest.mark.parametrize(
    'A_s_conf, f_yd, b_csx, b_csy, s, expected',
    [
        (500.0, 500.0, 250.0, 300.0, 150.0, 11.11),
    ],
)
def test_confinement_sigma_c2d_rectangular(
    A_s_conf: float,
    f_yd: float,
    b_csx: float,
    b_csy: float,
    s: float,
    expected: float,
):
    """Test test_confinement_sigma_c2d_rectangular."""
    result = _section_8_1_bending.confinement_sigma_c2d_rectangular(
        A_s_conf, f_yd, b_csx, b_csy, s
    )
    assert pytest.approx(result, 0.01) == expected


@pytest.mark.parametrize(
    'A_s_confx, A_s_confy, f_yd, b_csx, b_csy, s, expected',
    [
        ([100.0, 100.0], [100.0, 100.0], 500.0, 250.0, 300.0, 150.0, 2.22),
    ],
)
def test_confinement_sigma_c2d_multiple(
    A_s_confx: List[float],
    A_s_confy: List[float],
    f_yd: float,
    b_csx: float,
    b_csy: float,
    s: float,
    expected: float,
):
    """Test test_confinement_sigma_c2d_multiple."""
    result = _section_8_1_bending.confinement_sigma_c2d_multiple(
        A_s_confx, A_s_confy, f_yd, b_csx, b_csy, s
    )
    assert pytest.approx(result, 0.01) == expected


@pytest.mark.parametrize(
    'A_s_confx, A_s_confy, f_yd, b_csy, b_csx, s, expected',
    [
        ([100.0, 100.0], [100.0, 100.0], 500.0, 250.0, 300.0, 150.0, 2.22),
    ],
)
def test_confinement_sigma_c2d_compression_zones(
    A_s_confx: float,
    A_s_confy: float,
    f_yd: float,
    b_csy: float,
    b_csx: float,
    s: float,
    expected: float,
):
    """Test confinement_sigma_c2d_compression_zones."""
    result = _section_8_1_bending.confinement_sigma_c2d_compression_zones(
        A_s_confx, A_s_confy, f_yd, b_csy, b_csx, s
    )
    assert pytest.approx(result, 0.01) == expected


@pytest.mark.parametrize(
    'fcd, kconf_b, kconf_s, delta_fcd, expected',
    [
        (30, 0.5, 0.6, 5, 30 + 0.5 * 0.6 * 5),
        (40, 0.7, 0.8, 10, 40 + 0.7 * 0.8 * 10),
        (25, 1, 1, 15, 25 + 1 * 1 * 15),
    ],
)
def test_fcd_c(fcd, kconf_b, kconf_s, delta_fcd, expected):
    """Test fcd_c."""
    assert _section_8_1_bending.fcd_c(
        fcd, kconf_b, kconf_s, delta_fcd
    ) == pytest.approx(expected, rel=1e-5)


@pytest.mark.parametrize(
    'fcd, kconf_b, kconf_s, delta_fcd',
    [
        (-10, 0.5, 0.6, 5),
        (30, -0.5, 0.6, 5),
        (30, 0.5, -0.6, 5),
        (30, 0.5, 0.6, -5),
    ],
)
def test_fcd_c_raises_value_error(fcd, kconf_b, kconf_s, delta_fcd):
    """Test fcd_c_raises_value_error."""
    with pytest.raises(ValueError):
        _section_8_1_bending.fcd_c(fcd, kconf_b, kconf_s, delta_fcd)


@pytest.mark.parametrize(
    'bcs, b, expected',
    [
        (100, 200, (1 / 3) * (100 / 200) ** 2),
        (150, 300, (1 / 3) * (150 / 300) ** 2),
        (200, 400, (1 / 3) * (200 / 400) ** 2),
    ],
)
def test_kconf_b_square_single(bcs, b, expected):
    """Test kconf_b_square_single."""
    assert _section_8_1_bending.kconf_b_square_single(bcs, b) == pytest.approx(
        expected, rel=1e-5
    )


@pytest.mark.parametrize(
    'bcs, b',
    [
        (-100, 200),
        (100, -200),
    ],
)
def test_kconf_b_square_single_raises_value_error(bcs, b):
    """Test kconf_b_square_single_raises_value_error."""
    with pytest.raises(ValueError):
        _section_8_1_bending.kconf_b_square_single(bcs, b)


@pytest.mark.parametrize(
    's, bcs, expected',
    [
        (50, 200, max(1 - 50 / (2 * 200), 0) ** 2),
        (30, 150, max(1 - 30 / (2 * 150), 0) ** 2),
        (20, 100, max(1 - 20 / (2 * 100), 0) ** 2),
    ],
)
def test_kconf_s_square_single(s, bcs, expected):
    """Test kconf_s_square_single."""
    assert _section_8_1_bending.kconf_s_square_single(s, bcs) == pytest.approx(
        expected, rel=1e-5
    )


@pytest.mark.parametrize(
    's, bcs',
    [
        (-50, 200),
        (50, -200),
    ],
)
def test_kconf_s_square_single_raises_value_error(s, bcs):
    """Test kconf_s_square_single_raises_value_error."""
    with pytest.raises(ValueError):
        _section_8_1_bending.kconf_s_square_single(s, bcs)


@pytest.mark.parametrize(
    'bcs, b, expected',
    [
        (100, 200, (100 / 200) ** 2),
        (150, 300, (150 / 300) ** 2),
        (200, 400, (200 / 400) ** 2),
    ],
)
def test_kconf_b_circular(bcs, b, expected):
    """Test kconf_b_circular."""
    assert _section_8_1_bending.kconf_b_circular(bcs, b) == pytest.approx(
        expected, rel=1e-5
    )


@pytest.mark.parametrize(
    'bcs, b',
    [
        (-100, 200),
        (100, -200),
    ],
)
def test_kconf_b_circular_raises_value_error(bcs, b):
    """Test kconf_b_circular_raises_value_error."""
    with pytest.raises(ValueError):
        _section_8_1_bending.kconf_b_circular(bcs, b)


@pytest.mark.parametrize(
    'bcsx, bcsy, b_i, bx, by, expected',
    [
        (
            100,
            200,
            [50, 50],
            400,
            400,
            (100 * 200 - (1 / 6) * (50**2 + 50**2)) / (400 * 400),
        ),
        (
            150,
            250,
            [75, 75],
            500,
            500,
            (150 * 250 - (1 / 6) * (75**2 + 75**2)) / (500 * 500),
        ),
    ],
)
def test_kconf_b_multiple(bcsx, bcsy, b_i, bx, by, expected):
    """Test kconf_b_multiple."""
    assert _section_8_1_bending.kconf_b_multiple(
        bcsx, bcsy, b_i, bx, by
    ) == pytest.approx(expected, rel=1e-5)


@pytest.mark.parametrize(
    'bcsx, bcsy, b_i, bx, by',
    [
        (-100, 200, [50, 50], 400, 400),
        (100, -200, [50, 50], 400, 400),
        (100, 200, [-50, 50], 400, 400),
        (100, 200, [50, 50], -400, 400),
        (100, 200, [50, 50], 400, -400),
    ],
)
def test_kconf_b_multiple_raises_value_error(bcsx, bcsy, b_i, bx, by):
    """Test kconf_b_multiple_raises_value_error."""
    with pytest.raises(ValueError):
        _section_8_1_bending.kconf_b_multiple(bcsx, bcsy, b_i, bx, by)


@pytest.mark.parametrize(
    's, bcsx, bcsy, expected',
    [
        (
            50,
            200,
            200,
            max((1 - 50 / (2 * 200)), 0) * max((1 - 50 / (2 * 200)), 0),
        ),
        (
            30,
            150,
            150,
            max((1 - 30 / (2 * 150)), 0) * max((1 - 30 / (2 * 150)), 0),
        ),
    ],
)
def test_kconf_s_multiple(s, bcsx, bcsy, expected):
    """Test kconf_s_multiple."""
    assert _section_8_1_bending.kconf_s_multiple(
        s, bcsx, bcsy
    ) == pytest.approx(expected, rel=1e-5)


@pytest.mark.parametrize(
    's, bcsx, bcsy',
    [
        (-50, 200, 200),
        (50, -200, 200),
        (50, 200, -200),
    ],
)
def test_kconf_s_multiple_raises_value_error(s, bcsx, bcsy):
    """Test kconf_s_multiple_raises_value_error."""
    with pytest.raises(ValueError):
        _section_8_1_bending.kconf_s_multiple(s, bcsx, bcsy)


@pytest.mark.parametrize(
    'Ac_conf, Acc, b_i, expected',
    [
        (50000, 10000, [50, 50], (50000 - (1 / 6) * (50**2 + 50**2)) / 10000),
        (
            100000,
            20000,
            [75, 75],
            (100000 - (1 / 6) * (75**2 + 75**2)) / 20000,
        ),
    ],
)
def test_kconf_b_bending(Ac_conf, Acc, b_i, expected):
    """Test kconf_b_bending."""
    assert _section_8_1_bending.kconf_b_bending(
        Ac_conf, Acc, b_i
    ) == pytest.approx(expected, rel=1e-5)


@pytest.mark.parametrize(
    'Ac_conf, Acc, b_i',
    [
        (-50000, 10000, [50, 50]),
        (50000, -10000, [50, 50]),
        (50000, 10000, [-50, 50]),
    ],
)
def test_kconf_b_bending_raises_value_error(Ac_conf, Acc, b_i):
    """Test kconf_b_bending_raises_value_error."""
    with pytest.raises(ValueError):
        _section_8_1_bending.kconf_b_bending(Ac_conf, Acc, b_i)


@pytest.mark.parametrize(
    's, xcs, bcsx, bcsy, expected',
    [
        (
            50,
            100,
            200,
            200,
            max((1 - 50 / (4 * 100)), 0) * max((1 - 50 / (2 * 200)), 0),
        ),
        (
            30,
            75,
            150,
            150,
            max((1 - 30 / (4 * 75)), 0) * max((1 - 30 / (2 * 150)), 0),
        ),
    ],
)
def test_kconf_s_bending(s, xcs, bcsx, bcsy, expected):
    """Test kconf_s_bending."""
    assert _section_8_1_bending.kconf_s_bending(
        s, xcs, bcsx, bcsy
    ) == pytest.approx(expected, rel=1e-5)


@pytest.mark.parametrize(
    's, xcs, bcsx, bcsy',
    [
        (-50, 100, 200, 200),
        (50, -100, 200, 200),
        (50, 100, -200, 200),
        (50, 100, 200, -200),
    ],
)
def test_kconf_s_bending_raises_value_error(s, xcs, bcsx, bcsy):
    """Test kconf_s_bending_raises_value_error."""
    with pytest.raises(ValueError):
        _section_8_1_bending.kconf_s_bending(s, xcs, bcsx, bcsy)


@pytest.mark.parametrize(
    'epsc2, delta_fcd, fcd, expected',
    [
        (0.002, 5, 30, 0.002 * (1 + 5 * 5 / 30)),
        (0.003, 10, 40, 0.003 * (1 + 5 * 10 / 40)),
    ],
)
def test_epsc2_c(epsc2, delta_fcd, fcd, expected):
    """Test epsc2_c."""
    assert _section_8_1_bending.epsc2_c(
        epsc2, delta_fcd, fcd
    ) == pytest.approx(expected, rel=1e-5)


@pytest.mark.parametrize(
    'epsc2, delta_fcd, fcd',
    [
        (-0.002, 5, 30),
        (0.002, -5, 30),
        (0.002, 5, -30),
    ],
)
def test_epsc2_c_raises_value_error(epsc2, delta_fcd, fcd):
    """Test epsc2_c_raises_value_error."""
    with pytest.raises(ValueError):
        _section_8_1_bending.epsc2_c(epsc2, delta_fcd, fcd)


@pytest.mark.parametrize(
    'epscu, sigma_c2d, fcd, expected',
    [
        (0.003, 5, 30, 0.003 + 0.2 * 5 / 30),
        (0.004, 10, 40, 0.004 + 0.2 * 10 / 40),
    ],
)
def test_epscu_c(epscu, sigma_c2d, fcd, expected):
    """Test epscu_c."""
    assert _section_8_1_bending.epscu_c(
        epscu, sigma_c2d, fcd
    ) == pytest.approx(expected, rel=1e-5)


@pytest.mark.parametrize(
    'epscu, sigma_c2d, fcd',
    [
        (-0.003, 5, 30),
        (0.003, -5, 30),
        (0.003, 5, -30),
    ],
)
def test_epscu_c_raises_value_error(epscu, sigma_c2d, fcd):
    """Test epscu_c_raises_value_error."""
    with pytest.raises(ValueError):
        _section_8_1_bending.epscu_c(epscu, sigma_c2d, fcd)
