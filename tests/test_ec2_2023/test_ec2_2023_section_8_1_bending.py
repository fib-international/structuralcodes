"""EC2-2023 Section 8.1 Bending Test."""

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
            0.3**2.0 + 0.4**2.0,
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
