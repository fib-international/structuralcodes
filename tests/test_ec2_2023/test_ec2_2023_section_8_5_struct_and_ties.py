"""Test for functions from Section 8.5 of EN 1992-1-1:2023."""

import pytest

from structuralcodes.codes.ec2_2023 import _section_8_5_strut_and_ties


@pytest.mark.parametrize(
    'F_cd, b_c, t, expected',
    [
        (100, 200, 20, 25.0),
        (150, 350, 40, 10.714),
        (-150, 350, 40, 10.714),
    ],
)
def test_calculate_sigma_cd(F_cd, b_c, t, expected):
    """Test the calculation of compressive stress σ_cd."""
    assert _section_8_5_strut_and_ties.sigma_cd_strut(
        F_cd, b_c, t
    ) == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'F_cd, b_c, t',
    [
        (100, -200, 20),
        (200, 350, -40),
    ],
)
def test_calculate_sigma_cd_exceptions(F_cd, b_c, t):
    """Test exceptions in compressive stress σ_cd calculation."""
    with pytest.raises(ValueError):
        _section_8_5_strut_and_ties.sigma_cd_strut(F_cd, b_c, t)


@pytest.mark.parametrize(
    'input_value, expected',
    [
        (0, 0),
        (25, 0.471),
        (35, 0.642),
        (50, 0.791),
        (75, 0.888),
    ],
)
def test_calculate_nu(input_value, expected):
    """Test the calculation of strength reduction factor ν."""
    assert _section_8_5_strut_and_ties.nu_strut(input_value) == pytest.approx(
        expected, rel=1e-2
    )


@pytest.mark.parametrize(
    'input_value',
    [(-10), (180)],
)
def test_calculate_nu_exceptions(input_value):
    """Test exceptions in strength reduction factor ν calculation."""
    with pytest.raises(ValueError):
        _section_8_5_strut_and_ties.nu_strut(input_value)


@pytest.mark.parametrize(
    'input_value, expected',
    [
        (0.0, 1.0),
        (0.001, 0.9009),
        (0.01, 0.4762),
        (0.1, 0.0833),
    ],
)
def test_calculate_nu_refined(input_value, expected):
    """Test the calculation of refined strength reduction factor ν."""
    assert _section_8_5_strut_and_ties.nu_refined(
        input_value
    ) == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'input_value',
    [
        (-0.01),
    ],
)
def test_calculate_nu_refined_exceptions(input_value):
    """Test exceptions in refined strength reduction factor ν calculation."""
    with pytest.raises(ValueError):
        _section_8_5_strut_and_ties.nu_refined(input_value)


@pytest.mark.parametrize(
    'As, fyd, Ap, fpd, sigma_pd, expected',
    [
        (1000, 500, 500, 1600, 0, 1300),
        (1000, 500, 500, 1600, 100, 1250),
        (1500, 450, 300, 1450, 0, 1110),
        (1200, 550, 400, 1500, 50, 1240),
        (0, 500, 500, 1600, 0, 800),
        (1000, 0, 500, 1600, 0, 800),
        (1000, 500, 0, 1600, 0, 500),
        (1000, 500, 500, 0, 0, 500),
    ],
)
def test_calculate_tie_resistance(As, fyd, Ap, fpd, sigma_pd, expected):
    """Test the calculate_tie_resistance function with various inputs."""
    result = _section_8_5_strut_and_ties.FRd_tie(As, fyd, Ap, fpd, sigma_pd)
    assert pytest.approx(result, 0.01) == expected


@pytest.mark.parametrize(
    'Fd, a, b, H, near_edge, expected_F_td',
    [
        (1000, 200, 400, 600, False, 125),
        (1000, 200, 450, 600, False, 138.888),
        (1000, 200, 300, 600, True, 250),
        (800, 150, 300, 500, False, 100),
        (800, 150, 500, 300, True, 400),
        (800, 100, 200, 400, False, 100),
    ],
)
def test_calculate_transverse_reinforcement(
    Fd, a, b, H, near_edge, expected_F_td
):
    """Test the calculate_transverse_reinforcement
    function with various inputs.
    """
    F_td = _section_8_5_strut_and_ties.Ftd_conc(Fd, a, b, H, near_edge)
    assert pytest.approx(F_td, 0.01) == expected_F_td
