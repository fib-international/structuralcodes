"""Tests for functions from Section 5 of EN 1992-1-1:2023."""

import pytest

from structuralcodes.codes.ec2_2023 import _section_8_3_torsion


@pytest.mark.parametrize(
    'TEd, Ak, teff_i, expected',
    [
        (100, 50000, 200, 5),
        (200, 100000, 400, 2.5),
        (50, 25000, 100, 10),
    ],
)
def test_tao_t_i(TEd, Ak, teff_i, expected):
    """Test calculate_torsional_shear_stress with parameterized inputs."""
    result = _section_8_3_torsion.tau_t_i(TEd, Ak, teff_i)
    assert result == pytest.approx(expected, rel=1e-3)


@pytest.mark.parametrize(
    'tau_t_i, teff_i, zi, expected',
    [
        (10, 200, 300, 600),
        (5, 100, 150, 75),
        (8, 250, 200, 400),
    ],
)
def test_VEd_i(tau_t_i, teff_i, zi, expected):
    """Test calculate_shear_force with parameterized inputs."""
    result = _section_8_3_torsion.VEd_i(tau_t_i, teff_i, zi)
    assert result == pytest.approx(expected, rel=1e-3)


# Test for tau_t_rd_sw
@pytest.mark.parametrize(
    'Asw, fywd, teff, s, cot_theta, cot_theta_min, expected',
    [
        (100, 500, 50, 200, 1.5, 1.0, 5.0),
        (150, 400, 60, 150, 2.0, 1.5, 10.0),
        (200, 450, 80, 250, 2.5, 2.0, 9.0),
        (250, 600, 100, 300, 3.0, 2.5, 12.5),
        (300, 550, 90, 240, 1.0, 0.8, 9.5486),
    ],
)
def test_tau_t_rd_sw(Asw, fywd, teff, s, cot_theta, cot_theta_min, expected):
    """Test test_tau_t_rd_sw."""
    result = _section_8_3_torsion.tau_t_rd_sw(
        Asw, fywd, teff, s, cot_theta, cot_theta_min
    )
    assert pytest.approx(result, rel=1e-2) == expected


@pytest.mark.parametrize(
    'Asw, fywd, teff, s, cot_theta, cot_theta_min, expected_exception',
    [
        (-100, 500, 50, 200, 1.5, 1.0, ValueError),
        (100, -500, 50, 200, 1.5, 1.0, ValueError),
        (100, 500, -50, 200, 1.5, 1.0, ValueError),
    ],
)
def test_tau_t_rd_sw_exceptions(
    Asw, fywd, teff, s, cot_theta, cot_theta_min, expected_exception
):
    """Test _tau_t_rd_sw exceptions."""
    with pytest.raises(expected_exception):
        _section_8_3_torsion.tau_t_rd_sw(
            Asw, fywd, teff, s, cot_theta, cot_theta_min
        )


# Test for tau_t_rd_sl
@pytest.mark.parametrize(
    'Asl, fyd, teff, uk, cot_theta, cot_theta_min, expected',
    [
        ([100, 200], [500, 600], 50, 300, 1.5, 1.0, 11.334),
        ([150, 250], [400, 550], 60, 400, 2.0, 1.5, 5.486),
        ([200, 300], [450, 500], 80, 500, 2.5, 2.0, 3.0),
        ([250, 350], [600, 650], 100, 600, 3.0, 2.5, 2.5166),
        ([300, 400], [550, 700], 90, 700, 1.0, 0.8, 5.650),
    ],
)
def test_tau_t_rd_sl(Asl, fyd, teff, uk, cot_theta, cot_theta_min, expected):
    """Test tau_t_rd_sl."""
    result = _section_8_3_torsion.tau_t_rd_sl(
        Asl, fyd, teff, uk, cot_theta, cot_theta_min
    )
    assert pytest.approx(result, rel=1e-2) == expected


@pytest.mark.parametrize(
    'Asl, fyd, teff, uk, cot_theta, cot_theta_min, expected_exception',
    [
        ([100], [500, 600], 50, 300, 1.5, 1.0, ValueError),
        ([100, -200], [500, 600], 50, 300, 1.5, 1.0, ValueError),
        ([100, 200], [500, -600], 50, 300, 1.5, 1.0, ValueError),
        ([100, 200], [500, 600], -50, 300, 1.5, 1.0, ValueError),
        ([100, 200], [500, 600], 50, -300, 1.5, 1.0, ValueError),
    ],
)
def test_tau_t_rd_sl_exceptions(
    Asl, fyd, teff, uk, cot_theta, cot_theta_min, expected_exception
):
    """Test tau_t_rd_sl excptions."""
    with pytest.raises(expected_exception):
        _section_8_3_torsion.tau_t_rd_sl(
            Asl, fyd, teff, uk, cot_theta, cot_theta_min
        )


# Test for tau_t_rd_max
@pytest.mark.parametrize(
    'nu, fcd, cot_theta, cot_theta_min, expected',
    [
        (0.6, 30, 1.5, 1.0, 9.0),
        (0.7, 35, 2.0, 1.5, 11.307),
        (0.8, 40, 2.5, 2.0, 12.8),
        (0.9, 45, 3.0, 2.5, 13.965),
        (0.5, 25, 1.0, 0.8, 6.0975),
    ],
)
def test_tau_t_rd_max(nu, fcd, cot_theta, cot_theta_min, expected):
    """Test tau_t_rd_max."""
    result = _section_8_3_torsion.tau_t_rd_max(
        nu, fcd, cot_theta, cot_theta_min
    )
    assert pytest.approx(result, rel=1e-2) == expected


@pytest.mark.parametrize(
    'nu, fcd, cot_theta, cot_theta_min, expected_exception',
    [
        (-0.6, 30, 1.5, 1.0, ValueError),
        (0.6, -30, 1.5, 1.0, ValueError),
    ],
)
def test_tau_t_rd_max_exceptions(
    nu, fcd, cot_theta, cot_theta_min, expected_exception
):
    """Test tau_t_rd_max exception."""
    with pytest.raises(expected_exception):
        _section_8_3_torsion.tau_t_rd_max(nu, fcd, cot_theta, cot_theta_min)


@pytest.mark.parametrize(
    'tau_t_rd_sw_val, tau_t_rd_sl_val, tau_t_rd_max_val, expected',
    [
        (10.0, 15.0, 20.0, 10.0),
        (25.0, 20.0, 15.0, 15.0),
        (30.0, 35.0, 40.0, 30.0),
        (45.0, 50.0, 55.0, 45.0),
    ],
)
def test_tau_t_rd(
    tau_t_rd_sw_val, tau_t_rd_sl_val, tau_t_rd_max_val, expected
):
    """Test for tau_t_rd."""
    result = _section_8_3_torsion.tau_t_rd(
        tau_t_rd_sw_val, tau_t_rd_sl_val, tau_t_rd_max_val
    )
    assert result == expected


@pytest.mark.parametrize(
    'tau_t_rd_sw_val, tau_t_rd_sl_val, tau_t_rd_max_val, expected_exception',
    [
        (-10.0, 15.0, 20.0, ValueError),
        (10.0, -15.0, 20.0, ValueError),
        (10.0, 15.0, -20.0, ValueError),
    ],
)
def test_tau_t_rd_exceptions(
    tau_t_rd_sw_val, tau_t_rd_sl_val, tau_t_rd_max_val, expected_exception
):
    """Test tau_t_rd exceptions."""
    with pytest.raises(expected_exception):
        _section_8_3_torsion.tau_t_rd(
            tau_t_rd_sw_val, tau_t_rd_sl_val, tau_t_rd_max_val
        )
