"""Testing for functions from Section 8.4 of EN 1992-1-1:2023."""

import pytest

from structuralcodes.codes.ec2_2023 import _section_8_4_punching


@pytest.mark.parametrize(
    'dvx, dvy, expected',
    [
        (200.0, 200.0, 200.0),  # Case with equal values
        (250.0, 150.0, 200.0),  # Case with different values
        (0.0, 0.0, 0.0),  # Zero values
    ],
)
def test_shear_resisting_effective_depth(dvx, dvy, expected):
    """Test dv."""
    assert _section_8_4_punching.dv(dvx, dvy) == expected


def test_shear_resisting_effective_depth_negative():
    """Test dv errors."""
    with pytest.raises(ValueError):
        _section_8_4_punching.dv(-1, 100)

    with pytest.raises(ValueError):
        _section_8_4_punching.dv(100, -1)


@pytest.mark.parametrize(
    'support_type, ebx, eby, bb_min, bb_max, refined, expected',
    [
        ('internal columns', None, None, None, None, False, 1.15),
        ('edge columns', None, None, None, None, False, 1.4),
        ('corner columns', None, None, None, None, False, 1.5),
        ('ends of walls', None, None, None, None, False, 1.4),
        ('corners of walls', None, None, None, None, False, 1.2),
        ('internal columns', 100.0, 100.0, 400.0, 500.0, True, 1.3478),
        ('edge columns', 100.0, 100.0, 400.0, 500.0, True, 1.3689),
        ('corners of walls', 100.0, 100.0, 400.0, 500.0, True, 1.1328),
        ('internal columns', 10.0, 10.0, 400.0, 500.0, True, 1.05),
        ('edge columns', 10.0, 10.0, 400.0, 500.0, True, 1.05),
        ('corners of walls', 10.0, 10.0, 400.0, 500.0, True, 1.05),
    ],
)
def test_beta_e(support_type, ebx, eby, bb_min, bb_max, refined, expected):
    """Test the beta_e function with various inputs."""
    assert _section_8_4_punching.beta_e(
        support_type, ebx, eby, bb_min, bb_max, refined
    ) == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'support_type, ebx, eby, bb_min, bb_max, refined, error',
    [
        ('internal columns', 100.0, 100.0, -400.0, 500.0, True, ValueError),
        ('internal columns', 100.0, 100.0, 400.0, -500.0, True, ValueError),
        ('internal columns', 100.0, 100.0, 500.0, 400.0, True, ValueError),
        ('internal columns', None, 100.0, 400.0, 500.0, True, ValueError),
        ('internal columns', 100.0, None, 400.0, 500.0, True, ValueError),
    ],
)
def test_beta_e_invalid_inputs(
    support_type, ebx, eby, bb_min, bb_max, refined, error
):
    """Test that the beta_e function raises a ValueError for invalid inputs."""
    with pytest.raises(error):
        _section_8_4_punching.beta_e(
            support_type, ebx, eby, bb_min, bb_max, refined
        )


@pytest.mark.parametrize(
    'V_ed, beta_e, b_0_5, d_v, expected_tau_Ed',
    [
        (1000, 1.0, 500, 200, 10.0),
        (1500, 1.5, 600, 250, 15.0),
        (2000, 1.2, 700, 300, 11.42857),
        (0, 1.0, 500, 200, 0.0),
        (1000, 0.8, 400, 150, 13.33333),
    ],
)
def test_tau_Ed_punch(V_ed, beta_e, b_0_5, d_v, expected_tau_Ed):
    """Test the calculate_design_shear_stress function with various inputs."""
    assert (
        pytest.approx(
            _section_8_4_punching.tau_Ed_punch(V_ed, beta_e, b_0_5, d_v),
            rel=1e-3,
        )
        == expected_tau_Ed
    )


def test_calculate_design_shear_stress_invalid_inputs():
    """Test the calculate_design_shear_stress function for invalid inputs."""
    with pytest.raises(ValueError):
        _section_8_4_punching.tau_Ed_punch(1000, -1.0, 500, 200)
    with pytest.raises(ValueError):
        _section_8_4_punching.tau_Ed_punch(1000, 1.0, -500, 200)
    with pytest.raises(ValueError):
        _section_8_4_punching.tau_Ed_punch(1000, 1.0, 500, -200)


@pytest.mark.parametrize(
    'v_ed, d_v, expected_tau_Ed',
    [
        (200, 200, 1.0),
        (300, 250, 1.2),
        (150, 150, 1.0),
        (0, 200, 0.0),
        (250, 125, 2.0),
    ],
)
def test_tau_Ed_punch_d(v_ed, d_v, expected_tau_Ed):
    """Test the calculate_design_shear_stress_detailed
    function with various inputs.
    """
    assert (
        pytest.approx(
            _section_8_4_punching.tau_Ed_punch_d(v_ed, d_v), rel=1e-6
        )
        == expected_tau_Ed
    )


def test_tau_Ed_punch_d_invalid_inputs():
    """Test the calculate_design_shear_stress_detailed
    function for invalid inputs.
    """
    with pytest.raises(ValueError):
        _section_8_4_punching.tau_Ed_punch_d(200, -200)


@pytest.mark.parametrize(
    'gamma_v, rho_l_x, rho_l_y, f_ck, d_v, d_dg, b_0, b_0_5, expected',
    [
        (1.5, 0.01, 0.01, 30, 200, 16, 500, 550, 0.5355),
        (1.4, 0.015, 0.015, 40, 250, 20, 600, 650, 0.7229),
    ],
)
def test_tau_Rd_punch(
    gamma_v, rho_l_x, rho_l_y, f_ck, d_v, d_dg, b_0, b_0_5, expected
):
    """Test the calculation of punching shear
    resistance with different parameters.
    """
    result = _section_8_4_punching.tau_Rd_punch(
        gamma_v, rho_l_x, rho_l_y, f_ck, d_v, d_dg, b_0, b_0_5
    )
    assert result == pytest.approx(expected, rel=1e-3)


@pytest.mark.parametrize(
    'a_p_x, a_p_y, d_v, expected',
    [
        (1600, 1600, 200, 200),
        (1000, 1000, 300, 193.649),
        (2000, 1500, 250, 232.651),
    ],
)
def test_a_pd(a_p_x, a_p_y, d_v, expected):
    """Test the calculation of the modified effective
    depth with different parameters.
    """
    result = _section_8_4_punching.a_pd(a_p_x, a_p_y, d_v)
    assert result == pytest.approx(expected, rel=1e-3)


@pytest.mark.parametrize(
    'kpp_x, kpp_y, expected',
    [
        (2, 1, 1.4142),
        (3, 2, 2.449),
    ],
)
def test_kpp_xy(kpp_x, kpp_y, expected):
    """Test kpp_xy."""
    result = _section_8_4_punching.kpp_xy(kpp_x, kpp_y)
    assert result == pytest.approx(expected, rel=1e-3)


@pytest.mark.parametrize(
    'tau_Rd_c, tau_Ed, rho_w, f_ywd, d_v, phi_w, d_dg, k_pb, expected',
    [
        (0.5, 0.4, 0.01, 500, 200, 10, 20, 0.6, 5.0),
        (0.6, 0.5, 0.02, 450, 250, 12, 16, 0.7, 9.0),
    ],
)
def test_calculate_tau_Rd_cs(
    tau_Rd_c, tau_Ed, rho_w, f_ywd, d_v, phi_w, d_dg, k_pb, expected
):
    """Test the calculation of punching shear resistance."""
    result = _section_8_4_punching.tau_Rd_cs(
        tau_Rd_c, tau_Ed, rho_w, f_ywd, d_v, phi_w, d_dg, k_pb
    )
    assert result == pytest.approx(expected, rel=1e-3)


@pytest.mark.parametrize(
    'Asw, sr, st, expected',
    [
        (100, 150, 200, 0.00333),
        (150, 100, 250, 0.006),
    ],
)
def test_calculate_rho_w(Asw, sr, st, expected):
    """Test the calculation of vertical shear reinforcement ratio."""
    result = _section_8_4_punching.rho_w(Asw, sr, st)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'd_x, d_y, c_v, expected',
    [
        (300, 400, 50, 300),
        (500, 500, 40, 460),
    ],
)
def test_calculate_dv_out(d_x, d_y, c_v, expected):
    """Test the calculation of dv,out for shear reinforcement."""
    result = _section_8_4_punching.dv_out(d_x, d_y, c_v)
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'tau_Rd_c, b_0, d_v, reinforcement_type, expected',
    [
        (0.5, 300, 200, 'studs', 0.699),
        (0.6, 250, 150, 'links_and_stirrups', 0.729),
    ],
)
def test_calculate_tau_Rd_max(
    tau_Rd_c, b_0, d_v, reinforcement_type, expected
):
    """Test the calculation of maximum punching shear resistance."""
    result = _section_8_4_punching.tau_Rd_max_punch(
        tau_Rd_c, b_0, d_v, reinforcement_type
    )
    assert result == pytest.approx(expected, rel=1e-2)


@pytest.mark.parametrize(
    'b0_5, d_v, dv_out, eta_c, expected',
    [
        (2000, 250, 300, 1.2, 964.506),
        (1800, 300, 350, 1.1, 1092.933),
    ],
)
def test_calculate_b0_5_out(b0_5, d_v, dv_out, eta_c, expected):
    """Test the calculation of the outer control perimeter
    where shear reinforcement is not required.
    """
    result = _section_8_4_punching.b0_5_out(b0_5, d_v, dv_out, eta_c)
    assert result == pytest.approx(expected, rel=1e-2)
