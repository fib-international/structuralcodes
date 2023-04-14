"""Tests for the _concrete_creep_and_shrinkage module"""
import numpy as np

import pytest

from structuralcodes.codes.mc2010 import _concrete_creep_and_shrinkage


@pytest.mark.parametrize(
    'fcm, expected',
    [
        (10, None), # Raises ValueError
        (20, None), # Raises ValueError
        (21, None),
        (129, None),
        (130, None), # Raises ValueError
        (150, None) # Raises ValueError
    ],
)
def test_check_fcm(fcm, expected):
    """Test check_fcm function."""
    try:
        assert _concrete_creep_and_shrinkage._check_fcm(fcm) == expected
    except ValueError:
        with pytest.raises(ValueError) as exc_info:
            _concrete_creep_and_shrinkage._check_fcm(fcm) == expected
            assert str(exc_info.value).startswith(
                    "The specified mean compressive strength")


@pytest.mark.parametrize(
    'sigma, fcm, expected',
    [
        (5, 20, None),
        (10, 20, None),
        (15, 20, None) # Raises ValueError
    ],
)
def test_check_initial_stress(sigma, fcm, expected):
    """Test check_initial_stress function."""
    try:
        assert _concrete_creep_and_shrinkage._check_initial_stress(sigma, fcm) == expected
    except ValueError:
        with pytest.raises(ValueError) as exc_info:
            _concrete_creep_and_shrinkage._check_initial_stress(sigma, fcm) == expected
            assert str(exc_info.value).startswith(
                    "The stress level exceeds")


@pytest.mark.parametrize(
    't0, expected',
    [
        (0.5, None), # raises ValueError
        (1, None)
    ],
)
def test_check_age_at_loading(t0, expected):
    """Test check_age_at_loading function."""
    try:
        assert _concrete_creep_and_shrinkage._check_age_at_loading(t0) == expected
    except ValueError:
        with pytest.raises(ValueError) as exc_info:
            _concrete_creep_and_shrinkage._check_age_at_loading(t0) == expected
            assert str(exc_info.value).startswith(
                    "The load is applied too soon")


@pytest.mark.parametrize(
    'rh, expected',
    [
        (0.2, None), # raises ValueError
        (0.5, None),
        (0.99, None),
        (1.2, None), # raises ValueError
        (20, None), # raises ValueError
        (40, None),
        (99, None),
    ],
)
def test_check_RH(rh, expected):
    """Test check_RH function."""
    try:
        assert _concrete_creep_and_shrinkage._check_RH(rh) == expected
    except ValueError:
        with pytest.raises(ValueError) as exc_info:
            _concrete_creep_and_shrinkage._check_RH(rh) == expected
            assert str(exc_info.value).startswith(
                    "The specified relative humidity")


@pytest.mark.parametrize(
    'T, expected',
    [
        (1, None),
        (6, None),
        (35, None),
    ],
)
def test_check_env_temp(T, expected):
    """Test check_env_temp function."""
    assert _concrete_creep_and_shrinkage._check_env_temp(T) == expected


@pytest.mark.parametrize(
    'cem_class, expected',
    [
        ("32.5 N", None),
        ("32.5 R", None),
        ("42.5 N", None),
        ("42.5 R", None),
        ("52.5 N", None),
        ("52.5 R", None),
        ("", None), # raises ValueError
    ],
)
def test_check_cem_strength_class(cem_class, expected):
    """Test check_cem_strength_class function."""
    try:
        assert _concrete_creep_and_shrinkage._check_cem_strength_class(cem_class) == expected
    except ValueError:
        with pytest.raises(ValueError) as exc_info:
            _concrete_creep_and_shrinkage._check_cem_strength_class(cem_class) == expected
            assert str(exc_info.value).startswith(
                    "Unknown cem_class used.")


@pytest.mark.parametrize(
    'cem_class, expected',
    [
        ("32.5 R", {"alpha": 0, "alpha_as": 700, "alpha_ds1": 4, "alpha_ds2": 0.012}),
        ("42.5 N", {"alpha": 0, "alpha_as": 700, "alpha_ds1": 4, "alpha_ds2": 0.012}),
        ("42.5 R", {"alpha": 1, "alpha_as": 600, "alpha_ds1": 6, "alpha_ds2": 0.012}),
        ("52.5 N", {"alpha": 1, "alpha_as": 600, "alpha_ds1": 6, "alpha_ds2": 0.012}),
        ("52.5 R", {"alpha": 1, "alpha_as": 600, "alpha_ds1": 6, "alpha_ds2": 0.012}),
        ("32.5 N", {"alpha": -1, "alpha_as": 800, "alpha_ds1": 3, "alpha_ds2": 0.013})
    ],
)
def test_get_creep_shrinkage_coeffs(cem_class, expected):
    """Test get_creep_shrinkage_coeffs function."""
    assert _concrete_creep_and_shrinkage._get_creep_shrinkage_coeffs(cem_class) == expected


@pytest.mark.parametrize(
    't0, T_cur, expected',
    [
        (1, 20, 0.99812),
        (7, 20, 6.98687),
        (14, 20, 13.97374),
        (28, 20, 27.94749),
    ],
)
def test_calc_temp_corr_age(t0, T_cur, expected):
    """Test _calc_temp_corr_age function"""
    assert np.isclose(_concrete_creep_and_shrinkage._calc_temp_corr_age(t0, T_cur), expected, rtol=1e-5)


@pytest.mark.parametrize(
    't0, T_cur, TABULAR_VALUES, expected',
    [
        (7, 20, {"alpha": 0, "S": 0.25, "alpha_as": 700,
            "alpha_ds1": 4, "alpha_ds2": 0.012, "alpha_e": 1.2},
            6.98687
        ),
        (7, 20, {"alpha": -1, "S": 0.25, "alpha_as": 700,
            "alpha_ds1": 4, "alpha_ds2": 0.012, "alpha_e": 1.2},
            4.03566
        ),
        (7, 20, {"alpha": 1, "S": 0.25, "alpha_as": 700,
            "alpha_ds1": 4, "alpha_ds2": 0.012, "alpha_e": 1.2},
            12.09623
        ),
        (1, 20, {"alpha": -1, "S": 0.25, "alpha_as": 700,
            "alpha_ds1": 4, "alpha_ds2": 0.012, "alpha_e": 1.2},
            0.5
        )
    ],
)
def test_calc_modified_t0(t0, T_cur, TABULAR_VALUES, expected):
    """Test _calc_modified_tp function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage._calc_modified_t0(
            t0, T_cur, TABULAR_VALUES),
        expected, rtol=1e-5)


@pytest.mark.parametrize(
    'fcm, TABULAR_VALUES, expected',
    [
        (28, {"alpha": 0, "S": 0.25, "alpha_as": 700,
            "alpha_ds1": 4, "alpha_ds2": 0.012, "alpha_e": 1.2},
            4.71651e-4
        ),
        (60, {"alpha": 0, "S": 0.25, "alpha_as": 700,
            "alpha_ds1": 4, "alpha_ds2": 0.012, "alpha_e": 1.2},
            3.21256e-4
        )
    ],
)
def test_calc_notional_drying_shrinkage(fcm, TABULAR_VALUES, expected):
    """Test _calc_notional_drying_shrinkage function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage._calc_notional_drying_shrinkage(
            fcm, TABULAR_VALUES),
        expected, rtol=1e-5, atol=1e-10
    )


@pytest.mark.parametrize(
    'time, ts, notional_size, expected',
    [
        (1 + 1e-6, 1, 150, 3.56348e-5),
        (10, 1, 150, 0.106299),
        (100, 1, 150, 0.334178),
        (1000, 1, 150, 0.747793),
    ],
)
def test_calc_beta_ds(time, ts, notional_size, expected):
    """Test _calc_beta_ds function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage._calc_beta_ds(
            time, ts, notional_size),
        expected, rtol=1e-5, atol=1e-10
    )


@pytest.mark.parametrize(
    'fcm, expected',
    [
        (28, 1),
        (40, .98674),
        (60, .94753)
    ],
)
def test_calc_beta_s1(fcm, expected):
    """Test _calc_beta_s1 function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage._calc_beta_s1(fcm), expected,
        rtol=1e-5, atol=1e-10
    )


@pytest.mark.parametrize(
    'rh, beta_s1, expected',
    [
        (.99, 1, .25),
        (.9, .9, .25),
        (.5, 1, -1.35625),
        (.2, .9, -1.5376),  # Raises ValueError
        (99, 1, .25),
        (90, .9, .25),
        (50, 1, -1.35625),
        (20, .9, -1.5376),  # Raises ValueError
    ],
)
def test_calc_beta_RH(rh, beta_s1, expected):
    """Test _calc_beta_RH function."""
    try:
        assert np.isclose(
            _concrete_creep_and_shrinkage._calc_beta_RH(rh, beta_s1),
            expected, rtol=1e-5, atol=1e-10
        )
    except ValueError:
        with pytest.raises(ValueError) as exc_info:
            np.isclose(
                _concrete_creep_and_shrinkage._calc_beta_RH(
                    rh, beta_s1
                ), expected, rtol=1e-5, atol=1e-10
            )
            assert str(exc_info.value).startswith(
                    "The specified rh*beta_s1")


@pytest.mark.parametrize(
    'time, fcm, ts, notional_size, rh, cem_class, agg_type, expected',
    [
        (1+1e-6, 28, 1, 150, 99, '32.5 R', 'basalt', 4.201797e-9),
        (100, 28, 1, 150, 99, '32.5 R', 'basalt', 3.940385e-5),
        (1+1e-6, 28, 1, 150, 50, '32.5 R', 'basalt', -2.2794475e-8),
        (100, 28, 1, 150, 50, '32.5 R', 'basalt', -2.13766e-4)
    ],
)
def test_calc_drying_shrinkage(
    time, fcm, ts, notional_size, rh, cem_class, agg_type, expected):
    """Test calc_drying_shrinkage function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage.calc_drying_shrinkage(
            time, fcm, ts, notional_size, rh, cem_class, agg_type
        ), expected, rtol=1e-5, atol=1e-10
    )


@pytest.mark.parametrize(
    'fcm, TABULAR_VALUE, expected',
    [
        (28, {"alpha": 0, "S": 0.25, "alpha_as": 700,
            "alpha_ds1": 4, "alpha_ds2": 0.012, "alpha_e": 1.2},
            -3.997481e-5
        ),
        (40, {"alpha": 0, "S": 0.25, "alpha_as": 700,
            "alpha_ds1": 4, "alpha_ds2": 0.012, "alpha_e": 1.2},
            -7.083502e-5
        ),
        (60, {"alpha": 0, "S": 0.25, "alpha_as": 700,
            "alpha_ds1": 4, "alpha_ds2": 0.012, "alpha_e": 1.2},
            -1.237436e-4
        )
    ],
)
def test_calc_notional_autogenous_shrinkage(fcm, TABULAR_VALUE, expected):
    """Test _calc_notional_autogenous_shrinkage function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage._calc_notional_autogenous_shrinkage(
            fcm, TABULAR_VALUE
        ), expected, rtol=1e-5, atol=1e-10
    )


@pytest.mark.parametrize(
    'time, expected',
    [
        (1, .181269),
        (10, .468714),
        (100, .864665),
        (1000, .998208)
    ],
)
def test_calc_beta_au(time, expected):
    """Test _calc_beta_au function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage._calc_beta_au(time), expected,
        rtol=1e-5
    )


@pytest.mark.parametrize(
    'time, fcm, cem_class, agg_type, expected',
    [
        (1, 28, '32.5 R', 'basalt', -7.246194e-6),
        (10, 40, '32.5 R', 'basalt', -3.320137e-5),
        (100, 60, '32.5 R', 'basalt', -1.069968e-4)
    ],
)
def test_calc_autogenous_shrinkage(time, fcm, cem_class, agg_type, expected):
    """Test calc_autogenous_shrinkage function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage.calc_autogenous_shrinkage(
            time, fcm, cem_class, agg_type
        ), expected, rtol=1e-5, atol=1e-10
    )

# # Check total shrinkage with values given in table 5.1-13:
# # From fib MC2010: "In case where a lower level of accuracy is
# # sufficient, the values given in Table 5.1-13 can be accepted as
# # representative values for total shrinkage after 50 years ..."
# @pytest.mark.parametrize(
#     'time, fcm, ts, notional_size, rh, cem_class, agg_type, expected',
#     [
#         (50*365, 28, 1, 50, 50, '32.5 R', 'quartzite', -0.61),  # Is more than 10% off.
#         (50*365, 28, 1, 150, 50, '32.5 R', 'quartzite', -0.60),  # Is more than 10% off.
#         (50*365, 28, 1, 600, 50, '32.5 R', 'quartzite', -0.49),
#         (50*365, 28, 1, 50, 80, '32.5 R', 'quartzite', -0.38),
#         (50*365, 28, 1, 150, 80, '32.5 R', 'quartzite', -0.38),
#         (50*365, 28, 1, 600, 80, '32.5 R', 'quartzite', -0.31)
#     ],
# )

# def test_shrinkage(
#     time, fcm, ts, notional_size, rh, cem_class, agg_type, expected):
#     """Test full shrinkage calculation."""
#     assert np.isclose((
#         _concrete_creep_and_shrinkage.calc_autogenous_shrinkage(
#             time, fcm, cem_class, agg_type)
#         + _concrete_creep_and_shrinkage.calc_drying_shrinkage(
#             time, fcm, ts, notional_size, rh, cem_class, agg_type)
#         )*1e3, expected, rtol=1e-1
#     )


@pytest.mark.parametrize(
    'fcm, expected',
    [
        (28, .174687),
        (40, .136091),
        (60, .102463)
    ],
)
def test_calc_beta_bc_fcm(fcm, expected):
    """Test _calc_beta_bc_fcm function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage._calc_beta_bc_fcm(fcm),
        expected, rtol=1e-5
    )


@pytest.mark.parametrize(
    'time, t0, t0_adj, expected',
    [
        (14, 7, 6.98687, 4.88407),
        (100, 7, 6.98687, 7.46374),
        (100, 7, 12.09623, 6.37893)
    ],
)
def test_calc_beta_bc_t(time, t0, t0_adj, expected):
    """Test _calc_beta_bc_t function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage._calc_beta_bc_t(time, t0, t0_adj),
        expected, rtol=1e-5
    )


@pytest.mark.parametrize(
    'beta_bc_fcm, beta_bc_t, expected',
    [
        (.174687, 4.88407, .853184),
        (.136091, 7.46374, 1.015748),
        (.102463, 6.37893, .653604)
    ],
)
def test_calc_basic_creep_coefficient(beta_bc_fcm, beta_bc_t, expected):
    """Test calc_basic_creep_coefficient function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage.calc_basic_creep_coefficient(
            beta_bc_fcm, beta_bc_t
        ), expected, rtol=1e-5
    )


@pytest.mark.parametrize(
    'fcm, expected',
    [
        (28, 3.88040),
        (40, 2.35512),
        (60, 1.33501)
    ],
)
def test_calc_beta_dc_fcm(fcm, expected):
    """Test _calc_beta_dc_fcm function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage._calc_beta_dc_fcm(fcm),
        expected, rtol=1e-5
    )


@pytest.mark.parametrize(
    'rh, notional_size, expected',
    [
        (0.4, 150, 1.129243),
        (0.9, 150, .188207),
        (40, 150, 1.129243),
        (90, 150, .188207)
    ],
)
def test_calc_beta_dc_RH(rh, notional_size, expected):
    """Test _calc_beta_dc_RH function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage._calc_beta_dc_RH(rh, notional_size),
        expected, rtol=1e-5
    )


@pytest.mark.parametrize(
    't0_adj, expected',
    [
        (6.98687, .634832),
        (4.03566, .703308),
        (12.09623, .572613)
    ],
)
def test_calc_beta_dc_t0(t0_adj, expected):
    """Test _calc_beat_dc_t0 function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage._calc_beta_dc_t0(t0_adj),
        expected, rtol=1e-5
    )


@pytest.mark.parametrize(
    'fcm, expected',
    [
        (28, 1.118034),
        (40, .935414),
        (60, .763763)
    ],
)
def test_calc_alpha_fcm(fcm, expected):
    """Test _calc_alpha_fcm function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage._calc_alpha_fcm(fcm), expected,
        rtol=1e-5
    )


@pytest.mark.parametrize(
    'notional_size, alpha_fcm, expected',
    [
        (150, 1.118034, 504.5085),
        (1500, 1.118034, 1677.051)
    ],
)
def test_calc_beta_h(notional_size, alpha_fcm, expected):
    """Test _calc_beta_h function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage._calc_beta_h(
            notional_size, alpha_fcm
        ), expected, rtol=1e-5
    )


@pytest.mark.parametrize(
    't0_adj, expected',
    [
        (6.98687, 0.275929),
        (4.03566, 0.247387),
        (12.09623, 0.302450)
    ],
)
def test_calc_gamma_t0(t0_adj, expected):
    """Test _calc_gamma_t0 function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage._calc_gamma_t0(t0_adj),
        expected, rtol=1e-5
    )


@pytest.mark.parametrize(
    'time, t0, beta_h, gamma_t0, expected',
    [
        (14, 7, 504.5085, 0.275929, 0.30601),
        (14, 7, 1677.051, 0.247387, 0.25758),
        (100, 7, 504.5085, 0.302450, 0.56972)
    ],
)
def test_calc_beta_cd_t(time, t0, beta_h, gamma_t0, expected):
    """Test _calc_beta_cd_t function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage._calc_beta_dc_t(
            time, t0, beta_h, gamma_t0
        ), expected, rtol=1e-5
    )


@pytest.mark.parametrize(
    'beta_dc_fcm, beta_dc_RH, beta_dc_t0, beta_dc_t, expected',
    [
        (3.88040, 1.129243, .634832, 0.30601, 0.85125)
    ],
)
def test_calc_drying_creep_coefficient(
    beta_dc_fcm, beta_dc_RH, beta_dc_t0, beta_dc_t, expected):
    """Test calc_drying_creep_coefficient function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage.calc_drying_creep_coefficient(
            beta_dc_fcm, beta_dc_RH, beta_dc_t0, beta_dc_t
        ), expected, rtol=1e-5
    )


@pytest.mark.parametrize(
    'sigma, fcm, expected',
    [
        (0, 28, 0),
        (10, 28, .35714),
        (-10, 28, .35714),
        (20, 28, .71429),
        (-20, 28, .71429)
    ],
)
def test_calc_k_sigma(sigma, fcm, expected):
    """Test _calc_k_sigma function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage._calc_k_sigma(sigma, fcm),
        expected, rtol=1e-5
    )


@pytest.mark.parametrize(
    'time, t0, T_cur, fcm, rh, notional_size, sigma, agg_type,'\
    'cem_class, expected',
    [
        (14, 7, 20, 28, 40, 150, 10, 'basalt', '32.5 R', 1.704434),
        (14, 7, 20, 28, 40, 150, 14, 'basalt', '32.5 R', 1.980270),
        (100, 7, 20, 28, 40, 150, 10, 'basalt', '32.5 R', 2.968806)
    ],
)
def test_calc_creep_coefficient(
    time, t0, T_cur, fcm, rh, notional_size, sigma, agg_type, cem_class,
    expected):
    """Test calc_creep_coefficient function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage.calc_creep_coefficient(
            time, t0, T_cur, fcm, rh, notional_size, sigma, agg_type,
            cem_class
        ), expected, rtol=1e-5
    )


@pytest.mark.parametrize(
    't0, fcm, T_cur, phi, agg_type, cem_class, expected',
    [
        (7, 28, 20, 1.704434, 'basalt', '32.5 R', 7.8032605e-5),
        # (7, 28, 20, np.array([1.704434, 2.968806]), 'basalt', '32.5 R', np.array([7.8032605e-5, 1.1281e-4]))
    ],
)
def test_calc_creep_compliance(
    t0, fcm, T_cur, phi, agg_type, cem_class, expected):
    """Test calc_creep_compliance function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage.calc_creep_compliance(
            t0, fcm, T_cur, phi, agg_type, cem_class
        ), expected, rtol=1e-5, atol=1e-8
    )

