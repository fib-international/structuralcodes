"""Tests for the _concrete_creep_and_shrinkage module."""

import warnings

import numpy as np
import pytest

from structuralcodes.codes.mc2010 import _concrete_creep_and_shrinkage
from structuralcodes.codes.mc2010._concrete_creep_and_shrinkage import (
    _check_age_at_loading,
    _check_env_temp,
    _check_fcm,
    _check_initial_stress,
    _check_RH,
)


@pytest.mark.parametrize(
    'fcm, expected',
    [
        (10, None),  # Raises ValueError
        (20, None),  # Raises ValueError
        (21, None),
        (129, None),
        (130, None),  # Raises ValueError
        (150, None),  # Raises ValueError
    ],
)
def test_check_fcm(fcm, expected):
    """Test check_fcm function."""
    try:
        assert _check_fcm(fcm) == expected
    except ValueError:
        with pytest.raises(ValueError) as exc_info:
            assert _check_fcm(fcm) == expected
        assert str(exc_info.value).startswith(
            'The specified mean compressive strength'
        )


@pytest.mark.parametrize(
    'sigma, fcm, expected',
    [(5, 20, None), (15, 20, None)],  # Raises ValueError
)
def test_check_initial_stress(sigma, fcm, expected):
    """Test check_initial_stress function."""
    try:
        assert _check_initial_stress(sigma, fcm) == expected
    except ValueError:
        with pytest.raises(ValueError) as exc_info:
            assert _check_initial_stress(sigma, fcm) == expected
        assert str(exc_info.value).startswith('The stress level exceeds')


def test_warning_initial_stress():
    """Test check_initial_stress function."""
    with warnings.catch_warnings(record=True) as w:
        _check_initial_stress(10, 20)
        assert len(w) == 1
        assert 'Initial stress is too high' in str(w[-1].message)


@pytest.mark.parametrize(
    't0, expected',
    [(0.5, None), (1, None)],  # raises ValueError
)
def test_check_age_at_loading(t0, expected):
    """Test check_age_at_loading function."""
    try:
        assert _check_age_at_loading(t0) == expected
    except ValueError:
        with pytest.raises(ValueError) as exc_info:
            assert _check_age_at_loading(t0) == expected
        assert str(exc_info.value).startswith('The load is applied too soon')


@pytest.mark.parametrize(
    'rh, expected',
    [
        (0.2, None),  # raises ValueError
        (0.5, None),
        (0.99, None),
        (1.2, None),  # raises ValueError
        (20, None),  # raises ValueError
        (40, None),
        (99, None),
    ],
)
def test_check_RH(rh, expected):
    """Test check_RH function."""
    try:
        assert _check_RH(rh) == expected
    except ValueError:
        with pytest.raises(ValueError) as exc_info:
            assert _check_RH(rh) == expected
        assert str(exc_info.value).startswith(
            'The specified relative humidity'
        )


def test_warning_env_temp():
    """Test check_initial_stress function."""
    with warnings.catch_warnings(record=True) as w:
        _check_env_temp(1)
        assert len(w) == 1
        assert 'The given environmental temperature' in str(w[-1].message)


@pytest.mark.parametrize(
    't0, T_cur, dt, expected',
    [
        (1, 20, None, 0.99812),
        (7, 20, None, 6.98687),
        (14, [20, 25], [7, 7], 15.7723),
        (28, 20, [7, 7, 7, 7], None),  # raises ValueError.
        (28, [20, 25], [7, 7], None),  # raises ValueError.
        (28, [20, 25], [28, 28], None),  # raises ValueError.
    ],
)
def test_t_T(t0, T_cur, dt, expected):
    """Test t_T function."""
    try:
        assert np.isclose(
            _concrete_creep_and_shrinkage.t_T(t0, T_cur, dt),
            expected,
            rtol=1e-5,
        )
    except ValueError:
        with pytest.raises(ValueError) as exc_info:
            assert _concrete_creep_and_shrinkage.t_T(t0, T_cur, dt) == expected
        assert str(exc_info.value).startswith('Dimensions of T_cur') or str(
            exc_info.value
        ).startswith('Curing time')

    except TypeError:
        with pytest.raises(TypeError) as exc_info:
            assert _concrete_creep_and_shrinkage.t_T(t0, T_cur, dt) == expected
        assert str(exc_info.value) == 'T_cur has to be provided as list.'


@pytest.mark.parametrize(  # Complete test cases
    '_tT, cem_class, expected',
    [
        (7, '32.5 N', 4.04647),
        (7, '32.5 R', 7),
        (7, '42.5 N', 7),
        (7, '42.5 R', 12.10932),
        (7, '52.5 N', 12.10932),
        (7, '52.5 R', 12.10932),
        (0, '32.5 N', 0.5),
    ],
)
def test_t0_adj(_tT, cem_class, expected):
    """Test t0_adj function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage.t0_adj(_tT, cem_class),
        expected,
        rtol=1e-5,
    )


@pytest.mark.parametrize(  # Complete test cases
    'fcm, cem_class, expected',
    [
        (28, '32.5 N', 3.82190e-4),
        (28, '32.5 R', 4.71651e-4),  # original
        (28, '42.5 N', 4.71651e-4),
        (28, '42.5 R', 6.28868e-4),
        (28, '52.5 N', 6.28868e-4),
        (28, '52.5 R', 6.28868e-4),
        (60, '32.5 R', 3.21256e-4),  # original
    ],
)
def test_eps_cds0(fcm, cem_class, expected):
    """Test eps_cds0 function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage.eps_cds0(fcm, cem_class),
        expected,
        rtol=1e-5,
        atol=1e-10,
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
def test_beta_ds(time, ts, notional_size, expected):
    """Test beta_ds function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage.beta_ds(time, ts, notional_size),
        expected,
        rtol=1e-5,
        atol=1e-10,
    )


@pytest.mark.parametrize(
    'fcm, expected',
    [(28, 1), (40, 0.98674), (60, 0.94753)],
)
def test_beta_s1(fcm, expected):
    """Test beta_s1 function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage.beta_s1(fcm),
        expected,
        rtol=1e-5,
        atol=1e-10,
    )


@pytest.mark.parametrize(
    'rh, beta_s1, expected',
    [
        (0.99, 1, 0.25),
        (0.9, 0.9, 0.25),
        (0.5, 1, -1.35625),
        (0.2, 0.9, -1.5376),  # Raises ValueError
        (99, 1, 0.25),
        (90, 0.9, 0.25),
        (50, 1, -1.35625),
        (20, 0.9, -1.5376),  # Raises ValueError
    ],
)
def test_beta_RH(rh, beta_s1, expected):
    """Test beta_RH function."""
    try:
        assert np.isclose(
            _concrete_creep_and_shrinkage.beta_RH(rh, beta_s1),
            expected,
            rtol=1e-5,
            atol=1e-10,
        )
    except ValueError:
        with pytest.raises(ValueError) as exc_info:
            np.isclose(
                _concrete_creep_and_shrinkage.beta_RH(rh, beta_s1),
                expected,
                rtol=1e-5,
                atol=1e-10,
            )
            assert str(exc_info.value).startswith('The specified rh*beta_s1')


@pytest.mark.parametrize(
    '_eps_cds0, _beta_ds, _beta_rh, expected',
    [
        (4.71651e-4, 3.56348e-5, 0.25, 4.201797e-9),
        (3.21256e-4, 0.106299, -1.35625, -4.6315e-5),
    ],
)
def test_eps_cds(_eps_cds0, _beta_ds, _beta_rh, expected):
    """Test eps_cds function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage.eps_cds(_eps_cds0, _beta_ds, _beta_rh),
        expected,
        rtol=1e-5,
        atol=1e-10,
    )


@pytest.mark.parametrize(
    'fcm, cem_class, expected',
    [
        (28, '32.5 N', -4.568550e-5),
        (28, '32.5 R', -3.997481e-5),
        (28, '42.5 N', -3.997481e-5),
        (28, '42.5 R', -3.426412e-5),
        (28, '52.5 N', -3.426412e-5),
        (28, '52.5 R', -3.426412e-5),
    ],
)
def test_eps_cas0(fcm, cem_class, expected):
    """Test eps_cas0 function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage.eps_cas0(fcm, cem_class),
        expected,
        rtol=1e-5,
        atol=1e-10,
    )


@pytest.mark.parametrize(
    'time, expected',
    [
        (1, 0.181269),
        (10, 0.468714),
        (100, 0.864665),
        (1000, 0.998208),
        (
            np.array([1, 10, 100, 1000]),
            np.array([0.181269, 0.468714, 0.864665, 0.998208]),
        ),
    ],
)
def test_beta_as(time, expected):
    """Test beta_as function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage.beta_as(time), expected, rtol=1e-5
    ).all()  # all() is to test numpy array.


@pytest.mark.parametrize(
    '_eps_cas0, _beta_as, expected',
    [
        (-3.997481e-5, 0.181269, -7.246194e-6),
        (-4.568550e-5, 0.468714, -2.141343e-5),
        (-3.426412e-5, 0.998208, -3.420272e-5),
    ],
)
def test_eps_cas(_eps_cas0, _beta_as, expected):
    """Test eps_cas function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage.eps_cas(_eps_cas0, _beta_as),
        expected,
        rtol=1e-5,
        atol=1e-10,
    )


@pytest.mark.parametrize(
    'fcm, expected',
    [(28, 0.174687), (40, 0.136091), (60, 0.102463)],
)
def test_beta_bc_fcm(fcm, expected):
    """Test beta_bc_fcm function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage.beta_bc_fcm(fcm),
        expected,
        rtol=1e-5,
    )


@pytest.mark.parametrize(
    'time, t0, t0_adj, expected',
    [
        (14, 7, 6.98687, 4.88407),
        (100, 7, 6.98687, 7.46374),
        (100, 7, 12.09623, 6.37893),
        (
            np.array([14, 100, 1000]),
            7,
            6.98687,
            np.array([4.88407, 7.46374, 9.83135]),
        ),
    ],
)
def test_beta_bc_t(time, t0, t0_adj, expected):
    """Test beta_bc_t function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage.beta_bc_t(time, t0, t0_adj),
        expected,
        rtol=1e-5,
    ).all()  # all() is to test np.array.


@pytest.mark.parametrize(
    'beta_bc_fcm, beta_bc_t, expected',
    [
        (0.174687, 4.88407, 0.853184),
        (0.136091, 7.46374, 1.015748),
        (0.102463, 6.37893, 0.653604),
        (
            0.174687,
            np.array([4.88407, 7.46374, 9.83135]),
            np.array([0.85318, 1.30381, 1.71741]),
        ),
    ],
)
def test_phi_bc(beta_bc_fcm, beta_bc_t, expected):
    """Test phi_bc function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage.phi_bc(beta_bc_fcm, beta_bc_t),
        expected,
        rtol=1e-5,
    ).all()  # all() is to check np.array.


@pytest.mark.parametrize(
    'fcm, expected',
    [(28, 3.88040), (40, 2.35512), (60, 1.33501)],
)
def test_beta_dc_fcm(fcm, expected):
    """Test beta_dc_fcm function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage.beta_dc_fcm(fcm),
        expected,
        rtol=1e-5,
    )


@pytest.mark.parametrize(
    'rh, notional_size, expected',
    [
        (0.4, 150, 1.129243),
        (0.9, 150, 0.188207),
        (40, 150, 1.129243),
        (90, 150, 0.188207),
    ],
)
def test_beta_dc_RH(rh, notional_size, expected):
    """Test beta_dc_RH function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage.beta_dc_RH(rh, notional_size),
        expected,
        rtol=1e-5,
    )


@pytest.mark.parametrize(
    't0_adj, expected',
    [(6.98687, 0.634832), (4.03566, 0.703308), (12.09623, 0.572613)],
)
def test_beta_dc_t0(t0_adj, expected):
    """Test beta_dc_t0 function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage.beta_dc_t0(t0_adj),
        expected,
        rtol=1e-5,
    )


@pytest.mark.parametrize(
    'fcm, expected',
    [(28, 1.118034), (40, 0.935414), (60, 0.763763)],
)
def test_alpha_fcm(fcm, expected):
    """Test alpha_fcm function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage.alpha_fcm(fcm), expected, rtol=1e-5
    )


@pytest.mark.parametrize(
    'notional_size, alpha_fcm, expected',
    [(150, 1.118034, 504.5085), (1500, 1.118034, 1677.051)],
)
def test_beta_h(notional_size, alpha_fcm, expected):
    """Test beta_h function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage.beta_h(notional_size, alpha_fcm),
        expected,
        rtol=1e-5,
    )


@pytest.mark.parametrize(
    't0_adj, expected',
    [(6.98687, 0.275929), (4.03566, 0.247387), (12.09623, 0.302450)],
)
def test_gamma_t0(t0_adj, expected):
    """Test gamma_t0 function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage.gamma_t0(t0_adj),
        expected,
        rtol=1e-5,
    )


@pytest.mark.parametrize(
    'time, t0, beta_h, gamma_t0, expected',
    [
        (14, 7, 504.5085, 0.275929, 0.30601),
        (14, 7, 1677.051, 0.247387, 0.25758),
        (100, 7, 504.5085, 0.302450, 0.56972),
    ],
)
def test_beta_dc_t(time, t0, beta_h, gamma_t0, expected):
    """Test beta_dc_t function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage.beta_dc_t(time, t0, beta_h, gamma_t0),
        expected,
        rtol=1e-5,
    )


@pytest.mark.parametrize(
    '_beta_dc_fcm, _beta_dc_RH, _beta_dc_t0, _beta_dc_t, expected',
    [
        (3.88040, 1.129243, 0.634832, 0.30601, 0.85125),
        (2.35512, 0.188207, 0.703308, 0.25758, 0.080298),
    ],
)
def test_phi_dc(_beta_dc_fcm, _beta_dc_RH, _beta_dc_t0, _beta_dc_t, expected):
    """Test phi_dc function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage.phi_dc(
            _beta_dc_fcm, _beta_dc_RH, _beta_dc_t0, _beta_dc_t
        ),
        expected,
        rtol=1e-5,
    )


@pytest.mark.parametrize(
    'sigma, fcm, expected',
    [
        (0, 28, 0),
        (10, 28, 0.35714),
        (-10, 28, 0.35714),
    ],
)
def test_k_sigma(sigma, fcm, expected):
    """Test k_sigma function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage.k_sigma(sigma, fcm),
        expected,
        rtol=1e-5,
    )


@pytest.mark.parametrize(
    '_phi_bc, _phi_dc, _sigma, _fcm, expected',
    [
        (0.853184, 0.85125, 10, 28, 1.70443),
        (0.853184, 0.85125, -15, 28, 2.08924),  # Provides a UserWarning.
        (0.653604, 0.08030, -10, 28, 0.733904),
    ],
)
def test_phi(_phi_bc, _phi_dc, _sigma, _fcm, expected):
    """Test phi function."""
    if abs(_sigma) / _fcm > 0.4:
        # Asserts that UserWarning is raised
        with pytest.warns(UserWarning):
            phi = _concrete_creep_and_shrinkage.phi(
                _phi_bc, _phi_dc, _sigma, _fcm
            )
    else:
        phi = _concrete_creep_and_shrinkage.phi(_phi_bc, _phi_dc, _sigma, _fcm)
    assert np.isclose(
        phi,
        expected,
        rtol=1e-5,
    )


@pytest.mark.parametrize(
    '_E_ci_t0, _phi, _E_ci, expected',
    [
        (23835.0, 1.70443, 27088.3, 1.04876e-4),
        (23835.0, 2.08924, 27088.3, 1.19082e-4),
    ],
)
def test_calc_creep_compliance(_E_ci_t0, _phi, _E_ci, expected):
    """Test calc_J function."""
    assert np.isclose(
        _concrete_creep_and_shrinkage.calc_J(_E_ci_t0, _phi, _E_ci),
        expected,
        rtol=1e-5,
        atol=1e-8,
    )
