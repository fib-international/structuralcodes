"""Tests for EUROCODE 2-1-2:2004 creep and shrinkage."""

import numpy as np
import pytest

from structuralcodes.codes.ec2_2004 import _concrete_creep_and_shrinkage


@pytest.mark.parametrize(
    'test_h_0, test_f_cm, test_RH, test_c_class, test_t0, test_t, expected',
    [
        (138.5, 43, 50, 'R', 7, 18263, 2.567),
        (136.5, 38, 55, 'N', 7, 18263, 3.083),
        (136.5, 28, 55, 'N', 7, 18263, 3.748),
    ],
)
def test_creep_returns_expected_values(
    test_h_0, test_f_cm, test_RH, test_c_class, test_t0, test_t, expected
):
    """Test that phi returns expected values."""
    creep = _concrete_creep_and_shrinkage.phi(
        _concrete_creep_and_shrinkage.phi_0(
            _concrete_creep_and_shrinkage.phi_RH(
                test_h_0,
                test_f_cm,
                test_RH,
                _concrete_creep_and_shrinkage.alpha_1(test_f_cm),
                _concrete_creep_and_shrinkage.alpha_2(test_f_cm),
            ),
            _concrete_creep_and_shrinkage.beta_fcm(test_f_cm),
            _concrete_creep_and_shrinkage.beta_t0(
                _concrete_creep_and_shrinkage.t0_adj(
                    test_t0,
                    _concrete_creep_and_shrinkage.alpha_cement(test_c_class),
                )
            ),
        ),
        _concrete_creep_and_shrinkage.beta_c(
            _concrete_creep_and_shrinkage.t0_adj(
                test_t0,
                _concrete_creep_and_shrinkage.alpha_cement(test_c_class),
            ),
            test_t,
            _concrete_creep_and_shrinkage.beta_H(
                test_h_0,
                test_f_cm,
                test_RH,
                _concrete_creep_and_shrinkage.alpha_3(test_f_cm),
            ),
        ),
    )
    assert np.isclose(creep, expected, atol=1e-3)


def test_creep_array_input():
    """Test calculating creep number with array input."""
    # Arrange
    h_0 = 352
    f_cm = 51
    RH = 67
    cement_class = 's'
    t0 = 8
    t_final = 23 * 365
    num_time = 13
    t = np.linspace(0, t_final, num_time)

    # Act
    phi_RH = _concrete_creep_and_shrinkage.phi_RH(
        h_0,
        f_cm,
        RH,
        _concrete_creep_and_shrinkage.alpha_1(f_cm),
        _concrete_creep_and_shrinkage.alpha_2(f_cm),
    )
    beta_fcm = _concrete_creep_and_shrinkage.beta_fcm(f_cm)
    beta_t0 = _concrete_creep_and_shrinkage.beta_t0(
        _concrete_creep_and_shrinkage.t0_adj(
            t0,
            _concrete_creep_and_shrinkage.alpha_cement(cement_class),
        )
    )
    phi_0 = _concrete_creep_and_shrinkage.phi_0(phi_RH, beta_fcm, beta_t0)
    beta_c_array = _concrete_creep_and_shrinkage.beta_c(
        _concrete_creep_and_shrinkage.t0_adj(
            t0, _concrete_creep_and_shrinkage.alpha_cement(cement_class)
        ),
        t,
        _concrete_creep_and_shrinkage.beta_H(
            h_0,
            f_cm,
            RH,
            _concrete_creep_and_shrinkage.alpha_3(f_cm),
        ),
    )
    beta_c_float = _concrete_creep_and_shrinkage.beta_c(
        _concrete_creep_and_shrinkage.t0_adj(
            t0, _concrete_creep_and_shrinkage.alpha_cement(cement_class)
        ),
        t[-1],
        _concrete_creep_and_shrinkage.beta_H(
            h_0,
            f_cm,
            RH,
            _concrete_creep_and_shrinkage.alpha_3(f_cm),
        ),
    )

    creep = _concrete_creep_and_shrinkage.phi(phi_0, beta_c_array)
    creep_final = _concrete_creep_and_shrinkage.phi(phi_0, beta_c_float)

    # Assert
    assert len(creep) == num_time
    assert np.isclose(creep[-1], creep_final)


@pytest.mark.parametrize(
    'test_h_0, test_f_cm, test_c_class, test_RH, test_t_S, test_t, expected',
    [
        (138.5, 43, 'R ', 50, 28, 18263, 0.0006560211272528516),
        (138.5, 43, 'r', 50, 28, 18263, 0.0006560211272528516),
        (600.0, 43, 'N', 50, 28, 18263, 0.0003704813729755042),
        (600.0, 43, 'n', 50, 28, 18263, 0.0003704813729755042),
        (500.0, 43, 'N', 50, 28, 18263, 0.00037280025633482355),
        (136.8, 38, 'N ', 55, 28, 18263, 0.0004825589169449677),
        (100.0, 38, 'N', 55, 28, 18263, 0.000508432488714812),
        (90.0, 38, 'N', 55, 28, 18263, 0.0005085792190408255),
        (136.8, 28, 'N', 55, 28, 18263, 0.0005127088169780931),
    ],
)
def test_shrinkage_returns_expected_values(
    test_h_0, test_f_cm, test_c_class, test_RH, test_t_S, test_t, expected
):
    """Test that eps_cs returns expected values."""
    eps_ca = _concrete_creep_and_shrinkage.eps_ca(
        _concrete_creep_and_shrinkage.beta_as(test_t),
        _concrete_creep_and_shrinkage.eps_ca_inf(test_f_cm - 8),
    )
    eps_cd = _concrete_creep_and_shrinkage.eps_cd(
        _concrete_creep_and_shrinkage.beta_ds(test_t, test_t_S, test_h_0),
        _concrete_creep_and_shrinkage.k_h(test_h_0),
        _concrete_creep_and_shrinkage.eps_cd_0(
            _concrete_creep_and_shrinkage.alpha_ds1(test_c_class),
            _concrete_creep_and_shrinkage.alpha_ds2(test_c_class),
            test_f_cm,
            _concrete_creep_and_shrinkage.beta_RH(test_RH),
        ),
    )
    shrinkage = _concrete_creep_and_shrinkage.eps_cs(
        eps_ca=eps_ca, eps_cd=eps_cd
    )
    assert np.isclose(shrinkage, expected, atol=1e-6)


invalid_cement_classes = pytest.mark.parametrize(
    'cement_class',
    [
        'not a valid class',
        'rns',
    ],
)


@invalid_cement_classes
def test_creep_wrong_cement_class(cement_class):
    """Test that ValueError is raised if a wrong cement class is provided."""
    with pytest.raises(ValueError):
        _concrete_creep_and_shrinkage.alpha_cement(cement_class)


@invalid_cement_classes
def test_shrinkage_wrong_cement_class(cement_class):
    """Test that ValueError is raised if a wrong cement class is provided."""
    with pytest.raises(ValueError):
        _concrete_creep_and_shrinkage.alpha_ds1(cement_class)
    with pytest.raises(ValueError):
        _concrete_creep_and_shrinkage.alpha_ds2(cement_class)


def test_shrinkage_array_input():
    """Test calculating the shrinkage strain with an array as input."""
    # Arrange
    h_0 = 400
    f_cm = 53
    cement_class = 'n'
    RH = 70
    t_s = 14
    t_final = 35 * 365
    num_time = 10
    t = np.linspace(0, t_final, num_time)

    # Act
    eps_ca_inf = _concrete_creep_and_shrinkage.eps_ca_inf(f_cm - 8)
    beta_as_array = _concrete_creep_and_shrinkage.beta_as(t)
    beta_as_float = _concrete_creep_and_shrinkage.beta_as(t_final)
    eps_ca_array = _concrete_creep_and_shrinkage.eps_ca(
        beta_as_array, eps_ca_inf
    )
    eps_ca_float = _concrete_creep_and_shrinkage.eps_ca(
        beta_as_float, eps_ca_inf
    )
    beta_ds_array = _concrete_creep_and_shrinkage.beta_ds(t, t_s, h_0)
    beta_ds_float = _concrete_creep_and_shrinkage.beta_ds(t_final, t_s, h_0)
    alpha_ds1, alpha_ds2 = (
        _concrete_creep_and_shrinkage.alpha_ds1(cement_class),
        _concrete_creep_and_shrinkage.alpha_ds2(cement_class),
    )
    k_h = _concrete_creep_and_shrinkage.k_h(h_0)
    beta_RH = _concrete_creep_and_shrinkage.beta_RH(RH)
    eps_cd_0 = _concrete_creep_and_shrinkage.eps_cd_0(
        alpha_ds1,
        alpha_ds2,
        f_cm,
        beta_RH,
    )
    beta_ds_array = _concrete_creep_and_shrinkage.beta_ds(t, t_s, h_0)
    beta_ds_float = _concrete_creep_and_shrinkage.beta_ds(t_final, t_s, h_0)
    eps_cd_array = _concrete_creep_and_shrinkage.eps_cd(
        beta_ds_array, k_h, eps_cd_0
    )
    eps_cd_float = _concrete_creep_and_shrinkage.eps_cd(
        beta_ds_float, k_h, eps_cd_0
    )
    shrinkage = _concrete_creep_and_shrinkage.eps_cs(
        eps_ca=eps_ca_array, eps_cd=eps_cd_array
    )
    shrinkage_final = _concrete_creep_and_shrinkage.eps_cs(
        eps_ca=eps_ca_float, eps_cd=eps_cd_float
    )

    # Assert
    assert len(shrinkage) == num_time
    assert np.isclose(shrinkage[-1], shrinkage_final)


def test_maturity():
    """Test calculating the maturity of concrete."""
    # Arrange
    T = [19, 20, 20, 21, 22, 23]
    dt = [1, 2, 1, 1, 1, 2]

    # Act
    t_T = _concrete_creep_and_shrinkage.t_T(T=T, dt=dt)

    # Assert
    assert np.isclose(t_T, 8.37986781)


@pytest.mark.parametrize(
    'T, dt',
    (
        ([19, 20, 20, 21, 22, 23], [1, 2, 1, 1, 1]),
        ([19, 20, 20, 21, 22], [1, 2, 1, 1, 1, 2]),
    ),
)
def test_maturity_invalid_arrays(T, dt):
    """Test if maturity calculation raises an exception."""
    # Act and assert
    with pytest.raises(ValueError):
        _concrete_creep_and_shrinkage.t_T(T=T, dt=dt)


def test_calculate_h0():
    """Test calculating the effective section thickness."""
    # Arrange
    width = 250
    height = 350
    Ac = width * height
    u = 2 * (width + height)
    expected_h0 = 145.8333333

    # Act
    h_0 = _concrete_creep_and_shrinkage.h_0(Ac=Ac, u=u)

    # Assert
    assert np.isclose(h_0, expected_h0)
