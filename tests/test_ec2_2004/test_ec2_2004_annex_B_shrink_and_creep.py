"""Tests for EUROCODE 2-1-2:2004 Annex B and chapter 3."""

import numpy as np
import pytest

from structuralcodes.codes.ec2_2004 import annex_b_shrink_and_creep


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
    creep = annex_b_shrink_and_creep.phi(
        annex_b_shrink_and_creep.phi_0(
            annex_b_shrink_and_creep.phi_RH(
                test_h_0,
                test_f_cm,
                test_RH,
                annex_b_shrink_and_creep.alpha_1(test_f_cm),
                annex_b_shrink_and_creep.alpha_2(test_f_cm),
            ),
            annex_b_shrink_and_creep.beta_fcm(test_f_cm),
            annex_b_shrink_and_creep.beta_t0(
                annex_b_shrink_and_creep.t0_adj(
                    test_t0,
                    annex_b_shrink_and_creep.alpha_cement(test_c_class),
                )
            ),
        ),
        annex_b_shrink_and_creep.beta_c(
            annex_b_shrink_and_creep.t0_adj(
                test_t0, annex_b_shrink_and_creep.alpha_cement(test_c_class)
            ),
            test_t,
            annex_b_shrink_and_creep.beta_H(
                test_h_0,
                test_f_cm,
                test_RH,
                annex_b_shrink_and_creep.alpha_3(test_f_cm),
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
    phi_RH = annex_b_shrink_and_creep.phi_RH(
        h_0,
        f_cm,
        RH,
        annex_b_shrink_and_creep.alpha_1(f_cm),
        annex_b_shrink_and_creep.alpha_2(f_cm),
    )
    beta_fcm = annex_b_shrink_and_creep.beta_fcm(f_cm)
    beta_t0 = annex_b_shrink_and_creep.beta_t0(
        annex_b_shrink_and_creep.t0_adj(
            t0,
            annex_b_shrink_and_creep.alpha_cement(cement_class),
        )
    )
    phi_0 = annex_b_shrink_and_creep.phi_0(phi_RH, beta_fcm, beta_t0)
    beta_c_array = annex_b_shrink_and_creep.beta_c(
        annex_b_shrink_and_creep.t0_adj(
            t0, annex_b_shrink_and_creep.alpha_cement(cement_class)
        ),
        t,
        annex_b_shrink_and_creep.beta_H(
            h_0,
            f_cm,
            RH,
            annex_b_shrink_and_creep.alpha_3(f_cm),
        ),
    )
    beta_c_float = annex_b_shrink_and_creep.beta_c(
        annex_b_shrink_and_creep.t0_adj(
            t0, annex_b_shrink_and_creep.alpha_cement(cement_class)
        ),
        t[-1],
        annex_b_shrink_and_creep.beta_H(
            h_0,
            f_cm,
            RH,
            annex_b_shrink_and_creep.alpha_3(f_cm),
        ),
    )

    creep = annex_b_shrink_and_creep.phi(phi_0, beta_c_array)
    creep_final = annex_b_shrink_and_creep.phi(phi_0, beta_c_float)

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
    eps_ca = annex_b_shrink_and_creep.eps_ca(
        annex_b_shrink_and_creep.beta_as(test_t),
        annex_b_shrink_and_creep.eps_ca_inf(test_f_cm - 8),
    )
    eps_cd = annex_b_shrink_and_creep.eps_cd(
        annex_b_shrink_and_creep.beta_ds(test_t, test_t_S, test_h_0),
        annex_b_shrink_and_creep.k_h(test_h_0),
        annex_b_shrink_and_creep.eps_cd_0(
            annex_b_shrink_and_creep.alpha_ds1(test_c_class),
            annex_b_shrink_and_creep.alpha_ds2(test_c_class),
            test_f_cm,
            annex_b_shrink_and_creep.beta_RH(test_RH),
        ),
    )
    shrinkage = annex_b_shrink_and_creep.eps_cs(eps_ca=eps_ca, eps_cd=eps_cd)
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
        annex_b_shrink_and_creep.alpha_cement(cement_class)


@invalid_cement_classes
def test_shrinkage_wrong_cement_class(cement_class):
    """Test that ValueError is raised if a wrong cement class is provided."""
    with pytest.raises(ValueError):
        annex_b_shrink_and_creep.alpha_ds1(cement_class)
    with pytest.raises(ValueError):
        annex_b_shrink_and_creep.alpha_ds2(cement_class)


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
    eps_ca_inf = annex_b_shrink_and_creep.eps_ca_inf(f_cm - 8)
    beta_as_array = annex_b_shrink_and_creep.beta_as(t)
    beta_as_float = annex_b_shrink_and_creep.beta_as(t_final)
    eps_ca_array = annex_b_shrink_and_creep.eps_ca(beta_as_array, eps_ca_inf)
    eps_ca_float = annex_b_shrink_and_creep.eps_ca(beta_as_float, eps_ca_inf)
    beta_ds_array = annex_b_shrink_and_creep.beta_ds(t, t_s, h_0)
    beta_ds_float = annex_b_shrink_and_creep.beta_ds(t_final, t_s, h_0)
    alpha_ds1, alpha_ds2 = (
        annex_b_shrink_and_creep.alpha_ds1(cement_class),
        annex_b_shrink_and_creep.alpha_ds2(cement_class),
    )
    k_h = annex_b_shrink_and_creep.k_h(h_0)
    beta_RH = annex_b_shrink_and_creep.beta_RH(RH)
    eps_cd_0 = annex_b_shrink_and_creep.eps_cd_0(
        alpha_ds1,
        alpha_ds2,
        f_cm,
        beta_RH,
    )
    beta_ds_array = annex_b_shrink_and_creep.beta_ds(t, t_s, h_0)
    beta_ds_float = annex_b_shrink_and_creep.beta_ds(t_final, t_s, h_0)
    eps_cd_array = annex_b_shrink_and_creep.eps_cd(
        beta_ds_array, k_h, eps_cd_0
    )
    eps_cd_float = annex_b_shrink_and_creep.eps_cd(
        beta_ds_float, k_h, eps_cd_0
    )
    shrinkage = annex_b_shrink_and_creep.eps_cs(
        eps_ca=eps_ca_array, eps_cd=eps_cd_array
    )
    shrinkage_final = annex_b_shrink_and_creep.eps_cs(
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
    t_T = annex_b_shrink_and_creep.t_T(T=T, dt=dt)

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
        annex_b_shrink_and_creep.t_T(T=T, dt=dt)


def test_calculate_h0():
    """Test calculating the effective section thickness."""
    # Arrange
    width = 250
    height = 350
    Ac = width * height
    u = 2 * (width + height)
    expected_h0 = 145.8333333

    # Act
    h_0 = annex_b_shrink_and_creep.h_0(Ac=Ac, u=u)

    # Assert
    assert np.isclose(h_0, expected_h0)
