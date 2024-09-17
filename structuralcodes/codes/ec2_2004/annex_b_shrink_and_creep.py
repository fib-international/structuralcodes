"""Calculation routine for shrinkage and creep from EUROCODE 1992-1-1:2004
Annex B.
"""

import math
import typing as t

import numpy as np

ALPHA_CEMENT_DICT = {'R': 1, 'N': 0, 'S': -1}
ALPHA_DS_DICT = {
    'R': {'alpha_ds1': 6, 'alpha_ds2': 0.11},
    'N': {'alpha_ds1': 4, 'alpha_ds2': 0.12},
    'S': {'alpha_ds1': 3, 'alpha_ds2': 0.13},
}


def phi(
    h_0: float,
    f_cm: float,
    RH: int = 50,
    cement_class: t.Literal['R', 'N', 'S'] = 'R',
    t0: int = 7,
    t: int = 18263,
) -> float:
    """Calculates the creep number.

    EN 1992-1-1:2004, Eq. (B.1).

    Args:
        h_0 (float): The product of 2 times the cross-sectional area divided by
            the exposed circumference according to (B.6).
        f_cm (float): The mean concrete strength,
        RH (int): The relative humidity in percent, defaults to 50%.

    Keyword Args:
        cement_class (str): The cement class, defaults to 'R'. Possible values:
            'R', 'N', 'S',
        t0: The age of the concrete at the time (in days) of loading.
        t: The age of the concrete at the time (in days) of evaluation (50
            years default).

    Returns:
        float: The creep value for the load, phi(t, t0).

    Raises:
        ValueError: checks if the cement class equals R, N or S.
    """
    # The cement class will decide the value of alpha_cement, ds1 and ds2.
    # The values are given in (B.9) and (B.12)

    _cement_class = cement_class.upper().strip()
    alpha_cement = ALPHA_CEMENT_DICT.get(_cement_class)
    if alpha_cement is None:
        raise ValueError(f'cement_class={cement_class}, expected R, N or S')

    _beta_H = beta_H(h_0, f_cm, RH)
    _phi_RH = phi_RH(h_0, f_cm, RH)
    _beta_fcm = beta_fcm(f_cm)
    _t0_adj = t0_adj(t0, alpha_cement)
    beta_t0 = 1 / (0.1 + _t0_adj**0.20)  # (B.5)
    _beta_c = beta_c(t0, t, _beta_H)

    phi_0 = _phi_RH * _beta_fcm * beta_t0  # (B.2)
    return phi_0 * _beta_c


def beta_c(t0: float, t: float, beta_H: float) -> float:
    """Calculates the factor beta_c.

    EN 1992-1-1:2004, Eq. (B.7).

    Args:
        t0 (float): The concrete age in days a the time of loading.
        t (float): The concrete age at the evaluated time.
        beta_H: Parameter defined in (B.8).

    Returns:
        float: Parameter defined by Equation (B.7), beta_c.
    """
    return ((t - t0) / (beta_H + t - t0)) ** 0.3


def t0_adj(t0: float, alpha_cement: float) -> float:
    """Calculates the adjusted age of the concrete.

    EN 1992-1-1:2004, Eq. (B.9).

    Args:
        t0 (float): The concrete age in days at the time of loading.
        alpha_cement (float): Exponent derived from the sement type.

    Returns:
        float: The adjusted age of the concrete.
    """
    return max(t0 * (9 / (2 + t0**1.2) + 1) ** alpha_cement, 0.5)


def beta_fcm(f_cm: float) -> float:
    """Calculates beta_f_cm.

    EN 1992-1-1:2004, Eq. (B.4).

    Args:
        f_cm (float): The mean concrete strength.

    Returns:
        float: The factor defined in Equation (B.4).
    """
    return 16.8 / f_cm**0.5


def phi_RH(h_0: float, f_cm: float, RH: int) -> float:
    """Calculates phi_RH.

    EN 1992-1-1:2004, Eq. (B.3).

    Args:
        h_0 (float): The effective cross sectional thickness, Equation (B.6).
        f_cm (float): The mean concrete strength.
        RH (int): The relative humidity in percent.

    Returns:
        float: The calculation parameter (B.3).
    """
    # (B.8c) Alpha 1 to 2 is a constant where the only variable is f_cm
    alpha_1, alpha_2 = (
        (35 / f_cm) ** 0.7,
        (35 / f_cm) ** 0.2,
    )

    if f_cm <= 35:
        return 1 + (1 - RH / 100) / (0.1 * h_0 ** (1 / 3))
    return (1 + (1 - RH / 100) / (0.1 * h_0 ** (1 / 3)) * alpha_1) * alpha_2


def beta_H(h_0: float, f_cm: float, RH: int) -> float:
    """Calculates beta_H.

    EN 1992-1-1:2004, Eq. (B.8a and b).

    Args:
        h_0 (float): The effective cross sectional thickness, Equation (B.6).
        f_cm (float): The mean concrete strength.
        RH (int): The relative humidity in percent.

    Returns:
        float: The calculation parameter defined in (B.8).
    """
    # (B.8c) Alpha 3 is a constant where f_cm is the only  variable
    alpha_3 = (35 / f_cm) ** 0.5

    if f_cm <= 35:
        # (B.8a) and (B.3a) applies
        return min(1.5 * (1 + (0.012 * RH) ** 18) * h_0 + 250, 1500)
    # (B.8b) and (B.3b) applies
    return min(
        1.5 * (1 + (0.012 * RH) ** 18) * h_0 + 250 * alpha_3,
        1500 * alpha_3,
    )


def eps_cs(
    h_0: float,
    f_cm: float,
    cement_class: t.Literal['R', 'N', 'S'] = 'R',
    RH: int = 50,
    t_S: int = 28,
    t: int = 18263,
) -> float:
    """Calculates the shrinkage strain.

    EN 1992-1-1:2004, Eq. (3.8).

    Args:
        h_0 (float): The effective cross sectional thickness, Equation (B.6).
        f_cm (float): The mean concrete strength.
        cement_class (str): The cement class, defaults to 'R'. Possible values:
            'R', 'N', 'S'.

    Keyword Args:
        RH (int): The relative humidity in percent, defaults to 50.
        t_S (int): the number of days when shrinkage begins, default: 28 days.
        t (int): the concrete age at the time (in days) of evaluation, default:
            50 years.

    Returns:
        float: The shrinkage. Given as absolute, not in percent or ppm.

    Raises:
        ValueError: Checks if the cement class equals R, N or S.
    """
    _cement_class = cement_class.upper().strip()
    beta_ds = (t - t_S) / (t - t_S + 0.04 * h_0 ** (1 / 3))  # (3.10)
    beta_as = 1 - math.exp(-0.2 * t**0.5)  # (3.13)

    # k_h is defined in Table 3.3 under (3.9)
    if h_0 >= 500:
        k_h = 0.70
    elif h_0 <= 100:
        k_h = 1.0
    else:
        k_h = np.interp(h_0, [100, 200, 300, 500], [1.0, 0.85, 0.75, 0.7])

    eps_ca_infinite = 2.5 * (f_cm - 18) * 1e-6  # (3.12)
    eps_ca = beta_as * eps_ca_infinite  # (3.11)
    eps_cd = beta_ds * k_h * eps_cd_0(_cement_class, f_cm, RH)  # (3.9)
    return eps_cd + eps_ca  # (3.8)


def beta_RH(RH: int, RH_0: int = 100) -> float:
    """Calculates beta_RH.

    EN 1992-1-1:2004, Eq. (B.12).

    Args:
        RH (int): The relative humidity in percent.

    Keyword Args:
        RH_0 (int): The reference relative humidity, default: 100%.

    Returns:
        float: Calculation parameter from Equation (B.12).
    """
    return 1.55 * (1 - (RH / RH_0) ** 3)


def eps_cd_0(cement_class: str, f_cm: float, RH: int) -> float:
    """Calculates eps_cd_0.

    EN 1992-1-1:2004, Eq. (B.11).

    Args:
        cement_class (str): The cement class, defaults to 'R'. Possible values:
            'R', 'N', 'S'.
        f_cm (float): The mean concrete strength.
        RH (int): The relative humidity in percent.

    Returns:
        float: The nominal value for shrinkage.

    Raises:
        ValueError: Checks if the cement class equals R, N or S.
    """
    _cement_class = cement_class.upper().strip()
    alpha = ALPHA_DS_DICT.get(_cement_class)
    if alpha is None:
        raise ValueError(f'cement_class={cement_class}, expected R, N or S')
    alpha_ds1 = alpha['alpha_ds1']
    alpha_ds2 = alpha['alpha_ds2']

    return (
        0.85
        * ((220 + 110 * alpha_ds1) * math.exp(-alpha_ds2 * f_cm / 10))
        * 1e-6
        * beta_RH(RH)
    )  # (B.11)
