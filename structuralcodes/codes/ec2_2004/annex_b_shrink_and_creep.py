"""Calculation routine for shrinkage and creep from EUROCODE 1992-1-1:2004
Annex B.
"""

import typing as t

import numpy as np
import numpy.typing as npt

from structuralcodes.codes.mc2010._concrete_creep_and_shrinkage import (
    t_T as t_T_mc2010,
)

ALPHA_CEMENT_DICT = {'R': 1.0, 'N': 0.0, 'S': -1.0}
ALPHA_DS_DICT = {
    'R': {'alpha_ds1': 6, 'alpha_ds2': 0.11},
    'N': {'alpha_ds1': 4, 'alpha_ds2': 0.12},
    'S': {'alpha_ds1': 3, 'alpha_ds2': 0.13},
}


def phi(phi_0: float, beta_c: npt.ArrayLike) -> npt.ArrayLike:
    """Calculate the creep number.

    EN 1992-1-1:2004, Eq. (B.1).

    Args:
        phi_0 (float): The standardized creep number defined in Eq. B.2.
        beta_c (npt.ArrayLike): A factor taking into account the creep
            development as a function of time after loading defined in Eq.
            (B.7).

    Returns:
        float: The creep number.
    """
    return phi_0 * beta_c


def phi_0(phi_RH: float, beta_fcm: float, beta_t0: float) -> float:
    """Calculate the standardized creep number.

    EN 1992-1-1:2004, Eq. (B.2).

    Args:
        phi_RH (float): The effect of relative humidity defined in Eq. B.3.
        beta_fcm (float): The effect of the concrete strength defined in Eq.
            B.4.
        beta_t0 (float): The effect of the age at loading defined in Eq. B.5.

    Returns:
        float: The standardized creep number.
    """
    return phi_RH * beta_fcm * beta_t0


def phi_RH(
    h_0: float, fcm: float, RH: float, alpha_1: float, alpha_2: float
) -> float:
    """Calculates the effect of relative humidity on the standardized creep
    number.

    EN 1992-1-1:2004, Eq. (B.3).

    Args:
        h_0 (float): The effective cross sectional thickness, Equation (B.6).
        fcm (float): The mean concrete strength in MPa.
        RH (float): The relative humidity in percent.
        alpha_1 (float): A factor describing the effect of concrete strength
            defined in Eq. (B.8c).
        alpha_2 (float): A factor describing the effect of concrete strength
            defined in Eq. (B.8c).

    Returns:
        float: The calculation parameter (B.3).
    """
    if fcm <= 35:
        return 1 + (1 - RH / 100) / (0.1 * h_0 ** (1 / 3))
    return (1 + (1 - RH / 100) / (0.1 * h_0 ** (1 / 3)) * alpha_1) * alpha_2


def beta_fcm(fcm: float) -> float:
    """Calculates the effect of the concrete strength on the standardized creep
    number.

    EN 1992-1-1:2004, Eq. (B.4).

    Args:
        fcm (float): The mean concrete strength in MPa.

    Returns:
        float: The effect of concrete strength.
    """
    return 16.8 / fcm**0.5


def beta_t0(t0: float) -> float:
    """Calculates the effect of age at loading on the standardized creep
    number.

    EN 1992-1-1:2004, Eq. (B.5).

    Args:
        t0 (float): The age at loading in days.

    Returns:
        float: The effect of age at loading.
    """
    return 1 / (0.1 + t0**0.20)


def h_0(Ac: float, u: float) -> float:
    """Calculates the effective thickness of the cross section.

    EN 1992-1-1:2004, Eq. (B.6).

    Args:
        Ac (float): The cross section area.
        u (float): The part of the circumference of the cross section subject
            to drying.

    Returns:
        float: The effective thickness.
    """
    return 2 * Ac / u


def beta_c(t0: float, t: npt.ArrayLike, beta_H: float) -> float:
    """Calculates the factor that describes the creep development as a function
    of time after loading.

    EN 1992-1-1:2004, Eq. (B.7).

    Args:
        t0 (float): The concrete age in days a the time of loading.
        t (ArrayLike): The concrete age at the evaluated time.
        beta_H (float): Parameter defined in (B.8).

    Returns:
        float: Parameter defined by Equation (B.7), beta_c.
    """
    t_load = np.atleast_1d(t - t0)
    t_load[t_load < 0.0] = 0.0
    return (t_load / (beta_H + t_load)) ** 0.3


def beta_H(h_0: float, fcm: float, RH: float, alpha_3) -> float:
    """Calculates the effect of relative humidity and the effective thickness
    of the structural element.

    EN 1992-1-1:2004, Eq. (B.8a and b).

    Args:
        h_0 (float): The effective cross sectional thickness, Equation (B.6).
        fcm (float): The mean concrete strength in MPa.
        RH (float): The relative humidity in percent.
        alpha_3 (float): A factor describing the effect of concrete strength
            defined in Eq. B.8c.

    Returns:
        float: The effect of humidity and the effective thickness of the
        element.
    """
    if fcm <= 35:
        # (B.8a)
        return min(1.5 * (1 + (0.012 * RH) ** 18) * h_0 + 250, 1500)
    # (B.8b)
    return min(
        1.5 * (1 + (0.012 * RH) ** 18) * h_0 + 250 * alpha_3,
        1500 * alpha_3,
    )


def alpha_1(fcm: float) -> float:
    """A factor taking into account the effect of concrete strength.

    EN 1992-1-1:2004, Eq. (B.8c).

    Args:
        fcm (float): The mean concrete strength in MPa.
    """
    return (35 / fcm) ** 0.7


def alpha_2(fcm: float) -> float:
    """A factor taking into account the effect of concrete strength.

    EN 1992-1-1:2004, Eq. (B.8c).

    Args:
        fcm (float): The mean concrete strength in MPa.
    """
    return (35 / fcm) ** 0.2


def alpha_3(fcm: float) -> float:
    """A factor taking into account the effect of concrete strength.

    EN 1992-1-1:2004, Eq. (B.8c).

    Args:
        fcm (float): The mean concrete strength in MPa.
    """
    return (35 / fcm) ** 0.5


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


def t_T(
    t: npt.ArrayLike, T_cur: npt.ArrayLike, dt: npt.ArrayLike = None
) -> npt.ArrayLike:
    """Calculates the maturity of the concrete.

    EN 1992-1-1:2004, Eq. (B.10).

    Args:
        t (npt.ArrayLike): The age of the concrete in days.
        T_cur (npt.ArrayLike): The temperature of the environment during curing
            in degrees Celcius.

    Keyword Args:
        dt (npt.ArrayLike): Number of days at which T_cur prevails. Required
            when providing a list for T_cur.

    Returns:
        np.ndarray: The temperature corrected age of the concrete in days at
        loading.
    """
    return t_T_mc2010(t0=t, T_cur=T_cur, dt=dt)


def alpha_cement(cement_class: t.Literal['S', 'N', 'R']) -> float:
    """An exponent that depends on the cement type.

    Args:
        cement_class (str): The cement class, either 'S', 'N' or 'R'.

    Returns:
        float: The exponent that depends on the cement type.

    Raises:
        ValueError: If an invalid cement class is provided.
    """
    _alpha_cement = ALPHA_CEMENT_DICT.get(cement_class.upper())

    if _alpha_cement is None:
        raise ValueError(
            (
                f'"{cement_class}" is not a valid cement class. '
                'Use either S, N or R.'
            )
        )
    return _alpha_cement


def eps_cs(
    h_0: float,
    f_cm: float,
    cement_class: t.Literal['R', 'N', 'S'] = 'R',
    RH: float = 50,
    t_S: float = 28,
    t: npt.ArrayLike = 18263,
) -> npt.ArrayLike:
    """Calculates the shrinkage strain.

    EN 1992-1-1:2004, Eq. (3.8).

    Args:
        h_0 (float): The effective cross sectional thickness, Equation (B.6).
        f_cm (float): The mean concrete strength.
        cement_class (str): The cement class, defaults to 'R'. Possible values:
            'R', 'N', 'S'.

    Keyword Args:
        RH (float): The relative humidity in percent, defaults to 50.
        t_S (float): the number of days when shrinkage begins, default: 28
            days.
        t (ArrayLike): the concrete age at the time (in days) of evaluation,
            default: 50 years.

    Returns:
        float: The shrinkage. Given as absolute, not in percent or ppm.

    Raises:
        ValueError: Checks if the cement class equals R, N or S.
    """
    _cement_class = cement_class.upper().strip()
    t_drying = np.atleast_1d(t - t_S)
    t_drying[t_drying < 0.0] = 0.0
    beta_ds = t_drying / (t_drying + 0.04 * h_0 ** (3 / 2))  # (3.10)
    beta_as = 1 - np.exp(-0.2 * t**0.5)  # (3.13)

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


def beta_RH(RH: float, RH_0: float = 100) -> float:
    """Calculates beta_RH.

    EN 1992-1-1:2004, Eq. (B.12).

    Args:
        RH (float): The relative humidity in percent.

    Keyword Args:
        RH_0 (float): The reference relative humidity, default: 100%.

    Returns:
        float: Calculation parameter from Equation (B.12).
    """
    return 1.55 * (1 - (RH / RH_0) ** 3)


def eps_cd_0(cement_class: str, f_cm: float, RH: float) -> float:
    """Calculates eps_cd_0.

    EN 1992-1-1:2004, Eq. (B.11).

    Args:
        cement_class (str): The cement class, defaults to 'R'. Possible values:
            'R', 'N', 'S'.
        f_cm (float): The mean concrete strength.
        RH (float): The relative humidity in percent.

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
        * ((220 + 110 * alpha_ds1) * np.exp(-alpha_ds2 * f_cm / 10))
        * 1e-6
        * beta_RH(RH)
    )  # (B.11)
