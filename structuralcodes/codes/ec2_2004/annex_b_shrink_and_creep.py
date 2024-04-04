"""Calculation routine for shrinkage and creep from EUROCODE 1992-1-1:2004
Annex B.
"""

# Importing only required modules from the packages
from math import e

from numpy import interp


def annex_B_creep(
    h_0: float,
    f_cm: float,
    RH: int = 50,
    cement_class: str = 'R',
    t0: int = 7,
    t: int = 18263,
) -> float:
    """This function will calculate creep according to Annex B in
     EUROCODE 1992-1-1:2004.

    Args:
        h_0 (float): The product of 2 times the cross-sectional area divided
        by the exposed circumference according to (B.6),
        f_cm (float): The mean concrete strength,
        RH (int): The relative humidity in percent, defaults to 50%,
        cement_class (str): The cement class, defaults to 'R'.
            Possible values: 'R', 'N', 'S',
        t0: the age of the concrete at the time (in days) of loading
        t: the age of the concrete at the time (in days) of evaluation
            (50 years default)

    Returns:
        phi(t, t0): the creep value for the load

    Raises:
        ValueError: checks if the cement class equals R, N or S.
    """
    # The cement class will decide the value of alpha_cement, ds1 and ds2.
    # The values are given in (B.9) and (B.12)
    alpha_cement_dict = {'R': 1, 'N': 0, 'S': -1}
    alpha_cement = alpha_cement_dict.get(cement_class)
    if alpha_cement is None:
        raise ValueError(f'cement_class={cement_class}, expected R, N or S')

    beta_H = calculate_beta_H(h_0, f_cm, RH)
    phi_RH = calculate_phi_RH(h_0, f_cm, RH)
    beta_fcm = calculate_beta_fcm(f_cm)
    t0_adjusted = calculate_t0_adjusted(t0, alpha_cement)
    beta_t0 = 1 / (0.1 + t0_adjusted**0.20)  # (B.5)
    beta_c = calculate_beta_c(t0, t, beta_H)

    phi_0 = phi_RH * beta_fcm * beta_t0  # (B.2)
    return phi_0 * beta_c


def calculate_beta_c(t0: float, t: float, beta_H: float) -> float:
    """Calculates the factor beta_c as defined by Equation (B.7).

    Args:
        t0 (float) the concrete age in days a the time of loading
        t (float) the concrete age at the evaluated time
        beta_H: parameter defined in (B.8)

    Return:
        beta_c (float): parameter defined by Equation (B.7)
    """
    return ((t - t0) / (beta_H + t - t0)) ** 0.3


def calculate_t0_adjusted(t0: float, alpha_cement: float) -> float:
    """Calculates the adjusted age of the concrete according to Annex B (2)
    Equation (B.9).

    Args:
        t0 (float): the concrete age in days at the time of loading
        alpha_cement (float): exponent derived from the sement type

    Return:
        float: the adjusted age of the concrete
    """
    return max(t0 * (9 / (2 + t0**1.2) + 1) ** alpha_cement, 0.5)


def calculate_beta_fcm(f_cm: float) -> float:
    """Calculates beta_f_cm according to Equation (B.4).

    Args:
        f_cm (float): the mean concrete strength

    Return:
        float: the factor defined in Equation (B.4)
    """
    return 16.8 / f_cm**0.5


def calculate_phi_RH(h_0: float, f_cm: float, RH: int) -> float:
    """This function calculates phi_RH according to Annex B in EUROCODE
    1992-1-1:2004.

    Args:
        h_0 (float): the effective cross sectional thickness, Equation (B.6)
        f_cm (float): The mean concrete strength,
        RH (int): The relative humidity in percent

    Returns:
        phi_RH (float): the calculation parameter (B.3)
    """
    # (B.8c) Alpha 1 to 2 is a constant where the only variable is f_cm
    alpha_1, alpha_2 = (
        (35 / f_cm) ** 0.7,
        (35 / f_cm) ** 0.2,
    )

    if f_cm <= 35:
        phi_RH = 1 + (1 - RH / 100) / (0.1 * h_0 ** (1 / 3))
    else:
        phi_RH = (
            1 + (1 - RH / 100) / (0.1 * h_0 ** (1 / 3)) * alpha_1
        ) * alpha_2
    return phi_RH


def calculate_beta_H(h_0: float, f_cm: float, RH: int) -> float:
    """This function will calculate eps_cd_0 according to Annex B in EUROCODE
    1992-1-1:2004.

    Args:
        h_0 (float): The cement class, defaults to 'R'.
            Possible values: 'R', 'N', 'S',
        f_cm (float): The mean concrete strength,
        RH (int): Calculation parameter from equation (B.12).

    Returns:
        beta_H (float): the calculation parameters defined in (B.8)
    """
    # (B.8c) Alpha 3 is a constant where f_cm is the only  variable
    alpha_3 = (35 / f_cm) ** 0.5

    if f_cm <= 35:
        # (B.8a) and (B.3a) applies
        beta_H = min(1.5 * (1 + (0.012 * RH) ** 18) * h_0 + 250, 1500)
    else:
        # (B.8b) and (B.3b) applies
        beta_H = min(
            1.5 * (1 + (0.012 * RH) ** 18) * h_0 + 250 * alpha_3,
            1500 * alpha_3,
        )
    return beta_H


def annex_B_shrinkage(
    h_0: float,
    f_cm: float,
    cement_class: str,
    RH: int = 50,
    t_S: int = 28,
    t: int = 18263,
) -> float:
    """This function will calculate shrinkage according to Chapter 3 by using
     the shrinkage as calculated by Annex B in EUROCODE 1992-1-1:2004.

    Args:
        h_0 (float): The product of 2 times the cross-sectional area divided by
        the exposed circumference according to (B.6),
        f_cm (float): The mean concrete strength,
        cement_class (str): The cement class, defaults to 'R'.
            Possible values: 'R', 'N', 'S',
        RH (int): The relative humidity in percent, defaults to 50,
        t_S (int): the number of days when shrinkage begins, default: 28 days,
        t (int): the concrete age at the time (in days) of evaluation,
            default: 50 years

    Returns:
        float: the shrinkage. Given as absolute, not in percent or ppm.
    """
    beta_ds = (t - t_S) / (t - t_S + 0.04 * h_0 ** (1 / 3))  # (3.10)
    beta_as = 1 - e ** (-0.2 * t**0.5)  # (3.13)

    # k_h is defined in Table 3.3 under (3.9)
    if h_0 >= 500:
        k_h = 0.70
    elif h_0 <= 100:
        k_h = 1.0
    else:
        k_h = interp(h_0, [100, 200, 300, 500], [1.0, 0.85, 0.75, 0.7])

    eps_ca_infinite = 2.5 * (f_cm - 18) * 1e-6  # (3.12)
    eps_ca = beta_as * eps_ca_infinite  # (3.11)
    eps_cd = (
        beta_ds * k_h * calculate_eps_cd_0(cement_class, f_cm, RH)
    )  # (3.9)
    return eps_cd + eps_ca  # (3.8)


def calculate_beta_RH(RH: int, RH_0: int = 100) -> float:
    """This function calculates beta_RH according to Annex B.12
    Args:
        RH (int): The relative humidity
        RH_0 (int): 100%
    Returns:
        beta_RH (float): Calculation parameter from Equation (B.12).
    """
    return 1.55 * (1 - (RH / RH_0) ** 3)


def calculate_eps_cd_0(cement_class: str, f_cm: float, RH: int) -> float:
    """This function will calculate eps_cd_0 according to Annex B in
    EUROCODE 1992-1-1:2004.

    Args:
        cement_class (str): The cement class, defaults to 'R'.
            Possible values: 'R', 'N', 'S',
        f_cm (float): The mean concrete strength,
        RH (int): The relative humidity.

    Returns:
        float: the nominal value for shrinkage

    Raises:
        ValueError: checks if the cement class equals R, N or S.
    """
    alpha_values = {
        'R': {'alpha_ds1': 6, 'alpha_ds2': 0.11},
        'N': {'alpha_ds1': 4, 'alpha_ds2': 0.12},
        'S': {'alpha_ds1': 3, 'alpha_ds2': 0.13},
    }
    alpha = alpha_values.get(cement_class)
    alpha_ds1 = alpha['alpha_ds1']
    alpha_ds2 = alpha['alpha_ds2']

    return (
        0.85
        * ((220 + 110 * alpha_ds1) * e ** (-alpha_ds2 * f_cm / 10))
        * 1e-6
        * calculate_beta_RH(RH)
    )  # (B.11)
