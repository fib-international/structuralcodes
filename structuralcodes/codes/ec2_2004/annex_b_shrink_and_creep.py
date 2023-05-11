""" Calculation routine for shrinkage and creep from EUROCODE 1992-1-1:2004
Annex B """

# Importing only required modules from the packages
from math import e
from numpy import interp


def annex_B_creep(
    h_0: float,
    fck: float,
    RH: int,
    cement_class: str = 'R',
    t_dead_load: int = 7,
    t_live_load: int = 28,
) -> list[2]:
    """This function will calculate creep according to Annex B in EUROCODE 1992-1-1:2004.

    Args:
        h_0 (float): The product of 2 times the cross-sectional area divided by the
        exposed circumference according to (B.6),
        fck (float): The characteristic concrete strength,
        RH (int): The relative humidity, defaults to 50%,
        cement_class (str): The cement class, defaults to 'R'. Possible values: 'R', 'N', 'S',
        t_dead_load (int): the number of days before the dead load is applied, defaults to 7 days,
        t_live_load (int): the number of days before the live load is applied, defaults to 28 days.

    Returns:
        list[2]: a vector containing the two creep values for the dead and live load, respectively.

    Raises:
        ValueError: checks if the cement class equals R, N or S.
    """
    f_cm = fck + 8

    # The cement class will decide the value of the three constants alpha_cement, ds1 and ds2.
    # The values are given in (B.9) and (B.12)
    if cement_class == 'R':
        alpha_cement = 1
    elif cement_class == 'N':
        alpha_cement = 0
    elif cement_class == 'S':
        alpha_cement = -1
    else:
        raise ValueError(f'cement_class={cement_class}, expected R, N or S')

    beta_H, phi_RH = calculate_phiRH_betaH(h_0, f_cm, RH)
    beta_fcm = 16.8 / f_cm**0.5  # (B.4)

    # Collects and readies the time variables for the applied dead and live load
    t_forever = 18262  # OVE SLETTEN AS software defaults to 25 years, the support doesn't know why
    t0 = [t_dead_load, t_live_load]

    # This part calculates creep for each of the two load cases
    phi_t_t0 = [0, 0]  # just an initialization
    for i in range(2):
        t0_adjusted = max(
            t0[i] * (9 / (2 + t0[i] ** 1.2) + 1) ** alpha_cement, 0.5
        )  # (B.9)
        beta_t0 = 1 / (0.1 + t0_adjusted**0.20)  # (B.5)
        beta_c = (
            (t_forever - t0[i]) / (beta_H + t_forever - t0[i])
        ) ** 0.3  # (B.7)
        phi_0 = phi_RH * beta_fcm * beta_t0  # (B.2)
        phi_t_t0[i] = phi_0 * beta_c  # (B.1)
    return phi_t_t0


def calculate_phiRH_betaH(h_0, f_cm, RH):
    """This function will calculate eps_cd_0 according to Annex B in EUROCODE 1992-1-1:2004.

    Args:
        h_0 (float): The cement class, defaults to 'R'. Possible values: 'R', 'N', 'S',
        f_cm (float): The characteristic concrete strength,
        RH (float): Calculation parameter from equation (B.12).

    Returns:
        (float, float): the calculation parameters defined in (B.3) and (B.8)

    """
    # (B.8c) Alpha 1 to 3 is a constant where only f_cm is a variable
    alpha_1, alpha_2, alpha_3 = (
        (35 / f_cm) ** 0.7,
        (35 / f_cm) ** 0.2,
        (35 / f_cm) ** 0.5,
    )  # (B.8)

    if f_cm <= 35:
        # (B.8a) and (B.3a) applies
        beta_H = min(1.5 * (1 + (0.012 * RH) ** 18) * h_0 + 250, 1500)
        phi_RH = 1 + (1 - RH / 100) / (0.1 * h_0 ** (1 / 3))
    else:
        # (B.8b) and (B.3b) applies
        beta_H = min(
            1.5 * (1 + (0.012 * RH) ** 18) * h_0 + 250 * alpha_3,
            1500 * alpha_3,
        )
        phi_RH = (
            1 + (1 - RH / 100) / (0.1 * h_0 ** (1 / 3)) * alpha_1
        ) * alpha_2
    return beta_H, phi_RH


def annex_B_shrinkage(
    h_0: float,
    fck: float,
    RH: int,
    cement_class: str,
    t_S: int = 28,
) -> float:
    """This function will calculate shrinkage according to Chapter 3 by using
     the shrinkage as calculated by Annex B in EUROCODE 1992-1-1:2004.

    Args:
        h_0 (float): The product of 2 times the cross-sectional area divided by the
        exposed circumference according to (B.6),
        fck (float): The characteristic concrete strength,
        RH (int): The relative humidity, defaults to 50%,
        cement_class (str): The cement class, defaults to 'R'. Possible values: 'R', 'N', 'S',
        t_S (int): the number of days when shrinkage begins, defaults to 28 days.

    Returns:
        float: the shrinkage. Given as absolute, not in percent or ppm.
    """
    t_forever = 18262  # OVE SLETTEN AS software defaults to 25 years, the support doesn't know why

    beta_RH = 1.55 * (1 - (RH / 100) ** 3)  # (B.12)
    beta_ds = (t_forever - t_S) / (
        t_forever - t_S + 0.04 * h_0 ** (1 / 3)
    )  # (3.10)
    beta_as = 1 - e ** (-0.2 * t_forever**0.5)  # (3.13)

    # k_h is defined in Table 3.3 under (3.9)
    if h_0 >= 500:
        k_h = 0.70
    elif h_0 <= 100:
        k_h = 1.0
    else:
        k_h = interp(h_0, [100, 200, 300, 500], [1.0, 0.85, 0.75, 0.7])

    eps_ca_evig = 2.5 * (fck - 10) * 1e-6  # (3.12)
    eps_ca = beta_as * eps_ca_evig  # (3.11)
    eps_cd = (
        beta_ds * k_h * calculate_eps_cd_0(cement_class, fck, beta_RH)
    )  # (3.9)
    eps_cs = eps_cd + eps_ca  # (3.8)

    return eps_cs


def calculate_eps_cd_0(cement_class, fck, beta_RH):
    """This function will calculate eps_cd_0 according to Annex B in EUROCODE 1992-1-1:2004.

    Args:
        cement_class (str): The cement class, defaults to 'R'. Possible values: 'R', 'N', 'S',
        fck (float): The characteristic concrete strength,
        beta_RH (float): Calculation parameter from equation (B.12).

    Returns:
        float: the nominal value for shrinkage

    Raises:
        ValueError: checks if the cement class equals R, N or S.
    """
    f_cm = fck + 8
    if cement_class == 'R':
        alpha_ds1 = 6
        alpha_ds2 = 0.11
    elif cement_class == 'N':
        alpha_ds1 = 4
        alpha_ds2 = 0.12
    elif cement_class == 'S':
        alpha_ds1 = 3
        alpha_ds2 = 0.13
    else:
        raise ValueError(f'cement_class={cement_class}, expected R, N or S')

    eps_cd_0 = (
        0.85
        * ((220 + 110 * alpha_ds1) * e ** (-alpha_ds2 * f_cm / 10))
        * 1e-6
        * beta_RH
    )  # (B.11)
    return eps_cd_0
