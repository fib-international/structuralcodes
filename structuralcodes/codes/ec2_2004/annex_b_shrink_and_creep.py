"""Calculation routine for shrinkage and creep from EUROCODE 1992-1-1:2004
Annex B.
"""

import typing as t

import numpy as np
import numpy.typing as npt

ALPHA_CEMENT_DICT = {'R': 1.0, 'N': 0.0, 'S': -1.0}
ALPHA_DS_DICT = {
    'R': {'alpha_ds1': 6, 'alpha_ds2': 0.11},
    'N': {'alpha_ds1': 4, 'alpha_ds2': 0.12},
    'S': {'alpha_ds1': 3, 'alpha_ds2': 0.13},
}


def eps_cs(eps_cd: npt.ArrayLike, eps_ca: npt.ArrayLike) -> npt.ArrayLike:
    """Calculates the total shrinkage strain.

    EN 1992-1-1:2004, Eq. (3.8).

    Args:
        eps_cd (npt.ArrayLike): The drying shrinkage defined in Eq. (3.9).
        eps_ca (npt.ArrayLike): The autogenous shrinkage defined in Eq. (3.11).

    Returns:
        npt.ArrayLike: The total shrinkage strain.
    """
    return eps_cd + eps_ca


def eps_cd(
    beta_ds: npt.ArrayLike, k_h: float, eps_cd_0: float
) -> npt.ArrayLike:
    """Calculates the drying shrinkage.

    EN 1992-1-1:2004, Eq. (3.9).

    Args:
        beta_ds (npt.ArrayLike): A coefficient taking into account the time of
            drying defined in Eq. (3.10).
        k_h (float): A coefficient depending on the effective thickness of the
            section defined in Tab. 3.3.
        eps_cd_0 (float): The nominal value of drying shrinkage defined in Eq.
            (B.11).

    Returns:
        npt.ArrayLike: The drying shrinkage.

    Note:
        In EC2 (2004), the shrinkage strain is calculated as a positive number.
    """
    return beta_ds * k_h * eps_cd_0


def beta_ds(t: npt.ArrayLike, t_s: float, h_0: float):
    """Calculates the coefficient taking into account the time of drying.

    EN 1992-1-1:2004, Eq. (3.10).

    Args:
        t (npt.ArrayLike): The age of the concrete in days.
        t_s (float): The age of the concrete in days at the start of drying.
        h_0 (float): The effective cross section thichkness.

    Returns:
        npt.ArrayLike: The coefficient taking into account the time of drying.
    """
    t_drying = np.atleast_1d(t - t_s)
    t_drying[t_drying < 0.0] = 0.0
    return t_drying / (t_drying + 0.04 * h_0 ** (3 / 2))


def k_h(h_0: float) -> float:
    """Calculate the coefficient depending on the effective section thickness.

    EN 1992-1-1:2004, Tab. 3.3.

    Args:
        h_0 (float): The effective section thickness in mm.

    Returns:
        float: The coefficient depending on the effective section thickness.
    """
    if h_0 >= 500:
        k_h = 0.70
    elif h_0 <= 100:
        k_h = 1.0
    else:
        k_h = np.interp(h_0, [100, 200, 300, 500], [1.0, 0.85, 0.75, 0.7])
    return k_h


def eps_cd_0(
    alpha_ds1: float,
    alpha_ds2: float,
    fcm: float,
    beta_RH: float,
    fcm_0: float = 10,
) -> float:
    """Calculates the nominal value of drying shrinkage.

    EN 1992-1-1:2004, Eq. (B.11).

    Args:
        alpha_ds1 (float): A coefficient depending on the cement type, defined
            in EC2 (2004), Sec. B.2.
        alpha_ds2 (float): A coefficient depending on the cement type, defined
            in EC2 (2004), Sec. B.2.
        fcm (float): The mean compressive strength in MPa.
        beta_RH (float): A factor describing the effect of relative humidity,
            defined in Eq. (B.12).

    Keyword Args:
        fcm_0 (float): A reference strength in MPa, default 10 MPa.

    Returns:
        float: The nominal value of drying shrinkage.
    """
    return (
        0.85
        * ((220 + 110 * alpha_ds1) * np.exp(-alpha_ds2 * fcm / fcm_0))
        * 1e-6
        * beta_RH
    )


def alpha_ds1(cement_class: t.Literal['S', 'N', 'R']) -> float:
    """A coefficient depending on the cement class.

    EN 1992-1-1:2004, Sec. B.2.

    Args:
        cement_class (str): The cement class, either 'S', 'N' or 'R'.

    Returns:
        float: The exponent that depends on the cement type.

    Raises:
        ValueError: If an invalid cement class is provided.
    """
    _alpha_ds = _get_alpha_ds_dict(cement_class=cement_class)
    return _alpha_ds['alpha_ds1']


def alpha_ds2(cement_class: t.Literal['S', 'N', 'R']) -> float:
    """A coefficient depending on the cement class.

    EN 1992-1-1:2004, Sec. B.2.

    Args:
        cement_class (str): The cement class, either 'S', 'N' or 'R'.

    Returns:
        float: The exponent that depends on the cement type.

    Raises:
        ValueError: If an invalid cement class is provided.
    """
    _alpha_ds = _get_alpha_ds_dict(cement_class=cement_class)
    return _alpha_ds['alpha_ds2']


def _get_alpha_ds_dict(cement_class: str) -> t.Dict:
    """Internal function for getting a dictionary with values for aplha_ds1 and
    alpha_ds2.
    """
    _alpha_ds = ALPHA_DS_DICT.get(cement_class.upper().strip())

    if _alpha_ds is None:
        raise ValueError(
            (
                f'"{cement_class}" is not a valid cement class. '
                'Use either S, N or R.'
            )
        )
    return _alpha_ds


def beta_RH(RH: float, RH_0: float = 100) -> float:
    """Calculates the factor describing the effect of relative humidity.

    EN 1992-1-1:2004, Eq. (B.12).

    Args:
        RH (float): The relative humidity in percent.

    Keyword Args:
        RH_0 (float): The reference relative humidity, default: 100%.

    Returns:
        float: The factor taking into account the relative humidity.
    """
    return 1.55 * (1 - (RH / RH_0) ** 3)


def eps_ca(beta_as: npt.ArrayLike, eps_ca_inf: float) -> npt.ArrayLike:
    """Calculates the autogenous shrinkage.

    EN 1992-1-1:2004, Eq. (3.11).

    Args:
        beta_as (npt.ArrayLike): A factor describing the autogenous shrinkage
            development.
        eps_ca_inf (float): The final autogenous shrinkage.

    Returns:
        npt.ArrayLike: The autogenous shrinkage.

    Note:
        In EC2 (2004), the shrinkage strain is calculated as a positive number.
    """
    return beta_as * eps_ca_inf


def eps_ca_inf(fck: float) -> float:
    """Calculates the final autogenous shrinkage.

    EN 1992-1-1:2004, Eq. 3.12.

    Args:
        fck (float): The characteristic compressive strength in MPa.

    Returns:
        float: The final autogenous shrinkage.
    """
    return 2.5 * (fck - 10) * 1e-6


def beta_as(t: npt.ArrayLike) -> npt.ArrayLike:
    """Calculates the factor describing the development of autogenous
    shrinkage.

    EN 1992-1-1:2004, Eq. (3.13).

    Args:
        t (npt.ArrayLike): The age of the concrete in days.

    Returns:
        npt.ArrayLike: The factor describing the development of autogenous
        shrinkage.
    """
    return 1 - np.exp(-0.2 * t**0.5)


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


def t_T(T: npt.ArrayLike, dt: npt.ArrayLike) -> float:
    """Calculates the maturity of the concrete.

    EN 1992-1-1:2004, Eq. (B.10).

    Args:
        T (npt.ArrayLike): The curing temperature history in degrees Celcius.
        dt (npt.ArrayLike): The number of days with temperature T.

    Returns:
        float: The maturity of the concrete.

    Note:
        The two arrays T and dt should have the same length.
    """
    # Prepare the input
    T = np.atleast_1d(T)
    dt = np.atleast_1d(dt)

    # Check the shape of the input arrays
    if T.shape != dt.shape:
        raise ValueError(
            f'T {T.shape} and dt {dt.shape} should have the same shape.'
        )

    # Return the sum of the temperature adjusted time increments
    return float(np.sum(np.exp(13.65 - 4000 / (273 + T)) * dt))


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
