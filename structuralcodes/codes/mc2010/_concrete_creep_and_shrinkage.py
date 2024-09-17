"""Creep and shrinkage models MC2010 for concrete.

[1] fib Model Code for Concrete Structures 2010 (2013).

Todo:
 - Add the light-weight aggregate formulas
 - Include the additional creep and shrinkage terms for other temperatures
   (see 5.1.10.7 of [1])
"""

from __future__ import annotations  # To have clean hints of ArrayLike in docs

import typing as t
import warnings

import numpy as np
import numpy.typing as npt

ALPHA = {
    '32.5 R': 0,
    '42.5 N': 0,
    '42.5 R': 1,
    '52.5 N': 1,
    '52.5 R': 1,
    '32.5 N': -1,
}
ALPHA_BS = {
    '32.5 R': 700,
    '42.5 N': 700,
    '42.5 R': 600,
    '52.5 N': 600,
    '52.5 R': 600,
    '32.5 N': 800,
}
ALPHA_DS1 = {
    '32.5 R': 4,
    '42.5 N': 4,
    '42.5 R': 6,
    '52.5 N': 6,
    '52.5 R': 6,
    '32.5 N': 3,
}
ALPHA_DS2 = {
    '32.5 R': 0.012,
    '42.5 N': 0.012,
    '42.5 R': 0.012,
    '52.5 N': 0.012,
    '52.5 R': 0.012,
    '32.5 N': 0.013,
}


def _check_fcm(fcm: float) -> None:
    """Check if the mean compressive strength is in the range of applicability
    of the creep and shrinkage models.

    Defined in fib Model Code 2010 (2013), section 5.1.9.4.2.

    Args:
        fcm (float): The mean compressive strength of the concrete in MPa.

    Returns:
        Raises a ValueError if the mean compressive strength is outside of the
        specified range.
    """
    if fcm < 20 or fcm > 130:
        raise ValueError(
            'The specified mean compressive strength is '
            'outside of the range of applicability for the creep and '
            'shrinkage laws given in the fib Model Code 2010. The mean'
            ' compressive strength has to be within the range of '
            '20-130 MPa. Current compressive strength is'
            f': {fcm} MPa.'
        )


def _check_initial_stress(sigma: float, fcm: float) -> None:
    """Check if the initial compressive stress (based on load at t0) is within
    the range of applicability.

    Defined in fib Model Code 2010 (2013), section 5.1.9.4.2.

    Args:
        sigma (float): The compressive stress applied to the concrete at t0 in
            MPa.
        fcm (float): The mean compressive strength of the concrete in MPa. Note
            that its value is non-negative.

    Raises:
        Raises a warning if the initial compressive stress is greater than
            0.4*fcm.
        Raises a ValueError if the compressive stress is greater than 0.6*fcm.
    """
    if abs(sigma) > 0.6 * fcm:
        raise ValueError(
            'The stress level exceeds the range of application.'
            'Maximum allowable stress is 0.6*fcm. Current stress level '
            f'is {round(abs(sigma)/fcm, 3)}*fcm.'
        )
    if abs(sigma) > 0.4 * fcm:
        warnings.warn(
            'Initial stress is too high to consider the '
            'concrete as an aging linear visco-elastic material: '
            f'sigma = {round(abs(sigma)/fcm,3)}*fcm > 0.4*fcm. Nonlinear'
            ' creep calculations are performed according to subclause '
            '5.1.9.4.3 (d) of the fib Model Code 2010 to account for '
            'large compressive stresses.'
        )


def _check_age_at_loading(t0: float) -> None:
    """Check if the age of the concrete is greater than the minimum concrete
    age.

    Defined in fib Model Code 2010 (2013), section 5.1.9.4.2.

    Args:
        t0 (float): The age of the concrete in days at which the loading is
            applied.

    Raises:
        Raises a ValueError if the age of the concrete is too low.
    """
    if t0 < 1:
        raise ValueError(
            'The load is applied too soon to the concrete'
            ' in order to calculate the creep and shrinkage behaviour '
            'according to the fib Model Code 2010. The minimum age of the '
            'concrete is 1 day, whereas according to the input the load is'
            f' applied after {t0} day.'
        )


def _check_RH(rh: float) -> None:
    """Check if the given relative humidity is within the range of
    applicability.

    Defined in fib Model Code 2010 (2013), section 5.1.9.4.2.

    Args:
        rh (float): The relative humidity of the environment. Value can be
            provided as percentage (i.e. 40--100), or as ratio (i.e. 0.4--1).

    Raises:
        Raises a ValueError if the relative humidity is outside of the range of
            applicability.
    """
    if (rh < 0.4 or rh > 1) and (rh < 40 or rh > 100):
        raise ValueError(
            'The specified relative humidity is outside '
            'of the range of applicability to calculate the creep and '
            'shrinkage according to the fib Model Code 2010. The '
            'relative humidity has to be within the range of 0.4-1.0 or '
            f'40-100%. Currently rh={rh}.'
        )


def _check_env_temp(T: float) -> None:
    """Check if the given environmental temperature is within the range of
    applicability.

    Defined in fib Model Code 2010 (2013), section 5.1.9.4.2.

    Args:
        T (float): The environmental temperature in degrees Celcius.

    Raises:
        Raises a warning if the applied environmental temperature is outside of
            the range of applicability.
    """
    if T < 5 or T > 30:
        warnings.warn(
            'The given environmental temperature is outside'
            ' of the applicable range of 5-30 degrees Celcius for the'
            ' creep and shrinkage calculations according to the fib Model '
            f'Code 2010, T={T} degrees Celcius. Creep and shrinkage will'
            ' be calculated according to subclause 5.1.10 of the fib Model'
            ' Code 2010.'
        )


def t_T(
    t0: npt.ArrayLike, T_cur: npt.ArrayLike, dt: npt.ArrayLike = None
) -> np.ndarray:
    """Calculate the temperature corrected concrete age in days at t0.

    Defined in fib Model Code 2010 (2013). Eq. 5.1-85 (only for a single time
    value input, as required in Eq. 5.1-73).

    Args:
        t0 (np.typing.ArrayLike): The age of the concrete in days at which the
            loading is applied.
        T_cur (np.typing.ArrayLike): The temperature of the environment during
            curing in degrees Celcius.

    Keyword Args:
        dt (np.typing.ArrayLike): Number of days at which T_cur prevails.
            Required when providing a list for T_cur.

    Returns:
        np.ndarray: The temperature corrected age of the concrete in days at
        loading.
    """
    _check_age_at_loading(t0)
    if dt is None:
        dt = t0
    else:
        T_cur = np.asarray(T_cur)
        dt = np.asarray(dt)
        if T_cur.size != dt.size:
            raise ValueError('Dimensions of T_cur and dt do not match.')
        if np.sum(dt) != t0:
            raise ValueError(
                f'Curing time {np.sum(dt)} and time of loading {t0} do not'
                ' match.'
            )
    return np.sum(dt * np.exp(13.65 - (4000 / (273 + T_cur))))


def t0_adj(
    t_T: float,
    cem_class: t.Literal[
        '32.5 N',
        '32.5 R',
        '42.5 N',
        '42.5 R',
        '52.5 N',
        '52.5 R',
    ],
) -> float:
    """Calculate the modified age at loading (t0) to account for the effect of
    the type of cement and curing temperature on the degree of hydration and -
    in turn - on creep.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-73.

    Args:
        t_T (float): Temperature adjusted concrete age in days.
        cem_class (str): The cement strength class that is used. The choices
            are: '32.5 N', '32.5 R', '42.5 N', '42.5 R', '52.5 N', '52.5 R'.

    Returns:
        float: The temperature corrected age of the concrete in days at
        loading, accounting for the effect of the cement type. For slow
        hardening concrete, the creep coefficient is increased due to the lower
        modified age at loading.
    """
    return max(t_T * ((9 / (2 + t_T**1.2)) + 1) ** ALPHA[cem_class], 0.5)


def eps_cds0(
    fcm: float,
    cem_class: t.Literal[
        '32.5 N', '32.5 R', '42.5 N', '42.5 R', '52.5 N', '52.5 R'
    ],
) -> float:
    """Calculate the notional drying shrinkage.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-80.

    Args:
        fcm (float): The mean compressive strength of the concrete in MPa.
        cem_class (str): The cement strength class that is used. The choices
            are: '32.5 N', '32.5 R', '42.5 N', '42.5 R', '52.5 N', '52.5 R'.

    Returns:
        float: The notional drying shrinkage, no units.
    """
    return (
        (220 + 110 * ALPHA_DS1[cem_class])
        * np.exp(-ALPHA_DS2[cem_class] * fcm)
        * 1e-6
    )


def beta_ds(
    time: npt.ArrayLike, ts: float, notional_size: float
) -> np.ndarray:
    """Calculate the multiplication factor beta_ds.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-82.

    Args:
        time (numpy.typing.ArrayLike): The different times in days at which the
            shrinkage strain is determined.
        ts (float): Age of the concrete when exposed to the environment.
        notional_size (float): The notional size of the considered element in
            mm, defined as 2A/u.

    Returns:
        numpy.ndarray: Multiplication factor used for calculating the drying
        shrinkage as a function of time.
    """
    return np.sqrt((time - ts) / (0.035 * (notional_size) ** 2 + (time - ts)))


def beta_s1(fcm: float) -> float:
    """Calculate the correction factor beta_s1.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-83.

    Args:
        fcm (float): The mean compressive strength of the concrete in MPa.

    Returns:
        float: Multiplication factor used when calculating the drying
        shrinkage.
    """
    return min((35 / fcm) ** 0.1, 1.0)


def beta_RH(rh: float, beta_s1: float) -> float:
    """Calculate the multiplication factor beta_RH.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-81

    Args:
        rh (float): The relative humidity of the environment. Value can be
            provided as percentage (i.e. 40--100), or as ratio (i.e. 0.4--1).
        beta_s1 (float): Multiplication factor as calculated in the fib Model
            Code 2010 (2013), Eq. 5.1-83.

    Returns:
        float: Multiplication factor used when calculating the drying
        shrinkage.
    """
    _check_RH(rh)
    if rh > 1:
        rh = rh / 100

    if rh >= 0.99 * beta_s1:
        return 0.25
    if 0.4 * beta_s1 <= rh < 0.99 * beta_s1:
        return -1.55 * (1 - rh**3)
    raise ValueError(
        'The specified rh*beta_s1 is not in the range of application.'
    )


def eps_cds(
    eps_cds0: float,
    beta_ds: npt.ArrayLike,
    beta_rh: float,
) -> np.ndarray:
    """Calculate the drying shrinkage of the concrete element.

    Defined in fib Model Code 2010 (2013), Eqs. 5.1-77.

    Args:
        eps_cds0 (float):  The notional drying shrinkage, no units, as defined
            in fib Model Code 2010 (2013), Eq. 5.1-80.
        beta_ds (numpy.typing.ArrayLike): Multiplication factor used for
            calculating the drying shrinkage as a function of time, as defined
            in fib Model Code 2010 (2013), Eq. 5.1-82.
        beta_rh (float): Multiplication factor used when calculating the
            drying shrinkage.

    Returns:
        numpy.ndarray: The drying shrinkage strains for the given times, no
        units.
    """
    return eps_cds0 * beta_rh * beta_ds


def eps_cbs0(
    fcm: float,
    cem_class: t.Literal[
        '32.5 N', '32.5 R', '42.5 N', '42.5 R', '52.5 N', '52.5 R'
    ],
) -> float:
    """Calculate the notional basic shrinkage.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-78.

    Args:
        fcm (float): The mean compressive strength of the concrete in
            MPa.
        cem_class (str): The cement strength class that is used. The choices
            are: '32.5 N', '32.5 R', '42.5 N', '42.5 R', '52.5 N', '52.5 R'.

    Returns:
        float: The notional basic shrinkage, no units.
    """
    return -ALPHA_BS[cem_class] * ((0.1 * fcm) / (6 + 0.1 * fcm)) ** 2.5 * 1e-6


def beta_bs(time: npt.ArrayLike) -> np.ndarray:
    """Calculate multiplication factor beta_bs which is used to determine the
    basic shrinkage.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-79.

    Args:
        time (numpy.typing.ArrayLike): The different times in days at which the
            basic strain is determined.

    Returns:
        numpy.ndarray: Multiplication factor that is used to determine the
        basic shrinkage.
    """
    return 1 - np.exp(-0.2 * np.sqrt(time))


def eps_cbs(eps_cbs0: float, beta_bs: npt.ArrayLike) -> np.ndarray:
    """Calculate the basic shrinkage.

    Defined in fib Model Code 2010 (2013), Eqs. 5.1-76.

    Args:
        eps_cbs0 (float): Notional basic shrinkage, as defined in fib Model
            Code 2010 (2013), Eq. 5.1-78.
        beta_bs (numpy.typing.ArrayLike): Time function for basic shrinkage,
            as defined in fib Model Code 2010 (2013), Eq. 5.1-79.

    Returns:
        numpy.ndarray: The basic shrinkage strains for the given times.
    """
    return eps_cbs0 * beta_bs


def beta_bc_fcm(fcm: float) -> float:
    """Calculate multiplication factor that accounts for the effect of the
    compressive strength of the concrete to calculate the basic creep
    coefficient.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-65.

    Args:
        fcm (float): The mean compressive strength of the concrete in MPa.

    Returns:
        float: Multiplication factor beta_bc_fcm.
    """
    return 1.8 / fcm**0.7


def beta_bc_t(time: npt.ArrayLike, t0: float, t0_adj: float) -> np.ndarray:
    """Calculate multiplication factor that accounts for the effect of the age
    of the of the concrete to calculate the basic creep coefficient.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-66.

    Args:
        time (numpy.typing.ArrayLike): The different times in days at which the
            basic creep coefficient is determined.
        t0 (float): The age of the concrete in days at which the loading is
            applied.
        t0_adj (float): The temperature corrected age of the concrete when the
            loading is applied in days, as defined in fib Model Code 2010
            (2013). Eq. 5.1-85.

    Returns:
        numpy.ndarray: Multiplication factors beta_bc_t.
    """
    return np.log(((30 / t0_adj + 0.035) ** 2) * (time - t0) + 1)


def phi_bc(beta_bc_fcm: float, beta_bc_t: npt.ArrayLike) -> np.ndarray:
    """Calculate the basic creep coefficient.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-64.

    Args:
        beta_bc_fcm (float): Multiplication factor that accounts for the
            influence of the concrete strength of the creep behaviour, as
            defined in fib Model Code 2010 (2013), Eq. 5.1-65.
        beta_bc_t (numpy.typing.ArrayLike): Multiplication factor that
            accounts for the influence of the age of the concrete of the creep
            behaviour, as defined in fib Model Code 2010 (2013), Eq. 5.1-66.

    Returns:
        numpy.ndarray: The basic creep coefficient.
    """
    return beta_bc_fcm * beta_bc_t


def beta_dc_fcm(fcm: float) -> float:
    """Calculate multiplication factor that accounts for the effect of the
    strength of the concrete to calculate the drying creep coefficient.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-68.

    Args:
        fcm (float): The mean compressive strength of the concrete in MPa.

    Returns:
        float: Multiplication factor beta_dc_fcm.
    """
    return 412 / fcm**1.4


def beta_dc_RH(rh: float, notional_size: float) -> float:
    """Calculate multiplication factor that accounts for the effect of the
    relative humidity of the environment to calculate the drying creep
    coefficient.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-69.

    Args:
        rh (float): The relative humidity of the environment. Value can be
            provided as percentage (i.e. 40--100), or as ratio (i.e. 0.4--1).
        notional_size (float): The notional size of the considered element in
            mm, defined as 2A/u.

    Returns:
        float: Multiplication factor beta_RH.
    """
    _check_RH(rh)
    if rh > 1:
        rh = rh / 100
    return (1 - rh) / ((0.1 * notional_size / 100) ** (1 / 3))


def beta_dc_t0(t0_adj: float) -> float:
    """Calculate multiplication factor that accounts for the effect of the
    (temperature corrected) age of the concrete when loading is applied to
    calculate the drying creep coefficient.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-70.

    Args:
        t0_adj (float): The temperature corrected age of the concrete when the
            loading is applied in days, as defined in fib Model Code 2010
            (2013). Eq. 5.1-85.

    Returns:
        float: Multiplication factor beta_dc_t0.
    """
    return 1 / (0.1 + t0_adj**0.2)


def alpha_fcm(fcm: float) -> float:
    """Calculate multiplication factor that accounts for the effect of the
    strength of the concrete to calculate the drying creep coefficient.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-71d.

    Args:
        fcm (float): The mean compressive strength of the concrete in MPa.

    Returns:
        float: Multiplication factor alpha_fcm.
    """
    return np.sqrt(35 / fcm)


def beta_h(notional_size: float, alpha_fcm: float) -> float:
    """Calculate multiplication factor that accounts for the effect of the
    notional size of the to calculate the drying creep coefficient.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-71c

    Args:
        notional_size (float): The notional size of the considered element in
            mm, defined as 2A/h.
        alpha_fcm (float): Multiplication factor that accounts for the effect
            of the strength of the concrete on the drying creep coefficient.

    Returns:
        float: Multiplication factor beta_h.
    """
    return min(1.5 * notional_size + 250 * alpha_fcm, 1500 * alpha_fcm)


def gamma_t0(t0_adj: float) -> float:
    """Calculate exponent that accounts for the effect of the (temperature
    corrected) age of the concrete when loaded. Used to calculate the drying
    creep coefficient.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-71b.

    Args:
        t0_adj (float): The temperature corrected age of the concrete when the
            loading is applied in days, as defined in fib Model Code 2010
            (2013). Eq. 5.1-85.

    Returns:
        float: Exponent gamma_t0.
    """
    return 1 / (2.3 + 3.5 / np.sqrt(t0_adj))


def beta_dc_t(
    time: npt.ArrayLike, t0: float, beta_h: float, gamma_t0: float
) -> np.ndarray:
    """Calculate multiplication factor that accounts for the different
    considered values of time. Used to calculate the drying creep coefficient.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-71a.

    Args:
        time (numpy.typing.ArrayLike): The different times in days at which the
            drying creep coefficient is determined.
        t0 (float): The age of the concrete concrete when the loading is
            applied in days.
        beta_h (float): Multiplication factor that accounts for the effect of
            the notional size, as calculated by Eq. 5.1-71c
        gamma_t0 (float): Exponent that accounts for the effect of the
            (temperature corrected) age of the concrete when loaded, as
            calculated by Eq. 5.1-71b.

    Returns:
        numpy.ndarray: Multiplcation factor beta_dc_t for the considered values
        of time.
    """
    return ((time - t0) / (beta_h + (time - t0))) ** gamma_t0


def phi_dc(
    beta_dc_fcm: float,
    beta_dc_RH: float,
    beta_dc_t0: float,
    beta_dc_t: npt.ArrayLike,
) -> np.ndarray:
    """Calculate drying creep coefficient.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-67.

    Args:
        beta_dc_fcm (float): Multiplication factor that accounts for the
            effect of the strength of the concrete, as calculated by Eq.
            5.1-68.
        beta_dc_RH (float): Multiplication factor that accounts for the effect
            of the relative humidity of the environment, as calculated by Eq.
            5.1-69.
        beta_dc_t0 (float): multiplication factor that accounts for the effect
            of the (temperature corrected) age of the concrete when loading is
            applied, as calculated by Eq. 5.1-70.
        beta_dc_t (numpy.typing.ArrayLike): multiplication factor that
            accounts for the different considered values of time, as calculated
            by Eq. 5.1-71a.

    Returns:
        numpy.ndarray: Drying creep coeffcient.
    """
    return beta_dc_fcm * beta_dc_RH * beta_dc_t0 * beta_dc_t


def k_sigma(sigma: float, fcm: float) -> float:
    """Calculate the ratio between the applied stress and the mean concrete
    compressive strength.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-74.

    Args:
        sigma (float): The compressive stress applied to the concrete at ts in
            MPa.
        fcm (float): The mean compressive strength of the concrete in MPa.

    Returns:
        float: Absolute value of the ratio between the stress in the concrete
        and the mean concrete strength.
    """
    _check_initial_stress(sigma, fcm)
    return abs(sigma / fcm)


def phi(
    phi_bc: npt.ArrayLike, phi_dc: npt.ArrayLike, sigma: float, fcm: float
) -> np.ndarray:
    """Calculate the creep coefficient distinguishing between linear and
    non-linear creep for compressive stresses sigma <= 0.4fcm, and 0.4fcm <
    sigma <= 0.6fcm, respectively.

    Defined in fib Model Code 2010, Eqs. 5.1-63 and 5.1-74.

    Args:
        phi_bc (numpy.typing.ArrayLike): Basic creep coefficient, as defined
            in fib Model Code 2010 (2013), Eq. 5.1-64.
        phi_dc (numpy.typing.ArrayLike): Drying creep coefficient, as defined
            in fib Model Code 2010 (2013), Eq. 5.1-67.
        sigma (float): The compressive stress applied to the concrete at ts in
            MPa.
        fcm (float): The mean compressive strength of the concrete in MPa.

    Returns:
        numpy.ndarray: The creep coefficient.
    """
    # Calculate the creep coefficient (phi) (see Eq. 5.1-63)
    _phi = phi_bc + phi_dc
    # Include the effect of high stress if needed (see Eq. 5.1-74 of [1]):
    _k_sigma = k_sigma(sigma, fcm)
    if 0.4 < _k_sigma <= 0.6:
        return np.exp(1.5 * (_k_sigma - 0.4)) * _phi
    return _phi


def calc_J(E_ci_t0: float, phi: npt.ArrayLike, E_ci: float) -> np.ndarray:
    """Calculate the creep compliance function.

    Defined in fib Model Code 2010, Eq. 5.1-61.

    Args:
        E_ci_t0 (float): Modulus of elasticity at time of loading t0, as
            defined in fib Model Code 2010 (2013), Eq. 5.1-56.
        phi (numpy.typing.ArrayLike): Creep coefficient, as defined in fib
            Model Code 2010 (2013), Eq. 5.1-63.
        E_ci (float) Modulus of elasticity of the concrete at 28 days, as
            defined in fib Model Code 2010 (2013), Eq. 5.1-21.

    Returns:
        numpy.ndarray: The creep compliance function.
    """
    return (1 / E_ci_t0) + (phi / E_ci)
