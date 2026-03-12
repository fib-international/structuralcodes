"""Concrete material properties according to Tab. 3.1."""

from __future__ import annotations  # To have clean hints of ArrayLike in docs

import math
import typing as t

import numpy as np
from numpy.typing import ArrayLike

from structuralcodes.codes import mc2010

S_TIME_DEVELOPMENT_DICT = {
    'R': 0.20,
    'N': 0.25,
    'S': 0.38,
}  # As defined in Eq. (3.2)


def fcm(fck: float, delta_f: float = 8) -> float:
    """The mean compressive strength of concrete.

    EN 1992-1-1:2004, Table 3.1.

    Args:
        fck (float): The characteristic compressive strength of concrete in
            MPa.

    Keyword Args:
        delta_f (float): The difference between the mean and the
            characteristic strength.

    Returns:
        float: The mean compressive strength in MPa.
    """
    return mc2010.fcm(fck=abs(fck), delta_f=abs(delta_f))


def fctm(fck: float) -> float:
    """The mean tensile strength of concrete.

    EN 1992-1-1: 2004, Table 3.1.

    Args:
        fck (float): The characteristic compressive strength of concrete in
            MPa.

    Returns:
        float: The mean tensile strength in MPa.
    """
    return mc2010.fctm(fck=abs(fck))


def fctk_5(fctm: float) -> float:
    """The 5% fractile of the tensile strength of concrete.

    EN 1992-1-1: 2004, Table 3.1.

    Args:
        fctm (float): The mean tensile strength of concrete in MPa.

    Returns:
        float: The 5% fractile of the tensile strength in MPa.
    """
    return mc2010.fctkmin(fctm=abs(fctm))


def fctk_95(fctm: float) -> float:
    """The 95% fractile of the tensile strength of concrete.

    EN 1992-1-1: 2004, Table 3.1.

    Args:
        fctm (float): The mean tensile strength of concrete in MPa.

    Returns:
        float: The 95% fractile of the tensile strength in MPa.
    """
    return mc2010.fctkmax(fctm=abs(fctm))


def Ecm(fcm: float) -> float:
    """The secant modulus of concrete.

    EN 1992-1-1:2004, Table 3.1.

    Args:
        fcm (float): The mean compressive strength of concrete in MPa.

    Returns:
        float: The secant modulus of concrete in MPa.
    """
    return 22000.0 * math.pow(abs(fcm) / 10, 0.3)


def eps_c1(fcm: float) -> float:
    """The strain at maximum compressive stress of concrete (fcm) for the
    Sargin constitutive law.

    EN 1992-1-1:2004, Table 3.1.

    Args:
        fcm (float): The mean compressive strength of concrete in MPa.

    Returns:
        float: The strain at maximum compressive stress, absolute value, no
        unit.
    """
    return min(0.7 * math.pow(abs(fcm), 0.31), 2.8) / 1000


def eps_cu1(fck: float) -> float:
    """The ultimate strain for the Sargin constitutive law.

    EN 1992-1-1:2004, Table 3.1.

    Args:
        fck (float): The characteristic compressive strength of concrete in
            MPa.

    Returns:
        float: The ultimate strain, absolute value, no unit.
    """
    fck = abs(fck)
    return (
        3.5 / 1000
        if fck < 50
        else (2.8 + 27 * ((98 - fcm(fck)) / 100) ** 4) / 1000
    )


def k_sargin(
    Ecm: float,
    fcm: float,
    eps_c1: float,
) -> float:
    """Computation of k parameter for Sargin constitutive Law.

    EN 1992-1-1:2004, Eq. (3.14).

    Args:
        Ecm (float): the mean elastic modulus of concrete in MPa.
        fcm (float): the mean compressive strength in MPa.
        eps_c1 (float): the strain corresponding to peak stress.
    """
    return 1.05 * Ecm * abs(eps_c1) / fcm


def eps_c2(fck: float) -> float:
    """The strain at maximum compressive stress of concrete for the
    parabolic-rectangular law.

    EN 1992-1-1:2004, Table 3.1.

    Args:
        fck (float): The characteristic compressive strength of concrete in
            MPa.

    Returns:
        float: The strain at maximum compressive stress, absolute value, no
        unit.
    """
    fck = abs(fck)
    return (
        2.0 / 1000 if fck <= 50 else (2.0 + 0.085 * (fck - 50) ** 0.53) / 1000
    )


def eps_cu2(fck: float) -> float:
    """The ultimate strain of the parabolic-rectangular law.

    EN 1992-1-1:2004, Table 3.1.

    Args:
        fck (float): The characteristic compressive strength of concrete in
            MPa.

    Returns:
        float: The ultimate strain, absolute value, no unit.
    """
    fck = abs(fck)
    return (
        3.5 / 1000
        if fck <= 50
        else (2.6 + 35 * ((90 - fck) / 100) ** 4) / 1000
    )


def n_parabolic_rectangular(fck: float) -> float:
    """The exponent in the parabolic-rectangular law.

    EN 1992-1-1:2004, Table 3.1.

    Args:
        fck (float): The characteristic compressive strength of concrete in
            MPa.

    Returns:
        float: The exponent n, absolute value, no unit.
    """
    fck = abs(fck)
    return 2.0 if fck <= 50 else (1.4 + 23.4 * ((90 - fck) / 100) ** 4)


def eps_c3(fck: float) -> float:
    """The strain at maximum compressive stress of the bi-linear law.

    EN 1992-1-1:2004, Table 3.1.

    Args:
        fck (float): The characteristic compressive strength of concrete in
            MPa.

    Returns:
        float: The strain at maximum compressive stress, absolute value, no
        unit.
    """
    fck = abs(fck)
    return 1.75 / 1000 if fck <= 50 else (1.75 + 0.55 * (fck - 50) / 40) / 1000


def eps_cu3(fck: float) -> float:
    """The ultimate strain of the bi-linear law.

    EN 1992-1-1:2004, Table 3.1.

    Args:
        fck (float): The characteristic compressive strength of concrete in
            MPa.

    Returns:
        float: The ultimate strain, absolute value, no unit.
    """
    return eps_cu2(fck)


def fcd(fck: float, alpha_cc: float, gamma_c: float) -> float:
    """The design compressive strength of concrete.

    EN 1992-1-1:2004, Eq. (3.15).

    Args:
        fck (float): The characteristic compressive strength in MPa.
        alpha_cc (float): A factor for considering long-term effects on the
            strength, and effects that arise from the way the load is applied.
        gamma_c (float): The partial factor of concrete.

    Returns:
        float: The design compressive strength of concrete in MPa
    """
    return abs(alpha_cc) * abs(fck) / abs(gamma_c)


def beta_cc(t: ArrayLike, s: float) -> ArrayLike:
    """The time development function for compressive strength of concrete.

    EN 1992-1-1:2004, Eq. (3.2).

    Args:
        t (ArrayLike): The time in days to evaluate the development function
            for.
        s (float): The scale factor in the exponent for the time development
            function. s = 0.20 for class R, 0.25 for class N, and 0.38 for
            class N.

    Returns:
        ArrayLike: The value of the time development function.
    """
    return np.exp(s * (1 - np.sqrt(28 / t)))


def beta_ct(t: ArrayLike, s: float) -> ArrayLike:
    """The time development function for tensile strength of concrete.

    EN 1992-1-1:2004, part of Eq. (3.4).

    Args:
        t (ArrayLike): The time in days to evaluate the development function
            for.
        s (float): The scale factor in the exponent for the time development
            function. s = 0.20 for class R, 0.25 for class N, and 0.38 for
            class N.

    Returns:
        ArrayLike: The value of the time development function.
    """
    if np.isscalar(t):
        beta = beta_cc(t, s)
        if t < 28:
            return beta
        return np.pow(beta, 2 / 3)
    t = np.atleast_1d(t)
    beta = beta_cc(t, s)
    beta[t >= 28] = np.pow(beta[t >= 28], 2 / 3)
    return beta


def beta_E(t: ArrayLike, s: float) -> ArrayLike:
    """The time development function for Young's modulus of concrete.

    EN 1992-1-1:2004, part of Eq. (3.5).

    Args:
        t (ArrayLike): The time in days to evaluate the development function
            for.
        s (float): The scale factor in the exponent for the time development
            function. s = 0.20 for class R, 0.25 for class N, and 0.38 for
            class N.

    Returns:
        ArrayLike: The value of the time development function.
    """
    return np.pow(beta_cc(t, s), 0.3)


def s_time_development(cement_class: t.Literal['S', 'N', 'R']) -> float:
    """Return the scale factor for the exponent for the time development
    function.

    EN 1992-1-1:2004, Eq. (3.2).

    Args:
        cement_class (str): The cement class, either 'S', 'N' or 'R'.

    Returns:
        float: The scale factor that depends on the cement type.

    Raises:
        ValueError: If an invalid cement class is provided.

    """
    cement_class = (
        cement_class.upper().strip() if cement_class is not None else ''
    )
    s = S_TIME_DEVELOPMENT_DICT.get(cement_class)

    if s is None:
        raise ValueError(
            (
                f'"{cement_class}" is not a valid cement class. '
                'Use either S, N or R.'
            )
        )
    return s


def fcm_time(fcm: float, beta_cc: ArrayLike) -> ArrayLike:
    """Calculate the compressive strength as function of time.

    EN 1992-1-1:2004, Eq. (3.1).

    Args:
        fcm (float): The reference value for the compressive strength.
        beta_cc (ArrayLike): The value(s) of the time development function.

    Returns:
        ArrayLike: The calculated value(s) of the compressive strength.

    Note:
        The value of beta_cc should be calculated with the function beta_cc.
    """
    return fcm * beta_cc


def fctm_time(fctm: float, beta_cc: ArrayLike, alpha: ArrayLike) -> ArrayLike:
    """Calculate the tensile strength as function of time.

    EN 1992-1-1:2004, Eq. (3.4).

    Args:
        fctm (float): The reference value for the tensile strength.
        beta_cc (ArrayLike): The value(s) of the time development function.
        alpha (ArrayLike): An exponent for the time development function. It
            should be set to 1 for t < 28, and 2/3 else.

    Returns:
        ArrayLike: The calculated value(s) of the tensile strength.

    Note:
        The value of beta_cc should be calculated with the function beta_cc.
        Alternatively, the time development function for the tensile strength
        could be calculated directly with the function beta_ct.
    """
    return np.pow(beta_cc, alpha) * fctm


def Ecm_time(fcm: float, fcm_time: ArrayLike, Ecm: float) -> ArrayLike:
    """Calculate the Young's modulus as function of time.

    EN 1992-1-1:2004, Eq. (3.5).

    Args:
        fcm (float): The reference value for the compressive strength.
        fcm_time (float): The value(s) of the compressive strength at the
            point(s) in time.
        Ecm (float): The reference value for the Young's modulus.

    Returns:
        ArrayLike: The calculated value(s) of the Young's modulus.

    Note:
        The value of fcm_time should be calculated with the function fcm_time.
        Alternatively, the time development function for the Young's modulus
        could be calculated directly with the function beta_E.
    """
    return np.pow(fcm_time / fcm, 0.3) * Ecm
