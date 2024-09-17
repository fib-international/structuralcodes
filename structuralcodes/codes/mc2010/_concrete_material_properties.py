"""A collection of material properties for concrete."""

from __future__ import annotations  # To have clean hints of ArrayLike in docs

import math
import typing as t

import numpy as np
import numpy.typing as npt

# Values from Table 5.1-6.
ALPHA_E = {
    'basalt': 1.2,
    'quartzite': 1.0,
    'limestone': 0.9,
    'sandstone': 0.7,
}
# Values for normal strength concrete, from Table 5.1-9.
S_CEM = {
    '32.5 R': 0.25,
    '42.5 N': 0.25,
    '42.5 R': 0.2,
    '52.5 N': 0.2,
    '52.5 R': 0.2,
    '32.5 N': 0.38,
}


def fcm(fck: float, delta_f: float = 8.0) -> float:
    """Compute the mean concrete compressive strength from the characteristic
    strength.

    fib Model Code 2010, Eq. (5.1-1).

    Args:
        fck (float): The characteristic compressive strength in MPa.

    Keyword Args:
        delta_f (float): The difference between the mean and the characteristic
            strength.

    Returns:
        float: The mean compressive strength in MPa.
    """
    return abs(fck) + abs(delta_f)


def fctm(fck: float) -> float:
    """Compute the mean concrete tensile strength from the characteristic
    compressive strength.

    fib Model Code 2010, Eqs. (5.1-3a) and (5.1-3b).

    Args:
        fck (float): The characteristic compressive strength in MPa.

    Returns:
        float: The mean tensile strength in MPa.
    """
    if abs(fck) <= 50:
        return 0.3 * math.pow(abs(fck), 2 / 3)
    return 2.12 * math.log(1 + 0.1 * fcm(fck))


def fctkmin(fctm: float) -> float:
    """Compute the lower bound value of the characteristic tensile strength
    from the mean tensile strength.

    fib Model Code 2010, Eq. (5.1-4).

    Args:
        fctm (float): The mean tensile strength in MPa.

    Returns:
        float: Lower bound of the characteristic tensile strength in MPa.
    """
    return 0.7 * fctm


def fctkmax(fctm: float) -> float:
    """Compute the upper bound value of the characteristic tensile strength
    from the mean tensile strength.

    fib Model Code 2010, Eq. (5.1-5).

    Args:
        fctm (float): The mean tensile strength in MPa.

    Returns:
        float: Upper bound of the characteristic tensile strength in MPa.
    """
    return 1.3 * fctm


def Gf(fck: float) -> float:
    """Compute tensile fracture energy from characteristic compressive
    strength.

    fib Model Code 2010, Eq. (5.1-9).

    Args:
        fck (float): The characteristic compressive strength in MPa.

    Returns:
        float: The tensile fracture energy in N/m.
    """
    return 73 * fcm(fck) ** 0.18


def Eci(
    fcm: float,
    agg_type: t.Literal[
        'basalt', 'quartzite', 'limestone', 'sandstone'
    ] = 'quartzite',
    EC0: float = 21500,
) -> float:
    """Calculate the modulus of elasticity for normal weight concrete at 28
    days.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-21.

    Args:
        fcm (float): The mean value of the compressive strength of the
            concrete in MPa.

    Keyword Args:
        agg_type (str): Type of coarse grain aggregate used in the concrete.
            Choices are: 'basalt', 'quartzite', 'limestone', 'sandstone'.
        EC0 (float): Initial value of modulus of elasticity in MPa.

    Returns:
        float: The modulus of elasticity for normal weight concrete at 28 days
        in MPa.
    """
    return EC0 * ALPHA_E[agg_type.lower()] * (fcm / 10) ** (1 / 3)


def beta_cc(
    time: npt.ArrayLike,
    fcm: float,
    cem_class: t.Literal[
        '32.5 N', '32.5 R', '42.5 N', '42.5 R', '52.5 N', '52.5 R'
    ],
) -> np.ndarray:
    """Calculate multiplication factor beta_cc, used to determine the
    compressive strength at an arbitrary time.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-51.

    Args:
        time (numpy.typing.ArrayLike): The time in days at which the
            compressive strength is to be determined.
        fcm (float): The mean compressive strength of the concrete in MPa.
        cem_class (str): The cement strength class that is used. The choices
            are: '32.5 N', '32.5 R', '42.5 N', '42.5 R', '52.5 N', '52.5 R'.

    Returns:
        numpy.ndarray: Multiplication factor beta_cc.
    """
    if fcm > 60:
        return np.exp(0.2 * (1 - np.sqrt(28 / time)))
    return np.exp(S_CEM[cem_class.upper()] * (1 - np.sqrt(28 / time)))


def beta_e(beta_cc: npt.ArrayLike) -> np.ndarray:
    """Calculate multiplication factor beta_e, used to determine the modulus of
    elasticity at an arbitrary time.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-57.

    Args:
        beta_cc (numpy.typing.ArrayLike): Multiplication factor as defined in
            the fib Model Code 2010 (2013), Eq. 5.1-51.

    Returns:
        numpy.ndarray: Multiplication factor beta_e.
    """
    return np.sqrt(beta_cc)


def Eci_t(beta_e: npt.ArrayLike, Eci: float) -> np.ndarray:
    """Calculate the modulus of elasticity for normal weight concrete at time
    'time' (not 28 days).

    Defined in fib Model Code 2010 (2013), Eq. 5.1-56.

    Args:
        beta_e (numpy.typing.ArrayLike): Multiplication factor to determine
            the modulus of elasticity at an arbitrary time, as defined in fib
            Model Code 2010 (2013), Eq. 5.1-51.
        Eci (float): Modulus of elasticity of normal weight concrete at 28
            days, as defined in fib Model Code 2010 (2013), Eq. 5.1-21.

    Returns:
        numpy.ndarray: The modulus of elasticity for normal weight concrete at
        time 'time' (not 28 days) in MPa.
    """
    return beta_e * Eci


def fcd(fck: float, alpha_cc: float = 1.0, gamma_c: float = 1.5) -> float:
    """The design compressive strength of concrete.

    Defined in fib Model Code 2010 (2013), Eq. 7.2-11.

    Args:
        fck (float): The characteristic compressive strength in MPa.

    Keyword Args:
        alpha_cc (float): A factor for considering long-term effects on the
            strength, and effects that arise from the way the load is applied.
            Default value 1.0.
        gamma_c (float): The partial factor of concrete. Default value 1.5.

    Returns:
        float: The design compressive strength of concrete in MPa.
    """
    return abs(alpha_cc) * abs(fck) / abs(gamma_c)


def eps_c1(fck: float) -> float:
    """The strain at maximum compressive stress of concrete (fcm) for the
    Sargin constitutive law.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-26 and Table 5.1-8

    The value is computing using interpolation from Table 5.1-8

    Args:
        fck (float): The characteristic compressive strength of concrete in
            MPa.

    Returns:
        float: The strain at maximum compressive stress, absolute value, no
        unit.

    Raises:
        ValueError: if fck is less than 12 or greater than 120
    """
    grade = np.array(
        [12, 16, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 110, 120]
    )
    eps_c1 = np.array(
        [
            -1.9,
            -2.0,
            -2.1,
            -2.2,
            -2.3,
            -2.3,
            -2.4,
            -2.5,
            -2.6,
            -2.6,
            -2.7,
            -2.7,
            -2.8,
            -2.9,
            -3.0,
            -3.0,
            -3.0,
        ]
    )
    if fck < grade.min() or fck > grade.max():
        raise ValueError(
            f'fck must be between {grade.min()} MPa and {grade.max()} MPa.'
            ' fck = {fck} given.'
        )
    return np.interp(fck, grade, eps_c1) / 1000


def eps_clim(fck: float) -> float:
    """The ultimate strain for the Sargin constitutive law.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-26 and Table 5.1-8

    The value is computing using interpolation from Table 5.1-8

    Args:
        fck (float): The characteristic compressive strength of concrete in
            MPa.

    Returns:
        float: The ultimate strain, absolute value, no unit.
    """
    grade = np.array(
        [12, 16, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 110, 120]
    )
    eps_clim = np.array(
        [
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.5,
            -3.4,
            -3.4,
            -3.3,
            -3.2,
            -3.1,
            -3.0,
            -3.0,
            -3.0,
            -3.0,
        ]
    )
    if fck < grade.min() or fck > grade.max():
        raise ValueError(
            f'fck must be between {grade.min()} MPa and {grade.max()} MPa.'
            ' fck = {fck} given.'
        )
    return np.interp(fck, grade, eps_clim) / 1000


def k_sargin(fck: float) -> float:
    """The plasticity number k for Sargin constitutive Law.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-26 and Table 5.1-8

    Args:
        fck (float): The characteristic compressive strength of concrete in
            MPa.

    Returns:
        float: The plasticity number k, absolute value, no unit.
    """
    grade = np.array(
        [12, 16, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 110, 120]
    )
    k = np.array(
        [
            2.44,
            2.36,
            2.28,
            2.15,
            2.04,
            1.92,
            1.82,
            1.74,
            1.66,
            1.61,
            1.55,
            1.47,
            1.41,
            1.36,
            1.32,
            1.24,
            1.18,
        ]
    )
    if fck < grade.min() or fck > grade.max():
        raise ValueError(
            f'fck must be between {grade.min()} MPa and {grade.max()} MPa.'
            ' fck = {fck} given.'
        )
    return np.interp(fck, grade, k)


def eps_cu1(fck: float) -> float:
    """The nominal ultimate strain for the Sargin constitutive law.

    Defined in fib Model Code 2010 (2013), Table 7.2-1 and Eq. 7.2-10

    Args:
        fck (float): The characteristic compressive strength of concrete in
            MPa.

    Returns:
        float: The strain at maximum compressive stress, absolute value, no
        unit.
    """
    return eps_clim(fck)


def eps_c2(fck: float) -> float:
    """The strain at maximum compressive stress of concrete for the
    parabolic-rectangular law.

    Defined in fib Model Code 2010 (2013), Table 7.2-1

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

    Defined in fib Model Code 2010 (2013), Table 7.2-1

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

    Defined in fib Model Code 2010 (2013), Table 7.2-1

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

    Defined in fib Model Code 2010 (2013), Table 7.2-1

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

    Defined in fib Model Code 2010 (2013), Table 7.2-1

    Args:
        fck (float): The characteristic compressive strength of concrete in
            MPa.

    Returns:
        float: The ultimate strain, absolute value, no unit.
    """
    return eps_cu2(fck)
