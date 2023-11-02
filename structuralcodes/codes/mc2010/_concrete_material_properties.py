"""A collection of material properties for concrete"""
import math
import numpy as np

# Values from Table 5.1-6.
_alpha_E = {
    "basalt": 1.2,
    "quartzite": 1.0,
    "limestone": 0.9,
    "sandstone": 0.7,
}
# Values for normal strength concrete, from Table 5.1-9.
_s = {
    "32.5 R": 0.25,
    "42.5 N": 0.25,
    "42.5 R": 0.2,
    "52.5 N": 0.2,
    "52.5 R": 0.2,
    "32.5 N": 0.38,
}


def fcm(fck: float, delta_f: float = 8.0) -> float:
    """Compute the mean concrete compressive strength from the characteristic
    strength.

    fib Model Code 2010, Eq. (5.1-1)

    Args:
        fck (float): The characteristic compressive strength in MPa.

    Keyword Args:
        delta_f (float): The difference between the mean and the
        characteristic strength.

    Returns:
        float: The mean compressive strength in MPa.
    """
    return abs(fck) + abs(delta_f)


def fctm(fck: float) -> float:
    """Compute the mean concrete tensile strength from the characteristic
    compressive strength.

    fib Model Code 2010, Eqs. (5.1-3a) and (5.1-3b)

    Args:
        fck (float): The characteristic compressive strength in MPa.

    Returns:
        float: The mean tensile strength in MPa.
    """
    if abs(fck) <= 50:
        return 0.3 * math.pow(abs(fck), 2 / 3)
    return 2.12 * math.log(1 + 0.1 * fcm(fck))


def fctkmin(_fctm: float) -> float:
    """Compute the lower bound value of the characteristic tensile strength
    from the mean tensile strength.

    fib Model Code 2010, Eq. (5.1-4)

    Args:
        _fctm (float): The mean tensile strength in MPa.

    Returns:
        float: Lower bound of the characteristic tensile strength in MPa.
    """
    return 0.7 * _fctm


def fctkmax(_fctm: float) -> float:
    """Compute the upper bound value of the characteristic tensile strength
    from the mean tensile strength.

    fib Model Code 2010, Eq. (5.1-5)

    Args:
        _fctm (float): The mean tensile strength in MPa.

    Returns:
        float: Upper bound of the characteristic tensile strength in MPa.
    """
    return 1.3 * _fctm


def Gf(fck: float) -> float:
    """Compute tensile fracture energy from characteristic compressive
    strength.

    fib Model Code 2010, Eq. (5.1-9)

    Args:
        fck (float): The characteristic compressive strength in MPa.

    Returns:
        float: The tensile fracture energy in N/m.
    """
    return 73 * fcm(fck) ** 0.18


def E_ci(
    _fcm: float, agg_type: str = 'quartzite', EC0: float = 21500
) -> float:
    """Calculate the modulus of elasticity for normal weight concrete
        at 28 days.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-21.

    Args:
        _fcm (float): The mean value of the compressive strength of the
            concrete in MPa.

    Keyword args:
        agg_type (str): Type of coarse grain aggregate used in the concrete.
            Choices are: 'basalt', 'quartzite', 'limestone', 'sandstone'.
        EC0 (float): Initial value of modulus of elasticity in MPa.

    Returns:
        float: The modulus of elasticity for normal weight concrete
            at 28 days in MPa.
    """

    return EC0 * _alpha_E[agg_type.lower()] * (_fcm / 10) ** (1 / 3)


def beta_cc(time: float, _fcm: float, cem_class: str) -> float:
    """Calculate multiplication factor beta_cc, used to determine the
        modulus of elasticity at an arbitrary time.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-51.

    Args:
        time (float): The time in days at which the modulus of
            elasticity is to be determined.
        _fcm The mean compressive strength of the concrete in
            MPa.
        cem_class (str): The cement strength class that is used.
            The choices are:
                '32.5 N',
                '32.5 R', '42.5 N',
                '42.5 R', '52.5 N', '52.5 R'.

    Returns:
        float: Multiplication factor beta_cc.
    """
    if _fcm > 60:
        return np.exp(0.2 * (1 - np.sqrt(28 / time)))
    return np.exp(_s[cem_class.upper()] * (1 - np.sqrt(28 / time)))


def beta_e(_beta_cc: float) -> float:
    """Calculate multiplication factor beta_e, used to determine the
        modulus of elasticity at an arbitrary time.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-57.

    Args:
        _beta_cc (float): multiplication factor as defined in the fib Model
            Code 2010 (2013), Eq. 5.1-51.

    Returns:
        float: Multiplication factor beta_e.
    """
    return np.sqrt(_beta_cc)


def E_ci_t(_beta_e: float, _E_ci: float) -> float:
    """Calculate the modulus of elasticity for normal weight concrete
        at time 'time' (not 28 days).

    Defined in fib Model Code 2010 (2013), Eq. 5.1-56.

    Args:
        _beta_e (float): Multiplication factor to determine the modulus of
            elasticity at an arbitrary time, as defined in fib Model Code 2010
            (2013), Eq. 5.1-51.
        _E_ci (float): Modulus of elasticity of normal weight concrete at 28
            days, as defined in fib Model Code 2010 (2013), Eq. 5.1-21.

    Returns:
        float: The modulus of elasticity for normal weight concrete at
            time 'time' (not 28 days) in MPa.
    """
    return _beta_e * _E_ci
