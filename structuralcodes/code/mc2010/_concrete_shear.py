"""A collection of shear formulas for concrete"""
import typing as t
import math


def vrd(fck: float, z: float, bw: float, gamma_c: float, asw: float, sw: float, fywd: float, theta: float) -> float:
    """Compute the shear resistance of a web or slab.

    fib Model Code 2010, Eq. (7.3-11)

    Args:
        vrdc (float): Design shear resistance of concrete.
        vrds (float): Design shear resistance of shear reinforcement.

    Returns:
        float: Design shear resistance
    """

    return abs(vrdc(fck, z, bw, gamma_c)) + abs(vrds(asw, sw, z, bw, fywd, theta))


def vrdc(fck: float, z: float, bw: float, gamma_c: float = 1.5) -> float:
    """The design shear resistance of a web or a slab without
    shear reinforcement.

    fib Model Code 2010, Eq. (7.3-17)

    Args:
        vck (float): The characteristic compressive strength in MPa.
        z (float): the effective shear depth.
        gamma_c: Material factor.
        bw:

    Returns:
        float: Design shear resistance without shear reinforcement
    """
    kv = 180/(1000+1.25*z)
    if fck**0.5 > 8:
        fsqr = 8
    else:
        fsqr = fck**0.5

    return (kv*fsqr*z*bw)/gamma_c

def vrds(asw: float, sw: float, z: float, fywd: float, theta: float, alpha: t.optional[float] = 0.0) -> float:
    """The design shear resistance provided by stirrups

    fib Model Code 2010, Eq. (7.3-29)

    Args:
        asw (float):
        sw (float):
        gamma_c:
        bw:

    Returns:
        float:
    """

    if alpha == 0.0:
        return (asw/sw)*z*fywd*math.cot(theta)
    else:
        return (asw/sw)*z*fywd*(cot(theta) + cot(alpha)) * sin(alpha)
