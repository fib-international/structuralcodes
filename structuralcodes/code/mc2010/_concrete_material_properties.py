"""A collection of material properties for concrete"""
import math


def fcm(fck: float, delta_f: float = 8.0) -> float:
    """Compute the mean concrete compressive strength from the characteristic
    strength.

    fib Model Code 2010, 5.1-1

    Args:
        fck (float): The characteristic compressive strength in MPa.

    Kwargs:
        delta_f (float): The difference between the mean and the
        characteristic strength.

    Returns:
        float: The mean compressive strength in MPa.
    """
    return abs(fck) + abs(delta_f)


def fctm(fck: float) -> float:
    """Compute the mean concrete tensile strength from the characteristic
    compressive strength.

    fib Model Code 2010, 5.1-3a, 5.1-3b

    Args:
        fck (float): The characteristic compressive strength in MPa.

    Returns:
        float: The mean tensile strength in MPa.
    """
    if abs(fck) <= 50:
        return 0.3 * math.pow(abs(fck), 2 / 3)
    return 2.12 * math.log(1 + 0.1 * fcm(fck))


def fctkmin(fck: float) -> float:
    """Compute the lower bound value of the characteristic tensile strength
    from the characteristic compressive strength.

    fib Model Code 2010, 5.1-4

    Args:
        fck (float): The characteristic compressive strength in MPa.

    Returns:
        float: Lower bound of the characteristic tensile strength in MPa.
    """
    return 0.7 * fctm(fck)


def fctkmax(fck: float) -> float:
    """Compute the upper bound value of the characteristic tensile strength
    from the characteristic compressive strength.

    fib Model Code 2010, 5.1-5

    Args:
        fck (float): The characteristic compressive strength in MPa.

    Returns:
        float: Upper bound of the characteristic tensile strength in MPa.
    """
    return 1.3 * fctm(fck)


def Gf(fck: float) -> float:
    """Compute tensile fracture energy from characteristic compressive
    strength.

    fib Model Code 2010, 5.1-9

    Args:
        fck (float): The characteristic compressive strength in MPa.

    Returns:
        float: The tensile fracture energy in N/m.
    """
    return 73 * fcm(fck) ** 0.18


def fcd(
    fck: float,
    gammaC: float = 1.5,
    alfaC: float = 0.85,
    existing: bool = False,
    FC: float = 1.0,
) -> float:
    """Compute fcd ... to be completed"""
    if not existing:
        return fck * alfaC / gammaC

    return fcm(fck) / FC


# For Eci: for existing is fcm/10, for new is (fck+deltaf/10)
# -> add a flag to concrete saying if existing?