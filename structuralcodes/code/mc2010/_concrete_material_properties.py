"""A collection of material properties for concrete"""
import math


def fcm(fck: float, delta_f: float = 8.0) -> float:
    """Compute the mean concrete compressive strength from the characteristic
    strength.

    Args:
        fck (float): The characteristic compressive strength.

    Kwargs:
        delta_f (float): The difference between the mean and the
        characteristic strength.
    """
    return fck + delta_f


def fctm(fck: float) -> float:
    """Compute the mean concrete tensile strength from the characteristic
    compressive strength.
    """
    if fck <= 50:
        return 0.3 * math.pow(fck, 2 / 3)
    return 2.12 * math.log(1 + 0.1 * fcm(fck))
