"""Functions from Section 5 of FprEN 1992-1-1:2022"""

import math

from structuralcodes.codes import mc2010
#import structuralcodes.codes.mc2010 as mc2010


def fcm(fck: float, delta_f: float = 8.0) -> float:
    """Determines the mean strength of concrete from its characteristic
    value"""
    return mc2010.fcm(fck, delta_f)


def fctm(fck: float) -> float:
    """Compute the mean concrete tensile strength from the characteristic
    compressive strength

    FprEN 1992-1-1, Table 5.1""""""

    Args:
        fck (float): The characteristic compressive strength in MPa.

    Returns:
        float: The mean tensile strength in MPa.
    """
    if abs(fck) <= 50:
        return 0.3 * math.pow(abs(fck), 2 / 3)
    return 1.1 * math.pow(abs(fck), 1 / 3)