"""Functions from Section 5 of FprEN 1992-1-1:2022"""

import mc2010

import structuralcodes


def fcm(fck: float, delta_f: float = 8.0) -> float:
    """Determines the mean strength of concrete form its characteristic value"""
    return mc2010.fcm(fck, delta_f)
