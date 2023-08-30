"""Utilities for sections"""
from dataclasses import dataclass
from structuralcodes.core.base import ConstitutiveLaw


@dataclass(frozen=True)
class Fiber():
    """Class for a generic fiber
    x (float): position in local x axis
    y (float): position in local y axis
    A (float): area
    mat (ConstitutiveLaw): the constitutive law for the fiber"""

    x: float
    y: float
    A: float
    mat: ConstitutiveLaw
