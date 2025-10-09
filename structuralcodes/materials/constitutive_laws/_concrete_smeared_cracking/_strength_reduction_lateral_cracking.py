"""Strength reduction due to lateral cracking for concrete smeared cracking."""

from __future__ import annotations  # To have clean hints of ArrayLike in docs

from abc import ABC, abstractmethod

from numpy.typing import ArrayLike


class StrengthReductionLateralCracking(ABC):
    """Base class for models for accounting for strength reduction due to
    lateral cracking.
    """

    @abstractmethod
    def reduction(self, strains: ArrayLike) -> float:
        """Calculate the strength reduction based on a set of strains."""
        raise NotImplementedError


class GeneralVecchioCollins(StrengthReductionLateralCracking):
    """General model for accounting for strength reduction due to lateral
    cracking based on the work of Vecchio and Collins.
    """

    _c_1: float
    _c_2: float
    _minimum_value: float

    def __init__(self, c_1: float, c_2: float, minimum_value: float = 0.0):
        """Initialize the model."""
        self._c_1 = c_1
        self._c_2 = c_2
        self._minimum_value = minimum_value

    @property
    def c_1(self) -> float:
        """Return the coefficient c_1."""
        return self._c_1

    @property
    def c_2(self) -> float:
        """Return the coefficient c_2."""
        return self._c_2

    @property
    def minimum_value(self) -> float:
        """Return the minimum value."""
        return self._minimum_value

    def reduction(self, strains: ArrayLike) -> float:
        """Calculate the strength reduction based on a set of principal
        strains.

        beta = 1 / (c_1 + c_2 * eps_1)
        """
        beta = 1 / (self.c_1 + self.c_2 * max(strains))

        return max(min(beta, 1), self.minimum_value)
