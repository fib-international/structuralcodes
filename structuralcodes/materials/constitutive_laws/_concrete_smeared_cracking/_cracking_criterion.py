"""Cracking criteria for concrete smeared cracking."""

from __future__ import annotations  # To have clean hints of ArrayLike in docs

from abc import ABC, abstractmethod

from numpy.typing import ArrayLike


class CrackingCriterion(ABC):
    """An abstract based class for modelling a cracking criterion."""

    @abstractmethod
    def cracked(self, eps_p: ArrayLike) -> bool:
        """Check if the material is cracked based on principal strains."""
        raise NotImplementedError


class NoTension(CrackingCriterion):
    """Cracking criterion to represent a material with no tensile strength."""

    def cracked(self, eps_p: ArrayLike) -> bool:
        """Check if the material is cracked based on the principal strains."""
        return max(eps_p) > 0
