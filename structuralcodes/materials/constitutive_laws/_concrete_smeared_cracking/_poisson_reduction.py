"""Reduction of Poisson ratio in concrete smeared cracking."""

from __future__ import annotations  # To have clean hints of ArrayLike in docs

import typing as t
from abc import ABC, abstractmethod

import numpy as np
from numpy.typing import ArrayLike


class PoissonReduction(ABC):
    """An abstract base class for modelling reduction of the Poisson ratio due
    to cracking.
    """

    @abstractmethod
    def nu(
        self,
        cracked: bool,
        eps: ArrayLike,
        eps_p: ArrayLike,
        eps_pf: ArrayLike,
    ):
        raise NotImplementedError

    @abstractmethod
    def poisson_matrix(
        self,
        cracked: bool,
        eps: ArrayLike,
        eps_p: ArrayLike,
        eps_pf: ArrayLike,
    ):
        raise NotImplementedError


class ConstantPoissonReduction(PoissonReduction):
    """A material with a constant reduction in Poisson ratio due to
    cracking.
    """

    _initial_nu: float
    _cracked_nu: float
    _poisson_matrix_cracked: t.Optional[np.ndarray]
    _poisson_matrix_uncracked: t.Optional[np.ndarray]

    def __init__(self, initial_nu: float, cracked_nu: float = 0.0):
        """Initialize the model."""
        self._initial_nu = initial_nu
        self._cracked_nu = cracked_nu
        self._poisson_matrix_cracked = None
        self._poisson_matrix_uncracked = None

    def nu(self, cracked: bool, *args, **kwargs):
        """Return the Poisson ratio."""
        del args, kwargs
        return self._initial_nu if not cracked else self._cracked_nu

    def poisson_matrix(self, cracked: bool, *args, **kwargs) -> np.ndarray:
        """Return the 2x2 Poisson matrix."""
        del args, kwargs

        # Return Poisson matrix
        if cracked:
            if self._poisson_matrix_cracked is None:
                self._poisson_matrix_cracked = np.eye(2)
            return self._poisson_matrix_cracked

        if self._poisson_matrix_uncracked is None:
            self._poisson_matrix_uncracked = (
                1 / (1 - self._initial_nu**2)
            ) * np.array(
                [
                    [1, self._initial_nu],
                    [self._initial_nu, 1],
                ]
            )
        return self._poisson_matrix_uncracked
