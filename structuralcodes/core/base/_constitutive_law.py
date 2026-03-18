"""Abstract base class for constitutive laws."""

from __future__ import annotations  # To have clean hints of ArrayLike in docs

import abc
import typing as t

import numpy as np
from numpy.typing import ArrayLike


class ConstitutiveLaw(abc.ABC):
    """Abstract base class for constitutive laws."""

    __materials__: t.Tuple[str] = ()
    constitutive_law_counter: t.ClassVar[int] = 0

    def __init__(self, name: t.Optional[str] = None) -> None:
        """Initialize a ConstitutiveLaw object."""
        self.id = self.constitutive_law_counter
        self._name = name if name is not None else f'ConstitutiveLaw_{self.id}'
        self._increase_global_counter()

    @property
    def name(self):
        """Returns the name of the constitutive law."""
        return self._name

    @classmethod
    def _increase_global_counter(cls):
        cls.constitutive_law_counter += 1

    @abc.abstractmethod
    def get_stress(
        self, eps: t.Union[float, ArrayLike]
    ) -> t.Union[float, ArrayLike]:
        """Each constitutive law should provide a method to return the
        stress given the strain level.
        """

    @abc.abstractmethod
    def get_tangent(
        self, eps: t.Union[float, ArrayLike]
    ) -> t.Union[float, ArrayLike]:
        """Each constitutive law should provide a method to return the
        tangent at a given strain level.
        """

    @abc.abstractmethod
    def get_ultimate_strain(self) -> t.Tuple[float, float]:
        """Each constitutive law should provide a method to return the
        ultimate strain (negative and positive).
        """

    def preprocess_strains_with_limits(
        self, eps: t.Union[float, ArrayLike]
    ) -> t.Union[float, ArrayLike]:
        """Preprocess strain arrays setting those strains sufficiently
        near to ultimate strain limits to exactly ultimate strain limit.
        """
        eps = eps if np.isscalar(eps) else np.atleast_1d(eps)
        eps_min, eps_max = self.get_ultimate_strain()

        if np.isscalar(eps):
            if np.isclose(eps, eps_max, atol=1e-6):
                return eps_max
            if np.isclose(eps, eps_min, atol=1e-6):
                return eps_min
            return eps
        idxs = np.isclose(eps, np.zeros_like(eps) + eps_max, atol=1e-6)
        eps[idxs] = eps_max
        idxs = np.isclose(eps, np.zeros_like(eps) + eps_min, atol=1e-6)
        eps[idxs] = eps_min

        return eps

    def _discretize_law(self) -> ConstitutiveLaw:
        """Discretize the law as a piecewise linear function."""

        # Discretize the constitutive law in a "smart way"
        def find_x_lim(x, y):
            # Check if there are non-zero values for x > 0
            if np.any(y[0:] != 0):
                # Find the last non-zero index for x > 0
                non_zero_indices = np.nonzero(y[0:])[0]
                x_lim_index = 0 + non_zero_indices[-1]
                return x[x_lim_index]
            # All values are zero for x > 0
            return None

        eps_min, eps_max = self.get_ultimate_strain()
        eps_max = min(eps_max, 1)
        # Analise positive branch
        eps = np.linspace(0, eps_max, 10000)
        sig = self.get_stress(eps)
        sig[(sig < np.max(sig) * 1e-6)] = 0
        eps_lim = find_x_lim(eps, sig)
        # Now discretize the function in 10 steps for positive part
        eps_pos = (
            np.linspace(0, -eps_min, 1)
            if eps_lim is None
            else np.linspace(0, eps_lim, 10)
        )
        # Analise negative branch
        eps = np.linspace(0, eps_min, 10000)
        sig = -self.get_stress(eps)
        sig[(sig < np.max(sig) * 1e-6)] = 0
        eps_lim = find_x_lim(-eps, sig)
        # Now discretize the function in 10 steps for negative part
        eps_neg = (
            np.linspace(eps_min, 0, 1, endpoint=False)
            if eps_lim is None
            else np.linspace(-eps_lim, 0, 10, endpoint=False)
        )

        eps = np.concatenate((eps_neg, eps_pos))
        sig = self.get_stress(eps)
        from structuralcodes.materials.constitutive_laws import (  # noqa: PLC0415
            UserDefined,
        )

        return UserDefined(eps, sig)

    def __marin__(self, **kwargs):
        """Function for getting the strain limits and coefficients
        for marin integration.

        By default the law is discretized as a piecewise linear
        function. Then marin coefficients are computed based on this
        discretization.
        """
        piecewise_law = self._discretize_law()
        return piecewise_law.__marin__(**kwargs)

    def __marin_tangent__(self, **kwargs):
        """Function for getting the strain limits and coefficients
        for marin integration of tangent modulus.

        By default the law is discretized as a piecewise linear
        function. Then marin coefficients are computed based on this
        discretization.
        """
        piecewise_law = self._discretize_law()
        return piecewise_law.__marin_tangent__(**kwargs)

    def get_secant(
        self, eps: t.Union[float, ArrayLike]
    ) -> t.Union[float, ArrayLike]:
        """Method to return the secant at a given strain level."""
        # Adjust eps if it is not scalar
        eps = eps if np.isscalar(eps) else np.atleast_1d(eps)

        # Calculate secant for scalar eps
        if np.isscalar(eps):
            if eps != 0:
                sig = self.get_stress(eps)
                return sig / eps
            return self.get_tangent(eps)

        # Calculate secant for array eps
        secant = np.zeros_like(eps)
        strain_is_zero = eps == 0
        strain_is_nonzero = eps != 0
        secant[strain_is_zero] = self.get_tangent(eps[strain_is_zero])
        secant[strain_is_nonzero] = (
            self.get_stress(eps[strain_is_nonzero]) / eps[strain_is_nonzero]
        )
        return secant
