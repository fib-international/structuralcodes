"""Bilinear compression constitutive law."""

from __future__ import annotations  # To have clean hints of ArrayLike in docs

import typing as t

import numpy as np
from numpy.typing import ArrayLike

from ...core.base import ConstitutiveLaw


class BilinearCompression(ConstitutiveLaw):
    """Class for Bilinear Elastic-PerfectlyPlastic Constitutive Law for
    Concrete (only compression behavior).
    """

    __materials__: t.Tuple[str] = ('concrete',)

    def __init__(
        self,
        fc: float,
        eps_c: float,
        eps_cu: t.Optional[float] = None,
        name: t.Optional[str] = None,
    ) -> None:
        """Initialize a BilinearCompression Material.

        Arguments:
            fc (float): Compressive strength (negative number).
            eps_c (float): Strain at compressive strength (pure number).

        Keyword Arguments:
            eps_cu (float): Ultimate strain (pure number).
            name (str): A descriptive name for the constitutive law.
        """
        name = name if name is not None else 'BilinearCompressionLaw'
        super().__init__(name=name)
        self._fc = -abs(fc)
        self._eps_c = -abs(eps_c)
        self._eps_cu = -abs(eps_cu)
        self._E = self._fc / self._eps_c

    def get_stress(
        self, eps: t.Union[float, ArrayLike]
    ) -> t.Union[float, ArrayLike]:
        """Return the stress given strain."""
        eps = eps if np.isscalar(eps) else np.atleast_1d(eps)
        # Preprocess eps array in order
        eps = self.preprocess_strains_with_limits(eps=eps)
        # Compute stress
        # If it is a scalar
        if np.isscalar(eps):
            if eps > 0 or eps < self._eps_cu:
               return 0
            return max(self._E * eps, self._fc)
        
        # If it is an array
        sig = self._E * eps
        sig[sig < self._fc] = self._fc
        sig[eps > 0] = 0
        sig[eps < self._eps_cu] = 0
        return sig

    def get_tangent(
        self, eps: t.Union[float, ArrayLike]
    ) -> t.Union[float, ArrayLike]:
        """Return the tangent for given strain."""
        eps = eps if np.isscalar(eps) else np.atleast_1d(eps)
        # If it is a scalar
        if np.isscalar(eps):
            if self._eps_c < eps <= 0:
                return self._E
            return 0

        # If it is an array
        tangent = np.ones_like(eps) * self._E
        tangent[eps >= 0] = 0
        tangent[eps < self._eps_c] = 0

        return tangent

    def __marin__(
        self, strain: t.Tuple[float, float]
    ) -> t.Tuple[t.List[t.Tuple], t.List[t.Tuple]]:
        """Returns coefficients and strain limits for Marin integration in a
        simply formatted way.

        Arguments:
            strain (float, float): Tuple defining the strain profile: eps =
                strain[0] + strain[1]*y.

        Example:
            [(0, -0.002), (-0.002, -0.003)]
            [(a0, a1, a2), (a0)]
        """
        strains = []
        coeff = []
        if strain[1] == 0:
            # Uniform strain equal to strain[0]
            # understand in which branch we are
            strain[0] = self.preprocess_strains_with_limits(strain[0])
            if strain[0] > 0:
                # We are in tensile branch
                strains = None
                coeff.append((0.0,))
            elif strain[0] > self._eps_c:
                # We are in the linear branch
                strains = None
                a0 = self._E * strain[0]
                a1 = self._E * strain[1]
                coeff.append((a0, a1))
            elif strain[0] >= self._eps_cu:
                # We are in the constant branch
                strains = None
                coeff.append((self._fc,))
            else:
                # We are in a branch of non-resisting concrete
                # Too much compression
                strains = None
                coeff.append((0.0,))
        else:
            # linear part
            strains.append((self._eps_c, 0))
            a0 = self._E * strain[0]
            a1 = self._E * strain[1]
            coeff.append((a0, a1))
            # Constant part
            strains.append((self._eps_cu, self._eps_c))
            coeff.append((self._fc,))
        return strains, coeff

    def __marin_tangent__(
        self, strain: t.Tuple[float, float]
    ) -> t.Tuple[t.List[t.Tuple], t.List[t.Tuple]]:
        """Returns coefficients and strain limits for Marin integration of
        tangent in a simply formatted way.

        Arguments:
            strain (float, float): Tuple defining the strain profile: eps =
                strain[0] + strain[1]*y.

        Example:
            [(0, -0.002), (-0.002, -0.003)]
            [(a0, a1, a2), (a0)]
        """
        strains = []
        coeff = []
        if strain[1] == 0:
            # Uniform strain equal to strain[0]
            # understand in which branch we are
            strain[0] = self.preprocess_strains_with_limits(strain[0])
            if strain[0] > 0:
                # We are in tensile branch
                strains = None
                coeff.append((0.0,))
            elif strain[0] > self._eps_c:
                # We are in the linear branch
                strains = None
                a0 = self._E
                coeff.append((a0,))
            else:
                # We are in the constant branch or
                # We are in a branch of non-resisting concrete
                # Too much compression
                strains = None
                coeff.append((0.0,))
        else:
            # linear part
            strains.append((self._eps_c, 0))
            a0 = self._E
            coeff.append((a0,))
            # Constant part
            strains.append((self._eps_cu, self._eps_c))
            coeff.append((0.0,))
        return strains, coeff

    def get_ultimate_strain(
        self, yielding: bool = False
    ) -> t.Tuple[float, float]:
        """Return the ultimate strain (negative and positive)."""
        if yielding:
            return (self._eps_c, 100)
        return (self._eps_cu, 100)
