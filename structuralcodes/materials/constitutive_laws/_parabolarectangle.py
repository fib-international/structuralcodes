"""Collection of some standard constitutive laws."""

from __future__ import annotations  # To have clean hints of ArrayLike in docs

import typing as t

import numpy as np
from numpy.typing import ArrayLike

from ...core.base import ConstitutiveLaw


class ParabolaRectangle(ConstitutiveLaw):
    """Class for parabola rectangle constitutive law.

    The stresses and strains are assumed negative in compression and positive
    in tension.
    """

    __materials__: t.Tuple[str] = ('concrete',)

    def __init__(
        self,
        fc: float,
        eps_0: float = -0.002,
        eps_u: float = -0.0035,
        n: float = 2.0,
        name: t.Optional[str] = None,
    ) -> None:
        """Initialize a Parabola-Rectangle Material.

        Arguments:
            fc (float): The strength of concrete in compression.

        Keyword Arguments:
            eps_0 (float): Peak strain of concrete in compression. Default
                value = -0.002.
            eps_u (float): Ultimate strain of concrete in compression. Default
                value = -0.0035.
            n (float): Exponent for the pre-peak branch. Default value = 2.
            name (str): A name for the constitutive law.
        """
        name = name if name is not None else 'ParabolaRectangleLaw'
        super().__init__(name=name)
        self._fc = -abs(fc)
        self._eps_0 = -abs(eps_0)
        self._eps_u = -abs(eps_u)
        self._n = n

    def get_stress(self, eps: ArrayLike) -> t.Union[float, ArrayLike]:
        """Return the stress given strain."""
        eps = eps if np.isscalar(eps) else np.atleast_1d(eps)
        # Preprocess eps array in order
        eps = self.preprocess_strains_with_limits(eps=eps)
        # Compute stress
        # If it is a scalar
        if np.isscalar(eps):
            sig = 0
            if self._eps_0 <= eps <= 0:
                sig = self._fc * (1 - (1 - eps / self._eps_0) ** self._n)
            if self._eps_u <= eps < self._eps_0:
                sig = self._fc
            return sig
        # If it is an array
        sig = np.zeros_like(eps)
        # Parabolic branch
        sig[(eps <= 0) & (eps >= self._eps_0)] = self._fc * (
            1
            - (1 - (eps[(eps <= 0) & (eps >= self._eps_0)] / self._eps_0))
            ** self._n
        )
        # Rectangle branch
        sig[eps < self._eps_0] = self._fc
        # Zero elsewhere
        sig[eps < self._eps_u] = 0
        sig[eps > 0] = 0
        return sig

    def get_tangent(self, eps: ArrayLike) -> t.Union[float, ArrayLike]:
        """Return the tangent given strain."""
        eps = eps if np.isscalar(eps) else np.atleast_1d(eps)
        # If it is a scalar
        if np.isscalar(eps):
            tangent = 0
            if self._eps_0 <= eps <= 0:
                tangent = (
                    self._n
                    * self._fc
                    / self._eps_0
                    * (1 - (eps / self._eps_0)) ** (self._n - 1)
                )
            return tangent
        # If it is an array
        # parabolic branch
        tangent = np.zeros_like(eps)
        tangent[(eps <= 0) & (eps >= self._eps_0)] = (
            self._n
            * self._fc
            / self._eps_0
            * (1 - (eps[(eps <= 0) & (eps >= self._eps_0)] / self._eps_0))
            ** (self._n - 1)
        )
        # Elsewhere tangent is zero
        tangent[eps < self._eps_0] = 0.0
        tangent[eps > 0] = 0.0
        return tangent

    def __marin__(
        self, strain: t.Tuple[float, float]
    ) -> t.Tuple[t.List[float], t.List[float]]:
        """Returns coefficients and strain limits for Marin integration in a
        simply formatted way.

        Arguments:
            strain (float, float): Tuple defining the strain profile: eps =
                strain[0] + strain[1]*y.

        Example:
            [(0, -0.002), (-0.002, -0.003)]
            [(a0, a1, a2), (a0)]
        """
        if self._n != 2:
            # The constitutive law is not writtable as a polynomial,
            # Call the generic distretizing method
            return super().__marin__(strain=strain)

        strains = []
        coeff = []
        if strain[1] == 0:
            # Uniform strain equal to strain[0]
            # understand in which branch are we
            strain[0] = self.preprocess_strains_with_limits(strain[0])
            if strain[0] > 0:
                # We are in tensile branch
                strains = None
                coeff.append((0.0,))
            elif strain[0] > self._eps_0:
                # We are in the parabolic branch
                strains = None
                a0 = (
                    2
                    * self._fc
                    * strain[0]
                    / self._eps_0
                    * (1 - 0.5 * (strain[0] / self._eps_0))
                )
                a1 = (
                    2
                    * self._fc
                    / self._eps_0
                    * strain[1]
                    * (1 - strain[0] / self._eps_0)
                )
                coeff.append((a0, a1, 0.0))
            elif strain[0] >= self._eps_u:
                # We are in the constant branch
                strains = None
                coeff.append((self._fc,))
            else:
                # We are in a branch of non-resisting concrete
                # Too much compression
                strains = None
                coeff.append((0.0,))
        else:
            # Parabolic part
            strains.append((self._eps_0, 0))
            a0 = (
                2
                * self._fc
                * strain[0]
                / self._eps_0
                * (1 - 0.5 * (strain[0] / self._eps_0))
            )
            a1 = (
                2
                * self._fc
                / self._eps_0
                * strain[1]
                * (1 - strain[0] / self._eps_0)
            )
            a2 = -self._fc * strain[1] ** 2 / self._eps_0**2
            coeff.append((a0, a1, a2))
            # Constant part
            strains.append((self._eps_u, self._eps_0))
            coeff.append((self._fc,))
        return strains, coeff

    def get_ultimate_strain(
        self, yielding: bool = False
    ) -> t.Tuple[float, float]:
        """Return the ultimate strain (negative and positive)."""
        if yielding:
            return (self._eps_0, 100)
        return (self._eps_u, 100)
