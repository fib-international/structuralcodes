"""Collection of some standard constitutive laws."""

from __future__ import annotations  # To have clean hints of ArrayLike in docs

import typing as t

import numpy as np
from numpy.typing import ArrayLike

from ...core.base import ConstitutiveLaw


class Sargin(ConstitutiveLaw):
    """Class for Sargin constitutive law.

    The stresses and strains are assumed negative in compression and positive
    in tension.

    References:
    Sargin, M. (1971), "Stress-strain relationship for concrete and the
    analysis of structural concrete section, Study No. 4,
    Solid Mechanics Division, University of Waterloo, Ontario, Canada
    """

    __materials__: t.Tuple[str] = ('concrete',)

    def __init__(
        self,
        fc: float,
        eps_c1: float = -0.0023,
        eps_cu1: float = -0.0035,
        k: float = 2.04,
        name: t.Optional[str] = None,
    ) -> None:
        """Initialize a Sargin Material.

        Arguments:
            fc (float): The strength of concrete in compression

        Keyword Arguments:
            eps_c1 (float): Peak strain of concrete in compression. Default
                value = -0.0023.
            eps_u (float): Ultimate strain of concrete in compression. Default
                value = -0.0035.
            k (float): Plasticity number. Default value = 2.04.
            name (str): A name for the constitutive law.

        Raises:
            ValueError: If k is less or equal to 0.

        Note:
            If positive values are input for fc, eps_c1 and eps_cu1 are input,
            they will be assumed negative.
        """
        name = name if name is not None else 'SarginLaw'
        super().__init__(name=name)
        self._fc = -abs(fc)
        self._eps_c1 = -abs(eps_c1)
        self._eps_cu1 = -abs(eps_cu1)
        self._k = k

    def get_stress(self, eps: ArrayLike) -> t.Union[float, ArrayLike]:
        """Return the stress given the strain."""
        eps = eps if np.isscalar(eps) else np.atleast_1d(eps)
        # Preprocess eps array in order
        eps = self.preprocess_strains_with_limits(eps=eps)
        # Compute stress
        # Polynomial branch
        eta = eps / self._eps_c1

        sig = self._fc * (self._k * eta - eta**2) / (1 + (self._k - 2) * eta)

        # Elsewhere stress is 0.0
        if np.isscalar(eps):
            if eps < self._eps_cu1 or eps > 0:
                return 0.0
        else:
            sig[eps < self._eps_cu1] = 0.0
            sig[eps > 0] = 0.0

        return sig

    def get_tangent(self, eps: ArrayLike) -> t.Union[float, ArrayLike]:
        """Return the tangent given strain."""
        eps = eps if np.isscalar(eps) else np.atleast_1d(eps)
        # polynomial branch
        eta = eps / self._eps_c1

        tangent = (
            self._fc
            / self._eps_c1
            * ((2 - self._k) * eta**2 - 2 * eta + self._k)
            / (1 + (self._k - 2) * eta) ** 2
        )
        # Elsewhere tangent is zero
        if np.isscalar(eps):
            if eps < self._eps_cu1 or eps > 0:
                return 0
        else:
            tangent[eps < self._eps_cu1] = 0.0
            tangent[eps > 0] = 0.0

        return tangent

    def get_ultimate_strain(
        self, yielding: bool = False
    ) -> t.Tuple[float, float]:
        """Return the ultimate strain (negative and positive)."""
        if yielding:
            return (self._eps_c1, 100)
        return (self._eps_cu1, 100)
