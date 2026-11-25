"""Popovics constitutive law."""

from __future__ import annotations  # To have clean hints of ArrayLike in docs

import typing as t

import numpy as np
from numpy.typing import ArrayLike

from ...core.base import ConstitutiveLaw


class Popovics(ConstitutiveLaw):
    """Class for Popovics-Mander constitutive law.

    The stresses and strains are assumed negative in compression and positive
    in tension.

    If the relation Ec = 5000 * sqrt(fc) is used for elastic modulus, the
    constitutive law is identical to the one proposed by Mander et al. (1988).

    References:
    Popovics, S., 1973, “A Numerical Approach to the Complete Stress-Strain
    Curve of Concrete”, Cement and Concrete Research, 3(4), 583-599.

    Mander, J.B., Priestley, M.J.N., Park, R., 1988, "Theoretical Stress-Strain
    Model for Confined Concrete", Journal of Structural Engineering, 114(8),
    1804-1826.
    """

    __materials__: t.Tuple[str] = ('concrete',)

    def __init__(
        self,
        fc: float,
        eps_c: float = -0.002,
        eps_cu: float = -0.0035,
        Ec: t.Optional[float] = None,
        name: t.Optional[str] = None,
    ) -> None:
        """Initialize a Popovics Material.

        Arguments:
            fc (float): the strength of concrete in compression

        Keyword Arguments:
            eps_c (float): Peak strain of concrete in compression. Default
                value = -0.002.
            eps_cu (float): Ultimate strain of concrete in compression. Default
                value = -0.0035.
            E (optional float): Elastic modulus of concrete. If None, the
                equation Ec = 5000 * fc**0.5 proposed by Mander et al. (1988)
                is adopted (fc in MPa). Default value = None.
            name (str): A name for the constitutive law.

        Raises:
            ValueError: If E is less or equal to 0.

        Note:
            If positive values are input for fc, eps_c and eps_cu are input,
            they will be assumed negative.
        """
        super().__init__(name=name, base_name='PopovicsLaw')
        self._fc = -abs(fc)
        self._eps_c = -abs(eps_c)
        self._eps_cu = -abs(eps_cu)
        if Ec is None:
            # fc in MPa, relation of Mander et al. (1988)
            Ec = 5000 * abs(fc) ** 0.5
        if Ec <= 0:
            raise ValueError('Elastic modulus must be a positive number.')
        E_sec = self._fc / self._eps_c
        self._n = Ec / (Ec - E_sec)

    def get_stress(
        self, eps: t.Union[float, ArrayLike]
    ) -> t.Union[float, ArrayLike]:
        """Return the stress given the strain."""
        eps = eps if np.isscalar(eps) else np.atleast_1d(eps)
        # Preprocess eps array in order
        eps = self.preprocess_strains_with_limits(eps=eps)
        # The law is defined only in compression, for eps_cu <= eps <= 0, and
        # the stress is zero elsewhere. The expression is evaluated only in
        # this range: for eps > 0 it would be eta < 0, and eta**n with a
        # non-integer n is not defined (NaN and RuntimeWarning).
        if np.isscalar(eps):
            if eps < self._eps_cu or eps > 0:
                return 0.0
            eta = eps / self._eps_c
            return self._fc * eta * self._n / (self._n - 1 + eta**self._n)

        sig = np.zeros_like(eps, dtype=float)
        in_range = (eps >= self._eps_cu) & (eps <= 0)
        eta = eps[in_range] / self._eps_c
        sig[in_range] = self._fc * eta * self._n / (self._n - 1 + eta**self._n)
        return sig

    def get_tangent(
        self, eps: t.Union[float, ArrayLike]
    ) -> t.Union[float, ArrayLike]:
        """Return the tangent given strain."""
        eps = eps if np.isscalar(eps) else np.atleast_1d(eps)
        # Preprocess eps array in order
        eps = self.preprocess_strains_with_limits(eps=eps)

        def _tangent(eta):
            return (
                (1 - eta**self._n)
                / (self._n - 1 + eta**self._n) ** 2
                * self._n
                * (self._n - 1)
                * self._fc
                / self._eps_c
            )

        # The law is defined only in compression, for eps_cu <= eps <= 0, and
        # the tangent is zero elsewhere. The expression is evaluated only in
        # this range (see get_stress).
        if np.isscalar(eps):
            if eps < self._eps_cu or eps > 0:
                return 0
            return _tangent(eps / self._eps_c)

        tangent = np.zeros_like(eps, dtype=float)
        in_range = (eps >= self._eps_cu) & (eps <= 0)
        tangent[in_range] = _tangent(eps[in_range] / self._eps_c)
        return tangent

    def get_ultimate_strain(
        self, yielding: bool = False
    ) -> t.Tuple[float, float]:
        """Return the ultimate strain (negative and positive)."""
        if yielding:
            return (self._eps_c, 100)
        return (self._eps_cu, 100)
