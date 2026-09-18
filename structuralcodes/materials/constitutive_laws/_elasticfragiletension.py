"""Elastic-fragile tension constitutive law."""

from __future__ import annotations  # To have clean hints of ArrayLike in docs

import typing as t

import numpy as np
from numpy.typing import ArrayLike

from ...core.base import ConstitutiveLaw


class ElasticFragileTension(ConstitutiveLaw):
    """Class elastic-fragile constitutive law in tension only.

    The law is linear elastic in tension up to the tensile strength, where
    it cracks and the stress drops to zero. It carries nothing in
    compression, and is meant to be combined with a compressive law
    through :class:`Parallel`::

    Example :
        compression = ParabolaRectangle(fc=45)
        tension = ElasticFragileTension(fct=2.7, Ec=36000)
        law = Parallel(constitutive_laws=[compression, tension])
    """

    __materials__: t.Tuple[str] = ('concrete',)

    def __init__(
        self,
        fct: float,
        Ec: float,
        name: t.Optional[str] = None,
    ) -> None:
        """Initialize an elastic-fragile tension Material.

        Arguments:
            fct (float): The strength of concrete in tension.
            Ec (float): The modulus of elasticity in tension.

        Keyword Arguments:
            name (str): A name for the constitutive law.

        Raises:
            ValueError: If fct or Ec is not strictly positive.

        Note:
            If negative values are input for fct or Ec, they will be
            assumed positive.
        """
        super().__init__(name=name, base_name='ElasticFragileTensionLaw')
        if fct == 0:
            raise ValueError('Tensile strength should be different from 0')
        if Ec == 0:
            raise ValueError('Elastic modulus should be different from 0')
        self._fct = +abs(fct)
        self._Ec = +abs(Ec)
        self._eps_ct_u = self._fct / self._Ec

    @property
    def fct(self) -> float:
        """Returns the tensile strength."""
        return self._fct

    @property
    def Ec(self) -> float:
        """Returns the modulus of elasticity."""
        return self._Ec

    @property
    def eps_ct_u(self) -> float:
        """Returns the cracking strain, where the law fails."""
        return self._eps_ct_u

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
            if 0 < eps <= self._eps_ct_u:
                return self._Ec * eps
            return 0.0
        # If it is an array
        sig = np.zeros_like(eps)
        elastic = (eps > 0) & (eps <= self._eps_ct_u)
        sig[elastic] = self._Ec * eps[elastic]
        return sig

    def get_tangent(
        self, eps: t.Union[float, ArrayLike]
    ) -> t.Union[float, ArrayLike]:
        """Return the tangent given strain."""
        eps = eps if np.isscalar(eps) else np.atleast_1d(eps)
        # If it is a scalar
        if np.isscalar(eps):
            if 0 < eps <= self._eps_ct_u:
                return self._Ec
            return 0.0
        # If it is an array
        tangent = np.zeros_like(eps)
        tangent[(eps > 0) & (eps <= self._eps_ct_u)] = self._Ec
        return tangent

    def get_ultimate_strain(self, **kwargs) -> t.Tuple[float, float]:
        """Return the ultimate strain (negative and positive).

        The law is fragile: cracking is both the end of the elastic branch.
        In compression it resists nothing, and reports zero
        """
        del kwargs
        return (-np.inf, self._eps_ct_u)
