"""Elastic-fragile tension constitutive law."""

from __future__ import annotations  # To have clean hints of ArrayLike in docs

import typing as t

import numpy as np
from numpy.typing import ArrayLike

from ...core.base import ConstitutiveLaw


class ElasticFragileTension(ConstitutiveLaw):
    """Class for an elastic-fragile constitutive law for concrete in tension
    (only tension behavior).

    The law is linear elastic in tension until the tensile strength is
    reached, followed by an abrupt loss of strength. It carries no stress in
    compression, and is meant to be combined with a constitutive law for
    concrete in compression in a :class:`Parallel` constitutive law.
    """

    __materials__: t.Tuple[str] = ('concrete',)

    def __init__(
        self,
        fct: float,
        Ec: float,
        name: t.Optional[str] = None,
    ) -> None:
        """Initialize an ElasticFragileTension Material.

        Arguments:
            fct (float): The tensile strength.
            Ec (float): The modulus of elasticity in tension.

        Keyword Arguments:
            name (str): A descriptive name for the constitutive law.

        Raises:
            ValueError: If fct or Ec is zero.

        Note:
            If negative values are input for fct or Ec, they will be assumed
            positive.
        """
        super().__init__(name=name, base_name='ElasticFragileTensionLaw')
        if fct == 0:
            raise ValueError('Tensile strength should be different from 0')
        if Ec == 0:
            raise ValueError('Elastic modulus should be different from 0')
        self._fct = abs(fct)
        self._Ec = abs(Ec)
        self._eps_ct = self._fct / self._Ec

    @property
    def fct(self) -> float:
        """Returns the tensile strength."""
        return self._fct

    @property
    def Ec(self) -> float:
        """Returns the modulus of elasticity."""
        return self._Ec

    @property
    def eps_ct(self) -> float:
        """Returns the strain at tensile strength."""
        return self._eps_ct

    def get_stress(
        self, eps: t.Union[float, ArrayLike]
    ) -> t.Union[float, ArrayLike]:
        """Return the stress given strain."""
        eps = eps if np.isscalar(eps) else np.atleast_1d(eps)
        # If it is a scalar
        if np.isscalar(eps):
            if 0 < eps <= self._eps_ct:
                return self._Ec * eps
            return 0.0
        # If it is an array
        sig = np.zeros(eps.shape)
        elastic = (eps > 0) & (eps <= self._eps_ct)
        sig[elastic] = self._Ec * eps[elastic]
        return sig

    def get_tangent(
        self, eps: t.Union[float, ArrayLike]
    ) -> t.Union[float, ArrayLike]:
        """Return the tangent given strain.

        Note:
            The tangent is zero at zero strain. When the law is combined with
            a law for concrete in compression in a Parallel constitutive law,
            the tangent at zero strain is that of the law in compression.
        """
        eps = eps if np.isscalar(eps) else np.atleast_1d(eps)
        # If it is a scalar
        if np.isscalar(eps):
            if 0 < eps <= self._eps_ct:
                return self._Ec
            return 0.0
        # If it is an array
        tangent = np.zeros(eps.shape)
        tangent[(eps > 0) & (eps <= self._eps_ct)] = self._Ec
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
            [(0, 0.0001)]
            [(a0, a1)]
        """
        strains = []
        coeff = []
        if strain[1] == 0:
            # Uniform strain equal to strain[0]
            # understand in which branch we are
            strains = None
            if 0 < strain[0] <= self._eps_ct:
                # We are in the elastic branch
                coeff.append((self._Ec * strain[0],))
            else:
                # We are in compression or after cracking: no stress
                coeff.append((0.0,))
        else:
            # Only the elastic branch carries stress
            strains.append((0, self._eps_ct))
            a0 = self._Ec * strain[0]
            a1 = self._Ec * strain[1]
            coeff.append((a0, a1))
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
            [(0, 0.0001)]
            [(a0,)]
        """
        strains = []
        coeff = []
        if strain[1] == 0:
            # Uniform strain equal to strain[0]
            # understand in which branch we are
            strains = None
            if 0 < strain[0] <= self._eps_ct:
                # We are in the elastic branch
                coeff.append((self._Ec,))
            else:
                # We are in compression or after cracking: no stiffness
                coeff.append((0.0,))
        else:
            # Only the elastic branch has stiffness
            strains.append((0, self._eps_ct))
            coeff.append((self._Ec,))
        return strains, coeff

    def get_ultimate_strain(self, **kwargs) -> t.Tuple[float, float]:
        """Return the ultimate strain (negative and positive).

        Note:
            The law returns the limits of the strain range where it carries
            stress, i.e. zero and the strain at tensile strength. When the law
            is combined with a law for concrete in compression in a Parallel
            constitutive law, the ultimate strains of the Parallel law are
            those of the law in compression.
        """
        del kwargs
        return (0.0, self._eps_ct)
