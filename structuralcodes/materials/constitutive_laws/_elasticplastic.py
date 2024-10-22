"""Collection of some standard constitutive laws."""

from __future__ import annotations  # To have clean hints of ArrayLike in docs

import typing as t

import numpy as np
from numpy.typing import ArrayLike

from ...core.base import ConstitutiveLaw


class ElasticPlastic(ConstitutiveLaw):
    """Class for elastic-plastic Constitutive Law."""

    __materials__: t.Tuple[str] = (
        'steel',
        'rebars',
    )

    def __init__(
        self,
        E: float,
        fy: float,
        Eh: float = 0.0,
        eps_su: t.Optional[float] = None,
        name: t.Optional[str] = None,
    ) -> None:
        """Initialize an Elastic-Plastic Material.

        Arguments:
            E (float): The elastic modulus.
            fy (float): The yield strength.

        Keyword Arguments:
            Eh (float): The hardening modulus.
            eps_su (float): The ultimate strain.
            name (str): A descriptive name for the constitutive law.
        """
        name = name if name is not None else 'ElasticPlasticLaw'
        super().__init__(name=name)
        if E > 0:
            self._E = E
        else:
            raise ValueError('Elastic modulus E must be greater than zero')
        self._fy = fy
        self._Eh = Eh
        self._eps_su = eps_su
        self._eps_sy = fy / E

    def get_stress(self, eps: ArrayLike) -> t.Union[float, ArrayLike]:
        """Return the stress given strain."""
        eps = eps if np.isscalar(eps) else np.atleast_1d(eps)
        # Preprocess eps array in order
        eps = self.preprocess_strains_with_limits(eps=eps)
        # Compute stress
        sig = self._E * eps
        delta_sig = self._fy * (1 - self._Eh / self._E)
        if np.isscalar(sig):
            if sig < -self._fy:
                sig = eps * self._Eh - delta_sig
            if sig > self._fy:
                sig = eps * self._Eh + delta_sig
            if (self._eps_su is not None) and (
                eps > self._eps_su or eps < -self._eps_su
            ):
                sig = 0
            return sig
        sig[sig < -self._fy] = eps[sig < -self._fy] * self._Eh - delta_sig
        sig[sig > self._fy] = eps[sig > self._fy] * self._Eh + delta_sig
        if self._eps_su is not None:
            sig[eps > self._eps_su] = 0
            sig[eps < -self._eps_su] = 0  # pylint: disable=E1130
        return sig

    def get_tangent(self, eps: ArrayLike) -> t.Union[float, ArrayLike]:
        """Return the tangent for given strain."""
        if np.isscalar(eps):
            tangent = (
                self._E if -self._eps_sy <= eps <= self._eps_sy else self._Eh
            )
            if (self._eps_su is not None) and (
                eps > self._eps_su or eps < -self._eps_su
            ):
                tangent = 0
            return tangent

        eps = np.atleast_1d(eps)
        tangent = np.ones_like(eps) * self._E
        tangent[eps > self._eps_sy] = self._Eh
        tangent[eps < -self._eps_sy] = self._Eh
        if self._eps_su is not None:
            tangent[eps > self._eps_su] = 0
            tangent[eps < -self._eps_su] = 0  # pylint: disable=E1130

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
        eps_sy_n, eps_sy_p = self.get_ultimate_strain(yielding=True)
        eps_su_n, eps_su_p = self.get_ultimate_strain()
        if strain[1] == 0:
            # Uniform strain equal to strain[0]
            # Understand in which branch are we
            strain[0] = self.preprocess_strains_with_limits(strain[0])
            if strain[0] > eps_sy_p and strain[0] <= eps_su_p:
                # We are in the Hardening part positive
                strains = None
                a0 = self._Eh * strain[0] + self._fy * (1 - self._Eh / self._E)
                a1 = self._Eh * strain[1]
                coeff.append((a0, a1))
            elif strain[0] < eps_sy_n and strain[0] >= eps_su_n:
                # We are in the Hardening part negative
                strains = None
                a0 = self._Eh * strain[0] - self._fy * (1 - self._Eh / self._E)
                a1 = self._Eh * strain[1]
                coeff.append((a0, a1))
            elif abs(strain[0]) <= self._eps_sy:
                # We are in the elastic part
                strains = None
                a0 = self._E * strain[0]
                a1 = self._E * strain[1]
                coeff.append((a0, a1))
            else:
                strains = None
                coeff.append((0.0,))
        else:
            # Hardening part negative
            strains.append((eps_su_n, eps_sy_n))
            a0 = self._Eh * strain[0] - self._fy * (1 - self._Eh / self._E)
            a1 = self._Eh * strain[1]
            coeff.append((a0, a1))
            # Elastic part
            strains.append((eps_sy_n, eps_sy_p))
            a0 = self._E * strain[0]
            a1 = self._E * strain[1]
            coeff.append((a0, a1))
            # Hardening part positive
            strains.append((eps_sy_p, eps_su_p))
            a0 = self._Eh * strain[0] + self._fy * (1 - self._Eh / self._E)
            a1 = self._Eh * strain[1]
            coeff.append((a0, a1))
        return strains, coeff

    def get_ultimate_strain(
        self, yielding: bool = False
    ) -> t.Tuple[float, float]:
        """Return the ultimate strain (negative and positive)."""
        if yielding:
            return (-self._eps_sy, self._eps_sy)
        # If not specified eps
        if self._eps_su is None:
            return (-self._eps_sy * 2, self._eps_sy * 2)
        return (-self._eps_su, self._eps_su)
