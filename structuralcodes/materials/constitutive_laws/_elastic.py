"""Elastic constitutive law."""

from __future__ import annotations  # To have clean hints of ArrayLike in docs

import typing as t

import numpy as np
from numpy.typing import ArrayLike

from ...core.base import ConstitutiveLaw


class Elastic(ConstitutiveLaw):
    """Class for elastic constitutive law."""

    __materials__: t.Tuple[str] = (
        'concrete',
        'steel',
        'rebars',
    )

    def __init__(self, E: float, name: t.Optional[str] = None) -> None:
        """Initialize an Elastic Material.

        Arguments:
            E (float): The elastic modulus.

        Keyword Arguments:
            name (str): A descriptive name for the constitutive law.
        """
        name = name if name is not None else 'ElasticLaw'
        super().__init__(name=name)
        self._E = E
        self._eps_su = None

    def get_stress(
        self, eps: t.Union[float, ArrayLike]
    ) -> t.Union[float, ArrayLike]:
        """Return stress given strain."""
        eps = eps if np.isscalar(eps) else np.atleast_1d(eps)
        return self._E * eps

    def get_tangent(
        self, eps: t.Union[float, ArrayLike]
    ) -> t.Union[float, ArrayLike]:
        """Return the tangent."""
        if np.isscalar(eps):
            return self._E
        eps = np.atleast_1d(eps)
        return np.ones_like(eps) * self._E

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
        strains = None
        a0 = self._E * strain[0]
        a1 = self._E * strain[1]
        coeff = [(a0, a1)]
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
        strains = None
        a0 = self._E
        coeff = [(a0,)]
        return strains, coeff

    def get_ultimate_strain(self, **kwargs) -> t.Tuple[float, float]:
        """Return the ultimate strain (negative and positive)."""
        # There is no real strain limit, so set it to very large values
        # unlesse specified by the user differently
        del kwargs
        return self._eps_su or (-100, 100)

    def set_ultimate_strain(
        self, eps_su=t.Union[float, t.Tuple[float, float]]
    ) -> None:
        """Set ultimate strains for Elastic Material if needed.

        Arguments:
            eps_su (float or (float, float)): Defining ultimate strain. If a
                single value is provided the same is adopted for both negative
                and positive strains. If a tuple is provided, it should be
                given as (negative, positive).
        """
        if isinstance(eps_su, float):
            self._eps_su = (-abs(eps_su), abs(eps_su))
        elif isinstance(eps_su, tuple):
            if len(eps_su) < 2:
                raise ValueError(
                    'Two values need to be provided when setting the tuple'
                )
            eps_su_n = eps_su[0]
            eps_su_p = eps_su[1]
            if eps_su_p < eps_su_n:
                eps_su_p, eps_su_n = eps_su_n, eps_su_p
            if eps_su_p < 0:
                raise ValueError(
                    'Positive ultimate strain should be non-negative'
                )
            if eps_su_n > 0:
                raise ValueError(
                    'Negative utimate strain should be non-positive'
                )
            self._eps_su = (eps_su_n, eps_su_p)
        else:
            raise ValueError(
                'set_ultimate_strain requires a single value or a tuple \
                with  two values'
            )
