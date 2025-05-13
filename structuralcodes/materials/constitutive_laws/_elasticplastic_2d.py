"""Elastic-plastic constitutive law."""

from __future__ import annotations

import typing as t

import numpy as np
from numpy.typing import ArrayLike

from structuralcodes.materials.constitutive_laws._elasticplastic import (
    ElasticPlastic,
)


class ElasticPlastic2D(ElasticPlastic):
    """Class for elastic-plastic Constitutive Law in 2D."""

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
        """Initialize an Elastic-Plastic 2D Material.

        Arguments:
            E (float): The elastic modulus.
            fy (float): The yield strength.

        Keyword Arguments:
            Eh (float): The hardening modulus.
            eps_su (float): The ultimate strain.
            name (str): A descriptive name for the constitutive law.
        """
        name = name if name is not None else 'ElasticPlasticLaw2D'
        super().__init__(E=E, fy=fy, Eh=Eh, eps_su=eps_su, name=name)

    @property
    def E(self) -> float:
        """Return the elastic modulus."""
        return self._E

    @property
    def C_s(self) -> np.ndarray:
        """Return the 2D constitutive matrix."""
        return self.E * np.array(
            [
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0],
            ]
        )

    def get_stress(self, eps: ArrayLike) -> np.ndarray:
        """Return the stress given strain."""
        sig_s = super().get_stress(eps)
        return sig_s @ self.C_s / self.E

    def get_tangent(self, eps: ArrayLike) -> np.ndarray:
        """Compute the 3x3 tangent stiffness matrix C."""
        eps = np.atleast_1d(eps)
        return self.C_s
