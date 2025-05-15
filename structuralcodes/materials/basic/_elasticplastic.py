"""A material class with elastic plastic properties."""

import typing as t

from ...core.base import Material
from ..constitutive_laws import create_constitutive_law


class ElasticPlasticMaterial(Material):
    """A material class with elastic plastic properties."""

    _E: float
    _fy: float
    _Eh: float
    _eps_su: float

    def __init__(
        self,
        E: float,
        fy: float,
        density: float,
        Eh: float = 0,
        eps_su: t.Optional[float] = None,
        name: t.Optional[str] = None,
    ):
        """Initialize a material with an elastic plastic constitutive law.

        Arguments:
            E (float): The Young's modulus.
            fy (float): The yield stress.
            density (float): The density.
            Eh (float, optional): The hardening modulus, default value 0.
            eps_su (float, optional): The ultimate strain, default value None.
            name (str, optional): The name of the material, default value None.
        """
        super().__init__(density=density, name=name)
        self._E = E
        self._fy = fy
        self._Eh = Eh
        self._eps_su = eps_su

        self._constitutive_law = create_constitutive_law(
            'elasticplastic', self
        )

    @property
    def E(self) -> float:
        """Returns the Young's modulus."""
        return self._E

    @property
    def fy(self) -> float:
        """Returns the yield stress."""
        return self._fy

    @property
    def Eh(self) -> float:
        """Returns the hardening modulus."""
        return self._Eh

    @property
    def eps_su(self) -> float:
        """Returns the ultimate strain."""
        return self._eps_su

    def __elasticplastic__(self) -> dict:
        """Returns kwargs for ElasticPlastic constitutive law with strain
        hardening.
        """
        return {
            'E': self.E,
            'fy': self.fy,
            'Eh': self.Eh,
            'eps_su': self.eps_su,
        }
