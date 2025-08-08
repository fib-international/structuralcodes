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
        initial_strain: t.Optional[float] = None,
        initial_stress: t.Optional[float] = None,
        strain_compatibility: t.Optional[float] = None,
        name: t.Optional[str] = None,
    ):
        """Initialize a material with an elastic plastic constitutive law.

        Arguments:
            E (float): The Young's modulus.
            fy (float): The yield stress.
            density (float): The density.
            Eh (float, optional): The hardening modulus, default value 0.
            eps_su (float, optional): The ultimate strain, default value None.
            initial_strain (float, optional): The initial strain of the
                material, default value None.
            initial_stress (float, optional): The initial stress of the
                material, default value None.
            strain_compatibility (float, optional): If False, the strain
                compatibility is not enforced, default value None.
            name (str, optional): The name of the material, default value None.
        """
        super().__init__(
            density=density,
            initial_strain=initial_strain,
            initial_stress=initial_stress,
            strain_compatibility=strain_compatibility,
            name=name if name else 'ElasticPlasticMaterial',
        )
        self._E = E
        self._fy = fy
        self._Eh = Eh
        self._eps_su = eps_su

        self._constitutive_law = create_constitutive_law(
            'elasticplastic', self
        )
        self._apply_initial_strain()

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
