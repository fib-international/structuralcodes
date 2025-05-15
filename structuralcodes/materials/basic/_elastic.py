"""A material class with elastic properties."""

import typing as t

from ...core.base import Material
from ..constitutive_laws import create_constitutive_law


class ElasticMaterial(Material):
    """A material class with elastic properties."""

    _E: float

    def __init__(
        self,
        E: float,
        density: float,
        name: t.Optional[str] = None,
    ):
        """Initialize a material with an elastic plastic constitutive law.

        Arguments:
            E (float): The Young's modulus.
            density (float): The density.
            name (str, optional): The name of the material, default value None.
        """
        super().__init__(density=density, name=name)
        self._E = E
        self._constitutive_law = create_constitutive_law('elastic', self)

    @property
    def E(self) -> float:
        """Returns the Young's modulus."""
        return self._E

    def __elastic__(self) -> dict:
        """Returns kwargs for creating an elastic constitutive law."""
        return {'E': self.E}
