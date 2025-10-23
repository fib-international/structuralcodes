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
        initial_strain: t.Optional[float] = None,
        initial_stress: t.Optional[float] = None,
        strain_compatibility: t.Optional[float] = None,
        name: t.Optional[str] = None,
    ):
        """Initialize a material with an elastic plastic constitutive law.

        Arguments:
            E (float): The Young's modulus.
            density (float): The density.
            initial_strain (Optional[float]): Initial strain of the material.
            initial_stress (Optional[float]): Initial stress of the material.
            strain_compatibility (Optional[bool]): Only relevant if
                initial_strain or initial_stress are different from zero. If
                True, the material deforms with the geometry. If False, the
                stress in the material upon loading is kept constant
                corresponding to the initial strain.
            name (str, optional): The name of the material, default value None.
        """
        super().__init__(
            density=density,
            initial_strain=initial_strain,
            initial_stress=initial_stress,
            strain_compatibility=strain_compatibility,
            name=name if name else 'ElasticMaterial',
        )
        self._E = E
        self._constitutive_law = create_constitutive_law('elastic', self)
        self._apply_initial_strain()

    @property
    def E(self) -> float:
        """Returns the Young's modulus."""
        return self._E

    def __elastic__(self) -> dict:
        """Returns kwargs for creating an elastic constitutive law."""
        return {'E': self.E}

    @classmethod
    def from_material(cls, other_material: Material):
        """Create an elastic material based on another material."""
        # Create name of elastic material
        name = other_material.name
        if name is not None:
            name += '_elastic'

        return cls(
            E=other_material.constitutive_law.get_tangent(eps=0),
            density=other_material.density,
            name=name,
        )
