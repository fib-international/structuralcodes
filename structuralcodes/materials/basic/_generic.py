"""A generic material that could hold any type of constitutive law."""

import typing as t

from ...core.base import ConstitutiveLaw, Material


class GenericMaterial(Material):
    """A material class that accepts any constitutive law."""

    def __init__(
        self,
        density: float,
        constitutive_law: ConstitutiveLaw,
        initial_strain: t.Optional[float] = None,
        initial_stress: t.Optional[float] = None,
        name: t.Optional[str] = None,
    ):
        """Initialize a material with a constitutive law.

        Arguments:
            density (float): The density.
            constitutive_law (ConstitutiveLaw): The constitutive law of the
                material.
            initial_strain (float, optional): The initial strain of the
                material, default value None.
            initial_stress (float, optional): The initial stress of the
                material, default value None.
            name (str, optional): The name of the material, default value None.
        """
        super().__init__(
            density=density,
            initial_strain=initial_strain,
            initial_stress=initial_stress,
            name=name,
        )
        self._constitutive_law = constitutive_law
        self._apply_initial_strain()
