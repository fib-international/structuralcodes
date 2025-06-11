"""A generic material that could hold any type of constitutive law."""

import typing as t

from ...core.base import ConstitutiveLaw, Material


class GenericMaterial(Material):
    """A material class that accepts any constitutive law."""

    def __init__(
        self,
        density: float,
        constitutive_law: ConstitutiveLaw,
        name: t.Optional[str] = None,
    ):
        """Initialize a material with a constitutive law.

        Arguments:
            density (float): The density.
            constitutive_law (ConstitutiveLaw): The constitutive law of the
                material.
            name (str, optional): The name of the material, default value None.
        """
        super().__init__(density=density, name=name)
        self._constitutive_law = constitutive_law
