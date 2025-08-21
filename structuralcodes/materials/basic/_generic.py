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
        strain_compatibility: t.Optional[bool] = None,
        name: t.Optional[str] = None,
    ):
        """Initialize a material with a constitutive law.

        Arguments:
            density (float): The density.
            constitutive_law (ConstitutiveLaw): The constitutive law of the
                material.
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
            name=name if name else 'GenericMaterial',
        )
        self._constitutive_law = constitutive_law
        self._apply_initial_strain()
