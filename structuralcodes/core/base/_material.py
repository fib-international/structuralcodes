"""Abstract base class for materials."""

from __future__ import annotations  # To have clean hints of ArrayLike in docs

import abc
import typing as t

from ._constitutive_law import ConstitutiveLaw


class Material(abc.ABC):
    """Abstract base class for materials."""

    _constitutive_law = None
    _initial_strain: t.Optional[float] = None
    _initial_stress: t.Optional[float] = None
    _strain_compatibility: t.Optional[bool] = None

    def __init__(
        self,
        density: float,
        initial_strain: t.Optional[float] = None,
        initial_stress: t.Optional[float] = None,
        strain_compatibility: t.Optional[bool] = None,
        name: t.Optional[str] = None,
    ) -> None:
        """Initializes an instance of a new material.

        Args:
            density (float): Density of the material in kg/m3.

        Keyword Args:
            initial_strain (Optional[float]): Initial strain of the material.
            initial_stress (Optional[float]): Initial stress of the material.
            strain_compatibility (Optional[bool]): Only relevant if
                initial_strain or initial_stress are different from zero. If
                True, the material deforms with the geometry. If False, the
                stress in the material upon loading is kept constant
                corresponding to the initial strain.
            name (Optional[str]): descriptive name of the material

        Raise:
            ValueError: if both initial_strain and initial_stress are provided
        """
        self._density = abs(density)
        if initial_strain is not None and initial_stress is not None:
            raise ValueError(
                'Both initial_strain and initial_stress cannot be provided.'
            )
        self._initial_strain = initial_strain
        self._initial_stress = initial_stress
        self._strain_compatibility = strain_compatibility
        self._name = name if name is not None else 'Material'

    @property
    def constitutive_law(self) -> ConstitutiveLaw:
        """Returns the ConstitutiveLaw of the object."""
        return self._constitutive_law

    @property
    def name(self):
        """Returns the name of the material."""
        return self._name

    @property
    def density(self):
        """Returns the density of the material in kg/m3."""
        return self._density

    @property
    def initial_strain(self):
        """Returns the initial strain of the material."""
        return self._initial_strain

    @property
    def initial_stress(self):
        """Returns the initial stress of the material."""
        return self._initial_stress

    @property
    def strain_compatibility(self):
        """Returns the strain compatibility of the material.

        If true (default), the strain compatibility is enforced
        haveing the same strain as in all other materials of the
        section at the same point. If false, the strain compatibility
        is not enforced and the initial strain is applied to the section
        independently.
        """
        return self._strain_compatibility

    def _apply_initial_strain(self):
        """Wraps the current constitutive law to apply initial strain."""
        strain_compatibility = (
            self._strain_compatibility
            if self._strain_compatibility is not None
            else True
        )
        if self._initial_stress is not None:
            # Specified a stress, compute the strain from it
            self._initial_strain_from_stress()
        if self._initial_strain is not None:
            # Lazy import to avoid circular dependency
            from structuralcodes.materials.constitutive_laws import (  # noqa: PLC0415
                InitialStrain,
            )

            if self._initial_stress is None:
                # Compute the stress from the strain
                self._initial_stress = self._constitutive_law.get_stress(
                    self._initial_strain
                )

            self._constitutive_law = InitialStrain(
                self._constitutive_law,
                self._initial_strain,
                strain_compatibility,
            )

    def _initial_strain_from_stress(self):
        """Computes the initial strain from the initial stress.

        This function is called internally so it assumes that the
        initial stress is not None
        """
        # Iteratively compute the initial strain that gives the desired
        # initial stress. Note that the wrapped law can be nonlinear
        tol = 1e-12
        max_iter = 100
        target_stress = self._initial_stress
        strain = 0.0
        stress = self._constitutive_law.get_stress(strain)
        d_stress = target_stress - stress
        num_iter = 0
        while abs(d_stress) > tol and num_iter < max_iter:
            tangent = self._constitutive_law.get_tangent(strain)
            if tangent == 0:
                raise ValueError(
                    'Tangent modulus = 0 during initial strain computation.'
                )
            d_strain = d_stress / tangent
            strain += d_strain
            stress = self._constitutive_law.get_stress(strain)
            d_stress = target_stress - stress
            num_iter += 1

        if abs(d_stress) > tol:
            raise RuntimeError('Failed to converge for given initial stress.')

        self._initial_strain = strain
