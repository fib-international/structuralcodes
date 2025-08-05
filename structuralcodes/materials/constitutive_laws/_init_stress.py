"""Initial strain constitutive law."""

from __future__ import annotations  # To have clean hints of ArrayLike in docs

import typing as t

from numpy.typing import ArrayLike

from ...core.base import ConstitutiveLaw


class InitStress(ConstitutiveLaw):
    """Class for initial stress Constitutive Law."""

    __materials__: t.Tuple[str] = (
        'steel',
        'rebars',
        'concrete',
    )

    def __init__(
        self,
        constitutive_law: ConstitutiveLaw,
        initial_stress: float,
        name: t.Optional[str] = None,
    ) -> None:
        """Initialize an Initial Stress Material.

        This constitutive law is a wrapper for another constitutive law
        that assigns an initial stress.

        Arguments:
            constitutive_law (ConstitutiveLaw): Wrapped constitutive law.
            initial_stress (float): The initial stress to be applied.
        """
        name = name if name is not None else 'InitStressLaw'
        super().__init__(name=name)
        if not isinstance(constitutive_law, ConstitutiveLaw):
            raise TypeError(
                f'Expected a ConstitutiveLaw instance, '
                f'got {type(constitutive_law)}'
            )
        self._wrapped_law = constitutive_law

        # Iteratively compute the initial strain that gives the desired
        # initial stress. Note that the wrapped law can be nonlinear
        TOL = 1e-12
        MAX_ITER = 100
        target_stress = initial_stress
        strain = 0.0
        stress = self._wrapped_law.get_stress(strain)
        d_stress = target_stress - stress
        count = 0
        while abs(d_stress) > TOL and count < MAX_ITER:
            tangent = self._wrapped_law.get_tangent(strain)
            if tangent == 0:
                raise ValueError(
                    'Tangent modulus = 0 during initial strain computation.'
                )
            d_strain = d_stress / tangent
            strain += d_strain
            stress = self._wrapped_law.get_stress(strain)
            d_stress = target_stress - stress
            count += 1

        if abs(d_stress) > TOL:
            raise RuntimeError('Failed to converge for given initial stress.')

        self._initial_strain = strain

    def get_stress(
        self, eps: t.Union[float, ArrayLike]
    ) -> t.Union[float, ArrayLike]:
        """Return the stress given strain."""
        return self._wrapped_law.get_stress(eps + self._initial_strain)

    def get_tangent(
        self, eps: t.Union[float, ArrayLike]
    ) -> t.Union[float, ArrayLike]:
        """Return the tangent for given strain."""
        return self._wrapped_law.get_tangent(eps + self._initial_strain)

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
        return self._wrapped_law.__marin__(
            strain=[strain[0] + self._initial_strain, strain[1]]
        )

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
        return self._wrapped_law.__marin_tangent__(
            strain=[strain[0] + self._initial_strain, strain[1]]
        )

    def get_ultimate_strain(
        self, yielding: bool = False
    ) -> t.Tuple[float, float]:
        """Return the ultimate strain (negative and positive)."""
        ult_strain = self._wrapped_law.get_ultimate_strain(yielding=yielding)
        return (
            ult_strain[0] - self._initial_strain,
            ult_strain[1] - self._initial_strain,
        )
