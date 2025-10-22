"""The core class for representing concrete smeared cracking."""

from __future__ import annotations  # To have clean hints of ArrayLike in docs

import typing as t

import numpy as np
from numpy.typing import ArrayLike

from structuralcodes.core.base import ConstitutiveLaw

from ._cracking_criterion import CrackingCriterion, NoTension
from ._poisson_reduction import PoissonReduction
from ._strength_reduction_lateral_cracking import (
    StrengthReductionLateralCracking,
)
from ._utils import (
    calculate_principal_strains,
    establish_strain_transformation_matrix,
)


class ConcreteSmearedCracking:
    """A plane stress model representing smeared cracking in reinforced
    concrete.
    """

    _initial_modulus_compression: t.Optional[float] = None
    _initial_modulus_tension: t.Optional[float] = None
    _uniaxial_compression: ConstitutiveLaw
    _strength_reduction_lateral_cracking: StrengthReductionLateralCracking
    _poisson_reduction: PoissonReduction
    _cracking_criterion: CrackingCriterion

    def __init__(
        self,
        uniaxial_compression: ConstitutiveLaw,
        strength_reduction_lateral_cracking: StrengthReductionLateralCracking,
        poisson_reduction: PoissonReduction,
        cracking_criterion: t.Optional[CrackingCriterion] = None,
    ):
        """Initialize the model."""
        self._uniaxial_compression = uniaxial_compression
        self._strength_reduction_lateral_cracking = (
            strength_reduction_lateral_cracking
        )
        self._poisson_reduction = poisson_reduction
        self._cracking_criterion = cracking_criterion or NoTension()

    def get_stress(
        self, eps: ArrayLike, return_global: bool = True
    ) -> ArrayLike:
        """Calculate the stresses based on strains.

        Arguments:
            eps (ArrayLike): The strains epsilon_x, epsilon_y and gamma_xy.

        """
        # Calculate principal strains and principal strain direction
        eps_p, phi = calculate_principal_strains(strain=eps)

        # Establish strain transformation matrix
        T = establish_strain_transformation_matrix(phi=phi)

        # Compressive-strength reduction factor due to lateral tension.
        beta = self.strength_reduction_lateral_cracking.reduction(eps_p)

        # Include the influence of Poisson's ratio on the principal strains.
        eps_pf = self.calculate_effective_principal_strains(eps_p)

        # Compute the principal stresses from the uniaxial compression law.
        sig_p = self.uniaxial_compression.get_stress(eps_pf)

        # Apply the compressive-strength reduction factor
        sig_p[sig_p < 0] *= beta

        if not return_global:
            return np.array([*sig_p, 0])

        # Transform back to global coords
        return T.T @ np.array([*sig_p, 0])

    def get_secant(
        self, eps: ArrayLike, return_global: bool = True
    ) -> np.ndarray:
        # Calculate principal strains and principal strain direction
        eps_p, phi = calculate_principal_strains(strain=eps)

        # Establish strain transformation matrix
        T = establish_strain_transformation_matrix(phi=phi)

        # Compressive-strength reduction factor due to lateral tension.
        beta = self.strength_reduction_lateral_cracking.reduction(eps_p)

        # Poisson correction
        eps_pf = self.calculate_effective_principal_strains(eps_p)

        # Establish the diagonal of material stiffness matrix
        D = self.uniaxial_compression.get_secant(eps_pf)

        # Scale down compressive stiffness with strength reduction factor
        D[eps_pf < 0] *= beta

        # Correct for the poisson effect
        C = self.poisson_reduction.poisson_matrix(
            self.cracking_criterion.cracked(eps_p=eps_p)
        ) @ np.diag(D)

        # Ensure symmetry
        C = 0.5 * (C + C.T)

        # Shear modulus
        G_value = D.mean() / (
            2
            * (
                1
                + self.poisson_reduction.nu(
                    self.cracking_criterion.cracked(eps_p=eps_p)
                )
            )
        )
        G_value = max(G_value, 1e-4 * self.initial_modulus_compression)
        G_12 = np.array([[G_value]])  # Shape (1, 1)

        Z_21 = np.zeros((2, 1))  # Shape (2, 1)

        # 3x3 secant stiffness matrix
        Cp = np.block(
            [
                [C, Z_21],
                [Z_21.T, G_12],
            ]
        )
        if not return_global:
            return Cp

        return T.T @ Cp @ T

    def get_tangent():
        pass

    def calculate_effective_principal_strains(
        self, eps_p: ArrayLike
    ) -> ArrayLike:
        """Compute the effective principal strains to include the influence of
        Poisson's ratio. Taken from 'Nonlinear Analysis of Reinforced-Concrete
        Shells' by M. A. Polak and F. J. Vecchio (1993).
        """
        return (
            self.poisson_reduction.poisson_matrix(
                cracked=self.cracking_criterion.cracked(eps_p=eps_p)
            )
            @ eps_p
        )

    @property
    def uniaxial_compression(self) -> ConstitutiveLaw:
        """Return the constitutive law for uniaxial compression."""
        return self._uniaxial_compression

    @property
    def strength_reduction_lateral_cracking(
        self,
    ) -> StrengthReductionLateralCracking:
        """Return the model for strength reduction due to lateral cracking."""
        return self._strength_reduction_lateral_cracking

    @property
    def poisson_reduction(self) -> PoissonReduction:
        """Return the model for reduction of Poisson ratio due to cracking."""
        return self._poisson_reduction

    @property
    def cracking_criterion(self) -> CrackingCriterion:
        """Return the cracking criterion."""
        return self._cracking_criterion

    @property
    def initial_modulus_compression(self) -> float:
        """Return the initial modulus in uniaxial compression."""
        if self._initial_modulus_compression is None:
            self._initial_modulus_compression = (
                self.uniaxial_compression.get_tangent(0.0)
            )

        return self._initial_modulus_compression
