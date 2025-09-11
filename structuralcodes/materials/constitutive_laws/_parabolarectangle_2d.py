"""Parabola-Rectangle 2D constitutive law."""

from __future__ import annotations  # To have clean hints of ArrayLike in docs

import typing as t

import numpy as np
from numpy.typing import ArrayLike

from ._parabolarectangle import ParabolaRectangle


class ParabolaRectangle2D(ParabolaRectangle):
    """Class for parabola rectangle constitutive law in 2D.

    The stresses and strains are assumed negative in compression and positive
    in tension.
    """

    __materials__: t.Tuple[str] = ('concrete',)
    _poisson_matrix_cracked: t.Optional[ArrayLike] = None
    _poisson_matrix_uncracked: t.Optional[ArrayLike] = None

    def __init__(
        self,
        fc: float,
        eps_0: float = -0.002,
        eps_u: float = -0.0035,
        n: float = 2.0,
        nu: float = 0.2,
        c_1: float = 0.8,
        c_2: float = 100,
        name: t.Optional[str] = None,
    ) -> None:
        """Initialize a Parabola-Rectangle 2D Material.

        Arguments:
            fc (float): The strength of concrete in compression.

        Keyword Arguments:
            eps_0 (float): Peak strain of concrete in compression. Default
                value = -0.002.
            eps_u (float): Ultimate strain of concrete in compression. Default
                value = -0.0035.
            n (float): Exponent for the pre-peak branch. Default value = 2.
            name (str): A name for the constitutive law.
            nu (float): Poisson's ratio. Default value = 0.2.
            c_1 (float): First coefficient for the compressive-strength
                reduction factor due to lateral tension. Default value = 0.8.
            c_2 (float): Second coefficient for the compressive-strength
                reduction factor due to lateral tension. Default value = 100.
        """
        super().__init__(fc=fc, eps_0=eps_0, eps_u=eps_u, n=n, name=name)
        self._nu = nu
        self._c_1 = c_1
        self._c_2 = c_2

    @property
    def c_1(self) -> float:
        """Return the first coefficient for the compressive-strength reduction
        factor due to lateral tension. Default value = 0.8.
        """
        return self._c_1

    @property
    def c_2(self) -> float:
        """Return the second coefficient for the compressive-strength reduction
        factor due to lateral tension. Default value = 100.
        """
        return self._c_2

    def poisson_matrix(self, cracked: bool) -> np.ndarray:
        """Return the Poisson matrix."""
        if cracked:
            return self.poisson_matrix_cracked
        return self.poisson_matrix_uncracked

    def nu(self, cracked: bool) -> float:
        """Return the Poisson ratio."""
        if cracked:
            return 0.0
        return self._nu

    @property
    def poisson_matrix_uncracked(self) -> ArrayLike:
        """Return the uncracked Poisson matrix."""
        if self._poisson_matrix_uncracked is None:
            self._poisson_matrix_uncracked = (
                1 / (1 - self._nu**2)
            ) * np.array(
                [
                    [1, self._nu],
                    [self._nu, 1],
                ]
            )
        return self._poisson_matrix_uncracked

    @property
    def poisson_matrix_cracked(self) -> ArrayLike:
        """Return the cracked Poisson matrix."""
        if self._poisson_matrix_cracked is None:
            self._poisson_matrix_cracked = np.eye(2)
        return self._poisson_matrix_cracked

    def get_effective_principal_strains(self, eps_p: ArrayLike) -> np.ndarray:
        """Compute the effective principal strains to include the influence of
        Poisson's ratio. Taken from 'Nonlinear Analysis of Reinforced-Concrete
        Shells' by M. A. Polak and F. J. Vecchio (1993).
        """
        return self.poisson_matrix(self.check_cracked(eps_p=eps_p)) @ eps_p

    def check_cracked(self, eps_p: ArrayLike) -> bool:
        """Check if the concrete is cracked. Returns True if any principal
        strain (eps_p1, eps_p2) is greater than zero.
        """
        return np.any(eps_p > 0)

    def strength_reduction_lateral_cracking(self, eps_p: ArrayLike) -> float:
        """Return the compressive-strength reduction factor due to lateral
        tension. This relation comes from Vecchio & Collins (1986), "The
        Modified Compression-Field Theory for Reinforced Concrete Elements
        Subjected to Shear.
        """
        if not self.check_cracked(eps_p):
            return 1

        beta = 1 / (self.c_1 + self.c_2 * max(eps_p))

        return min(beta, 1)

    def get_stress(self, strain: ArrayLike) -> np.ndarray:
        """Return a 2D stress vector [sigma_x, sigma_y, tau_xy]
        given a 2D strain vector [eps_x, epx_y, gamma_xy].
        """
        # Calculate principal strains and principal strain direction
        eps_p, phi = calculate_principal_strains(strain=strain)

        # Establish strain transformation matrix
        T = establish_strain_transformation_matrix(phi=phi)

        # Compressive-strength reduction factor due to lateral tension.
        beta = self.strength_reduction_lateral_cracking(eps_p)

        # Include the influence of Poisson's ratio on the principal strains.
        eps_pf = self.get_effective_principal_strains(eps_p)

        # Compute the principal stresses from 1D parabola-rectangle law.
        sig_p = super().get_stress(eps_pf)

        # Apply the compressive-strength reduction factor
        sig_p[sig_p < 0] *= beta

        # Transform back to global coords
        return T.T @ np.array([*sig_p, 0])

    def get_secant(self, strain: ArrayLike) -> np.ndarray:
        """Compute the 3x3 secant stiffness matrix C."""
        # Calculate principal strains and principal strain direction
        eps_p, phi = calculate_principal_strains(strain=strain)

        # Establish strain transformation matrix
        T = establish_strain_transformation_matrix(phi=phi)

        # Compressive-strength reduction factor due to lateral tension.
        beta = self.strength_reduction_lateral_cracking(eps_p)

        # Poisson correction
        eps_pf = self.get_effective_principal_strains(eps_p)

        sig_p = super().get_stress(eps_pf)
        sig_p[sig_p < 0] *= beta

        # Establish the diagonal of material stiffness matrix
        # First set the diagonal terms equal to the initial stiffness
        # Then calculate the secant stiffness for components with nonzero
        # strains
        D = self._fc * self._n / self._eps_0 * np.ones_like(sig_p)
        nonzero_strains = np.abs(eps_pf) > 0
        D[nonzero_strains] = sig_p[nonzero_strains] / eps_pf[nonzero_strains]

        # Correct for the poisson effect
        C = self.poisson_matrix(self.check_cracked(eps_p=eps_p)) @ np.diag(D)

        # Ensure symmetry
        C = 0.5 * (C + C.T)

        # Shear modulus
        G_12 = np.array(
            [[D.mean() / (2 * (1 + self.nu(self.check_cracked(eps_p=eps_p))))]]
        )  # Shape (1, 1)

        Z_21 = np.zeros((2, 1))  # Shape (2, 1)

        # 3x3 secant stiffness matrix
        Cp = np.block(
            [
                [C, Z_21],
                [Z_21.T, G_12],
            ]
        )

        return T.T @ Cp @ T


def calculate_principal_strains(
    strain: ArrayLike,
) -> t.Tuple[ArrayLike, float]:
    """Calculate the principal strains in 2D."""
    eps = np.atleast_1d(strain)
    # Use arctan2 to obtain the angle in the correct quadrant
    phi = 0.5 * np.arctan2(eps[2], eps[0] - eps[1])
    eps_p = np.array(
        [
            (eps[0] + eps[1]) / 2
            + (((eps[0] - eps[1]) / 2) ** 2 + 0.25 * eps[2] ** 2) ** 0.5,
            (eps[0] + eps[1]) / 2
            - (((eps[0] - eps[1]) / 2) ** 2 + 0.25 * eps[2] ** 2) ** 0.5,
        ]
    )

    return eps_p, phi


def establish_strain_transformation_matrix(phi: float) -> ArrayLike:
    """Establish the strain transformation matrix for a given angle."""
    c = np.cos(phi)
    s = np.sin(phi)
    return np.array(
        [
            [c * c, s * s, s * c],
            [s * s, c * c, -s * c],
            [-2 * s * c, 2 * s * c, c * c - s * s],
        ]
    )
