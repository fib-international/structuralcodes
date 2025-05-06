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

    def __init__(
        self,
        fc: float,
        eps_0: float = -0.002,
        eps_u: float = -0.0035,
        n: float = 2.0,
        nu: float = 0.2,
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
        """
        super().__init__(fc=fc, eps_0=eps_0, eps_u=eps_u, n=n, name=name)
        self._nu = nu

    @property
    def E_c(self) -> float:
        """Return the elastic modulus of concrete."""
        return 2 * self._fc / self._eps_0

    @property
    def f_cr(self) -> float:
        """Return the cracking strength of concrete."""
        return -0.33 * np.sqrt(-self._fc)

    @property
    def eps_cr(self) -> float:
        """Return the cracking strain of concrete."""
        return self.f_cr / self.E_c

    def transform(self, strain: ArrayLike) -> np.ndarray:
        """Transform the strain vector to principal directions."""
        eps = np.atleast_1d(strain)
        # Use arctan2 to obtain the angle in the correct quadrant
        theta = 0.5 * np.arctan2(2 * eps[2], eps[0] - eps[1])
        c = np.cos(theta)
        s = np.sin(theta)
        T = np.array(
            [
                [c * c, s * s, s * c],
                [s * s, c * c, -s * c],
                [-2 * s * c, 2 * s * c, c * c - s * s],
            ]
        )
        return T, T @ eps

    def get_effective_principal_strains(self, strain: ArrayLike) -> np.ndarray:
        """Compute the effective principal strains. Taken from 'Nonlinear
        Analysis of Reinforced-Concrete Shells' by M. A. Polak and F. J.
        Vecchio [eq. 5].
        """
        _, eps_p = self.transform(strain)

        # Effective principal strains accounting for Poisson's ratio
        nu = self._nu
        eps_1f = ((1 - nu) * eps_p[0] + nu * eps_p[1] + nu * eps_p[2]) / (
            (1 + nu) * (1 - 2 * nu)
        )
        eps_2f = ((1 - nu) * eps_p[1] + nu * eps_p[0] + nu * eps_p[2]) / (
            (1 + nu) * (1 - 2 * nu)
        )
        eps_3f = ((1 - nu) * eps_p[2] + nu * eps_p[0] + nu * eps_p[1]) / (
            (1 + nu) * (1 - 2 * nu)
        )

        eps_pf = np.array([eps_1f, eps_2f, eps_3f])

        # Check if the material is cracked
        cracked = np.any(eps_pf < self.eps_cr)

        # If cracked, Poisson's ratio (nu) is neglected
        if cracked:
            eps_pf = eps_p
            print("Material is cracked, Poisson's ratio is neglected.")
        else:
            # If not cracked, use the effective principal strains (eps_pf)
            print(
                'Material is not cracked, using effective principal strains.'
            )

        return eps_pf

    def get_stress(self, strain: ArrayLike) -> np.ndarray:
        """Return a 2D stress vector [sigma_x, sigma_y, tau_xy]
        given a 2D strain vector [eps_x, epx_y, gamma_xy].
        """
        T, _ = self.transform(strain)
        eps_pf = self.get_effective_principal_strains(strain)

        sig_p = super().get_stress(eps_pf)

        # Rotate back to global coords
        return T.T @ sig_p

    def get_secant(self, strain: ArrayLike) -> np.ndarray:
        """Compute the 3x3 secant stiffness matrix C."""
        T, _ = self.transform(strain)
        eps_pf = self.get_effective_principal_strains(strain)

        sig_p = super().get_stress(eps_pf)

        E11 = sig_p[0] / eps_pf[0]
        E22 = sig_p[1] / eps_pf[1]
        E12 = (E11 + E22) / 2

        Cp = np.array(
            [
                [E11, 0.0, 0.0],
                [0.0, E22, 0.0],
                [0.0, 0.0, 0.5 * E12],
            ]
        )

        return T.T @ Cp @ T
