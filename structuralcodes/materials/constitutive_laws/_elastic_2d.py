import typing as t

import numpy as np
from numpy.typing import ArrayLike

from ._elastic import Elastic


class Elastic2D(Elastic):
    """Class for elastic constitutive law for 2D operations."""

    def __init__(
        self,
        E: float,
        nu: float,
        name: t.Optional[str] = None,
    ) -> None:
        """Initialize an Elastic2D Material.

        Arguments:
            E (float): The elastic modulus.
            nu (float): Poisson's ratio.
            name (str, optional): A descriptive name for the constitutive law.
        """
        super().__init__(E=E, name=name)
        self._nu = nu

        self._C_matrix: t.Optional[ArrayLike] = None

    @property
    def E(self) -> float:
        """Return the elastic modulus."""
        return self._E

    @property
    def nu(self) -> float:
        """Return Poisson's ratio.

        Note:
        Including this property allows future checks or transformations
        (e.g., ensuring 0 ≤ nu < 0.5).
        """
        return self._nu

    @property
    def C_matrix(self) -> np.ndarray:
        """Return the 2D constitutive matrix."""
        if self._C_matrix is None:
            E = self.E
            nu = self.nu

            self._C_matrix = (
                E
                / (1 - nu**2)
                * np.array(
                    [
                        [1.0, nu, 0.0],
                        [nu, 1.0, 0.0],
                        [0.0, 0.0, (1 - nu) / 2.0],
                    ]
                )
            )
        return self._C_matrix

    def get_stress(self, eps: ArrayLike) -> np.ndarray:
        """Return a 2D stress vector [sigma_x, sigma_y, tau_xy]
        given a 2D strain vector [eps_x, epx_y, gamma_xy].
        """
        eps = np.atleast_1d(eps)
        try:
            return self.C_matrix @ eps
        # Sjekk at det er en ValueError
        except ValueError as e:
            raise ValueError(
                'The input strain vector must have a length of 3.'
            ) from e

    def get_tangent(self, *args, **kwargs) -> np.ndarray:
        """Return the 2D tangent stiffness matrix."""
        del args, kwargs
        return self.C_matrix
