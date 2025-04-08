import typing as t

import numpy as np
from _elastic import Elastic
from numpy.typing import ArrayLike


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

    @property
    def nu(self) -> float:
        """Return Poisson's ratio.

        Note:
        Including this property allows future checks or transformations
        (e.g., ensuring 0 ≤ nu < 0.5).
        """
        return self._nu

    def get_stress(self, eps: t.Union[ArrayLike, np.ndarray]) -> np.ndarray:
        """Return a 2D stress vector [sigma_x, sigma_y, tau_xy]
        given a 2D strain vector [eps_x, epx_y, gamma_xy].
        """
        E = self._E
        nu = self.nu

        eps = np.array(eps, dtype=float).reshape(3)

        factor = E / (1 - nu**2)
        C = factor * np.array(
            [
                [1.0, nu, 0.0],
                [nu, 1.0, 0.0],
                [0.0, 0.0, (1 - nu) / 2.0],
            ]
        )
        return eps @ C

    def get_tangent(self) -> np.ndarray:
        """Return the 2D tangent stiffness matrix."""
        E = self._E
        nu = self.nu

        factor = E / (1 - nu**2)
        C = factor * np.array(
            [
                [1.0, nu, 0.0],
                [nu, 1.0, 0.0],
                [0.0, 0.0, (1 - nu) / 2.0],
            ]
        )
        return C  # noqa: RET504
