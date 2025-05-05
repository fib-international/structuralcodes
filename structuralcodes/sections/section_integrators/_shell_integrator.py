"""A concrete implementation of a shell integrator using fiber
discretization.
"""

import math
import typing as t

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ...geometry._shell_geometry import ShellGeometry
from ._section_integrator import SectionIntegrator


class ShellFiberIntegrator(SectionIntegrator):
    """A concrete implementation of a fiber integrator for shell sections."""

    def prepare_input(
        self,
        geo: ShellGeometry,
        strain: ArrayLike,
        integrate: t.Literal['stress', 'modulus'] = 'stress',
        **kwargs,
    ) -> t.Tuple[t.List[t.Tuple[np.ndarray, np.ndarray]], None]:
        """Prepare general input to the integration of stress or material
        modulus in the shell section.

        Calculate the stress resultants or secant section stiffness based on
        strains at discrete fibers along the shell thickness.

        Arguments:
            geo (ShellGeometry): The shell geometry object.
            strain (ArrayLike): The strains and curvatures of the shell section
                [eps_x, eps_y, gamma_xy, kappa_x, kappa_y, kappa_xy].
            integrate (str): A string indicating the quantity to integrate over
                the shell section. It can be 'stress' or 'modulus'. When
                'stress' is selected, the return value will be the generalised
                stress resultants: membrane forces and bending moments Nx, Ny,
                Nxy, Mx, My, Mxy. When 'modulus' is selected, the return value
                will be the section stiffness matrix K (6x6), which relates
                generalised strains to stress resultants (default is 'stress').

        Keyword Arguments:
            z_coords (ArrayLike): The z-coordinates of the layers in the shell.
            mesh_size (float): fraction of the total shell thickness for each
                layer ([0,1]). Default is 0.01.


        Returns:
            Tuple: (prepared_input, z_coords)
        """
        z_coords = kwargs.get('z_coords', None)

        t_total = geo.thickness

        if z_coords is None:
            mesh_size = kwargs.get('mesh_size', 0.01)
            if not (0 < mesh_size <= 1):
                raise ValueError('mesh_size must be [0,1].')
            n_layers = max(1, math.ceil(1 / mesh_size))
            dz = t_total / n_layers
            z_coords = np.linspace(
                -t_total / 2 + dz / 2, t_total / 2 - dz / 2, n_layers
            )

        material = geo.material

        prepared_input = []

        IA = []

        for z in z_coords:
            fiber_strain = strain[:3] + z * strain[3:]

            if integrate == 'stress':
                integrand = material.get_stress(fiber_strain)
            elif integrate == 'modulus':
                integrand = material.get_secant(fiber_strain)
            else:
                raise ValueError(f'Unknown integrate type: {integrate}')

            IA.append(integrand * dz)

        prepared_input = [(np.array(z_coords), np.array(IA))]

        return prepared_input, z_coords

    def integrate_stress(
        self,
        prepared_input: t.List[t.Tuple[np.ndarray, np.ndarray]],
    ) -> t.Tuple[float, float, float, float, float, float]:
        """Integrate stresses over the shell thickness.

        Arguments:
            prepared_input (List): The prepared input from .prepare_input().

        Returns:
            Tuple(float, float, float, float, float, float):
            The stress resultants Nx, Ny, Nxy, Mx, My, Mxy
        """
        z, stress_resultants = prepared_input[0]

        Nx = np.sum(stress_resultants[:, 0])
        Ny = np.sum(stress_resultants[:, 1])
        Nxy = np.sum(stress_resultants[:, 2])

        Mx = np.sum(stress_resultants[:, 0] * z)
        My = np.sum(stress_resultants[:, 1] * z)
        Mxy = np.sum(stress_resultants[:, 2] * z)

        return Nx, Ny, Nxy, Mx, My, Mxy

    def integrate_modulus(
        self,
        prepared_input: t.List[t.Tuple[np.ndarray, np.ndarray]],
    ) -> NDArray:
        """Integrate material modulus over shell thickness to obtain shell
        section stiffness.

        Arguments:
            prepared_input (List): The prepared input from .prepare_input().

        Returns:
            NDArray: Section stiffness matrix (6x6).
        """
        z, MA = prepared_input[0]

        A = np.zeros((3, 3))
        B = np.zeros((3, 3))
        D = np.zeros((3, 3))

        for C_dz, z_i in zip(MA, z):
            A += C_dz
            B += z_i * C_dz
            D += z_i**2 * C_dz

        return np.block([[A, B], [B, D]])

    def integrate_strain_response_on_geometry(
        self,
        geo: ShellGeometry,
        strain: ArrayLike,
        integrate: t.Literal['stress', 'modulus'] = 'stress',
        **kwargs,
    ) -> t.Union[
        t.Tuple[float, float, float, float, float, float], np.ndarray
    ]:
        """Integrate strain profile through the shell thickness.

        This method evaluates the response of a shell section subjected to a
        generalised strain state, either by integrating the stresses to obtain
        resultant forces and moments, or by integrating the secant stiffness
        to obtain the section stiffness matrix.

        Arguments:
            geo (ShellGeometry): The shell geometry including thickness and
            material.
            strain (ArrayLike): Generalised strain vector of the form:
                [eps_x, eps_y, gamma_xy, kappa_x, kappa_y, kappa_xy]
            integrate (str): Type of integration to perform. Must be one of:
                - 'stress': returns the stress resultants (Nx, Ny, Nxy, Mx, My,
                Mxy)
                - 'modulus': returns the section stiffness matrix (6x6)

        Returns:
            Union[Tuple[float, float, float, float, float, float], np.ndarray]:
                - If integrate = 'stress': returns 6 stress resultants.
                - If integrate = 'modulus': returns the 6x6 stiffness matrix.

        Example:
            result = integrate_strain_response_on_geometry(geo, strain,
            integrate='modulus')
            'result will be the 6x6 stiffness matrix of the shell.

        Raises:
            ValueError: If `integrate` is not 'stress' or 'modulus'.
        """
        # Prepare the general input based on the geometry and the input strains
        prepared_input, _ = self.prepare_input(
            geo=geo, strain=strain, integrate=integrate, **kwargs
        )

        # Return the calculated response
        if integrate == 'stress':
            return self.integrate_stress(prepared_input)
        if integrate == 'modulus':
            return self.integrate_modulus(prepared_input)
        raise ValueError(f'Unknown integrate type: {integrate}')
