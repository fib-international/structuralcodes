"""A concrete implementation of a shell integrator."""

import typing as t

import numpy as np
from _section_integrator import SectionIntegrator
from numpy.typing import ArrayLike, NDArray

from ...geometry._shell_geometry import ShellGeometry


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

        Calculate the stress resultants or tangent section stiffness based on
        strains at discrete fibers along the shell thickness.

        Arguments:
            geo (ShellGeometry): The shell geometry object.
            strain (ArrayLike): The strains and curvatures of the shell
                section, given in the format (eps_x, eps_y, gamma_xy,
                kappa_x, kappa_y, kappa_xy).
            integrate (str): A string indicating the quantity to integrate over
                the shell section. It can be 'stress' or 'modulus'. When
                'stress' is selected, the return value will be the generalised
                stress resultants: membrane forces and bending moments Nx, Ny,
                Nxy, Mx, My, Mxy. When 'modulus' is selected, the return value
                will be the section stiffness matrix K (6x6), which relates
                generalised strains to stress resultants (default is 'stress').

        Keyword Arguments:
            layer_thickness (float): The thickness of each layer in
            the shell section in mm (default is 0.1 * (total shell thickness)).

        Returns:
        Tuple: (prepared input, triangulation_data=None)
        """
        eps_x, eps_y, gamma_xy, kappa_x, kappa_y, kappa_xy = strain

        t_total = geo.thickness
        t_layer = kwargs.get('layer_thickness', 0.1 * t_total)
        n_layers = max(1, int(round(t_total / t_layer)))
        dz = t_total / n_layers
        z_coords = np.linspace(
            -t_total / 2 + dz / 2, t_total / 2 - dz / 2, n_layers
        )

        material = geo.material

        prepared_input = []

        z_values = []
        IA = []

        for z in z_coords:
            fiber_strain = np.array(
                [
                    eps_x + z * kappa_x,
                    eps_y + z * kappa_y,
                    gamma_xy + z * kappa_xy,
                ]
            )

            if integrate == 'stress':
                integrand = material.get_stress(fiber_strain)
            elif integrate == 'modulus':
                integrand = material.get_tangent()
            else:
                raise ValueError(f'Unknown integrate type: {integrate}')

            z_values.append(z)
            IA.append(integrand.squeeze() * dz)

        prepared_input = [(np.array(z_values), np.array(IA))]

        return prepared_input, None

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

        for i in range(len(z)):
            C_dz = MA[i]
            z_i = z[i]

            A += C_dz
            B += z_i * C_dz
            D += z_i**2 * C_dz

        stiffness = np.block([[A, B], [B, D]])

        return stiffness  # noqa: RET504

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
        resultant forces and moments, or by integrating the tangent stiffness
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
