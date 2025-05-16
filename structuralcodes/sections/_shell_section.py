import typing as t

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ..core.base import Section, SectionCalculator
from ..geometry import ShellGeometry
from .section_integrators import ShellFiberIntegrator


class ShellSection(Section):
    """This is the implementation of the shell class section.

    Attributes:
        geometry ShellGeometry: The geometry of
            the section.
        name (str): The name of the section.
    """

    def __init__(
        self,
        geometry: ShellGeometry,
        name: t.Optional[str] = None,
        **kwargs,
    ) -> None:
        """Initialize a ShellSection.

        Note:
            The ShellSection uses a ShellSectionCalculator for all
            calculations. The ShellSectionCalculator uses a
            ShellFiberIntegrator for integrating over the thickness.
        """
        if name is None:
            name = 'ShellSection'
        super().__init__(name)
        self._geometry = geometry
        self.section_calculator = ShellSectionCalculator(
            section=self, **kwargs
        )

    @property
    def geometry(self):
        """Return the geometry of the section."""
        return self._geometry


class ShellSectionCalculator(SectionCalculator):
    """A calculator for shell sections."""

    section: ShellSection

    def __init__(
        self,
        section: ShellSection,
        **kwargs,
    ) -> None:
        """Initialize the ShellSectionCalculator.

        Arguments:
            section (ShellSection): The section object.
        """
        super().__init__(section=section)
        self.integrator = ShellFiberIntegrator()
        self.mesh_size = kwargs.get('mesh_size', 0.01)
        self.z_coords = None

    def _calculate_gross_section_properties(self):
        pass

    def calculate_bending_strength(self):
        pass

    def calculate_moment_curvature(self):
        pass

    def integrate_strain_profile(
        self,
        strain: ArrayLike,
        integrate: t.Literal['stress', 'modulus'] = 'stress',
        **kwargs,
    ) -> t.Union[t.Tuple[float, float, float, float, float, float], NDArray]:
        """Integrate a strain profile returning stress resultants or tangent
        section stiffness matrix.

        Arguments:
            strain (ArrayLike): Represents the deformation plane. The strain
                should have six entries representing respectively: eps_x,
                eps_y, gamma_xy, kappa_x, kappa_y, kappa_xy.
            integrate (str): a string indicating the quantity to integrate over
                the section. It can be 'stress' or 'modulus'. When 'stress'
                is selected, the return value will be the stress resultants Nx,
                Ny, Nxy, Mx, My and Mxy, while if 'modulus' is selected, the
                return will be the tangent section stiffness matrix (default is
                'stress').

        Returns:
            Union(Tuple(float, float, float, float, float, float),NDArray):
            Nx, Ny, Nxy, Mx, My and Mxy My when `integrate='stress'`, or a
            numpy array representing the stiffness matrix then
            `integrate='modulus'`.

        Examples:
            result = self.integrate_strain_profile(strain,integrate='modulus')
            # `result` will be the tangent stiffness matrix (a 6x6 numpy array)

            result = self.integrate_strain_profile(strain)
            # `result` will be a tuple containing section forces (Nx, Ny, Nxy,
               Mx, My, Mxy)

        Raises:
            ValueError: If a unkown value is passed to the `integrate`
            parameter.
        """
        return self.integrator.integrate_strain_response_on_geometry(
            geo=self.section.geometry,
            strain=strain,
            integrate=integrate,
            **kwargs,
        )

    def calculate_strain_profile(
        self,
        nx: float,
        ny: float,
        nxy: float,
        mx: float,
        my: float,
        mxy: float,
    ) -> t.List[float]:
        """Computes the strain profile for a given set of stress resultants.

        The strain profile is computed by solving:

            eps = K^(-1) * R

        where:
        - K is the integrated concrete stiffness (6x6) over the thickness plus
            the sum of the reinforcement contributions (6x6),
            computed from each reinforcement layer.

        Parameters:
        nx, ny, nxy, mx, my, mxy (float): The stress resultants.

        Returns:
        NDArray[float]: 6 floats: eps_x, eps_y, gamma_xy, kappa_x, kappa_y and
        kappa_xy.
        """
        stress = nx, ny, nxy, mx, my, mxy

        geo = self.section.geometry
        integrator = self.integrator

        strain = np.zeros(6)
        prepared_input, _ = integrator.prepare_input(
            geo=geo, strain=strain, integrate='modulus'
        )
        K = integrator.integrate_modulus(prepared_input)

        # Solve for the strain profile using the pseudoinverse
        return np.linalg.pinv(K) @ stress
