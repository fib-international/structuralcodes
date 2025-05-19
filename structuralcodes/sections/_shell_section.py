import typing as t

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.linalg import lu_factor, lu_solve

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
    mesh_size: float

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
        self.layers: t.Optional[t.Tuple] = None

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
        result = self.integrator.integrate_strain_response_on_geometry(
            geo=self.section.geometry,
            strain=strain,
            integrate=integrate,
            mesh_size=self.mesh_size,
            layers=self.layers,
        )

        # Save layers for future use
        if self.layers is None and result[-1] is not None:
            self.layers = result[-1]

        # Return the results without layers
        if len(result) == 2:
            # Return only the stiffness matrix
            return result[0]

        # Return only the stress resultants
        return result[:-1]

    def calculate_strain_profile(
        self,
        nx: float,
        ny: float,
        nxy: float,
        mx: float,
        my: float,
        mxy: float,
        initial: bool = False,
        max_iter: int = 10,
        tol: float = 1e-6,
    ) -> t.List[float]:
        """Computes the strain profile for a given set of stress resultants.

        Args:
            nx (float): Membrane force in x-direction.
            ny (float): Membrane force in y-direction.
            nxy (float): Membrane shear.
            mx (float): Plate moment giving stresses in x-direction.
            my (float): Plate moment giving stresses in y-direction.
            mxy (float): Plate moment giving shear stresses in the xy-plane.

        Returns:
            list(float): The strain response eps_x, eps_y, gamma_xy, chi_x,
            chi_y and chi_xy.
        """
        # Get the gometry
        geom = self.section.geometry

        # Collect loads in a numpy array
        loads = np.array([nx, ny, nxy, mx, my, mxy])

        # Compute initial tangent stiffness matrix
        stiffness, layers = (
            self.integrator.integrate_strain_response_on_geometry(
                geom,
                [0, 0, 0, 0, 0, 0],
                integrate='modulus',
                mesh_size=self.mesh_size,
                layers=self.layers,
            )
        )

        # Save layers if needed
        if self.layers is None and layers is not None:
            self.layers = layers

        # Calculate strain plane with Newton Rhapson Iterative method
        num_iter = 0
        strain = np.zeros(6)

        # Factorize once the stiffness matrix if using initial
        if initial:
            # LU factorization
            lu, piv = lu_factor(stiffness)

        # Do Newton loops
        while True:
            # Check if number of iterations exceeds the maximum
            if num_iter > max_iter:
                break

            # Calculate response and residuals
            response = np.array(self.integrate_strain_profile(strain=strain))
            residual = loads - response

            if initial:
                # Solve using the decomposed matrix
                delta_strain = lu_solve((lu, piv), residual)
            else:
                # Calculate the current stiffness
                stiffness, _ = (
                    self.integrator.integrate_strain_response_on_geometry(
                        geom,
                        strain,
                        integrate='modulus',
                        layers=self.layers,
                    )
                )

                # Solve using the current stiffness
                delta_strain = np.linalg.solve(stiffness, residual)

            # Update the strain
            strain += delta_strain

            num_iter += 1

            # Check for convergence:
            if np.linalg.norm(delta_strain) < tol:
                break

        if num_iter >= max_iter:
            raise StopIteration('Maximum number of iterations reached.')

        return strain.tolist()
