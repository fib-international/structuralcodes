""""Generic class section implemenetation."""

from __future__ import annotations

import typing as t

from math import sin, cos
import numpy as np

from structuralcodes.sections._geometry import (
    PointGeometry,
    SurfaceGeometry,
    CompoundGeometry,
)
from structuralcodes.core.base import Section, SectionCalculator
import structuralcodes.core._section_results as s_res
from ._section_integrators import integrator_factory


class GenericSection(Section):
    """This is the implementation of the generic class section."""

    def __init__(
        self,
        geometry: SurfaceGeometry | CompoundGeometry,
        name: t.Optional[str] = None,
        **kwargs,
    ) -> None:
        if name is None:
            name = 'GenericSection'
        super().__init__(name)
        self.geometry = geometry
        self.section_analyzer = GenericSectionCalculator(self, **kwargs)
        self.gross_properties = (
            self.section_analyzer._calculate_gross_section_properties()
        )


class GenericSectionCalculator(SectionCalculator):
    """Calculator class implementing analysis algortims for
    Code checks.
    """

    def __init__(
        self, sec: GenericSection, integrator: str = 'Marin', **kwargs
    ) -> None:
        """Initialize the GenericSectionCalculator
        Input:
        section (SectionCalculator): the section object.
        """
        super().__init__(section=sec)
        # Select the integrator if specified
        self.integrator = integrator_factory(integrator)()
        # Mesh size used for Fibre integrator
        self.mesh_size = kwargs.get('mesh_size', 0.001)
        # triangulated_data used for Fibre integrator
        self.triangulated_data = None

    def _calculate_gross_section_properties(self) -> s_res.GrossProperties:
        """Calculates the gross section properties of the GenericSection
        This function is private and called when the section is created
        It stores the result into the result object.

        Returns:
        gross_section_properties (GrossSection)
        """
        # It will use the algorithms for generic sections
        gp = s_res.GrossProperties()
        # Here compute area, centroid and other stuff from gross properties
        # 1. Computation of area
        gp.area = self.section.geometry.area
        return gp

    def get_balanced_failure_strain(
        self, geom: CompoundGeometry, yielding: bool = False
    ) -> [float, float, float]:
        """Returns the strain profile corresponding to balanced failure.
        This is found from all ultimate strains for all materials, checking
        the minimum value of curvature.
        It returns a tuple with:
        1. Value of y coordinate for negative failure
        2. Value of y coordinate for positive failure
        3. Strain profile as a list with three values: axial strain, curvature
           y*, curvature z* (assumed zero since in the rotated frame y*z* it is
           a case of uniaxial bending).
        """
        chi_min = 1e10
        for g in geom.geometries + geom.point_geometries:
            for other_g in geom.geometries + geom.point_geometries:
                if g != other_g:
                    eps_p = g.material.get_ultimate_strain(yielding)[0]
                    if isinstance(g, SurfaceGeometry):
                        y_p = g.polygon.bounds[1]
                    elif isinstance(g, PointGeometry):
                        y_p = g._point.coords[0][1]
                    eps_n = other_g.material.get_ultimate_strain(yielding)[1]
                    if isinstance(other_g, SurfaceGeometry):
                        y_n = other_g.polygon.bounds[3]
                    elif isinstance(other_g, PointGeometry):
                        y_n = other_g._point.coords[0][1]
                    if y_p >= y_n:
                        continue
                    chi = -(eps_p - eps_n) / (y_p - y_n)
                    # print(y_p,eps_p,y_n,eps_n,chi)
                    if chi < chi_min:
                        chi_min = chi
                        eps_0 = eps_n + chi_min * y_n
                        y_n_min = y_n
                        y_p_min = y_p
        y_p, y_n = y_p_min, y_n_min
        strain = [eps_0, chi_min, 0]
        return (y_n, y_p, strain)

    def find_equilibrium_fixed_pivot(
        self, geom: CompoundGeometry, n: float, yielding: bool = False
    ) -> [float, float, float]:
        """Find the equilibrium changing curvature fixed a pivot.
        The algorithm uses bisection algorithm between curvature
        of balanced failure and 0. Selected the pivot point as
        the top or the bottom one, the neutral axis is lowered or
        raised respectively.

        Arguments:
        geom: (CompoundGeometry) A geometry in the rotated reference system
        n: (float) value of external axial force needed to be equilibrated

        Returns:
        strain: list of 3 floats. Axial strain at (0,0) curvatures of y* and z*
            axes. Note that being uniaxial bending, curvature along z* is 0.0
        """
        # Number of maximum iteration for the bisection algorithm
        ITMAX = 100
        # 1. Start with a balanced failure: this is found from all ultimate
        # strains for all materials, checking the minimum curvature value
        y_n, y_p, strain = self.get_balanced_failure_strain(geom, yielding)
        eps_p = strain[0] - strain[1] * y_p
        eps_n = strain[0] - strain[1] * y_n
        # Integrate this strain profile corresponding to balanced failure
        (
            n_int,
            _,
            _,
            tri,
        ) = self.integrator.integrate_strain_response_on_geometry(geom, strain)
        self.triangulated_data = tri
        # Check if there is equilibrium with this strain distribution
        chi_a = strain[1]
        dn_a = n_int - n
        chi_b = 1e-13
        if n_int < n:
            pivot = y_p
            strain_pivot = eps_p
        else:
            # Too much tension, lower NA
            pivot = y_n
            strain_pivot = eps_n
        eps_0 = strain_pivot + chi_b * pivot
        n_int, _, _, _ = self.integrator.integrate_strain_response_on_geometry(
            geom, [eps_0, chi_b, 0], tri=self.triangulated_data
        )
        dn_b = n_int - n
        it = 0
        while (abs(dn_a - dn_b) > 1e-2) and (it < ITMAX):
            chi_c = (chi_a + chi_b) / 2.0
            eps_0 = strain_pivot + chi_c * pivot
            (
                n_int,
                _,
                _,
                _,
            ) = self.integrator.integrate_strain_response_on_geometry(
                geom, [eps_0, chi_c, 0], tri=self.triangulated_data
            )
            dn_c = n_int - n
            if dn_c * dn_a < 0:
                chi_b = chi_c
                dn_b = dn_c
            else:
                chi_a = chi_c
                dn_a = dn_c
            it += 1
        if it >= ITMAX:
            s = f'Last iteration reached a unbalance of {dn_c}'
            raise ValueError(f'Maximum number of iterations reached.\n{s}')
        # Found equilibrium
        # save the triangulation data
        if self.triangulated_data is None:
            self.triangulated_data = tri
        # Return the strain distribution
        return [eps_0, chi_c, 0]

    def _prefind_range_curvature_equilibrium(
        self,
        geom: CompoundGeometry,
        n: float,
        curv: float,
        x_a: float,
        dn_a: float,
    ):
        # With the following algorithm we find quickly a position of NA that
        # guarantees that the is at least one zero in the function dn vs. curv
        # in order to apply bisection algorithm
        ITMAX = 100
        sign = -1 if dn_a > 0 else 1
        found = False
        it = 0
        while not found and it < ITMAX:
            x_b = x_a + sign * (10 * it**2)
            (
                n_int,
                _,
                _,
                _,
            ) = self.integrator.integrate_strain_response_on_geometry(
                geom, [curv * x_b, curv, 0], tri=self.triangulated_data
            )
            dn_b = n_int - n
            if dn_a * dn_b < 0:
                found = True
            it += 1
        if it >= ITMAX and not found:
            s = f'Last iteration reached a unbalance of: \
                dn_a = {dn_a} dn_b = {dn_b})'
            raise ValueError(f'Maximum number of iterations reached.\n{s}')
        return (x_b, dn_b)

    def find_equilibrium_fixed_cruvature(
        self, geom: CompoundGeometry, n: float, curv: float, x: float
    ):
        """Find strain profile with equilibrium with fixed curvature.
        Given curvature and external axial force, find the strain profile
        that makes internal and external axial force in equilibrium.

        Arguments:
        geom: (CompounGeometry) the geometry
        n: (float) the external axial load
        curv: (float) the value of curvature
        x: (float) a first attempt for neutral axis position
        """
        # Useful for Moment Curvature Analysis
        # Number of maximum iteration for the bisection algorithm
        ITMAX = 100
        # Start from previous position of N.A.
        eps_0 = curv * x
        # find internal axial force by integration
        (
            n_int,
            _,
            _,
            tri,
        ) = self.integrator.integrate_strain_response_on_geometry(
            geom, [eps_0, curv, 0], tri=self.triangulated_data
        )
        if self.triangulated_data is None:
            self.triangulated_data = tri
        dn_a = n_int - n
        x_a = x
        x_b, dn_b = self._prefind_range_curvature_equilibrium(
            geom, n, curv, x_a, dn_a
        )
        # Found a range within there is the solution, apply bisection
        it = 0
        while (abs(dn_a - dn_b) > 1e-2) and (it < ITMAX):
            x_c = (x_a + x_b) / 2
            (
                n_int,
                _,
                _,
                _,
            ) = self.integrator.integrate_strain_response_on_geometry(
                geom, [curv * x_c, curv, 0], tri=self.triangulated_data
            )
            dn_c = n_int - n
            if dn_a * dn_c < 0:
                dn_b = dn_c
                x_b = x_c
            else:
                dn_a = dn_c
                x_a = x_c
            it += 1
        if it >= ITMAX:
            s = f'Last iteration reached a unbalance of: \
                dn_c = {dn_c}'
            raise ValueError(f'Maximum number of iterations reached.\n{s}')
        return (x_c, [curv * x_c, curv, 0])

    def calculate_bending_strength(
        self, theta=0, n=0
    ) -> s_res.UltimateBendingMomentResults:
        """Calculates the bending strength for given inclination of n.a.
        and axial load.

        Arguments:
        theta (float, default = 0): inclination of n.a. respect to y axis
        n (float, default = 0): axial load applied to the section
        (+: tension, -: compression)

        Return:
        ultimate_bending_moment_result (UltimateBendingMomentResult)
        """
        # Compute the bending strength with the bisection algorithm
        # Rotate the section of angle theta
        rotated_geom = self.section.geometry.rotate(theta)
        # Find the strain distribution corresponding to failure and equilibrium
        # with external axial force
        strain = self.find_equilibrium_fixed_pivot(rotated_geom, n)
        # Compute the internal forces with this strain distribution
        N, Mx, My, _ = self.integrator.integrate_strain_response_on_geometry(
            geo=rotated_geom, strain=strain, tri=self.triangulated_data
        )

        # Rotate back to section CRS
        T = np.array([[cos(theta), sin(theta)], [-sin(theta), cos(theta)]])
        M = T @ np.array([[Mx], [My]])

        # Create result object
        res = s_res.UltimateBendingMomentResults()
        res.theta = theta
        res.n = N
        res.chi_x = strain[1]
        res.chi_y = strain[2]
        res.eps_a = strain[0]
        res.m_x = M[0, 0]
        res.m_y = M[1, 0]

        return res

    def calculate_moment_curvature(
        self, theta=0, n=0
    ) -> s_res.MomentCurvatureResults:
        """Calculates the moment-curvature relation for given inclination of
        n.a. and axial load.

        Arguments:
        theta (float, default = 0): inclination of n.a. respect to y axis
        n (float, default = 0): axial load applied to the section
            (+: tension, -: compression)

        Return:
        moment_curvature_result (MomentCurvatureResults)
        """
        # Create an empty response object
        res = s_res.MomentCurvatureResults()
        res.n = n
        # Rotate the section of angle theta
        rotated_geom = self.section.geometry.rotate(theta)
        # Find ultimate curvature from the strain distribution corresponding
        # to failure and equilibrium with external axial force
        strain = self.find_equilibrium_fixed_pivot(rotated_geom, n)
        chi_ultimate = strain[1]
        # Find the yielding curvature
        strain = self.find_equilibrium_fixed_pivot(rotated_geom, n, True)
        chi_yield = strain[1]
        # Define the array of curvatures
        chi = np.concatenate(
            (
                np.linspace(1e-13, chi_yield, 10, endpoint=False),
                np.linspace(chi_yield, chi_ultimate, 100),
            )
        )
        # prepare results
        eps_a = np.zeros_like(chi)
        my = np.zeros_like(chi)
        mz = np.zeros_like(chi)
        # Previous position of neutral axes
        x = 0
        # For each value of curvature
        for i, curv in enumerate(chi):
            # find the new position of neutral axis for mantaining equilibrium
            # store the information in the results object for the current
            # value of curvature
            x, strain = self.find_equilibrium_fixed_cruvature(
                rotated_geom, n, curv, x
            )
            (
                _,
                Mx,
                My,
                _,
            ) = self.integrator.integrate_strain_response_on_geometry(
                geo=rotated_geom, strain=strain, tri=self.triangulated_data
            )
            # Rotate back to section CRS
            T = np.array([[cos(theta), sin(theta)], [-sin(theta), cos(theta)]])
            M = T @ np.array([[Mx], [My]])
            eps_a[i] = strain[0]
            my[i] = M[0, 0]
            mz[i] = M[1, 0]
        res.chi = chi
        res.eps_axial = eps_a
        res.my = my
        res.mz = mz

        return res


# Use examples:

# sec = Section(...xxx...)
# res = sec.section_analyzer.ultmateBendingStrength(n = -1000)
# stress = sec.section_analyzer.compute_stress_distribution(res)
