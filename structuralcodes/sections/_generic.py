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
import structuralcodes.sections._section_results as s_res
from ._section_integrators import integrator_factory


class GenericSection(Section):
    """This is the implementation of the generic class section."""

    def __init__(
        self,
        geometry: SurfaceGeometry | CompoundGeometry,
        name: t.Optional[str] = None,
    ) -> None:
        if name is None:
            name = 'GenericSection'
        super().__init__(name)
        self.geometry = geometry
        self.section_analyzer = GenericSectionCalculator(self)
        self.gross_properties = (
            self.section_analyzer._calculate_gross_section_properties()
        )


class GenericSectionCalculator(SectionCalculator):
    """Calculator class implementing analysis algortims for
    Code checks.
    """

    def __init__(self, sec: GenericSection) -> None:
        """Initialize the GenericSectionCalculator
        Input:
        section (SectionCalculator): the section object.
        """
        super().__init__(section=sec)

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

    def calculate_bending_strength(
        self, theta=0, n=0, **kwargs
    ) -> s_res.UltimateBendingMomentResult:
        """Calculates the bending strength for given inclination of n.a.
        and axial load.

        Input:
        theta (float, default = 0): inclination of n.a. respect to y axis
        n (float, default = 0): axial load applied to the section
        (+: tension, -: compression)

        Return:
        ultimate_bending_moment_result (UltimateBendingMomentResult)
        """
        # Select the integrator if specified
        integr = integrator_factory(kwargs.get('integrator', 'Marin'))()

        ITMAX = 100
        # Compute the bending strength with the bisection algorithm
        # 1. Rotate the section of angle theta
        rotated_geom = self.section.geometry.rotate(theta)
        # 2. Start with a balanced failure: this is found from all ultimate
        # strains for all materials, checking the minimum curvature value
        chi_min = 1e10
        for g in rotated_geom.geometries + rotated_geom.point_geometries:
            for other_g in (
                rotated_geom.geometries + rotated_geom.point_geometries
            ):
                if g != other_g:
                    eps_p = g.material.get_ultimate_strain()[0]
                    if isinstance(g, SurfaceGeometry):
                        y_p = g.polygon.bounds[1]
                    elif isinstance(g, PointGeometry):
                        y_p = g._point.coords[0][1]
                    eps_n = other_g.material.get_ultimate_strain()[1]
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
                        eps_p_min = eps_p
                        eps_n_min = eps_n

        strain = [eps_0, chi_min, 0]
        # print('strain:',strain)
        y_p, y_n = y_p_min, y_n_min
        eps_p, eps_n = eps_p_min, eps_n_min
        # Integrate this strain profile
        n_int, _, _, tri = integr.integrate_strain_response_on_geometry(
            rotated_geom, strain, **kwargs
        )
        # print('n_int=',n_int)
        # 3. Check if we have equilibrium
        chi_a = chi_min
        dn_a = n_int - n
        it = 1
        chi_c = 0
        if n_int < n:
            # Too much compression, raise NA
            chi_b = 1e-13
            eps_0 = eps_p + chi_b * y_p
            # print('Too much compression')
            n_int, _, _, _ = integr.integrate_strain_response_on_geometry(
                rotated_geom, [eps_0, chi_b, 0], tri=tri
            )
            dn_b = n_int - n
            while (abs(dn_a - dn_b) > 1e-2) and (it < ITMAX):
                chi_c = (chi_a + chi_b) / 2.0
                eps_0 = eps_p + chi_c * y_p
                n_int, _, _, _ = integr.integrate_strain_response_on_geometry(
                    rotated_geom, [eps_0, chi_c, 0], tri=tri
                )
                dn_c = n_int - n
                if dn_c * dn_a < 0:
                    chi_b = chi_c
                    dn_b = dn_c
                else:
                    chi_a = chi_c
                    dn_a = dn_c
                it += 1
        elif n_int > n:
            # Too much tension, lower NA
            chi_b = 1e-13
            eps_0 = eps_n + chi_b * y_n
            # print('Too much tension')
            n_int, _, _, _ = integr.integrate_strain_response_on_geometry(
                rotated_geom, [eps_0, chi_b, 0], tri=tri
            )
            dn_b = n_int - n
            while (abs(dn_a - dn_b) > 1e-2) and (it < ITMAX):
                chi_c = (chi_a + chi_b) / 2.0
                eps_0 = eps_n + chi_c * y_n
                n_int, _, _, _ = integr.integrate_strain_response_on_geometry(
                    rotated_geom, [eps_0, chi_c, 0], tri=tri
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
            raise ValueError(f'Maximum number of iterations reached\n{s}')

        # Found equilibrium
        # print(f'Found equilibrium (after {it} iterations)')
        strain = [eps_0, chi_c, 0]
        N, Mx, My, _ = integr.integrate_strain_response_on_geometry(
            geo=rotated_geom, strain=strain, tri=tri
        )
        # Check if n and N are sufficiently near
        # print('check: ',n,N)

        # Rotate back to section CRS
        T = np.array([[cos(theta), sin(theta)], [-sin(theta), cos(theta)]])
        M = T @ np.array([[Mx], [My]])

        # Create result object
        res = s_res.UltimateBendingMomentResult()
        res.theta = theta
        res.n = N
        res.chi_x = chi_c
        res.chi_y = 0
        res.eps_a = eps_0
        res.m_x = M[0, 0]
        res.m_y = M[1, 0]

        return res

    def calculate_moment_curvature(
        self, theta=0, n=0
    ) -> s_res.MomentCurvatureResults:
        """Calculates the moment-curvature relation for given inclination of
        n.a. and axial load.

        Input:
        theta (float, default = 0): inclination of n.a. respect to y axis
        n (float, default = 0): axial load applied to the section
            (+: tension, -: compression)
        chi_incr (float, default = 1e-8): the curvature increment for the
            analysis

        Return:
        moment_curvature_result (MomentCurvatureResults)
        """
        # For now it returns an empty response. The proper algorithms for
        # generic section will be here
        res = s_res.MomentCurvatureResults()
        raise NotImplementedError
        return res


# Use examples:

# sec = Section(...xxx...)
# res = sec.section_analyzer.ultmateBendingStrength(n = -1000)
# stress = sec.section_analyzer.compute_stress_distribution(res)
