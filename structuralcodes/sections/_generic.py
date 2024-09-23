"""Generic class section implemenetation."""

from __future__ import annotations

import math
import typing as t
import warnings
from math import cos, sin

import numpy as np
from numpy.typing import ArrayLike
from shapely import MultiPolygon
from shapely.ops import unary_union

import structuralcodes.core._section_results as s_res
from structuralcodes.core.base import Section, SectionCalculator
from structuralcodes.geometry import (
    CompoundGeometry,
    PointGeometry,
    SurfaceGeometry,
)
from structuralcodes.materials.constitutive_laws import Elastic

from .section_integrators import integrator_factory


class GenericSection(Section):
    """This is the implementation of the generic class section.

    The section is a 2D geometry where Y axis is horizontal while
    Z axis is vertical.
    The moments and curvatures around Y and Z axes are assumed
    positive according to RHR.

    Attributes:
        geometry: (SurfaceGeometry or CompoundGeometry)
            the geometry of the section
        name: (str) the name of the section
        section_calculator: (GenericSectionCalculator)
            the object responsible for performing different
            calculations on the section (e.g. bending strength,
            moment curvature, etc.)
        gross_properties: (GrossProperties)
            the results of gross properties computation
    """

    def __init__(
        self,
        geometry: t.Union[SurfaceGeometry, CompoundGeometry],
        name: t.Optional[str] = None,
        integrator: t.Literal['marin', 'fiber'] = 'marin',
        **kwargs,
    ) -> None:
        if name is None:
            name = 'GenericSection'
        super().__init__(name)
        # Since only CompoundGeometry has the attribute geometries,
        # if a SurfaceGeometry is input, we create a CompoundGeometry
        # with only that geometry contained. After that all algorithms
        # work as usal.
        if isinstance(geometry, SurfaceGeometry):
            geometry = CompoundGeometry([geometry])
        self.geometry = geometry
        self.section_calculator = GenericSectionCalculator(
            sec=self, integrator=integrator, **kwargs
        )
        self._gross_properties = None

    @property
    def gross_properties(self) -> s_res.GrossProperties:
        """Return the gross properties of the section."""
        if self._gross_properties is None:
            self._gross_properties = (
                self.section_calculator._calculate_gross_section_properties()
            )
        return self._gross_properties


class GenericSectionCalculator(SectionCalculator):
    """Calculator class implementing analysis algortims for
    Code checks.
    """

    def __init__(
        self, sec: GenericSection, integrator: str = 'Marin', **kwargs
    ) -> None:
        """Initialize the GenericSectionCalculator.

        Args:
            section (SectionCalculator): the section object.
            integrator (str): the integrator to be used for computations
                            (default = 'marin')

        Note:
            when using 'fiber' integrator the kwarg 'mesh_size' can be used
            to specify a dimensionless number (between 0 and 1) specifying
            the size of the resulting mesh.
        """
        super().__init__(section=sec)
        # Select the integrator if specified
        self.integrator = integrator_factory(integrator)()
        # Mesh size used for Fibre integrator
        self.mesh_size = kwargs.get('mesh_size', 0.01)
        # triangulated_data used for Fibre integrator
        self.triangulated_data = None
        # Maximum and minimum axial load
        self._n_max = None
        self._n_min = None

    def _calculate_gross_section_properties(self) -> s_res.GrossProperties:
        """Calculates the gross section properties of the GenericSection.

        This function is private and called when the section is created
        It stores the result into the result object.

        Returns:
        gross_section_properties (GrossSection)
        """
        # It will use the algorithms for generic sections
        gp = s_res.GrossProperties()

        # Computation of perimeter using shapely
        polygon = unary_union(
            [geo.polygon for geo in self.section.geometry.geometries]
        )
        if isinstance(polygon, MultiPolygon):
            gp.perimeter = 0.0
            warnings.warn(
                'Perimiter computation for a multi polygon is not defined.'
            )

        gp.perimeter = polygon.exterior.length

        # Computation of area: this is taken directly from shapely
        gp.area = self.section.geometry.area
        # Computation of surface area, reinforcement area, EA (axial rigidity)
        # and mass: Morten -> problem with units! how do we deal with it?
        for geo in self.section.geometry.geometries:
            gp.ea += geo.area * geo.material.get_tangent(eps=0)[0]
            if geo.density is not None:
                # this assumes area in mm2 and density in kg/m3
                gp.mass += geo.area * geo.density * 1e-9

        for geo in self.section.geometry.point_geometries:
            gp.ea += geo.area * geo.material.get_tangent(eps=0)[0]
            gp.area_reinforcement += geo.area
            if geo.density is not None:
                # this assumes area in mm2 and density in kg/m3
                gp.mass += geo.area * geo.density * 1e-9

        # Computation of area moments
        #
        # Implementation idea:
        # Using integrator: we need to compute the following integrals:
        # E Sy = integr(E*z*dA)
        # E Sz = integr(E*y*dA)
        # E Iyy = integr(E*z*z*dA)
        # E Izz = integr(E*y*y*dA)
        # E Iyz = integr(E*y*z*dA)
        #
        # The first can be imagined as computing axial force
        # by integration of a E*z stress;
        # since eps = eps_a + ky * z - kz * y
        # E*z is the stress corresponding to an elastic material
        # with stiffness E and strain equal to z (i.e. eps_a and kz = 0,
        # ky = 1 )
        #
        # With the same idea we can integrate the other quantities.

        def compute_area_moments(geometry, material=None):
            # create a new dummy geometry from the original one
            # with dummy material
            geometry = geometry.from_geometry(
                geo=geometry, new_material=material
            )
            # Integrate a dummy strain profile for getting first and
            # second moment respect y axis and product moment
            (
                sy,
                iyy,
                iyz,
                tri,
            ) = self.integrator.integrate_strain_response_on_geometry(
                geometry,
                [0, 1, 0],
                mesh_size=self.mesh_size,
            )
            # Change sign due to moment sign convention
            iyz *= -1
            # Integrate a dummy strain profile for getting first
            # and second moment respect z axis and product moment
            (
                sz,
                izy,
                izz,
                _,
            ) = self.integrator.integrate_strain_response_on_geometry(
                geometry,
                [0, 0, -1],
                tri=tri,
                mesh_size=self.mesh_size,
            )
            # Change sign due to moment sign convention
            izz *= -1
            if abs(abs(izy) - abs(iyz)) > 10:
                error_str = 'Something went wrong with computation of '
                error_str += f'moments of area: iyz = {iyz}, izy = {izy}.\n'
                error_str += 'They should be equal but are not!'
                raise RuntimeError(error_str)

            return sy, sz, iyy, izz, iyz

        # Create a dummy material for integration of area moments
        # This is used for J, S etc, not for E_J E_S etc
        dummy_mat = Elastic(E=1)
        # Computation of moments of area (material-independet)
        # Note: this could be un-meaningfull when many materials
        # are combined
        gp.sy, gp.sz, gp.iyy, gp.izz, gp.iyz = compute_area_moments(
            geometry=self.section.geometry, material=dummy_mat
        )

        # Computation of moments of area times E
        gp.e_sy, gp.e_sz, gp.e_iyy, gp.e_izz, gp.e_iyz = compute_area_moments(
            geometry=self.section.geometry
        )

        # Compute Centroid coordinates
        gp.cy = gp.e_sz / gp.ea
        gp.cz = gp.e_sy / gp.ea

        # Compute of moments of area relative to yz centroidal axes
        translated_geometry = self.section.geometry.translate(
            dx=-gp.cy, dy=-gp.cz
        )
        _, _, gp.iyy_c, gp.izz_c, gp.iyz_c = compute_area_moments(
            geometry=translated_geometry, material=dummy_mat
        )

        # Computation of moments of area times E
        _, _, gp.e_iyy_c, gp.e_izz_c, gp.e_iyz_c = compute_area_moments(
            geometry=translated_geometry
        )

        # Compute principal axes of inertia and principal inertia
        def find_principal_axes_moments(iyy, izz, iyz):
            eigres = np.linalg.eig(np.array([[iyy, iyz], [iyz, izz]]))
            max_idx = np.argmax(eigres[0])
            min_idx = 0 if max_idx == 1 else 1
            i11 = eigres[0][max_idx]
            i22 = eigres[0][min_idx]
            theta = np.arccos(np.dot(np.array([1, 0]), eigres[1][:, max_idx]))
            return i11, i22, theta

        gp.i11, gp.i22, gp.theta = find_principal_axes_moments(
            gp.iyy_c, gp.izz_c, gp.iyz_c
        )
        gp.e_i11, gp.e_i22, gp.e_theta = find_principal_axes_moments(
            gp.e_iyy_c, gp.e_izz_c, gp.e_iyz_c
        )

        return gp

    def get_balanced_failure_strain(
        self, geom: CompoundGeometry, yielding: bool = False
    ) -> t.Tuple[float, float, float]:
        """Returns the strain profile corresponding to balanced failure.

        This is found from all ultimate strains for all materials, checking
        the minimum value of curvature.
        It returns a tuple with:
        1. Value of y coordinate for negative failure
        2. Value of y coordinate for positive failure
        3. Strain profile as a list with three values: axial strain, curvature
           y*, curvature z* (assumed zero since in the rotated frame y*z* it is
           a case of uniaxial bending).

        Args:
            geom: the compund geometry
            yielding: (bool, default = False) consider yielding instead of
                ultimate strain

        Returns:
            Returns a tuple (float, float, float) that define the strain
            distribution
        """
        chi_min = 1e10
        for g in geom.geometries + geom.point_geometries:
            for other_g in geom.geometries + geom.point_geometries:
                # if g != other_g:
                eps_p = g.material.get_ultimate_strain(yielding=yielding)[0]
                if isinstance(g, SurfaceGeometry):
                    y_p = g.polygon.bounds[1]
                elif isinstance(g, PointGeometry):
                    y_p = g._point.coords[0][1]
                eps_n = other_g.material.get_ultimate_strain(
                    yielding=yielding
                )[1]
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
        # In standard CRS negative curvature stretches bottom fiber
        strain = [eps_0, -chi_min, 0]
        return (y_n, y_p, strain)

    def find_equilibrium_fixed_pivot(
        self, geom: CompoundGeometry, n: float, yielding: bool = False
    ) -> t.Tuple[float, float, float]:
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
        eps_p = strain[0] + strain[1] * y_p
        eps_n = strain[0] + strain[1] * y_n
        # Integrate this strain profile corresponding to balanced failure
        (
            n_int,
            _,
            _,
            tri,
        ) = self.integrator.integrate_strain_response_on_geometry(
            geom, strain, tri=self.triangulated_data, mesh_size=self.mesh_size
        )
        # Check if there is equilibrium with this strain distribution
        chi_a = strain[1]
        dn_a = n_int - n
        # It may occur that dn_a is already almost zero (in equilibrium)
        if abs(dn_a) <= 1e-2:
            # return the equilibrium position
            return [strain[0], chi_a, 0]
        chi_b = -1e-13
        if n_int < n:
            # Too much compression, raise NA
            pivot = y_p
            strain_pivot = eps_p
        else:
            # Too much tension, lower NA
            pivot = y_n
            strain_pivot = eps_n
        eps_0 = strain_pivot - chi_b * pivot
        n_int, _, _, _ = self.integrator.integrate_strain_response_on_geometry(
            geom, [eps_0, chi_b, 0], tri=self.triangulated_data
        )
        dn_b = n_int - n
        it = 0
        while (abs(dn_a - dn_b) > 1e-2) and (it < ITMAX):
            chi_c = (chi_a + chi_b) / 2.0
            eps_0 = strain_pivot - chi_c * pivot
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
        eps_0_a: float,
        dn_a: float,
    ):
        """Perfind range where the curvature equilibrium is located.

        This algorithms quickly finds a position of NA that guaranteed the
        existence of at least one zero in the function dn vs. curv in order
        to apply the bisection algorithm.
        """
        ITMAX = 50
        sign = -1 if dn_a > 0 else 1
        found = False
        it = 0
        delta = 1e-3
        while not found and it < ITMAX:
            eps_0_b = eps_0_a + sign * delta * (it + 1)
            (
                n_int,
                _,
                _,
                _,
            ) = self.integrator.integrate_strain_response_on_geometry(
                geom, [eps_0_b, curv, 0], tri=self.triangulated_data
            )
            dn_b = n_int - n
            if dn_a * dn_b < 0:
                found = True
            elif abs(dn_b) > abs(dn_a):
                # we are driving aay from the solution, probably due
                # to failure of a material
                delta /= 2
                it -= 1
            it += 1
        if it >= ITMAX and not found:
            s = f'Last iteration reached a unbalance of: \
                dn_a = {dn_a} dn_b = {dn_b})'
            raise ValueError(f'Maximum number of iterations reached.\n{s}')
        return (eps_0_b, dn_b)

    def find_equilibrium_fixed_curvature(
        self, geom: CompoundGeometry, n: float, curv: float, eps_0: float
    ):
        """Find strain profile with equilibrium with fixed curvature.
        Given curvature and external axial force, find the strain profile
        that makes internal and external axial force in equilibrium.

        Arguments:
        geom: (CompounGeometry) the geometry
        n: (float) the external axial load
        curv: (float) the value of curvature
        eps_0: (float) a first attempt for neutral axis position
        """
        # Useful for Moment Curvature Analysis
        # Number of maximum iteration for the bisection algorithm
        ITMAX = 100
        # Start from previous position of N.A.
        eps_0_a = eps_0
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
        # It may occur that dn_a is already almost zero (in eqiulibrium)
        if abs(dn_a) <= 1e-2:
            # return the equilibrium position
            return [eps_0_a, curv, 0]
        eps_0_b, dn_b = self._prefind_range_curvature_equilibrium(
            geom, n, curv, eps_0_a, dn_a
        )
        # Found a range within there is the solution, apply bisection
        it = 0
        while (abs(dn_a - dn_b) > 1e-2) and (it < ITMAX):
            eps_0_c = (eps_0_a + eps_0_b) / 2
            (
                n_int,
                _,
                _,
                _,
            ) = self.integrator.integrate_strain_response_on_geometry(
                geom, [eps_0_c, curv, 0], tri=self.triangulated_data
            )
            dn_c = n_int - n
            if dn_a * dn_c < 0:
                dn_b = dn_c
                eps_0_b = eps_0_c
            else:
                dn_a = dn_c
                eps_0_a = eps_0_c
            it += 1
        if it >= ITMAX:
            s = f'Last iteration reached a unbalance of: \
                dn_c = {dn_c}'
            raise ValueError(f'Maximum number of iterations reached.\n{s}')
        return [eps_0_c, curv, 0]

    def calculate_limit_axial_load(self):
        """Compute maximum and minimum axial load.

        Returns:
            (float, float) Minimum and Maximum axial load.
        """
        # Find balanced failure to get strain limits
        y_n, y_p, strain = self.get_balanced_failure_strain(
            geom=self.section.geometry, yielding=False
        )
        eps_p = strain[0] + strain[1] * y_p
        eps_n = strain[0] + strain[1] * y_n

        n_min, _, _, tri = (
            self.integrator.integrate_strain_response_on_geometry(
                self.section.geometry,
                [eps_n, 0, 0],
                tri=self.triangulated_data,
                mesh_size=self.mesh_size,
            )
        )
        n_max, _, _, _ = self.integrator.integrate_strain_response_on_geometry(
            self.section.geometry, [eps_p, 0, 0], tri=tri
        )

        if self.triangulated_data is None:
            self.triangulated_data = tri
        return n_min, n_max

    @property
    def n_min(self) -> float:
        """Return minimum axial load."""
        if self._n_min is None:
            self._n_min, self._n_max = self.calculate_limit_axial_load()
        return self._n_min

    @property
    def n_max(self) -> float:
        """Return maximum axial load."""
        if self._n_max is None:
            self._n_min, self._n_max = self.calculate_limit_axial_load()
        return self._n_max

    def check_axial_load(self, n: float):
        """Check if axial load n is within section limits."""
        if n < self.n_min or n > self.n_max:
            error_str = f'Axial load {n} cannot be taken by section.\n'
            error_str += f'n_min = {self.n_min} / n_max = {self.n_max}'
            raise ValueError(error_str)

    def _rotate_triangulated_data(self, theta: float):
        """Rotate triangulated data of angle theta."""
        rotated_triangulated_data = []
        for tr in self.triangulated_data:
            T = np.array([[cos(theta), -sin(theta)], [sin(theta), cos(theta)]])
            coords = np.vstack((tr[0], tr[1]))
            coords_r = T @ coords
            rotated_triangulated_data.append(
                (coords_r[0, :], coords_r[1, :], tr[2], tr[3])
            )
        self.triangulated_data = rotated_triangulated_data

    def integrate_strain_profile(
        self, strain: ArrayLike
    ) -> t.Tuple[float, float, float]:
        """Integrate a strain profile returning internal forces.

        Args:
            strain: (Array like) represents the deformation plane
                    the strain should have three entries representing
                    respectively:
                    axial strain (At 0,0 coordinates)
                    curv_y
                    curv_z

        Returns:
            (float, float, float) N, My and Mz
        """
        N, My, Mz, _ = self.integrator.integrate_strain_response_on_geometry(
            geo=self.section.geometry,
            strain=strain,
            tri=self.triangulated_data,
        )
        return N, My, Mz

    def calculate_bending_strength(
        self, theta=0, n=0
    ) -> s_res.UltimateBendingMomentResults:
        """Calculates the bending strength for given inclination of n.a.
        and axial load.

        Arguments:
        theta (float, default = 0): inclination of n.a. respect to section
            y axis
        n (float, default = 0): axial load applied to the section
            (+: tension, -: compression)

        Return:
        ultimate_bending_moment_result (UltimateBendingMomentResult)
        """
        # Compute the bending strength with the bisection algorithm
        # Rotate the section of angle theta
        rotated_geom = self.section.geometry.rotate(-theta)
        if self.triangulated_data is not None:
            # Rotate also triangulated data!
            self._rotate_triangulated_data(-theta)

        # Check if the section can carry the axial load
        self.check_axial_load(n=n)
        # Find the strain distribution corresponding to failure and equilibrium
        # with external axial force
        strain = self.find_equilibrium_fixed_pivot(rotated_geom, n)
        # Compute the internal forces with this strain distribution
        N, My, Mz, _ = self.integrator.integrate_strain_response_on_geometry(
            geo=rotated_geom, strain=strain, tri=self.triangulated_data
        )

        # Rotate back to section CRS TODO Check
        T = np.array([[cos(theta), -sin(theta)], [sin(theta), cos(theta)]])
        M = T @ np.array([[My], [Mz]])
        if self.triangulated_data is not None:
            # Rotate back also triangulated data!
            self._rotate_triangulated_data(theta)

        # Create result object
        res = s_res.UltimateBendingMomentResults()
        res.theta = theta
        res.n = N
        res.chi_y = strain[1]
        res.chi_z = strain[2]
        res.eps_a = strain[0]
        res.m_y = M[0, 0]
        res.m_z = M[1, 0]

        return res

    def calculate_moment_curvature(
        self,
        theta: float = 0.0,
        n: float = 0.0,
        chi_first: float = 1e-8,
        num_pre_yield: int = 10,
        num_post_yield: int = 10,
        chi: t.Optional[ArrayLike] = None,
    ) -> s_res.MomentCurvatureResults:
        """Calculates the moment-curvature relation for given inclination of
        n.a. and axial load.

        Arguments:
        theta (float, default = 0): inclination of n.a. respect to y axis
        n (float, default = 0): axial load applied to the section (+: tension,
            -: compression)
        chi_first (float, default = 1e-8): the first value of the curvature
        num_pre_yield (int, default = 10): Number of points before yielding.
            Note that the yield curvature will be at the num_pre_yield-th point
            in the result array.
        num_post_yield (int, default = 10): Number of points after yielding
        chi (Optional[ArrayLike], default = None): An ArrayLike with curvatures
            to calculate the moment response for. If chi is None, the array is
            constructed from chi_first, num_pre_yield and num_post_yield. If
            chi is not None, chi_first, num_pre_yield and num_post_yield are
            disregarded, and the provided chi is used directly in the
            calculations.

        Return:
        moment_curvature_result (MomentCurvatureResults)
        """
        # Create an empty response object
        res = s_res.MomentCurvatureResults()
        res.n = n
        # Rotate the section of angle theta
        rotated_geom = self.section.geometry.rotate(-theta)
        if self.triangulated_data is not None:
            # Rotate also triangulated data!
            self._rotate_triangulated_data(-theta)

        # Check if the section can carry the axial load
        self.check_axial_load(n=n)

        if chi is None:
            # Find ultimate curvature from the strain distribution
            # corresponding to failure and equilibrium with external axial
            # force
            strain = self.find_equilibrium_fixed_pivot(rotated_geom, n)
            chi_ultimate = strain[1]
            # Find the yielding curvature
            strain = self.find_equilibrium_fixed_pivot(
                rotated_geom, n, yielding=True
            )
            chi_yield = strain[1]
            if chi_ultimate * chi_yield < 0:
                # They cannot have opposite signs!
                raise ValueError(
                    'curvature at yield and ultimate cannot have opposite '
                    'signs!'
                )

            # Make sure the sign of the first curvature matches the sign of the
            # yield curvature
            chi_first *= -1.0 if chi_first * chi_yield < 0 else 1.0

            # The first curvature should be less than the yield curvature
            if abs(chi_first) >= abs(chi_yield):
                chi_first = chi_yield / num_pre_yield

            # Define the array of curvatures
            if abs(chi_ultimate) <= abs(chi_yield) + 1e-8:
                # We don't want a plastic branch in the analysis
                # this is done to speed up analysis
                chi = np.linspace(chi_first, chi_yield, num_pre_yield)
            else:
                chi = np.concatenate(
                    (
                        np.linspace(
                            chi_first,
                            chi_yield,
                            num_pre_yield - 1,
                            endpoint=False,
                        ),
                        np.linspace(
                            chi_yield, chi_ultimate, num_post_yield + 1
                        ),
                    )
                )

        # prepare results
        eps_a = np.zeros_like(chi)
        my = np.zeros_like(chi)
        mz = np.zeros_like(chi)
        chi_y = np.zeros_like(chi)
        chi_z = np.zeros_like(chi)
        # Previous position of strain at (0,0)
        strain = [0, 0, 0]
        # For each value of curvature
        for i, curv in enumerate(chi):
            # find the new position of neutral axis for mantaining equilibrium
            # store the information in the results object for the current
            # value of curvature
            strain = self.find_equilibrium_fixed_curvature(
                rotated_geom, n, curv, strain[0]
            )
            (
                _,
                My,
                Mz,
                _,
            ) = self.integrator.integrate_strain_response_on_geometry(
                geo=rotated_geom, strain=strain, tri=self.triangulated_data
            )
            # Rotate back to section CRS
            T = np.array([[cos(theta), -sin(theta)], [sin(theta), cos(theta)]])
            M = T @ np.array([[My], [Mz]])
            eps_a[i] = strain[0]
            my[i] = M[0, 0]
            mz[i] = M[1, 0]
            chi_mat = T @ np.array([[curv], [0]])
            chi_y[i] = chi_mat[0, 0]
            chi_z[i] = chi_mat[1, 0]

        if self.triangulated_data is not None:
            # Rotate back also triangulated data!
            self._rotate_triangulated_data(theta)
        res.chi_y = chi_y
        res.chi_z = chi_z
        res.eps_axial = eps_a
        res.m_y = my
        res.m_z = mz

        return res

    def calculate_nm_interaction_domain(
        self, theta: float = 0, num_axial: int = 100
    ) -> s_res.NMInteractionDomain:
        """Calculate the NM interaction domain.

        Arguments:
        theta (float, default = 0): inclination of n.a. respect to y axis
        n_axial (int, default = 100): number of discretization for axial load

        Return:
        (NmInteractionDomain)
        """
        # Prepare the results
        res = s_res.NMInteractionDomain()
        res.theta = theta
        # with this version really n_axial is not used unless we perform
        # some resampling after computation?
        res.num_axial = num_axial

        # Get ultimate strain profiles for theta angle
        strains = self._compute_ultimate_strain_profiles(theta=theta)

        # integrate all strain profiles
        forces = np.zeros_like(strains)
        for i, strain in enumerate(strains):
            N, My, Mz, tri = (
                self.integrator.integrate_strain_response_on_geometry(
                    geo=self.section.geometry,
                    strain=strain,
                    tri=self.triangulated_data,
                )
            )
            if self.triangulated_data is None:
                self.triangulated_data = tri
            forces[i, 0] = N
            forces[i, 1] = My
            forces[i, 2] = Mz

        # Save to results
        res.strains = strains
        res.forces = forces
        res.m_z = forces[:, 2]
        res.m_y = forces[:, 1]
        res.n = forces[:, 0]

        return res

    def _compute_ultimate_strain_profiles(self, theta: float = 0):
        """Return an array of ultimate strain profiles.

        Args:
        theta: (float) the angle of neutral axis
        """
        rotated_geom = self.section.geometry.rotate(-theta)

        # Find yield failure
        y_n, y_p, strain = self.get_balanced_failure_strain(
            geom=rotated_geom, yielding=True
        )
        eps_p_y = strain[0] + strain[1] * y_p
        eps_n_y = strain[0] + strain[1] * y_n
        # Find balanced failure: this defines the transition
        # between fields 2 and 3
        y_n, y_p, strain = self.get_balanced_failure_strain(
            geom=rotated_geom, yielding=False
        )
        eps_p_b = strain[0] + strain[1] * y_p
        eps_n_b = strain[0] + strain[1] * y_n

        # get h of the rotated geometry
        _, _, min_y, _ = rotated_geom.calculate_extents()
        h = y_n - min_y

        # For generation of fields 1 and 2 pivot on positive strain
        # Field 1: pivot on positive strain
        n = 3
        eps_n = np.linspace(eps_p_b, 0, n, endpoint=False)
        eps_p = np.zeros_like(eps_n) + eps_p_b
        # Field 2: pivot on positive strain
        n = 10  # 10
        eps_n = np.append(eps_n, np.linspace(0, eps_n_b, n, endpoint=False))
        eps_p = np.append(eps_p, np.zeros(n) + eps_p_b)
        # For fields 3-4-5 pivot on negative strain
        # Field 3: pivot on negative strain
        n = 40  # 40
        eps_n = np.append(eps_n, np.zeros(n) + eps_n_b)
        eps_p = np.append(
            eps_p, np.linspace(eps_p_b, eps_p_y, n, endpoint=False)
        )
        # Field 4: pivot on negative strain
        n = 10  # 6
        eps_n = np.append(eps_n, np.zeros(n) + eps_n_b)
        eps_p = np.append(eps_p, np.linspace(eps_p_y, 0, n, endpoint=False))
        # Field 5: pivot on negative strain
        n = 40  # 8
        eps_p_lim = eps_n_b * (h - (y_n - y_p)) / h
        eps_n = np.append(eps_n, np.zeros(n) + eps_n_b)
        eps_p = np.append(eps_p, np.linspace(0, eps_p_lim, n, endpoint=False))
        # Field 6: pivot on eps_n_y point or eps_n_b
        n = 40  # 8
        # If reinforced concrete section pivot on eps_n_y (default -0.002)
        # otherwise pivot on eps_n_b (in top chord)
        if self.section.geometry.reinforced_concrete:
            z_pivot = y_n - (1 - eps_n_y / eps_n_b) * h
            eps_p_6 = np.linspace(eps_p_lim, eps_n_y, n, endpoint=True)
            eps_n_6 = (
                -(eps_n_y - eps_p_6) * (z_pivot - y_n) / (z_pivot - y_p)
                + eps_n_y
            )
        else:
            eps_n_6 = np.zeros(n) + eps_n_b
            eps_p_6 = np.linspace(eps_p_lim, eps_n_b, n, endpoint=False)
        eps_n = np.append(eps_n, eps_n_6)
        eps_p = np.append(eps_p, eps_p_6)

        # rotate them
        kappa_y = (eps_n - eps_p) / (y_n - y_p)
        eps_a = eps_n - kappa_y * y_n
        kappa_z = np.zeros_like(kappa_y)

        # rotate back components to work in section CRS
        T = np.array([[cos(theta), -sin(theta)], [sin(theta), cos(theta)]])
        components = np.vstack((kappa_y, kappa_z))
        rotated_components = T @ components
        return np.column_stack((eps_a, rotated_components.T))

    def calculate_nmm_interaction_domain(
        self, num_theta: int = 32, num_axial: int = 20
    ) -> s_res.NMMInteractionDomain:
        """Calculates the NMM interaction domain."""
        res = s_res.NMMInteractionDomain()
        res.num_theta = num_theta
        # TODO: if we want to resample and structure better data we can
        # use n_axial as the number of discretizations in axial direction
        res.num_axial = num_axial

        # cycle for all n_thetas
        thetas = np.linspace(0, np.pi * 2, num_theta)
        # Initialize an empty array with the correct shape
        strains = np.empty((0, 3))
        for theta in thetas:
            # Get ultimate strain profiles for theta angle
            strain = self._compute_ultimate_strain_profiles(theta=theta)
            strains = np.vstack((strains, strain))

        # integrate all strain profiles
        forces = np.zeros_like(strains)
        for i, strain in enumerate(strains):
            N, My, Mz, tri = (
                self.integrator.integrate_strain_response_on_geometry(
                    geo=self.section.geometry,
                    strain=strain,
                    tri=self.triangulated_data,
                )
            )
            if self.triangulated_data is None:
                self.triangulated_data = tri
            forces[i, 0] = N
            forces[i, 1] = My
            forces[i, 2] = Mz

        # Save to results
        res.strains = strains
        res.forces = forces

        return res

    def calculate_mm_interaction_domain(
        self, n: float = 0, num_theta: int = 32
    ) -> s_res.MMInteractionDomain:
        """Calculate the My-Mz interaction domain.

        Arguments:
        n (float, default = 0): axial force
        n_theta (int, default = 32): number of discretization for theta

        Return:
        (MMInteractionDomain)
        """
        # Prepare the results
        res = s_res.MMInteractionDomain()
        res.num_theta = num_theta
        res.n = n
        # Create array of thetas
        res.theta = np.linspace(0, np.pi * 2, num_theta)
        # Initialize the result's arrays
        res.m_y = np.zeros_like(res.theta)
        res.m_z = np.zeros_like(res.theta)
        # Compute strength for given angle of NA
        for i, th in enumerate(res.theta):
            res_bend_strength = self.calculate_bending_strength(theta=th, n=n)
            res.m_y[i] = res_bend_strength.m_y
            res.m_z[i] = res_bend_strength.m_z

        return res

    def calculate_strain_profile(self, n_ed, my_ed, mz_ed, num_chis=5):
        """Get the strain plane for a given axial force and biaxial bending.

        Args:
            n_ed (float): Axial load
            my_ed (float): Bending moment around y-axis
            mz_ed (float): Bending moment around z-axis
            num_chis (int) [optional]: number of points in M-curvature in
            each iteration of the neutral axe angle.

        Returns:
            eps_a (float): Strain at (0,0)
            chiy (float): Curvature of the section around y-axis
            chiz (float): Curvature of the section around z-axis
        """

        def interpolate(x1, x2, y1, y2, x):
            if x2 - x1 == 0:
                return y1
            y = y1 + (y2 - y1) * (x - x1) / (x2 - x1)
            return y

        def angle(x, y):
            """Obtain the angle of the vector with respect to (1,0) in
            counterclockwise direction.
            """
            angle = math.atan2(y, x)

            # check angle between [0, 2π]
            if angle < 0:
                angle += 2 * math.pi
            return angle

        def check_points_in_2D_diagram(
            boundary_x, boundary_y, forces, debug=False
        ):  # TODO Should this method be moved to the section_results classes?
            """Checks whether given points are inside a boundary defined by
            capacity_pts, and calculates efficiency for rays from the origin
            to each point.

            Parameters:
            boundary_x,boundary_y : array-like, shape (n_points, 1)
                The points defining the capacity boundary.
            forces : array-like, shape (n_forces, 2)
                The points to test.

            Returns:
            results : list of dict
                A list containing dictionaries with the following keys:
                - 'point': The point tested.
                - 'is_inside': True if the point is inside the boundary
                - 'efficiency': The efficiency for the ray from the origin to
                the intersection point.
                - 'intersection_point': The intersection point on the boundary
            """
            from shapely import LineString, Point, Polygon

            # Create boundary points in shape (n_points, 2)
            capacity_pts = np.column_stack((boundary_x, boundary_y))

            # Create the polygon from capacity points
            capacity_polygon = Polygon(capacity_pts)

            # Prepare the results list
            results = []

            # Convert forces to a NumPy array
            forces = np.array(forces)

            # Define the origin
            origin = np.array([0.0, 0.0])

            # Iterate over each point in forces
            for point in forces:
                # Create a line from the origin to the point
                factored_point = np.array([point[0] * 1e10, point[1] * 1e10])
                ray_line = LineString([origin, factored_point])

                # Check if the point is inside the polygon
                is_inside = capacity_polygon.contains(Point(point))

                # Find the intersection of the ray with the polygon boundary
                intersection = ray_line.intersection(capacity_polygon.boundary)

                if intersection.is_empty:
                    # No intersection found
                    efficiency = None
                    intersection_point = None
                else:
                    # There is an intersection
                    # Intersection could be a Point or MultiPoint
                    if isinstance(intersection, Point):
                        intersection_point = np.array(
                            [intersection.x, intersection.y]
                        )
                    elif isinstance(intersection, LineString):
                        # The ray lies along an edge; take the point closest
                        # to the origin
                        coords = np.array(intersection.coords)
                        distances = np.linalg.norm(coords - origin, axis=1)
                        min_index = np.argmin(distances)
                        intersection_point = coords[min_index]
                    elif intersection.geom_type == 'MultiPoint':
                        # Choose the closest intersection point to the origin
                        points = np.array(
                            [[pt.x, pt.y] for pt in intersection.geoms]
                        )
                        distances = np.linalg.norm(points - origin, axis=1)
                        min_index = np.argmin(distances)
                        intersection_point = points[min_index]
                    else:
                        # Unexpected geometry type
                        efficiency = None
                        intersection_point = None

                    if intersection_point is not None:
                        # Calculate distances
                        distance_origin_to_intersection = np.linalg.norm(
                            intersection_point - origin
                        )
                        distance_origin_to_point = np.linalg.norm(
                            point - origin
                        )

                        # Calculate efficiency
                        if distance_origin_to_intersection != 0:
                            efficiency = (
                                distance_origin_to_point
                                / distance_origin_to_intersection
                            )
                        else:
                            efficiency = np.inf  # Avoid division by zero

                        # Determine if the point is inside based on efficiency
                        is_inside = efficiency <= 1

                # Append the result for this point
                results.append(
                    {
                        'point': point,
                        'is_inside': is_inside,
                        'efficiency': efficiency,
                        'intersection_point': intersection_point,
                    }
                )

            if debug:
                # Print the results
                for res in results:
                    point = res['point']
                    is_inside = res['is_inside']
                    efficiency = res['efficiency']
                    print(
                        f"Point {point} is "
                        f"{'inside' if is_inside else 'outside'} the boundary."
                    )
                    if efficiency is not None:
                        print(f'Efficiency: {efficiency:.4f}')
                    if res['intersection_point'] is not None:
                        print(
                            f"Intersection Point: {res['intersection_point']}"
                        )
                    print('---')

            return results

        # Step 1: Check if the section can carry the forces
        mm_domain = self.calculate_mm_interaction_domain(n=n_ed, num_theta=16)
        res = check_points_in_2D_diagram(
            mm_domain.m_y, mm_domain.m_z, [[my_ed, mz_ed]]
        )
        if not res[0]['is_inside']:
            raise ValueError(
                f"Forces cannot be taken by section -> "
                f"efficiency = {res[0]['efficiency']:.3f} > 1"
            )

        # Step 2: Obtain the boundaries of theta (angle of moments)
        # corresponding to the quadrants of the alpha angle
        res = self.calculate_bending_strength(math.pi, n_ed)
        theta_0 = angle(res.m_y, res.m_z)

        res = self.calculate_bending_strength(3 * math.pi / 2, n_ed)
        theta_1 = angle(res.m_y, res.m_z)

        res = self.calculate_bending_strength(0, n_ed)
        theta_2 = angle(res.m_y, res.m_z)

        res = self.calculate_bending_strength(math.pi / 2, n_ed)
        theta_3 = angle(res.m_y, res.m_z)

        # Angle of moments: theta
        theta = angle(my_ed, mz_ed)

        # get first attemp of alfa for iterations and quadrant of aplication
        if (theta_0 <= theta < theta_1) or (0 <= theta < theta_1):
            alfa_1 = 0
            alfa_2 = math.pi / 2
            # alfa = interpolate(theta_0, theta_1, alfa_1, alfa_2, theta)
            if theta_0 <= theta_1:  # theta_0>0
                alfa = interpolate(theta_0, theta_1, alfa_1, alfa_2, theta)
            else:  # theta_0<0
                alfa = interpolate(
                    theta_0 - 2 * math.pi, theta_1, alfa_1, alfa_2, theta
                )
        elif theta_1 <= theta < theta_2:
            alfa_1 = math.pi / 2
            alfa_2 = math.pi
            alfa = interpolate(theta_1, theta_2, alfa_1, alfa_2, theta)
        elif theta_2 <= theta < theta_3:
            alfa_1 = math.pi
            alfa_2 = 3 * math.pi / 2
            alfa = interpolate(theta_0, theta_1, alfa_1, alfa_2, theta)
        else:
            alfa_1 = 3 * math.pi / 2
            alfa_2 = 2 * math.pi
            alfa = interpolate(theta_0, theta_1, alfa_1, alfa_2, theta)

        # Step 3: Iterative process to refine alfa
        ITMAX = 50
        iter = 0
        m_ed = (my_ed**2 + mz_ed**2) ** 0.5
        _chi_pre, _chi_post = None, None
        _theta = -9e9
        while ((abs(_theta - theta)) > 1e-4) and iter < ITMAX:
            if (alfa_2 - alfa_1) > math.pi / 18:  # more than 10° -> simplifies
                mc_r = self.calculate_moment_curvature(
                    alfa + math.pi,
                    n_ed,
                    num_pre_yield=num_chis,
                    num_post_yield=num_chis,
                )
            else:  # close to solutuion
                mc_r = self.calculate_moment_curvature(
                    alfa + math.pi,
                    n_ed,
                    chi=np.linspace(-_chi_pre, -_chi_post, num_chis),
                )
            _M = np.array((mc_r.m_y**2 + mc_r.m_z**2) ** 0.5)

            # Interpolation of curvatures and strain
            chiy = np.interp(m_ed, _M, mc_r.chi_y)
            chiz = np.interp(m_ed, _M, mc_r.chi_z)
            My = np.interp(m_ed, _M, mc_r.m_y)
            Mz = np.interp(m_ed, _M, mc_r.m_z)
            eps_a = np.interp(m_ed, _M, mc_r.eps_axial)

            # region faster algorithm for (alfa_2 - alfa_1)<10º
            _chi = np.array((mc_r.chi_y**2 + mc_r.chi_z**2) ** 0.5)
            _chi_current = np.interp(m_ed, _M, _chi)
            delta_chi = (_chi[-1] - _chi[0]) * 0.20
            _chi_pre = _chi_current - delta_chi
            _chi_post = _chi_current + delta_chi
            # endregion
            print(
                f'My {round(My/1e6)} - '
                f'Mz {round(Mz/1e6)} - '
                f'alfa1 {round(math.degrees(alfa_1),1)} - '
                f'alfa2 {round(math.degrees(alfa_2),1)} - '
            )

            _theta = angle(My, Mz)
            if _theta > theta:
                alfa_2 = alfa
            else:
                alfa_1 = alfa

            alfa = 0.5 * (alfa_1 + alfa_2)
            iter += 1
        if iter == ITMAX:
            s = f'Last iteration reached: \
                My = {My} Mz = {Mz})'
            raise ValueError(f'Maximum number of iterations reached.\n{s}')
        else:
            return eps_a, chiy, chiz
