"""Generic class section implemenetation."""

from __future__ import annotations

import typing as t
from math import cos, sin

import numpy as np

import structuralcodes.core._section_results as s_res
from structuralcodes.core.base import Section, SectionCalculator
from structuralcodes.geometry import (
    CompoundGeometry,
    PointGeometry,
    SurfaceGeometry,
)

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
        section_analyzer: (GenericSectionCalculator)
            the object responsible for performing different
            computations on the section (e.g. bending strength,
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
        self.geometry = geometry
        self.section_analyzer = GenericSectionCalculator(
            sec=self, integrator=integrator, **kwargs
        )
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
        self.mesh_size = kwargs.get('mesh_size', 0.001)
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
        # Here compute area, centroid and other stuff from gross properties
        # 1. Computation of area
        gp.area = self.section.geometry.area
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
        ITMAX = 20
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
            )
        )
        n_max, _, _, _ = self.integrator.integrate_strain_response_on_geometry(
            self.section.geometry, [eps_p, 0, 0], tri=self.triangulated_data
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
            # Morten, I would like to raise this error!:
            # raise ValueError(error_str)

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
        rotated_geom = self.section.geometry.rotate(-theta)
        if self.triangulated_data is not None:
            # Rotate also triangulated data!
            self._rotate_triangulated_data(-theta)
        # Find ultimate curvature from the strain distribution corresponding
        # to failure and equilibrium with external axial force
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
                'curvature at yield and ultimate cannot have opposite signs!'
            )
        # Define the array of curvatures
        if abs(chi_ultimate) <= abs(chi_yield) + 1e-8:
            # We don't want a plastic branch in the analysis
            # this is done to speed up analysis
            chi = np.linspace(chi_yield / 10, chi_yield, 10)
        else:
            chi = np.concatenate(
                (
                    np.linspace(chi_yield / 10, chi_yield, 10, endpoint=False),
                    np.linspace(chi_yield, chi_ultimate, 100),
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

    def calculate_interaction_domain(
        self, theta: float = 0, n_axial: int = 100
    ) -> s_res.NmInteractionDomain:
        """Calculate the NM interaction domain.

        Arguments:
        theta (float, default = 0): inclination of n.a. respect to y axis
        n_axial (int, default = 100): number of discretization for axial load

        Return:
        (NmInteractionDomain)
        """
        # Prepare the results
        res = s_res.NmInteractionDomain()
        res.theta = theta
        res.n_axial = n_axial
        # Create an array of axial loads
        res.n = np.linspace(self.n_min, self.n_max, n_axial)
        # Initialize the result's arrays
        res.m_y = np.zeros_like(res.n)
        res.m_z = np.zeros_like(res.n)
        # Compute strength for given angle of NA
        for i, n_i in enumerate(res.n):
            res_bend_strength = self.calculate_bending_strength(
                theta=theta, n=n_i
            )
            res.m_y[i] = res_bend_strength.m_y
            res.m_z[i] = res_bend_strength.m_z

        return res

    def calculate_interaction_domain_2(
        self, theta: float = 0, n_axial: int = 100
    ) -> s_res.NmInteractionDomain:
        """Calculate the NM interaction domain (attempt to make it faster).

        Arguments:
        theta (float, default = 0): inclination of n.a. respect to y axis
        n_axial (int, default = 100): number of discretization for axial load

        Return:
        (NmInteractionDomain)
        """
        # Prepare the results
        res = s_res.NmInteractionDomain()
        res.theta = theta
        # with this version really n_axial is not used unless we perform
        # some resampling after computation?
        res.n_axial = n_axial

        # Get ultimate strain profiles for theta angle
        strains = self._compute_ultimate_strain_profiles(theta=theta)

        # integrate all strain profiles
        forces = np.zeros_like(strains)
        for i in range(strains.shape[0]):
            N, My, Mz, tri = (
                self.integrator.integrate_strain_response_on_geometry(
                    geo=self.section.geometry,
                    strain=strains[i, :],
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
        # Find balanced failure: this defines the transition
        y_n, y_p, strain = self.get_balanced_failure_strain(
            geom=rotated_geom, yielding=True
        )
        eps_p_y = strain[0] + strain[1] * y_p
        eps_n_y = strain[0] + strain[1] * y_n
        # between fields 2 and 3
        y_n, y_p, strain = self.get_balanced_failure_strain(
            geom=rotated_geom, yielding=False
        )
        _, _, min_y, _ = rotated_geom.calculate_extents()
        h = y_n - min_y
        eps_p_b = strain[0] + strain[1] * y_p
        eps_n_b = strain[0] + strain[1] * y_n

        eps_p = []
        eps_n = []
        # For generation of fields 1 and 2 pivot on steel
        # Field 1: pivot on steel
        n = 3
        eps_n = np.linspace(eps_p_b, 0, n, endpoint=False)
        eps_p = np.zeros_like(eps_n) + eps_p_b
        # Field 2: pivot on steel
        n = 10  # 10
        eps_n = np.append(eps_n, np.linspace(0, eps_n_b, n, endpoint=False))
        eps_p = np.append(eps_p, np.zeros(n) + eps_p_b)
        # Field 3: pivot on concrete
        n = 40  # 40
        eps_n = np.append(eps_n, np.zeros(n) + eps_n_b)
        eps_p = np.append(
            eps_p, np.linspace(eps_p_b, eps_p_y, n, endpoint=False)
        )
        # Field 4: pivot on concrete
        n = 6  # 6
        eps_n = np.append(eps_n, np.zeros(n) + eps_n_b)
        eps_p = np.append(eps_p, np.linspace(eps_p_y, 0, n, endpoint=False))
        # Field 5: pivot on concrete
        n = 8  # 8
        eps_p_lim = eps_n_b * (h - (y_n - y_p)) / h
        eps_n = np.append(eps_n, np.zeros(n) + eps_n_b)
        eps_p = np.append(eps_p, np.linspace(0, eps_p_lim, n, endpoint=False))
        # Field 6: pivot on eps_n_y point! Morten:
        # How to make this general? This is concrete!
        n = 8  # 8
        z_pivot = y_n - (1 - eps_n_y / eps_n_b) * h
        eps_p_6 = np.linspace(eps_p_lim, eps_n_y, n, endpoint=True)
        eps_n_6 = (
            -(eps_n_y - eps_p_6) * (z_pivot - y_n) / (z_pivot - y_p) + eps_n_y
        )
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
        self, n_theta: int = 32, n_axial: int = 20
    ) -> s_res.NmmInteractionDomain:
        """Calculates the NMM interaction domain."""
        res = s_res.NmmInteractionDomain()
        res.n_theta = n_theta
        # TODO: if we want to resample and structure better data we can
        # use n_axial as the number of discretizations in axial direction
        res.n_axial = n_axial

        # cycle for all n_thetas
        thetas = np.linspace(0, np.pi * 2, n_theta)
        # Initialize an empty array with the correct shape
        strains = np.empty((0, 3))
        for theta in thetas:
            # Get ultimate strain profiles for theta angle
            strain = self._compute_ultimate_strain_profiles(theta=theta)
            strains = np.vstack((strains, strain))

        # integrate all strain profiles
        forces = np.zeros_like(strains)
        for i in range(strains.shape[0]):
            N, My, Mz, tri = (
                self.integrator.integrate_strain_response_on_geometry(
                    geo=self.section.geometry,
                    strain=strains[i, :],
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
        self, n: float = 0, n_theta: int = 32
    ) -> s_res.MmInteractionDomain:
        """Calculate the My-Mz interaction domain.

        Arguments:
        n (float, default = 0): axial force
        n_theta (int, default = 32): number of discretization for theta

        Return:
        (MmInteractionDomain)
        """
        # Prepare the results
        res = s_res.NmInteractionDomain()
        res.n_theta = n_theta
        res.n = n
        # Create array of thetas
        res.theta = np.linspace(0, np.pi * 2, n_theta)
        # Initialize the result's arrays
        res.m_y = np.zeros_like(res.theta)
        res.m_z = np.zeros_like(res.theta)
        # Compute strength for given angle of NA
        for i, th in enumerate(res.theta):
            res_bend_strength = self.calculate_bending_strength(theta=th, n=n)
            res.m_y[i] = res_bend_strength.m_y
            res.m_z[i] = res_bend_strength.m_z

        return res


# Use examples:

# sec = Section(...xxx...)
# res = sec.section_analyzer.ultmateBendingStrength(n = -1000)
# stress = sec.section_analyzer.compute_stress_distribution(res)
