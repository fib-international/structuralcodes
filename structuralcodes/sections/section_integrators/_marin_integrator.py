"""The Marin section integrator."""

from __future__ import annotations  # To have clean hints of ArrayLike in docs

import typing as t
from math import atan2, cos, sin

import numpy as np
from numpy.typing import ArrayLike, NDArray
from shapely import MultiLineString, MultiPolygon, Polygon
from shapely.geometry.polygon import orient

from structuralcodes.geometry import CompoundGeometry, create_line_point_angle

from ._marin_integration import marin_integration
from ._section_integrator import SectionIntegrator


class MarinIntegrator(SectionIntegrator):
    """Section integrator based on the Marin algorithm."""

    def _rotate_geometry(
        self, geo: CompoundGeometry, strain: ArrayLike
    ) -> t.Tuple[float, CompoundGeometry, ArrayLike]:
        """This function returns the geometry, strain and angle in a rotated
        CRS for which bending is uniaxial.
        """
        # Rotate section in order to have neutral axis horizontal
        angle = -atan2(strain[2], strain[1])

        rotated_geom = geo.rotate(angle)
        # Determine strain in this rotated CRS
        strain_rotated = [strain[0], (strain[2] ** 2 + strain[1] ** 2) ** 0.5]

        return angle, rotated_geom, strain_rotated

    def _get_input_polygon(
        self, polygon: Polygon, coeffs: ArrayLike, input: list
    ):
        """Appends to input list the coordinates and coefficient of polygon."""
        # Let's be sure to orient in the right way
        if polygon.is_empty:
            return
        polygon = orient(polygon, 1)
        if not polygon.exterior.is_ccw:
            raise ValueError(
                'The exterior of a polygon should have vertices \
                            ordered ccw'
            )
        # Manage exterior part
        x, y = polygon.exterior.coords.xy
        x = np.array(x)
        y = np.array(y)
        input.append((0, np.array(x), np.array(y), np.array(coeffs)))
        # Manage holes
        for i in polygon.interiors:
            if i.is_ccw:
                raise ValueError('A inner hole should have cw coordinates')
            x, y = i.coords.xy
            input.append((0, np.array(x), np.array(y), np.array(coeffs)))

    def _process_surface_geometries(
        self, geo: CompoundGeometry, strain: ArrayLike, input: t.List
    ):
        """Process Surface geometries filling the input data for each one."""
        # For each SurfaceGeometry on the CompoundGeometry:
        for g in geo.geometries:
            # Get coefficients and strain limits from constitutive law
            if hasattr(g.material, '__marin__'):
                strains, coeffs = g.material.__marin__(strain=strain)
            else:
                raise AttributeError(
                    f'The material object {g.material} of geometry {g} does \
                    not have implement the __marin__ function. \
                    Please implement the function or use another integrator, \
                    like '
                    'Fibre'
                    ''
                )

            # Subdivide the polygon at the different strain limits
            if strains is None:
                self._get_input_polygon(g.polygon, coeffs[0], input)
            else:
                for p in range(len(strains)):
                    # Create the two lines for selecting the needed part
                    y0 = -(strain[0] - strains[p][0]) / strain[1]
                    y1 = -(strain[0] - strains[p][1]) / strain[1]
                    if y0 > y1:
                        y0, y1 = y1, y0
                    bbox = g.polygon.bounds
                    line0 = create_line_point_angle((0, y0), 0, bbox)
                    line1 = create_line_point_angle((0, y1), 0, bbox)
                    lines = MultiLineString((line0, line1))
                    # intersection
                    result = g.split_two_lines(lines=lines)

                    if isinstance(result, Polygon):
                        # If the result is a single polygon
                        self._get_input_polygon(result, coeffs[p], input)
                    elif isinstance(result, MultiPolygon):
                        # If the result is a MultiPolygon
                        for polygon in result.geoms:
                            self._get_input_polygon(polygon, coeffs[p], input)

    def _process_point_geometries(
        self, geo: CompoundGeometry, strain: ArrayLike, input: t.List
    ):
        """Process Point geometries filling the input data."""
        # Tentative proposal for managing reinforcement (PointGeometry)
        x = []
        y = []
        F = []
        for pg in geo.point_geometries:
            xp, yp = pg._point.coords.xy
            xp = xp[0]
            yp = yp[0]
            A = pg.area
            strain_ = strain[0] + strain[1] * yp
            x.append(xp)
            y.append(yp)
            F.append(pg.material.get_stress(strain_) * A)
        input.append((1, np.array(x), np.array(y), np.array(F)))

    def prepare_input_stress(
        self, geo: CompoundGeometry, strain: ArrayLike
    ) -> t.Tuple[float, t.Tuple[np.ndarray, np.ndarray, np.ndarray]]:
        """Prepare general input to the stress integration.

        Calculate the stresses based on strains in a set of points.

        Keyword Arguments:
            geo (CompoundGeometry): The geometry of the section.
            strain (ArrayLike): The strains and curvatures of the section,
                given in the format (ea, ky, kz) which are i) strain at 0,0,
                ii) curvature y axis, iii) curvature z axis.
            mesh_size: Percentage of area (number from 0 to 1) max for triangle
                elements.

        Returns:
            Tuple(float, Tuple(ndarray, ndarray, ndarray)): The prepared input
            represented as the angle of rotation computed (needed for rotating
            back the resultants) and as a tuple with 3 ndarrys collecting
            respectively y, z and stress coefficients for each sub-part.
        """
        # The method should therefore return a tuple that collects the y, z,
        # and stress coefficients for each part.
        prepared_input = []
        # Rotate section in order to have neutral axis horizontal
        angle, rotated_geom, strain_rotated = self._rotate_geometry(
            geo, strain
        )

        # Process all the surface geometries splitting them and appending
        # coefficients and strain limits to prepared_input
        self._process_surface_geometries(
            rotated_geom, strain_rotated, prepared_input
        )

        # Process all the point geometries (i.e. reinforcement) splitting them
        # and appending coefficients and strain limits to prepared_input
        self._process_point_geometries(
            rotated_geom, strain_rotated, prepared_input
        )

        return angle, prepared_input

    def integrate_stress(
        self,
        angle: float,
        prepared_input: t.List[
            t.Tuple[int, np.ndarray, np.ndarray, np.ndarray]
        ],
    ) -> t.Tuple[float, float, float]:
        """Integrate stresses over the geometry.

        Arguments:
            prepared_input (List): The prepared input from .prepare_input().

        Returns:
            Tuple(float, float, float): The stress resultants N, Mx and My.
        """
        # Set the stress resultants to zero
        N, Mx, My = 0.0, 0.0, 0.0

        # Loop through all parts of the section and add contributions
        for i, y, z, stress_coeff in prepared_input:
            if i == 0:
                # Find integration order from shape of stress coeff array
                n = stress_coeff.shape[0]

                # Calculate area moments
                area_moments_N = np.array(
                    [marin_integration(y, z, 0, k) for k in range(n)]
                )
                area_moments_Mx = np.array(
                    [marin_integration(y, z, 0, k + 1) for k in range(n)]
                )
                area_moments_My = np.array(
                    [marin_integration(y, z, 1, k) for k in range(n)]
                )

                # Calculate contributions to stress resultants
                N += sum(stress_coeff * area_moments_N)
                Mx += sum(stress_coeff * area_moments_Mx)
                My += sum(stress_coeff * area_moments_My)
            elif i == 1:
                # Reinforcement
                N += sum(stress_coeff)
                Mx += sum(stress_coeff * z)
                My += sum(stress_coeff * y)

        # Rotate back to section CRS
        T = np.array([[cos(-angle), -sin(-angle)], [sin(-angle), cos(-angle)]])
        M = T @ np.array([[Mx], [-My]])

        return N, M[0, 0], M[1, 0]

    def integrate_strain_response_on_geometry_stress(
        self, geo: CompoundGeometry, strain: ArrayLike, **kwargs
    ):
        """Integrate the strees from strain response with the Marin algorithm.

        Arguments:
            geo (CompoundGeometry): The geometry of the section.
            strain (ArrayLike): The strains and curvatures of the section,
                given in the format (ea, ky, kz) which are i) strain at 0,0,
                ii) curvature y axis, iii) curvature z axis.

        Returns:
            Tuple(Tuple(float, float, float), Dict): The stress resultants N,
            Mx and My and the triangulation data.
        """
        del kwargs
        # Prepare the general input based on the geometry and the input strains
        angle, prepared_input = self.prepare_input_stress(geo, strain)

        # Return the calculated response
        return *self.integrate_stress(angle, prepared_input), None

    def prepare_input_tangent(
        self, geo: CompoundGeometry, strain: ArrayLike
    ) -> t.Tuple[float, t.Tuple[np.ndarray, np.ndarray, np.ndarray]]:
        """Prepare general input to the tangent modulus integration.

        Keyword Arguments:
            geo (CompoundGeometry): The geometry of the section.
            strain (ArrayLike): The strains and curvatures of the section,
                given in the format (ea, ky, kz) which are i) strain at 0,0,
                ii) curvature y axis, iii) curvature z axis.
            mesh_size: Percentage of area (number from 0 to 1) max for triangle
                elements.

        Returns:
            Tuple(float, Tuple(ndarray, ndarray, ndarray)): The prepared input
            represented as the angle of rotation computed (needed for rotating
            back the resultants) and as a tuple with 3 ndarrys collecting
            respectively y, z and stress coefficients for each sub-part.
        """
        del geo, strain
        raise NotImplementedError

    def integrate_tangent(
        self,
        angle: float,
        prepared_input: t.List[
            t.Tuple[int, np.ndarray, np.ndarray, np.ndarray]
        ],
    ) -> NDArray:
        """Integrate tangent modulus over the geometry.

        Arguments:
            prepared_input (List): The prepared input from .prepare_input().

        Returns:
            NDArray: The tangent modulus resultant as section stiffness.
        """
        del angle, prepared_input
        raise NotImplementedError

    def integrate_strain_response_on_geometry_tangent(
        self, geo: CompoundGeometry, strain: ArrayLike, **kwargs
    ):
        """Integrate the tangent from strain response with the Marin algorithm.

        Arguments:
            geo (CompoundGeometry): The geometry of the section.
            strain (ArrayLike): The strains and curvatures of the section,
                given in the format (ea, ky, kz) which are i) strain at 0,0,
                ii) curvature y axis, iii) curvature z axis.

        Returns:
            NDArray: The section tangent stiffness.
        """
        del geo, strain, kwargs
        raise NotImplementedError
