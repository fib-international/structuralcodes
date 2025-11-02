"""The Marin section integrator."""

from __future__ import annotations  # To have clean hints of ArrayLike in docs

import typing as t
from math import atan2, cos, sin

import numpy as np
from numpy.typing import ArrayLike, NDArray
from shapely import MultiLineString, MultiPolygon, Polygon
from shapely.geometry.polygon import orient

from structuralcodes.core._marin_integration import marin_integration
from structuralcodes.core.base import ConstitutiveLaw
from structuralcodes.geometry import (
    CompoundGeometry,
    SurfaceGeometry,
    create_line_point_angle,
)

from ._section_integrator import SectionIntegrator


class MarinIntegrator(SectionIntegrator):
    """Section integrator based on the Marin algorithm."""

    def _rotate_geometry(
        self, geo: CompoundGeometry, strain: ArrayLike, **kwargs
    ) -> t.Tuple[float, CompoundGeometry, ArrayLike]:
        """This function returns the geometry, strain and angle in a rotated
        CRS for which bending is uniaxial.
        """
        # Rotate section in order to have neutral axis horizontal
        angle = -atan2(strain[2], strain[1])

        rotated_geom = geo.rotate(angle)
        # If integration_data is present, rotate it too
        integration_data = kwargs.get('integration_data')
        if integration_data is not None:
            rotated_integration_data = []
            for reinf_data in integration_data:
                T = np.array(
                    [
                        [np.cos(angle), -np.sin(angle)],
                        [np.sin(angle), np.cos(angle)],
                    ]
                )
                coords = np.vstack((reinf_data[0], reinf_data[1]))
                rotated_coords = T @ coords
                rotated_integration_data.append(
                    (
                        rotated_coords[0, :],
                        rotated_coords[1, :],
                        reinf_data[2],
                        reinf_data[3],
                    )
                )
        else:
            rotated_integration_data = None
        # Determine strain in this rotated CRS
        strain_rotated = [strain[0], (strain[2] ** 2 + strain[1] ** 2) ** 0.5]

        return angle, rotated_geom, strain_rotated, rotated_integration_data

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

    def _get_coefficcients(
        self,
        geo: SurfaceGeometry,
        strain: ArrayLike,
        integrate: t.Literal['stress', 'modulus'] = 'stress',
    ) -> t.Tuple[t.List[t.Tuple], t.List[t.Tuple]]:
        """Get Marin coefficients."""
        if integrate == 'stress':
            if hasattr(geo.material.constitutive_law, '__marin__'):
                strains, coeffs = geo.material.constitutive_law.__marin__(
                    strain=strain
                )
            else:
                raise AttributeError(
                    'The constitutive law object '
                    f'{geo.material.constitutive_law} of geometry {geo} does '
                    'not have implement the __marin__ function. Please '
                    'implement the function or use another integrator, like '
                    'Fiber.'
                )
        elif integrate == 'modulus':
            if hasattr(geo.material.constitutive_law, '__marin_tangent__'):
                strains, coeffs = (
                    geo.material.constitutive_law.__marin_tangent__(
                        strain=strain
                    )
                )
            else:
                raise AttributeError(
                    'The constitutive law object '
                    f'{geo.material.constitutive_law} of geometry {geo} does '
                    'not have implement the __marin_tangent__ function. Please'
                    ' implement the function or use another integrator, like '
                    'Fibre.'
                )
        else:
            raise ValueError(f'Unknown integrate type: {integrate}')

        return strains, coeffs

    def _process_surface_geometries(
        self,
        geo: CompoundGeometry,
        strain: ArrayLike,
        input: t.List,
        integrate: t.Literal['stress', 'modulus'] = 'stress',
    ):
        """Process Surface geometries filling the input data for each one."""
        # For each SurfaceGeometry on the CompoundGeometry:
        for g in geo.geometries:
            # Get coefficients and strain limits from constitutive law
            strains, coeffs = self._get_coefficcients(g, strain, integrate)

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

    def _create_integration_data_point_geometries(
        self, geo: CompoundGeometry
    ) -> t.List[t.Tuple[np.ndarray, np.ndarray, np.ndarray, ConstitutiveLaw]]:
        """Return integration data structure for storing point geometries
        in numpy arrays.
        """
        # Initialize data structure
        integration_data = []
        # Map reinf_data to group point geometries by constitutive law
        reinf_data = {}
        # Preprocess geometries having the same material
        for pg in geo.point_geometries:
            x, y = pg._point.coords.xy
            x = x[0]
            y = y[0]
            area = pg.area
            constitutive_law = pg.material.constitutive_law
            if reinf_data.get(constitutive_law) is None:
                reinf_data[constitutive_law] = [
                    np.array(x),
                    np.array(y),
                    np.array(area),
                ]
            else:
                reinf_data[constitutive_law][0] = np.hstack(
                    (reinf_data[constitutive_law][0], x)
                )
                reinf_data[constitutive_law][1] = np.hstack(
                    (reinf_data[constitutive_law][1], y)
                )
                reinf_data[constitutive_law][2] = np.hstack(
                    (reinf_data[constitutive_law][2], area)
                )
        for constitutive_law, value in reinf_data.items():
            integration_data.append(
                (value[0], value[1], value[2], constitutive_law)
            )
        return integration_data

    def _process_point_geometries(
        self,
        geo: CompoundGeometry,
        strain: ArrayLike,
        input: t.List,
        integrate: t.Literal['stress', 'modulus'] = 'stress',
        **kwargs,
    ):
        """Process Point geometries filling the input data."""
        # See if there is a integration_data stored
        integration_data = kwargs.get('integration_data')
        if integration_data is None:
            # No integration data is present, process it now and store it
            # for future use
            integration_data = self._create_integration_data_point_geometries(
                geo=geo
            )
        # Now use the integration data (either just created or passed as
        # argument) to fill the input
        y = []
        z = []
        IA = []
        for reinf_data in integration_data:
            # All have the same constitutive law
            strains = strain[0] + strain[1] * reinf_data[1]
            # Compute integration quantity (stress or modulus)
            if integrate == 'stress':
                integrand = reinf_data[3].get_stress(strains)
            elif integrate == 'modulus':
                integrand = reinf_data[3].get_tangent(strains)
            else:
                raise ValueError(f'Unknown integrate type: {integrate}')
            y.append(reinf_data[0])
            z.append(reinf_data[1])
            IA.append(integrand * reinf_data[2])

        input.append((1, np.array(y), np.array(z), np.array(IA)))
        return integration_data

    def prepare_input(
        self,
        geo: CompoundGeometry,
        strain: ArrayLike,
        integrate: t.Literal['stress', 'modulus'] = 'stress',
        **kwargs,
    ) -> t.Tuple[float, t.Tuple[np.ndarray, np.ndarray, np.ndarray]]:
        """Prepare general input to the integration of stress or material
        modulus in the section.

        Calculate the stress resultants or tangent section stiffness based on
        strains in a set of points.

        Keyword Arguments:
            geo (CompoundGeometry): The geometry of the section.
            strain (ArrayLike): The strains and curvatures of the section,
                given in the format (ea, ky, kz) which are i) strain at 0,0,
                ii) curvature y axis, iii) curvature z axis.
            integrate (str): a string indicating the quantity to integrate over
                the section. It can be 'stress' or 'modulus'. When 'stress'
                is selected, the return value will be the stress resultants N,
                My, Mz, while if 'modulus' is selected, the return will be the
                section stiffness matrix (default is 'stress').

        Returns:
            Tuple(float, Tuple(ndarray, ndarray, ndarray)): The prepared input
            represented as the angle of rotation computed (needed for rotating
            back the resultants) and as a tuple with 3 ndarrys collecting
            respectively y, z and stress coefficients for each sub-part.

        Raises:
            ValueError: If a unkown value is passed to the `integrate`
            parameter.
        """
        # Get the integration data if passed
        integration_data = kwargs.get('integration_data')
        # The method should therefore return a tuple that collects the y, z,
        # and stress coefficients for each part.
        prepared_input = []
        # Rotate section in order to have neutral axis horizontal
        angle, rotated_geom, strain_rotated, rotated_integration_data = (
            self._rotate_geometry(
                geo, strain, integration_data=integration_data
            )
        )

        # Process all the surface geometries splitting them and appending
        # coefficients and strain limits to prepared_input
        self._process_surface_geometries(
            rotated_geom, strain_rotated, prepared_input, integrate
        )

        # Process all the point geometries (i.e. reinforcement) splitting them
        # and appending coefficients and strain limits to prepared_input
        integration_data = self._process_point_geometries(
            rotated_geom,
            strain_rotated,
            prepared_input,
            integrate,
            integration_data=rotated_integration_data,
        )

        return angle, prepared_input, integration_data

    def integrate_stress(
        self,
        angle: float,
        prepared_input: t.List[
            t.Tuple[int, np.ndarray, np.ndarray, np.ndarray]
        ],
        **kwargs,
    ) -> t.Tuple[float, float, float]:
        """Integrate stresses over the geometry.

        Arguments:
            prepared_input (List): The prepared input from .prepare_input().

        Returns:
            Tuple(float, float, float): The stress resultants N, My and Mz.
        """
        # Set the stress resultants to zero
        N, My, Mz = 0.0, 0.0, 0.0

        # Loop through all parts of the section and add contributions
        for i, y, z, stress_coeff in prepared_input:
            if i == 0:
                # Find integration order from shape of stress coeff array
                n = stress_coeff.shape[0]

                # Calculate area moments
                area_moments_N = np.array(
                    [marin_integration(y, z, 0, k) for k in range(n)]
                )
                area_moments_My = np.array(
                    [marin_integration(y, z, 0, k + 1) for k in range(n)]
                )
                area_moments_Mz = np.array(
                    [marin_integration(y, z, 1, k) for k in range(n)]
                )

                # Calculate contributions to stress resultants
                N += sum(stress_coeff * area_moments_N)
                My += sum(stress_coeff * area_moments_My)
                Mz -= sum(stress_coeff * area_moments_Mz)
            elif i == 1:
                # Reinforcement
                N += np.sum(stress_coeff)
                My += np.sum(stress_coeff * z)
                Mz -= np.sum(stress_coeff * y)

        # Rotate back to section CRS
        T = np.array([[cos(-angle), -sin(-angle)], [sin(-angle), cos(-angle)]])
        M = T @ np.array([[My], [Mz]])

        # Rotate back the integration_data
        integration_data = kwargs.get('integration_data')
        if integration_data is not None:
            rotated_integration_data = []
            for reinf_data in integration_data:
                coords = np.vstack((reinf_data[0], reinf_data[1]))
                rotated_coords = T @ coords
                rotated_integration_data.append(
                    (
                        rotated_coords[0, :],
                        rotated_coords[1, :],
                        reinf_data[2],
                        reinf_data[3],
                    )
                )
        else:
            rotated_integration_data = None

        return N, M[0, 0], M[1, 0], rotated_integration_data

    def integrate_modulus(
        self,
        angle: float,
        prepared_input: t.List[
            t.Tuple[int, np.ndarray, np.ndarray, np.ndarray]
        ],
        **kwargs,
    ) -> NDArray[np.float64]:
        """Integrate material modulus over the geometry.

        Arguments:
            prepared_input (List): The prepared input from .prepare_input().

        Returns:
            ndarray: The section stiffness matrix as a (3, 3) ndarray.
        """
        # Create the stiffness matrix
        stiffness = np.zeros((3, 3), dtype=np.float64)

        # Create a map of indices
        # (i,j), m, offset, sign
        indices = [
            ((0, 0), 0, 0, 1),
            ((0, 1), 0, 1, 1),
            ((0, 2), 1, 0, -1),
            ((1, 1), 0, 2, 1),
            ((1, 2), 1, 1, -1),
            ((2, 2), 2, 0, 1),
        ]

        # Loop through all parts of the section and add contributions
        for id, y, z, modulus_coeff in prepared_input:
            if id == 0:
                # Find integration order from shape of stress coeff array
                n = modulus_coeff.shape[0]

                # Calculate needed area moments
                for (i, j), m, offset, sign in indices:
                    area_moments = np.array(
                        [
                            marin_integration(y, z, m, k + offset)
                            for k in range(n)
                        ]
                    )
                    stiffness[i, j] += sign * sum(modulus_coeff * area_moments)
            elif id == 1:
                # Reinforcement
                stiffness[0, 0] += np.sum(modulus_coeff)
                stiffness[0, 1] += np.sum(modulus_coeff * z)
                stiffness[0, 2] -= np.sum(modulus_coeff * y)
                stiffness[1, 1] += np.sum(modulus_coeff * z * z)
                stiffness[1, 2] -= np.sum(modulus_coeff * y * z)
                stiffness[2, 2] += np.sum(modulus_coeff * y * y)

        # Apply for simmetry
        stiffness[1, 0] = stiffness[0, 1]
        stiffness[2, 0] = stiffness[0, 2]
        stiffness[2, 1] = stiffness[1, 2]

        # Rotate back to section CRS
        T = np.array([[cos(-angle), -sin(-angle)], [sin(-angle), cos(-angle)]])
        # Rotate back the integration_data if present
        integration_data = kwargs.get('integration_data')
        if integration_data is not None:
            T = np.array(
                [[cos(-angle), -sin(-angle)], [sin(-angle), cos(-angle)]]
            )
            rotated_integration_data = []
            for reinf_data in integration_data:
                coords = np.vstack((reinf_data[0], reinf_data[1]))
                rotated_coords = T @ coords
                rotated_integration_data.append(
                    (
                        rotated_coords[0, :],
                        rotated_coords[1, :],
                        reinf_data[2],
                        reinf_data[3],
                    )
                )
        else:
            rotated_integration_data = None

        T = np.array(
            [
                [1, 0, 0],
                [0, cos(-angle), sin(-angle)],
                [0, -sin(-angle), cos(-angle)],
            ]
        )

        return T.T @ stiffness @ T, rotated_integration_data

    def integrate_strain_response_on_geometry(
        self,
        geo: CompoundGeometry,
        strain: ArrayLike,
        integrate: t.Literal['stress', 'modulus'] = 'stress',
        **kwargs,
    ):
        """Integrate strees or material modulus in the section with the marin
        algorithm.

        Calculate the stress resultants or tangent section stiffness based on
        strains in a set of points.

        Arguments:
            geo (CompoundGeometry): The geometry of the section.
            strain (ArrayLike): The strains and curvatures of the section,
                given in the format (ea, ky, kz) which are i) strain at 0,0,
                ii) curvature y axis, iii) curvature z axis.
            integrate (str): a string indicating the quantity to integrate over
                the section. It can be 'stress' or 'modulus'. When 'stress'
                is selected, the return value will be the stress resultants N,
                My, Mz, while if 'modulus' is selected, the return will be the

        Returns:
            Tuple(Union(Tuple(float, float, float), np.ndarray), None): The
            first element is either a tuple of floats (for the stress
            resultants (N, My, Mz) when `integrate='stress'`, or a numpy
            array representing the stiffness matrix then `integrate='modulus'`.

        Example:
            result, _ = integrate_strain_response_on_geometry(geo, strain,
            integrate='tanent')
            # `result` will be the stiffness matrix (a 3x3 numpy array) if
            # `integrate='modulus'`

        Raises:
            ValueError: If a unkown value is passed to the `integrate`
            parameter.
        """
        # Prepare the general input based on the geometry and the input strains
        angle, prepared_input, integration_data = self.prepare_input(
            geo, strain, integrate, **kwargs
        )

        # Return the calculated response
        if integrate == 'stress':
            return self.integrate_stress(
                angle, prepared_input, integration_data=integration_data
            )
        if integrate == 'modulus':
            return self.integrate_modulus(
                angle, prepared_input, integration_data=integration_data
            )
        raise ValueError(f'Unknown integrate type: {integrate}')
