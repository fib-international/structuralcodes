"""The fiber section integrator."""

from __future__ import annotations  # To have clean hints of ArrayLike in docs

import typing as t

import numpy as np
import triangle
from numpy.typing import ArrayLike, NDArray
from shapely import Polygon

from structuralcodes.core.base import ConstitutiveLaw
from structuralcodes.geometry import CompoundGeometry, SurfaceGeometry

from ._section_integrator import SectionIntegrator


class FiberIntegrator(SectionIntegrator):
    """Section integrator based on the Marin algorithm."""

    def prepare_triangulation(self, geo: SurfaceGeometry) -> t.Dict:
        """Prepare data for triangulating it with triangle.

        Arguments:
            geo (SurfaceGeometry): The geometry to triangulate.

        Returns:
            Dict: The triangulation data.
        """
        # Create the tri dictionary
        tri: dict[str:ArrayLike] = {}
        # 1. External boundary process
        # 1a. Get vertices, skipping the last one
        vertices = np.column_stack(geo.polygon.exterior.xy)[:-1, :]
        n_vertices = vertices.shape[0]
        # 1b. Create segments
        node_i = np.arange(n_vertices)
        node_j = np.roll(node_i, -1)
        segments = np.column_stack((node_i, node_j))

        # 2. Process holes
        holes = []
        for interior in geo.polygon.interiors:
            # 2a. Get vertices, skipping the last one
            vertices_int = np.column_stack(interior.xy)[:-1, :]
            n_vertices_int = vertices_int.shape[0]
            # 2b. Create segments
            node_i = np.arange(n_vertices_int) + n_vertices
            node_j = np.roll(node_i, -1)
            segments_int = np.column_stack((node_i, node_j))
            c = Polygon(interior)
            holes.append([c.centroid.x, c.centroid.y])
            # Append to the global arrays
            vertices = np.vstack((vertices, vertices_int))
            segments = np.vstack((segments, segments_int))
            n_vertices += n_vertices_int
        # Return the dictionary with data for triangulate
        tri['vertices'] = vertices
        tri['segments'] = segments
        if len(holes) > 0:
            tri['holes'] = holes
        return tri

    def triangulate(
        self, geo: CompoundGeometry, mesh_size: float
    ) -> t.List[t.Tuple[np.ndarray, np.ndarray, np.ndarray, ConstitutiveLaw]]:
        """Triangulate the geometry discretizing it into fibers.

        Arguments:
            geo (CompoundGeometry): The geometry of the section.
            mesh_size (float): Percentage of area (number from 0 to 1) max for
                triangle elements.

        Raises:
            ValueError if mesh_size is <= 0.0 or if mesh_size is > 1.0.
        """
        # Check mesh_size being between 0 and 1.
        if mesh_size <= 0 or mesh_size > 1:
            raise ValueError('mesh_size is a number from 0 to 1')

        triangulated_data = []

        # For the surface geometries
        for g in geo.geometries:
            # prepare data structure for triangle module
            tri = self.prepare_triangulation(g)
            # define the maximum area of the triangles
            max_area = g.area * mesh_size
            # triangulate the geometry getting back the mesh
            mesh = triangle.triangulate(tri, f'pq{30:.1f}Aa{max_area}o1')
            mat = g.material
            # Get x and y coordinates (centroid) and area for each fiber
            x = []
            y = []
            area = []
            for tr in mesh['triangles']:
                # get centroid of triangle
                xc = (
                    mesh['vertices'][tr[0]][0]
                    + mesh['vertices'][tr[1]][0]
                    + mesh['vertices'][tr[2]][0]
                )
                xc /= 3.0
                x.append(xc)
                yc = (
                    mesh['vertices'][tr[0]][1]
                    + mesh['vertices'][tr[1]][1]
                    + mesh['vertices'][tr[2]][1]
                )
                yc /= 3.0
                y.append(yc)
                # compute area
                a = (
                    mesh['vertices'][tr[0]][0] * mesh['vertices'][tr[1]][1]
                    - mesh['vertices'][tr[0]][1] * mesh['vertices'][tr[1]][0]
                )
                a += (
                    mesh['vertices'][tr[1]][0] * mesh['vertices'][tr[2]][1]
                    - mesh['vertices'][tr[1]][1] * mesh['vertices'][tr[2]][0]
                )
                a += (
                    mesh['vertices'][tr[2]][0] * mesh['vertices'][tr[0]][1]
                    - mesh['vertices'][tr[2]][1] * mesh['vertices'][tr[0]][0]
                )
                a = abs(a) * 0.5
                area.append(a)
                # pointer to the material

            # return back the triangulation data
            triangulated_data.append(
                (np.array(x), np.array(y), np.array(area), mat)
            )
        # For the reinforcement
        # Tentative proposal for managing reinforcement (PointGeometry)
        reinf_data = {}
        # Preprocess geometries having the same material
        for pg in geo.point_geometries:
            x, y = pg._point.coords.xy
            x = x[0]
            y = y[0]
            area = pg.area
            mat = pg.material
            if reinf_data.get(mat) is None:
                reinf_data[mat] = [
                    np.array(x),
                    np.array(y),
                    np.array(area),
                ]
            else:
                reinf_data[mat][0] = np.hstack((reinf_data[mat][0], x))
                reinf_data[mat][1] = np.hstack((reinf_data[mat][1], y))
                reinf_data[mat][2] = np.hstack((reinf_data[mat][2], area))
        for mat, value in reinf_data.items():
            triangulated_data.append((value[0], value[1], value[2], mat))

        return triangulated_data

    def prepare_input(
        self,
        geo: CompoundGeometry,
        strain: ArrayLike,
        integrate: t.Literal['stress', 'modulus'] = 'stress',
        **kwargs,
    ) -> t.Tuple[t.Tuple[np.ndarray, np.ndarray, np.ndarray]]:
        """Prepare general input to the integration of stress or material
        modulus in the section.

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
                tangent section stiffness matrix (default is 'stress').

        Keyword Arguments:
            mesh_size (float): Percentage of area (number from 0 to 1) max for
                triangle elements.

        Returns:
            Tuple(List): The prepared input representing a list with
            x-coordinates, y-coordinates and force for each fiber and a
            list containing the triangulation data that can be stored and
            used later to avoid repetition of triangulation.

        Raises:
            ValueError: If a unkown value is passed to the `integrate`
            parameter.
        """
        # This method should:
        #  - discretize the section in a number of fibers (mesh_size)
        #  - prepare input as y coordiates z coordinates forces

        prepared_input = []

        triangulated_data = kwargs.get('tri')
        if triangulated_data is None:
            # No triangulation is provided, triangulate the section
            # Fiber integrator for generic section uses delaunay triangulation
            # for discretizing in fibers

            mesh_size = kwargs.get('mesh_size', 0.01)
            triangulated_data = self.triangulate(geo, mesh_size)

        y = []
        z = []
        IA = []  # integrand (stress or tangent) * area
        for tr in triangulated_data:
            # All have the same material
            strains = strain[0] - strain[2] * tr[0] + strain[1] * tr[1]
            # compute stresses in all materials
            if integrate == 'stress':
                integrand = tr[3].get_stress(strains)
            elif integrate == 'modulus':
                integrand = tr[3].get_tangent(strains)
            else:
                raise ValueError(f'Unknown integrate type: {integrate}')
            y.append(tr[0])
            z.append(tr[1])
            IA.append(integrand * tr[2])
        prepared_input = [(np.hstack(y), np.hstack(z), np.hstack(IA))]

        return prepared_input, triangulated_data

    def integrate_stress(
        self,
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
        # Integration over all fibers
        y, z, F = prepared_input[0]

        N = np.sum(F)
        My = np.sum(F * z)
        Mz = np.sum(-F * y)

        return N, My, Mz

    def integrate_modulus(
        self,
        prepared_input: t.List[
            t.Tuple[int, np.ndarray, np.ndarray, np.ndarray]
        ],
    ) -> NDArray:
        """Integrate material modulus over the geometry to obtain section
        stiffness.

        Arguments:
            prepared_input (List): The prepared input from .prepare_input().

        Returns:
            NDArray[float]: The stiffness matrix with shape (3,3).
        """
        # Integration over all fibers
        y, z, MA = prepared_input[0]

        stiffness = np.zeros((3, 3))
        stiffness[0, 0] = np.sum(MA)
        stiffness[0, 1] = stiffness[1, 0] = np.sum(MA * z)
        stiffness[0, 2] = stiffness[2, 0] = np.sum(-MA * y)
        stiffness[1, 1] = np.sum(z * z * MA)
        stiffness[1, 2] = stiffness[2, 1] = np.sum(-y * z * MA)
        stiffness[2, 2] = np.sum(y * y * MA)

        return stiffness

    def integrate_strain_response_on_geometry(
        self,
        geo: CompoundGeometry,
        strain: ArrayLike,
        integrate: t.Literal['stress', 'modulus'] = 'stress',
        **kwargs,
    ) -> t.Tuple[t.Union[t.Tuple[float, float, float], np.ndarray], t.List]:
        """Integrate stress or material modulus in the section with the fiber
        algorithm.

        Arguments:
            geo (CompoundGeometry): The geometry of the section.
            strain (ArrayLike): The strains and curvatures of the section,
                given in the format (ea, ky, kz) which are i) strain at 0,0,
                ii) curvature y axis, iii) curvature z axis.
            integrate (str): a string indicating the quantity to integrate over
                the section. It can be 'stress' or 'modulus'. When 'stress'
                is selected, the return value will be the stress resultants N,
                My, Mz, while if 'modulus' is selected, the return will be the
                section stiffness matrix (default is 'stress').
            mesh_size: Percentage of area (number from 0 to 1) max for triangle
                elements.

        Returns:
            Tuple(Union(Tuple(float, float, float), np.ndarray), List): The
            first element is either a tuple of floats (for the stress
            resultants (N, My, Mz) when `integrate='stress'`, or a numpy
            array representing the stiffness matrix then `integrate='modulus'`.
            The second element is a list with data from triangulation.

        Example:
            result, tri = integrate_strain_response_on_geometry(geo, strain,
            integrate='tangent')
            `result` will be the stiffness matrix (a 3x3 numpy array) if


        Raises:
            ValueError: If a unkown value is passed to the `integrate`
            parameter.
        """
        # Prepare the general input based on the geometry and the input strains
        prepared_input, triangulated_data = self.prepare_input(
            geo, strain, integrate, **kwargs
        )

        # Return the calculated response
        if integrate == 'stress':
            return *self.integrate_stress(prepared_input), triangulated_data
        if integrate == 'modulus':
            return self.integrate_modulus(prepared_input), triangulated_data
        raise ValueError(f'Unknown integrate type: {integrate}')
