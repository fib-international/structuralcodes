"""The fiber section integrator."""

from __future__ import annotations  # To have clean hints of ArrayLike in docs

import typing as t

import numpy as np
import triangle
from numpy.typing import ArrayLike
from shapely import Polygon

from structuralcodes.geometry import CompoundGeometry, SurfaceGeometry

from ._section_integrator import SectionIntegrator


class FiberIntegrator(SectionIntegrator):
    """Section integrator based on the Marin algorithm."""

    def prepare_triangulation(self, geo: SurfaceGeometry) -> t.Dict:
        """Triangulate a SurfaceGeometry object.

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

    def prepare_input(
        self, geo: CompoundGeometry, strain: ArrayLike, **kwargs
    ) -> t.Tuple[t.Tuple[np.ndarray, np.ndarray, np.ndarray]]:
        """Prepare general input to the integration.

        Calculate the stresses based on strains in a set of points.

        Keyword Arguments:
            geo (CompoundGeometry): The geometry of the section.
            strain (ArrayLike): The strains and curvatures of the section,
                given in the format (ea, ky, kz) which are i) strain at 0,0,
                ii) curvature y axis, iii) curvature z axis.
            mesh_size: Percentage of area (number from 0 to 1) max for triangle
                elements.

        Returns:
            Tuple(List, Dict): The prepared input representing a list with
            x-coordinates, y-coordinates and force for each fiber and a
            dictionary containing the triangulation data that can be stored and
            used later to avoid repetition of triangulation.
        """
        # This method should:
        #  - discretize the section in a number of fibers (mesh_size)
        #  - prepare input as y coordiates z coordinates forces

        prepared_input = []

        triangulated_data = kwargs.get('tri', None)
        if triangulated_data is None:
            # No triangulation is provided, triangulate the section
            # Fiber integrator for generic section uses delaunay triangulation
            # for discretizing in fibers
            triangulated_data = []
            mesh_size = kwargs.get('mesh_size', 0.01)
            if mesh_size <= 0 or mesh_size > 1:
                raise ValueError('mesh_size is a number from 0 to 1')
            # For the surface geometries
            for g in geo.geometries:
                # prepare data structure for triangle module
                tri = self.prepare_triangulation(g)
                # define the maximum area of the triangles
                max_area = g.area * mesh_size
                # triangulate the geometry getting back the mesh
                mesh = triangle.triangulate(
                    tri, f'pq{30:.1f}Aa{max_area:.1f}o1'
                )
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
                        - mesh['vertices'][tr[0]][1]
                        * mesh['vertices'][tr[1]][0]
                    )
                    a += (
                        mesh['vertices'][tr[1]][0] * mesh['vertices'][tr[2]][1]
                        - mesh['vertices'][tr[1]][1]
                        * mesh['vertices'][tr[2]][0]
                    )
                    a += (
                        mesh['vertices'][tr[2]][0] * mesh['vertices'][tr[0]][1]
                        - mesh['vertices'][tr[2]][1]
                        * mesh['vertices'][tr[0]][0]
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

        x = []
        y = []
        F = []
        for tr in triangulated_data:
            # All have the same material
            strains = strain[0] - strain[2] * tr[0] + strain[1] * tr[1]
            # compute stresses in all materials
            stresses = tr[3].get_stress(strains)
            x.append(tr[0])
            y.append(tr[1])
            F.append(stresses * tr[2])
        prepared_input = [(np.hstack(x), np.hstack(y), np.hstack(F))]

        return prepared_input, triangulated_data

    def integrate(
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
        x, y, F = prepared_input[0]

        N = np.sum(F)
        Mx = np.sum(F * y)
        My = np.sum(-F * x)

        return N, Mx, My

    def integrate_strain_response_on_geometry(
        self, geo: CompoundGeometry, strain: ArrayLike, **kwargs
    ):
        """Integrate the strain response with the fiber algorithm.

        Arguments:
            geo (CompoundGeometry): The geometry of the section.
            strain (ArrayLike): The strains and curvatures of the section,
                given in the format (ea, ky, kz) which are i) strain at 0,0,
                ii) curvature y axis, iii) curvature z axis.
            mesh_size: Percentage of area (number from 0 to 1) max for triangle
                elements.

        Returns:
            Tuple(Tuple(float, float, float), Dict): The stress resultants N,
            Mx and My and the triangulation data.
        """
        # Prepare the general input based on the geometry and the input strains
        prepared_input, triangulated_data = self.prepare_input(
            geo, strain, **kwargs
        )

        # Return the calculated response
        return *self.integrate(prepared_input), triangulated_data
