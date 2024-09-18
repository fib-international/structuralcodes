# import sys

import numpy as np
import trimesh
from shapely import Point
from shapely.geometry import LineString, Point, Polygon

# sys.path.append('../')


def get_stress_point(section, y, z, eps_a, chi_y, chi_z):
    """Get the stress in a given point (y,z) inside the cross section.

    Args:
        y,z : coordinates of point inside the c.s.
        eps_a :  strain at (0,0)
        chi_y,chi_z :  curvatures of the cross secction (1/mm)

    Returns:
        stress in y,z
    """
    pt = Point(y, z)
    strain = []
    strain.append(eps_a + z * chi_y + y * chi_z)
    for i, g in enumerate(section.geometry.point_geometries):
        center = g._point
        r = g._diameter / 2
        poly = center.buffer(r)
        if poly.contains(pt) or poly.touches(pt):
            return g.material.get_stress(strain)[0]
    for i, g in enumerate(section.geometry.geometries):
        poly = g.polygon
        if poly.contains(pt) or poly.touches(pt):
            return g.material.get_stress(strain)[0]
    return None


def check_points_in_N_My_Mz(capacity_pts, forces):
    """Checks whether given points are inside a mesh defined by capacity_pts (s_res.NMMInteractionDomain.forces),
    and calculates efficiency for rays from the origin to each point.

    Parameters:
    capacity_pts : array-like, shape (n_points, 3) <- (s_res.NMMInteractionDomain.forces)
        The points defining the capacity mesh.
    forces : array-like, shape (n_forces, 3)
        The points to test.

    Returns:
    results : list of dict
        A list containing dictionaries with the following keys:
        - 'point': The point tested.
        - 'is_inside': True if the point is inside the mesh, False otherwise.
        - 'efficiency': The efficiency for the ray from the origin to the intersection point.
    mesh : trimesh.Trimesh object
        The mesh created from capacity_pts.
    """
    # Adjust capacity_pts according to efficiency
    capacity_pts = np.array(capacity_pts)

    # Create the convex hull and mesh from capacity_pts
    hull = trimesh.convex.convex_hull(capacity_pts)
    mesh = trimesh.Trimesh(vertices=hull.vertices, faces=hull.faces)

    # Initialize the ray tracing module
    try:
        ray_tracer = trimesh.ray.ray_pyembree.RayMeshIntersector(mesh)
    except ImportError:
        # Fall back to the built-in ray tracing method if pyembree is not installed
        ray_tracer = mesh.ray

    # Prepare the results list
    results = []

    # Convert forces to a NumPy array
    forces = np.array(forces)

    # Iterate over each point in forces
    for point in forces:
        # Determine if the point is inside the mesh
        # is_inside = mesh.contains([point])[0]

        # Define the origin
        origin = np.array([0.0, 0.0, 0.0])

        # Calculate the ray direction from the origin to the point
        ray_direction = point - origin
        ray_direction_normalized = ray_direction / np.linalg.norm(
            ray_direction
        )

        # Perform the ray intersection with the mesh
        locations, index_ray, index_tri = ray_tracer.intersects_location(
            ray_origins=[origin],
            ray_directions=[ray_direction_normalized],
            multiple_hits=False,
        )

        if len(locations) > 0:
            # There is an intersection
            intersection_point = locations[0]
            distance_origin_to_intersection = np.linalg.norm(
                intersection_point - origin
            )
            distance_origin_to_point = np.linalg.norm(point - origin)
            efficiency = (
                distance_origin_to_point / distance_origin_to_intersection
            )
        else:
            # No intersection found
            efficiency = None

        # Determine if the point is inside the mesh
        is_inside = efficiency <= 1

        # Append the result for this point
        results.append(
            {
                'point': point,
                'is_inside': is_inside,
                'efficiency': efficiency,
                'intersection_point': intersection_point
                if len(locations) > 0
                else None,
            }
        )

    return results, mesh


def check_points_in_2D_diagram(boundary_x, boundary_y, forces, debug=False):
    """Checks whether given points are inside a boundary defined by capacity_pts,
    and calculates efficiency for rays from the origin to each point.

    Parameters:
    boundary_x,boundary_y : array-like, shape (n_points, 1)
        The points defining the capacity boundary.
    forces : array-like, shape (n_forces, 2)
        The points to test.

    Returns:
    results : list of dict
        A list containing dictionaries with the following keys:
        - 'point': The point tested.
        - 'is_inside': True if the point is inside the boundary, False otherwise.
        - 'efficiency': The efficiency for the ray from the origin to the intersection point.
        - 'intersection_point': The intersection point on the boundary, if any.
    """
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
                intersection_point = np.array([intersection.x, intersection.y])
            elif isinstance(intersection, LineString):
                # The ray lies along an edge; take the point closest to the origin
                coords = np.array(intersection.coords)
                distances = np.linalg.norm(coords - origin, axis=1)
                min_index = np.argmin(distances)
                intersection_point = coords[min_index]
            elif intersection.geom_type == 'MultiPoint':
                # Choose the closest intersection point to the origin
                points = np.array([[pt.x, pt.y] for pt in intersection.geoms])
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
                distance_origin_to_point = np.linalg.norm(point - origin)

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
                f"Point {point} is {'inside' if is_inside else 'outside'} the boundary."
            )
            if efficiency is not None:
                print(f'Efficiency: {efficiency:.4f}')
            if res['intersection_point'] is not None:
                print(f"Intersection Point: {res['intersection_point']}")
            print('---')

    return results
