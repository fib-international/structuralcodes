import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.append('../')

import math

from matplotlib.lines import Line2D
from matplotlib.patches import Arc
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from results_methods import get_stress_point, get_y_coordinate_inside_poly
from shapely import Polygon

from structuralcodes.materials.constitutive_laws import (
    ParabolaRectangle,
    Sargin,
)


def draw_section(section, title='', reduction_reinf=None):
    """Draw the cross section in matplotlib.

    Args:
        reduction_reinf : float to reduce the diameter in homogenised
        section
    """
    fig, ax = plt.subplots(figsize=(4, 4))
    for i, g in enumerate(section.geometry.geometries):
        poly = g.polygon
        x, y = poly.exterior.xy
        ax.fill(x, y, color='lightcoral')
    for i, g in enumerate(section.geometry.point_geometries):
        center = g._point
        if reduction_reinf is not None:
            r = g._diameter / 2 / reduction_reinf
            color = 'navy'
        else:
            r = g._diameter / 2
            color = 'darkred'
        poly = center.buffer(r)
        x, y = poly.exterior.xy
        ax.fill(x, y, color=color)
    # ax.invert_yaxis()  # Invert Y axis

    # same scale in X and Y
    ax.set_aspect('equal', 'box')
    ax.set_title(title if title else 'Section')

    # plt.grid(True)
    plt.show()


def draw_section_response3D(
    section, eps_a, chi_y, chi_z, lim_Sneg=None, lim_Spos=None, title=None
):
    """Draw the stress and strains of a cross section in 3D.

    Args:
        eps_a : strain at (0,0)
        chi_y : curvature in Y axis
        chi_z : curvature in Z axis
        lim_Sneg, lim_Spos: limits of stress in the plot
    """
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(121, projection='3d')
    ax2 = fig.add_subplot(122, projection='3d')

    if title is not None:
        fig.suptitle(title)

    ax.set_title('Section response - stress')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_xlabel('Stress [MPa]')

    ax2.set_title('Section response - strain')
    ax2.set_ylabel('Y')
    ax2.set_zlabel('Z')
    ax2.set_xlabel('Strain [mm/m]')

    for g in section.geometry.geometries:
        poly = g.polygon
        y_min, z_min, y_max, z_max = poly.bounds

        exterior_coords = np.array(poly.exterior.coords)
        distances = np.sqrt(
            np.sum(np.diff(exterior_coords, axis=0) ** 2, axis=1)
        )
        cumulative_distances = np.insert(np.cumsum(distances), 0, 0)
        total_length = cumulative_distances[-1]
        interpolated_distances = np.linspace(0, total_length, 50)

        border_points = np.empty((50, 2))
        for i, d in enumerate(interpolated_distances):
            index = np.searchsorted(cumulative_distances, d) - 1
            t = (d - cumulative_distances[index]) / distances[index]
            border_points[i] = (1 - t) * exterior_coords[
                index
            ] + t * exterior_coords[index + 1]

        stress = [
            get_stress_point(
                section,
                border_points[i, 0],
                border_points[i, 1],
                eps_a,
                chi_y,
                chi_z,
            )
            for i in range(len(border_points[:, 0]))
        ]

        stress = np.array(stress).reshape(-1, 1)
        stress = np.where(stress == None, 0, stress)
        vertices = np.hstack((stress, border_points))

        # Crear la sombra ed la seccion
        poly_section = Poly3DCollection(
            [
                np.hstack(
                    (np.zeros((border_points.shape[0], 1)), border_points)
                )
            ],
            facecolors='gray',
            edgecolors='b',
            alpha=0.5,
        )
        ax.add_collection3d(poly_section)

        # Crear el polígono 3D para tensiones
        poly_stress = Poly3DCollection(
            [vertices], facecolors='cyan', edgecolors='r', alpha=0.5
        )
        ax.add_collection3d(poly_stress)

        # Actualizar los límites de los ejes basados en los datos
        ax.auto_scale_xyz(vertices[:, 0], vertices[:, 1], vertices[:, 2])

        # Aquí iría un código similar para las deformaciones (strains)

    # Mostrar las figuras fuera del bucle
    plt.show()


def draw_section_response(
    section,
    eps_a,
    chi_y,
    chi_z,
    lim_Sneg=None,
    lim_Spos=None,
    title=None,
):
    """Draw the stress and strains of a cross section.

    Args:
        section: Cross-sectional geometry
        eps_a: Strain at the origin (0,0)
        chi_y: Curvature around the Y-axis
        chi_z: Curvature around the Z-axis
        lim_Sneg: Negative stress limit for plotting
        lim_Spos: Positive stress limit for plotting
        title: Title for the plot
    """

    def angle(x, y):
        """Calculate the angle of the vector (x, y) with respect to (1, 0) in the counterclockwise direction."""
        angle = math.atan2(y, x)
        # Ensure angle is between [0, 2π]
        if angle < 0:
            angle += 2 * math.pi
        return angle  # Angle in radians

    # Orient the section with the neutral axis horizontal
    alfa = angle(chi_y, chi_z)
    if alfa < math.pi / 2:  # I
        rotated_angle = -alfa
    elif alfa < math.pi:  # II
        rotated_angle = math.pi - alfa
    elif alfa < 3 / 2 * math.pi:  # III
        rotated_angle = -alfa + math.pi
    else:  # iV
        rotated_angle = 2 * math.pi - alfa

    sign = np.sign(chi_y) if abs(chi_y) > 0 else 1
    rotated_geom = section.geometry.rotate(rotated_angle)
    chi = (chi_y**2 + chi_z**2) ** 0.5 * sign

    # Print neutral axis position if curvature is non-zero
    if abs(chi) > 0:
        z_neutral_axis = -eps_a / chi
        print(f'z neutral axis = {round(z_neutral_axis, 2)}')
    elif eps_a < 0:
        z_neutral_axis = rotated_geom.geometries[0].polygon.bounds[1]
    else:
        z_neutral_axis = rotated_geom.geometries[0].polygon.bounds[3]

    # region General plot configuration
    # Create the figure and axes with the new criterion
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
    if title is not None:
        fig.suptitle(title)
    # Remove individual borders of axes
    for ax in [ax1, ax2, ax3]:
        for spine in ax.spines.values():
            spine.set_visible(False)  # Remove individual borders
    # Add a common border (spine) around ax1, ax2, and ax3.
    plt.subplots_adjust(
        top=0.85
    )  # Adjust the subplot to make space for a shared title or border
    fig.patch.set_linewidth(2)  # Set thickness of the common border
    fig.patch.set_edgecolor('black')  # Set color of the common border

    # Set equal aspect ratio for ax1
    ax1.set_aspect('equal', 'box')
    ax1.set_title(
        f'Section representation - $\\alpha$={math.degrees(rotated_angle):.1f}°'
    )
    # Set up ax2 for stress plot
    ax2.set_title('Section response - stress')
    ax2.set_xlabel('Stress [MPa]')
    ax2.set_ylabel('z')
    if lim_Sneg is not None and lim_Spos is not None:
        ax2.set_xlim(lim_Sneg, lim_Spos)
    else:
        lim_Sneg, lim_Spos = -1e11, 1e11
    ax2.invert_xaxis()

    # Set up ax3 for strain plot
    ax3.set_title('Section response - strain')
    ax3.set_xlabel('Strain [mm/m]')
    ax3.set_ylabel('z')
    ax3.invert_xaxis()
    # endregion

    # region Plot the section on ax1 with compression and tension zones
    for g in rotated_geom.geometries:
        poly = g.polygon
        minx, miny, maxx, maxy = poly.bounds

        # Create the compressed polygon based on chi
        if chi > 0:
            compressed_polygon = Polygon(
                [
                    (minx - 1, miny - 1),
                    (maxx + 1, miny - 1),
                    (maxx + 1, z_neutral_axis),
                    (minx - 1, z_neutral_axis),
                ]
            )
        else:
            compressed_polygon = Polygon(
                [
                    (minx - 1, z_neutral_axis),
                    (maxx + 1, z_neutral_axis),
                    (maxx + 1, maxy + 1),
                    (minx - 1, maxy + 1),
                ]
            )

        # Get compressed and tension parts of the polygon
        compressed_part = poly.intersection(compressed_polygon)
        tension_part = poly.difference(compressed_part)

        # Plot compressed part in light red
        if not compressed_part.is_empty:
            if compressed_part.geom_type == 'Polygon':
                x, y = compressed_part.exterior.xy
                ax1.fill(x, y, color='lightcoral')
            elif compressed_part.geom_type == 'MultiPolygon':
                for part in compressed_part.geoms:
                    x, y = part.exterior.xy
                    ax1.fill(x, y, color='lightcoral')

        # Plot tension part in light blue
        if not tension_part.is_empty:
            if tension_part.geom_type == 'Polygon':
                x, y = tension_part.exterior.xy
                ax1.fill(x, y, color='lightblue')
            elif tension_part.geom_type == 'MultiPolygon':
                for part in tension_part.geoms:
                    x, y = part.exterior.xy
                    ax1.fill(x, y, color='lightblue')

    # Plot point geometries with colors based on compression or tension
    for g in rotated_geom.point_geometries:
        center = g._point
        r = g._diameter / 2

        # Determine color based on position relative to neutral axis
        if chi > 0:
            if center.y <= z_neutral_axis:
                color = 'darkred'  # Compressed zone
            else:
                color = 'darkblue'  # Tension zone
        elif center.y >= z_neutral_axis:
            color = 'darkred'  # Compressed zone
        else:
            color = 'darkblue'  # Tension zone

        poly = center.buffer(r)
        x, y = poly.exterior.xy
        ax1.fill(x, y, color=color)

    # Annotate the neutral axis on ax1
    ax1.vlines(
        x=minx + 0.1 * (maxx - minx),
        ymin=min(0, z_neutral_axis),
        ymax=max(0, z_neutral_axis),
        color='black',
        linestyle='--',
        linewidth=1,
        label=f'z n.a = {z_neutral_axis:.1f}',
    )
    # Add a text label for the neutral axis
    ax1.text(
        x=minx + 0.1 * (maxx - minx),
        y=z_neutral_axis,
        s=f'z n.a =  {z_neutral_axis:.1f}',
        fontsize='small',
        color='black',
        ha='left',
        va='bottom',
    )

    # plot global axes
    module_axes = 0.25 * min((maxx - minx), (maxy - miny))
    _x = math.cos(rotated_angle) * module_axes
    _y = math.sin(rotated_angle) * module_axes
    ax1.plot([0, _x], [0, _y], color='DarkGray')
    ax1.text(
        x=_x,
        y=_y,
        s='X',
        fontsize='small',
        color='black',
        ha='center',
        va='center',
    )
    ax1.plot([0, -_y], [0, _x], color='DarkGray')
    ax1.text(
        x=-_y,
        y=_x,
        s='Y',
        fontsize='small',
        color='black',
        ha='center',
        va='center',
    )
    ax1.plot([0, module_axes], [0, 0], color='DarkGray')

    # rotated angle annotation
    arc_radius = module_axes * 0.25
    arc = Arc(
        (0, 0),
        arc_radius * 2,
        arc_radius * 2,
        angle=rotated_angle,
        color='black',
    )
    arc = Arc(
        (0, 0),
        width=arc_radius * 2,
        height=arc_radius * 2,
        theta1=0,
        theta2=math.degrees(rotated_angle),
        color='black',
    )

    ax1.add_patch(arc)
    ax1.text(
        arc_radius,
        0,
        r'$\alpha$',
        fontsize='small',
        ha='left',
        va='bottom',
        color='black',
    )

    # endregion

    # Initialize dictionaries for labels in ax2 and ax3
    labels_stress = {}
    labels_strain = {}

    # region Plot stress and strain for each geometry in the section
    for g in rotated_geom.geometries:
        poly = g.polygon
        x_min, y_min, x_max, y_max = poly.bounds
        z_range = np.linspace(y_min, y_max, 500)

        # Calculate stress along z_range
        stress = [
            get_stress_point(
                rotated_geom,
                get_y_coordinate_inside_poly(poly, z),
                z,
                eps_a,
                chi,
                0,
            )
            for z in z_range
        ]

        # Plot stress on ax2
        ax2.plot(stress, z_range, color='gray')

        # Fill positive stress areas
        ax2.fill_betweenx(
            z_range,
            0,
            stress,
            where=(np.array(stress) >= 0),
            facecolor='lightblue',
            interpolate=True,
        )

        # Fill negative stress areas
        ax2.fill_betweenx(
            z_range,
            0,
            stress,
            where=(np.array(stress) < 0),
            facecolor='lightcoral',
            interpolate=True,
        )

        # Add labels for stress
        labels_stress[(stress[0], z_range[0])] = round(stress[0], 2)
        labels_stress[(stress[-1], z_range[-1])] = round(stress[-1], 2)

        # Calculate strain along z_range (converted to mm/m)
        eps = [(eps_a + chi * z) * 1000 for z in z_range]

        # Plot strain on ax3
        ax3.plot(eps, z_range, color='gray')

        # Fill positive strain areas
        ax3.fill_betweenx(
            z_range,
            0,
            eps,
            where=(np.array(eps) >= 0),
            facecolor='lightblue',
            interpolate=True,
        )

        # Fill negative strain areas
        ax3.fill_betweenx(
            z_range,
            0,
            eps,
            where=(np.array(eps) < 0),
            facecolor='lightcoral',
            interpolate=True,
        )

        # Add labels for strain
        labels_strain[(eps[0], z_range[0])] = round(eps[0], 2)
        labels_strain[(eps[-1], z_range[-1])] = round(eps[-1], 2)
    # endregion

    # region Plot stress and strain for point geometries (e.g., reinforcement bars)
    for g in rotated_geom.point_geometries:
        center = g._point
        r = g._diameter / 2
        poly = center.buffer(r)
        x_min, y_min, x_max, y_max = poly.bounds
        z_range = [
            y_min + 1e-5,
            y_max - 1e-5,
        ]  # Small offset to avoid duplicates

        # Calculate stress at two points
        stress = [
            get_stress_point(
                rotated_geom, 0.5 * (x_min + x_max), z, eps_a, chi, 0
            )
            for z in z_range
        ]

        # Plot stress on ax2
        ax2.plot(stress, z_range, color='gray')

        # Fill positive stress areas
        ax2.fill_betweenx(
            z_range,
            0,
            stress,
            where=(np.array(stress) >= 0),
            facecolor='lightblue',
            interpolate=True,
        )

        # Fill negative stress areas
        ax2.fill_betweenx(
            z_range,
            0,
            stress,
            where=(np.array(stress) < 0),
            facecolor='lightcoral',
            interpolate=True,
        )

        # Average stress for labeling
        avg_stress = 0.5 * (stress[0] + stress[1])
        labels_stress[
            (
                min(max(avg_stress, lim_Sneg), lim_Spos),
                0.5 * (z_range[0] + z_range[1]),
            )
        ] = round(avg_stress, 1)

        # Calculate strain at two points
        eps = [(eps_a + chi * z) * 1000 for z in z_range]

        # Plot strain on ax3
        ax3.plot(eps, z_range, color='gray')

        # Fill positive strain areas
        ax3.fill_betweenx(
            z_range,
            0,
            eps,
            where=(np.array(eps) >= 0),
            facecolor='lightblue',
            interpolate=True,
        )

        # Fill negative strain areas
        ax3.fill_betweenx(
            z_range,
            0,
            eps,
            where=(np.array(eps) < 0),
            facecolor='lightcoral',
            interpolate=True,
        )

        # Average strain for labeling
        avg_eps = 0.5 * (eps[0] + eps[1])
        labels_strain[(avg_eps, 0.5 * (z_range[0] + z_range[1]))] = round(
            avg_eps, 2
        )

    # endregion

    # Add stress labels to ax2
    for coords, text in labels_stress.items():
        ax2.text(
            *coords, text, fontsize=10, color='black', ha='left', va='center'
        )

    # Add strain labels to ax3
    for coords, text in labels_strain.items():
        ax3.text(
            *coords, text, fontsize=10, color='black', ha='left', va='center'
        )

    # Display the plot
    plt.show()


def draw_constitutive_law(ec_const, lim_strain_neg=None, lim_strain_pos=None):
    if isinstance(ec_const, ParabolaRectangle):
        strain_p = 0
        strain_n = ec_const._eps_u
    elif isinstance(ec_const, Sargin):
        strain_p = 0
        strain_n = ec_const._eps_cu1
    else:
        strain_p, strain_n = ec_const.get_ultimate_strain()
    strain = np.linspace(strain_n, strain_p, 100)
    stress = ec_const.get_stress(strain)

    fig, ax = plt.subplots(figsize=(4, 4))
    if (lim_strain_neg is not None) and (lim_strain_pos is not None):
        ax.set_xlim(lim_strain_neg, lim_strain_pos)
    else:
        lim_strain_neg, lim_strain_pos = -1e11, 1e11

    ax.plot(strain * 1000, stress, color='green')
    ax.invert_yaxis()
    ax.invert_xaxis()
    ax.set_title(f'stress-strain curve {ec_const.name}')
    ax.set_xlabel('Strain [mm/m]')
    ax.set_ylabel('Stress [MPa]')
    ax.axhline(0, color='gray', linewidth=1)
    ax.axvline(0, color='gray', linewidth=1)

    # ax.set_xticks([strain[0], strain[-1]])
    # ax.set_xticklabels([round(strain[0], 3), round(strain[-1], 3)])
    # ax.set_yticks([stress[0], stress[-1]])
    # ax.set_yticklabels([round(stress[0], 2), round(stress[-1], 2)])

    x_legend = (
        np.linspace(
            max(strain[0], lim_strain_neg), min(lim_strain_pos, strain[-1]), 5
        )
        * 1000
    )
    y_legend = np.linspace(stress[0], stress[-1], 5)
    ax.set_xticks(x_legend)
    ax.set_yticks(y_legend)
    ax.set_xticklabels(np.around(x_legend, decimals=3))
    ax.set_yticklabels(np.around(y_legend, decimals=2))

    plt.show()


def draw_2D_diagram(
    x_boundary,
    y_boundary,
    test_points,
    figsize=(10, 10),
    nticks_x=10,
    nticks_y=10,
    title='2D diagram',
    title_x='X-axis',
    title_y='Y-axis',
    color='red',
    debug=False,
    scale_x=1,
    scale_y=1,
):
    """Draw the 2D diagram of a cross section result representation.

    Args:
        x_boundary : array of x value for the capacity
        y_boundary : array of y value for the capacity
        test_points : Points to check if inside/outside boundary
        figsize : Size of the plot
        nticks_x,nticks_y : number of ticks in the axes
        title :  title of the plot
        title_x,title_y : title of the axis
        color: color of boundary line
        debug: print values of boundary
        scale_x, scale_y: for change units of boundary points and test_points
    """
    # Apply scale (units)
    x_boundary *= scale_x
    y_boundary *= scale_y

    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(x_boundary, y_boundary, color=color)
    ax.set_title(title)
    ax.set_xlabel(title_x)
    ax.set_ylabel(title_y)
    ax.axhline(0, color='gray', linewidth=1)
    ax.axvline(0, color='gray', linewidth=1)
    x_legend = np.linspace(np.min(x_boundary), np.max(x_boundary), nticks_x)
    y_legend = np.linspace(np.min(y_boundary), np.max(y_boundary), nticks_y)
    ax.set_xticks(x_legend)
    ax.set_yticks(y_legend)
    ax.set_xticklabels(np.around(x_legend, decimals=1))
    ax.set_yticklabels(np.around(y_legend, decimals=1))

    for x_i, y_i in zip(x_boundary, y_boundary):
        if debug:
            print(f'X = {x_i:.2f}\tY = {y_i:.2f}')
        # boundary points
        ax.plot(x_i, y_i, 'o', color=color, markersize=2)

    if test_points:
        for res in test_points:
            point = res['point']
            point[0] *= scale_x
            point[1] *= scale_y
            intersection_point = res['intersection_point']
            if intersection_point is not None:
                intersection_point[0] *= scale_x
                intersection_point[1] *= scale_y
            is_inside = res['is_inside']
            # point
            color_check = 'blue' if is_inside else 'red'
            ax.plot(point[0], point[1], 'o', color=color_check, markersize=5)
            if is_inside:
                ax.plot(
                    [0, point[0]],
                    [0, point[1]],
                    linestyle='-',
                    color='gray',
                    linewidth=0.5,
                )
            elif intersection_point is not None:
                ax.plot(
                    intersection_point[0],
                    intersection_point[1],
                    'o',
                    color='orange',
                    markersize=5,
                )
                ax.plot(
                    [0, intersection_point[0]],
                    [0, intersection_point[1]],
                    linestyle='-',
                    linewidth=0.5,
                    color='gray',
                )
                ax.plot(
                    [intersection_point[0], point[0]],
                    [intersection_point[1], point[1]],
                    linestyle='--',
                    linewidth=0.5,
                    color='red',
                )
    # Calculate the maximum efficiency
    if test_points:
        efficiency_points = [
            (res['efficiency'], res['point'])
            for res in test_points
            if res['efficiency'] is not None
        ]
        if efficiency_points:
            max_efficiency, max_eff_point = max(
                efficiency_points, key=lambda x: x[0]
            )
            max_ef_label = f'Max efficiency: {max_efficiency:.2f} at ({max_eff_point[0]:.0f}, {max_eff_point[1]:.0f})'
        else:
            max_ef_label = 'No efficiency'

        # Add the maximum efficiency to the legend
        custom_line = Line2D(
            [0],
            [0],
            color='white',
            marker='',
            linestyle='',
            label=max_ef_label,
        )

        # Adjust legend to avoid duplicates
        handles, labels = ax.get_legend_handles_labels()
        handles.append(custom_line)
        labels.append(max_ef_label)
        ax.legend(
            handles, labels, fontsize='small', loc='upper right', frameon=False
        )
    plt.grid(True, color='lightcyan')
    plt.show()


def draw_N_My_Mz_diagram(
    mesh,
    results=None,
    N_scale=1e-3,
    My_scale=1e-6,
    Mz_scale=1e-6,
    simplify_mesh: bool = False,
):
    """Visualizes the mesh and points using Matplotlib.

    Parameters:
    mesh : trimesh.Trimesh object
        The mesh to visualize (capacity of the section).
    results : list of dict
        The results from check_points_in_N_My_Mz function.
        If None, just plot the mesh capacity.
    N_scale, My_scale, Mz_scale : float, optional
        unit change for the N, My, and Mz axes, respectively.
    """
    # Simplify the mesh to improve performance
    if simplify_mesh:
        target_face_count = int(
            len(mesh.faces) * 0.3
        )  # Reduce to 10% of original faces
        simplified_mesh = mesh.simplify_quadratic_decimation(target_face_count)
    else:
        simplified_mesh = mesh

    # Prepare the figure and axes
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Plot the simplified mesh
    ax.plot_trisurf(
        simplified_mesh.vertices[:, 0] * N_scale,
        simplified_mesh.vertices[:, 1] * My_scale,
        simplified_mesh.vertices[:, 2] * Mz_scale,
        triangles=simplified_mesh.faces,
        color='lightblue',
        alpha=0.5,  # transparency
        edgecolor='gray',
        linewidth=0.1,  # width of mesh edges lines
    )

    # Plot the origin
    origin = np.array([0.0, 0.0, 0.0])
    ax.scatter(
        origin[0],
        origin[1],
        origin[2],
        color='green',
        s=50,
        label='Origin',
        marker='+',
    )

    # Initialize lists for legend handling
    plotted_labels = set()

    # Plot each point and its ray
    if results:
        for res in results:
            point = res['point'].copy()
            point[0] *= N_scale
            point[1] *= My_scale
            point[2] *= Mz_scale
            is_inside = res['is_inside']
            efficiency = res['efficiency']
            intersection_point = res['intersection_point']

            if intersection_point is not None:
                intersection_point = intersection_point.copy()
                intersection_point[0] *= N_scale
                intersection_point[1] *= My_scale
                intersection_point[2] *= Mz_scale

            # Color coding for inside/outside points
            color = 'blue' if is_inside else 'red'
            label = 'Point inside' if is_inside else 'Point outside'

            # Plot the point
            if label not in plotted_labels:
                ax.scatter(
                    point[0],
                    point[1],
                    point[2],
                    color=color,
                    s=20,
                    label=label,
                )
                plotted_labels.add(label)
            else:
                ax.scatter(
                    point[0],
                    point[1],
                    point[2],
                    color=color,
                    s=20,
                )

            # Plot the ray from origin to point
            ax.plot(
                [origin[0], point[0]],
                [origin[1], point[1]],
                [origin[2], point[2]],
                color='purple',
                linestyle='--',
            )

            # If there's an intersection, plot it and the ray to it
            if efficiency is not None and intersection_point is not None:
                if 'Intersection Point' not in plotted_labels:
                    ax.scatter(
                        intersection_point[0],
                        intersection_point[1],
                        intersection_point[2],
                        color='orange',
                        s=20,
                        label='Intersection Point',
                    )
                    plotted_labels.add('Intersection Point')
                else:
                    ax.scatter(
                        intersection_point[0],
                        intersection_point[1],
                        intersection_point[2],
                        color='orange',
                        s=20,
                    )
                # Plot the ray from origin to intersection point
                ax.plot(
                    [origin[0], intersection_point[0]],
                    [origin[1], intersection_point[1]],
                    [origin[2], intersection_point[2]],
                    color='cyan',
                )

    # Set axis labels and title
    ax.set_xlabel('N [kN]')
    ax.set_ylabel('My [mkN]')
    ax.set_zlabel('Mz [mkN]')
    # ax.set_title('Section capacity N-My-Mz')
    ax.set_title('Potato diagram N-My-Mz')

    # Calculate the maximum efficiency
    if results:
        efficiency_points = [
            (res['efficiency'], res['point'])
            for res in results
            if res['efficiency'] is not None
        ]
        if efficiency_points:
            max_efficiency, max_eff_point = max(
                efficiency_points, key=lambda x: x[0]
            )
            max_ef_label = f'Max efficiency: {max_efficiency:.2f} at \nN={max_eff_point[0]*N_scale:.0f}, My={max_eff_point[1]*My_scale:.0f}, Mz={max_eff_point[2]*Mz_scale:.0f}'
        else:
            max_ef_label = 'No efficiency'

        # Add the maximum efficiency to the legend
        custom_line = Line2D(
            [0],
            [0],
            color='white',
            marker='',
            linestyle='',
            label=max_ef_label,
        )

        # Adjust legend to avoid duplicates
        handles, labels = ax.get_legend_handles_labels()
        handles.append(custom_line)
        labels.append(max_ef_label)
        ax.legend(
            handles, labels, fontsize='small', loc='upper right', frameon=False
        )

    plt.show()


def draw_M_curvature_diagram(res, figsize=(10, 10), nticks_x=10, nticks_y=10):
    """Draw the N-M diagram of a cross section.

    Args:
        res : MomentCurvatureResults results
        figsize : Size of the plot
        nticks_x,nticks_y : number of ticks in the axes
    """
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_title(f'M-X   N={res.n/1e3} kN')
    ax.set_xlabel('X (1/km)')
    ax.set_ylabel('My (mkN)')
    ax.axhline(0, color='gray', linewidth=1)
    ax.axvline(0, color='gray', linewidth=1)
    xmin, xmax, ymin, ymax = 0, 0, 0, 0
    ax.plot(res.chi_y * 1e6, res.m_y / 1e6, color='red')
    xmin = min(xmin, np.min(res.chi_y * 1e6))
    ymin = min(ymin, np.min(res.m_y / 1e6))
    xmax = max(xmax, np.max(res.chi_y * 1e6))
    ymax = max(ymax, np.max(res.m_y / 1e6))
    x_legend = np.linspace(xmin, xmax, nticks_x)
    y_legend = np.linspace(ymin, ymax, nticks_y)
    ax.set_xticks(x_legend)
    ax.set_yticks(y_legend)
    ax.set_xticklabels(np.around(x_legend, decimals=1))
    ax.set_yticklabels(np.around(y_legend, decimals=1))
    # plt.legend()
    plt.show()
