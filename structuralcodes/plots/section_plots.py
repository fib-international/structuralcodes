import sys

import matplotlib.pyplot as plt
import numpy as np
from shapely import Point

sys.path.append('../')
from structuralcodes.materials.constitutive_laws import (
    ParabolaRectangle,
    Sargin,
)


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
    ax.invert_yaxis()  # Invert Y axis

    # same scale in X and Y
    ax.set_aspect('equal', 'box')
    ax.set_title(title if title else 'Section')

    # plt.grid(True)
    plt.show()


from mpl_toolkits.mplot3d.art3d import Poly3DCollection


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
    ax.set_xlabel('Stress [MPa]')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    ax2.set_title('Section response - strain')
    ax2.set_xlabel('Strain [mm/m]')
    ax2.set_ylabel('Y')
    ax2.set_zlabel('Z')

    if (lim_Sneg is not None) and (lim_Spos is not None):
        ax.set_zlim(lim_Sneg, lim_Spos)
        ax2.set_zlim(lim_Sneg, lim_Spos)
    else:
        lim_Sneg, lim_Spos = -1e11, 1e11

    labels_stress = {}
    labels_strain = {}

    for g in section.geometry.geometries:
        poly = g.polygon
        y_min, z_min, y_max, z_max = poly.bounds
        y_range = np.linspace(y_min, y_max, 10)
        x, y = poly.exterior.xy
        for y in y_range:
            stress = [
                get_stress_point(section, x[i], y[i], eps_a, chi_y, chi_z)
                for i in range(len(x))
            ]
            eps = [
                (eps_a + chi_y * y[i] + chi_z * z) * 1000
                for i in range(len(x))
            ]

            ax.plot(x, y, stress, color='gray')
            ax2.plot(x, y, eps, color='gray')

            poly_verts = [list(zip(x, y, stress))]
            poly_verts2 = [list(zip(x, y, eps))]

            poly3d = Poly3DCollection(poly_verts, alpha=0.5, edgecolor='r')
            poly3d2 = Poly3DCollection(poly_verts2, alpha=0.5, edgecolor='r')

            ax.add_collection3d(poly3d)
            ax2.add_collection3d(poly3d2)

    for i, g in enumerate(section.geometry.point_geometries):
        center = g._point
        r = g._diameter / 2
        poly = center.buffer(r)
        x_min, y_min, x_max, y_max = poly.bounds
        z_range = np.linspace(y_min, y_max, 3)

        for z in z_range:
            x, y = poly.exterior.xy
            stress = [
                get_stress_point(section, x[i], y[i], eps_a, chi_y, chi_z)
                for i in range(len(x))
            ]
            eps = [
                (eps_a + chi_y * y[i] + chi_z * z) * 1000
                for i in range(len(x))
            ]

            ax.plot(x, y, stress, color='gray')
            ax2.plot(x, y, eps, color='gray')

            poly_verts = [list(zip(x, y, stress))]
            poly_verts2 = [list(zip(x, y, eps))]

            poly3d = Poly3DCollection(poly_verts, alpha=0.5, edgecolor='r')
            poly3d2 = Poly3DCollection(poly_verts2, alpha=0.5, edgecolor='r')

            ax.add_collection3d(poly3d)
            ax2.add_collection3d(poly3d2)

    plt.show()


def draw_section_response(
    section, eps_a, chi_y, lim_Sneg=None, lim_Spos=None, title=None
):
    """Draw the stres and strains of a cross section.

    Args:
        eps_a : strain at (0,0)
        chi_y : curvature in Y axis
        lim_Sneg, lim_Spos: limits of stress in the plot
    """
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10, 5))
    if title is not None:
        fig.suptitle(title)
    ax.set_title('Section response - stress')
    ax.set_xlabel('Stress [MPa]')
    ax.set_ylabel('z')
    ax.invert_yaxis()  # Invert axis
    if (lim_Sneg is not None) and (lim_Spos is not None):
        ax.set_xlim(lim_Sneg, lim_Spos)
    else:
        lim_Sneg, lim_Spos = -1e11, 1e11
    ax.invert_xaxis()

    labels_stress = {}
    labels_strain = {}

    for i, g in enumerate(section.geometry.geometries):
        poly = g.polygon
        x_min, y_min, x_max, y_max = poly.bounds
        z_range = np.linspace(y_min, y_max, 500)
        stress = [
            get_stress_point(
                section, 0.5 * (x_min + x_max), z, eps_a, chi_y, 0
            )
            for z in z_range
        ]
        ax.plot(stress, z_range, color='gray')
        ax.fill_betweenx(
            z_range,
            0,
            stress,
            where=(np.array(stress) >= 0),
            facecolor='lightblue',
            interpolate=True,
        )
        ax.fill_betweenx(
            z_range,
            0,
            stress,
            where=(np.array(stress) < 0),
            facecolor='lightcoral',
            interpolate=True,
        )
        labels_stress[(stress[0], z_range[0])] = round(stress[0], 2)
        labels_stress[(stress[-1], z_range[-1])] = round(stress[-1], 2)

        # strain
        eps = [(eps_a + chi_y * z) * 1000 for z in z_range]
        ax2.plot(eps, z_range, color='gray')
        ax2.fill_betweenx(
            z_range,
            0,
            eps,
            where=(np.array(stress) >= 0),
            facecolor='lightblue',
            interpolate=True,
        )
        ax2.fill_betweenx(
            z_range,
            0,
            eps,
            where=(np.array(stress) < 0),
            facecolor='lightcoral',
            interpolate=True,
        )
        labels_strain[(eps[0], z_range[0])] = round(eps[0], 2)
        labels_strain[(eps[-1], z_range[-1])] = round(eps[-1], 2)

    for i, g in enumerate(section.geometry.point_geometries):
        center = g._point
        r = g._diameter / 2
        poly = center.buffer(r)
        x_min, y_min, x_max, y_max = poly.bounds
        z_range = [y_min, y_max]

        # stress
        stress = [
            get_stress_point(
                section, 0.5 * (x_min + x_max), z, eps_a, chi_y, 0
            )
            for z in z_range
        ]
        ax.plot(stress, z_range, color='gray')
        ax.fill_betweenx(
            z_range,
            0,
            stress,
            where=(np.array(stress) >= 0),
            facecolor='lightblue',
            interpolate=True,
        )
        ax.fill_betweenx(
            z_range,
            0,
            stress,
            where=(np.array(stress) < 0),
            facecolor='lightcoral',
            interpolate=True,
        )
        labels_stress[
            (
                min(max(stress[0], lim_Sneg), lim_Spos),
                0.5 * (z_range[0] + z_range[1]),
            )
        ] = round(0.5 * (stress[0] + stress[1]), 1)

        # strain
        eps = [(eps_a + chi_y * z) * 1000 for z in z_range]
        ax2.plot(eps, z_range, color='gray')
        ax2.fill_betweenx(
            z_range,
            0,
            eps,
            where=(np.array(stress) >= 0),
            facecolor='lightblue',
            interpolate=True,
        )
        ax2.fill_betweenx(
            z_range,
            0,
            eps,
            where=(np.array(stress) < 0),
            facecolor='lightcoral',
            interpolate=True,
        )

        labels_strain[(eps[0], 0.5 * (z_range[0] + z_range[1]))] = round(
            0.5 * (eps[0] + eps[1]), 2
        )

    for coords, text in labels_stress.items():
        ax.text(
            *coords, text, fontsize=10, color='black', ha='left', va='center'
        )
    for coords, text in labels_strain.items():
        ax2.text(
            *coords, text, fontsize=10, color='black', ha='left', va='center'
        )

    # strains
    ax2.set_title('Section response - strain')
    ax2.set_xlabel('strain [mm/m]')
    ax2.set_ylabel('z')
    ax2.invert_yaxis()  # Invert axis
    ax2.invert_xaxis()
    plt.show()
    if abs(chi_y) > 0:
        print(f'z neutral axis = {round(-eps_a/chi_y,2)}')


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


def draw_My_Mz_diagram(res, n=0, figsize=(10, 10), nticks_x=10, nticks_y=10):
    """Draw the My-Mz diagram of a cross section.

    Args:
        res : MMInteractionDomain results
        figsize : Size of the plot
        nticks_x,nticks_y : number of ticks in the axes
    """
    m_y = res.m_y / 1e6
    m_z = res.m_z / 1e6
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(m_y, m_z, color='green')
    ax.set_title(f'My-Mz   N={n/1e3} kN')
    ax.set_xlabel('My (mkN)')
    ax.set_ylabel('Mz (mkN)')
    ax.axhline(0, color='gray', linewidth=1)
    ax.axvline(0, color='gray', linewidth=1)

    x_legend = np.linspace(np.min(m_y), np.max(m_y), nticks_x)
    y_legend = np.linspace(np.min(m_z), np.max(m_z), nticks_y)
    ax.set_xticks(x_legend)
    ax.set_yticks(y_legend)
    ax.set_xticklabels(np.around(x_legend, decimals=1))
    ax.set_yticklabels(np.around(y_legend, decimals=1))
    plt.show()


def draw_N_M_diagram(res, figsize=(10, 10), nticks_x=10, nticks_y=10):
    """Draw the N-M diagram of a cross section.

    Args:
        res : NMInteractionDomain results
        figsize : Size of the plot
        nticks_x,nticks_y : number of ticks in the axes
    """
    n = res.n / 1e3
    m_y = res.m_y / 1e6
    m_z = res.m_z / 1e6
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(n, m_y, color='red')
    ax.set_title('N-My')
    ax.set_xlabel('N (kN)')
    ax.set_ylabel('My (mkN)')
    ax.axhline(0, color='gray', linewidth=1)
    ax.axvline(0, color='gray', linewidth=1)
    x_legend = np.linspace(np.min(n), np.max(n), nticks_x)
    y_legend = np.linspace(np.min(m_y), np.max(m_y), nticks_y)
    ax.set_xticks(x_legend)
    ax.set_yticks(y_legend)
    ax.set_xticklabels(np.around(x_legend, decimals=1))
    ax.set_yticklabels(np.around(y_legend, decimals=1))
    plt.show()

    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(n, m_z, color='red')
    ax.set_title('N-MZ')
    ax.set_xlabel('N (kN)')
    ax.set_ylabel('Mz (mkN)')
    ax.axhline(0, color='gray', linewidth=1)
    ax.axvline(0, color='gray', linewidth=1)
    x_legend = np.linspace(np.min(n), np.max(n), nticks_x)
    y_legend = np.linspace(np.min(m_z), np.max(m_z), nticks_y)
    ax.set_xticks(x_legend)
    ax.set_yticks(y_legend)
    ax.set_xticklabels(np.around(x_legend, decimals=1))
    ax.set_yticklabels(np.around(y_legend, decimals=1))
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
