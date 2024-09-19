import os
import sys

sys.path.append(
    os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
)


import math

import numpy as np
import section_plots
from results_methods import check_points_in_2D_diagram, check_points_in_N_My_Mz
from shapely import Polygon

from structuralcodes import codes, materials
from structuralcodes.geometry import SurfaceGeometry
from structuralcodes.sections._generic import GenericSection
from structuralcodes.sections._reinforcement import add_reinforcement_line

# from structuralcodes.plots.section_plots import draw_section_response,draw_section,get_stress_point

# ················· CREATE THE SECTION ·································
codes.set_design_code(design_code='ec2_2004')
concrete = materials.concrete.create_concrete(fck=25)
reinforcemnet = materials.reinforcement.create_reinforcement(
    fyk=500,
    Es=200000,
    density=7850,
    ftk=500,
    epsuk=0.07,
)
poly = Polygon(((0, 0), (350, 0), (350, 500), (0, 500)))
geo = SurfaceGeometry(poly, concrete)
geo = add_reinforcement_line(
    geo, (50, 450), (300, 450), 12, reinforcemnet, n=3
)
geo = add_reinforcement_line(geo, (50, 50), (300, 50), 20, reinforcemnet, n=6)
sec = GenericSection(
    geo,
)
# sec.geometry = sec.geometry.translate(-sec.gross_properties.cy, -sec.gross_properties.cz)
sec.geometry = sec.geometry.translate(-175, -250)
# ················· CREATE THE SECTION ·································

# CONSTITUTIVE LAW
if False:
    section_plots.draw_constitutive_law(sec.geometry.geometries[0].material)

# SECTION
if False:
    section_plots.draw_section(sec, 'Section 1')

# RESPONSE TO STRAIN PLANE (UNIAXIAL)
if True:
    res = sec.section_calculator.calculate_bending_strength(0, 0)
    section_plots.draw_section_response(
        sec,
        res.eps_a,
        res.chi_y,
        res.chi_z,
        lim_Sneg=-50,
        lim_Spos=50,
        title='Section 1',
    )

# RESPONSE TO STRAIN PLANE (BIAXIAL)
if True:
    res = sec.section_calculator.calculate_bending_strength(0, 0)
    epsa = -0.1
    chisy = [2, 0, -2]
    # chisz = [-3, -1, 0, 1, 3]
    chisz = [2, 0, -2]
    for chiz in chisz:
        for chiy in chisy:
            _N, _MY, _MZ = sec.section_calculator.integrate_strain_profile(
                (epsa * 1e-3, chiy * 1e-6, chiz * 1e-6)
            )

            section_plots.draw_section_response(
                sec,
                epsa * 1e-3,
                chiy * 1e-6,
                chiz * 1e-6,
                lim_Sneg=-50,
                lim_Spos=50,
                title=rf'Sect1,    $\epsilon_a$ [mm/m]={epsa},    $\chi_y$ [km$^{{-1}}$]={chiy},    $\chi_z$ [km$^{{-1}}$]={chiz}     ->     N={_N*1e-3:.0f}, My={_MY*1e-6:.0f}, Mz={_MZ*1e-6:.0f}',
            )

# MOMENT-CURVATURE
if False:
    # junta rama 0-pi con rama pi-2pi
    res1 = sec.section_calculator.calculate_moment_curvature(theta=0, n=0)
    res2 = sec.section_calculator.calculate_moment_curvature(
        theta=math.pi, n=0
    )
    res1.chi_y = np.concatenate((res1.chi_y, res2.chi_y))
    res1.chi_z = np.concatenate((res1.chi_z, res2.chi_z))
    res1.m_y = np.concatenate((res1.m_y, res2.m_y))
    res1.m_z = np.concatenate((res1.m_z, res2.m_z))
    # Ordena las dos ramas
    indices = np.argsort(res1.chi_y)
    res1.chi_y = res1.chi_y[indices]
    res1.chi_z = res1.chi_z[indices]
    res1.m_y = res1.m_y[indices]
    res1.m_z = res1.m_z[indices]
    section_plots.draw_M_curvature_diagram(res1)

# N-My ; N-Mz
if False:
    theta = math.pi / 4
    # both branches [0-pi) & [pi-2pi)
    bound_1 = sec.section_calculator.calculate_nm_interaction_domain(theta)
    bound_2 = sec.section_calculator.calculate_nm_interaction_domain(
        theta + math.pi
    )
    # [::-1] revert second array
    bound_1.n = np.concatenate((bound_1.n, bound_2.n[::-1]))
    bound_1.m_y = np.concatenate((bound_1.m_y, bound_2.m_y[::-1]))
    bound_1.m_z = np.concatenate((bound_1.m_z, bound_2.m_z[::-1]))

    # ...............  N-My ................
    # check tested points
    forces_n_my = [[-1000 * 1e3, 300 * 1e6], [-500 * 1e3, 100 * 1e6]]
    res = check_points_in_2D_diagram(
        bound_1.n, bound_1.m_y, forces_n_my, debug=True
    )

    # section_plots.draw_N_M_diagram(res1)
    section_plots.draw_2D_diagram(
        bound_1.n,
        bound_1.m_y,
        res,
        title='N-My',
        title_x='N [kN]',
        title_y='My [mkN]',
        color='coral',
        scale_x=1e-3,
        scale_y=1e-6,
    )

    # ...............  N-Mz ................
    # check tested points
    forces_n_mz = [[-5 * 1e3, 10 * 1e6], [-1 * 1e3, 4 * 1e6]]
    res = check_points_in_2D_diagram(
        bound_1.n, bound_1.m_z, forces_n_mz, debug=True
    )

    # section_plots.draw_N_M_diagram(res1)
    section_plots.draw_2D_diagram(
        bound_1.n,
        bound_1.m_z,
        res,
        title='N-Mz',
        title_x='N [kN]',
        title_y='Mz [mkN]',
        color='blue',
        scale_x=1e-3,
        scale_y=1e-6,
    )

# My-Mz; N=cte
if False:
    n_ed = -300 * 1e3
    bound_1 = sec.section_calculator.calculate_mm_interaction_domain(n_ed)
    # check tested points
    # forces_my_mz = [[-60 * 1e6, 47 * 1e6], [-500 * 1e3, 100 * 1e6]]
    my_ed = np.random.randint(-350 * 1e6, 300 * 1e6, 30)
    mz_ed = np.random.randint(-50 * 1e6, 50 * 1e6, 30)
    forces_my_mz = np.column_stack((my_ed, mz_ed))
    res = check_points_in_2D_diagram(
        bound_1.m_y, bound_1.m_z, forces_my_mz, debug=True
    )

    # section_plots.draw_N_M_diagram(res1)
    section_plots.draw_2D_diagram(
        bound_1.m_y,
        bound_1.m_z,
        res,
        title=f'My-Mz, N={round(n_ed*1e-3)}',
        title_x='My [mkN]',
        title_y='Mz [mkN]',
        color='coral',
        scale_x=1e-6,
        scale_y=1e-6,
    )


# N-My-Mz
if False:
    # Sample forces to test (list of points) [[N,My,Mz]]
    """forces = [
        [-500 * 1e3, -300 * 1e6, 0],
        [0 * 1e3, -200 * 1e6, 0],
        [-100 * 1e3, 0, 50 * 1e6],
        [-100 * 1e3, -500 * 1e6, 50 * 1e6],
    ]"""
    n_ed = np.random.randint(-2000 * 1e3, 400 * 1e3, 30)
    my_ed = np.random.randint(-350 * 1e6, 300 * 1e6, 30)
    mz_ed = np.random.randint(-50 * 1e6, 50 * 1e6, 30)
    forces = np.column_stack((n_ed, my_ed, mz_ed))
    # calculate domain
    res = sec.section_calculator.calculate_nmm_interaction_domain(
        num_theta=10, num_axial=10
    )
    # check tested points
    results, mesh = check_points_in_N_My_Mz(
        capacity_pts=res.forces, forces=forces
    )
    # draw
    section_plots.draw_N_My_Mz_diagram(mesh, results)
