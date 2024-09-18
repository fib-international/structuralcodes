import os
import sys

# Agregar la carpeta A al PYTHONPATH
sys.path.append(
    os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
)


import math

import numpy as np
import section_plots
from results_methods import check_points_in_N_My_Mz
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
if False:
    res = sec.section_calculator.calculate_bending_strength(0, 0)
    section_plots.draw_section_response(
        sec, res.eps_a, res.chi_y, title='Section 1'
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
    theta = math.pi / 10
    # junta rama 0-pi con rama pi-2pi
    res1 = sec.section_calculator.calculate_nm_interaction_domain(theta)
    res2 = sec.section_calculator.calculate_nm_interaction_domain(
        theta + math.pi
    )
    # [::-1] para revertir el segundo array
    res1.n = np.concatenate((res1.n, res2.n[::-1]))
    res1.m_y = np.concatenate((res1.m_y, res2.m_y[::-1]))
    res1.m_z = np.concatenate((res1.m_z, res2.m_z[::-1]))
    section_plots.draw_N_M_diagram(res1)

# My-Mz; N=cte
if False:
    n_ed = -300 * 1e3
    res1 = sec.section_calculator.calculate_mm_interaction_domain(n_ed)
    section_plots.draw_My_Mz_diagram(res1)

# N-My-Mz
if True:
    # Sample forces to test (list of points) [[N,My,Mz]]
    """forces = [
        [-500 * 1e3, -300 * 1e6, 0],
        [0 * 1e3, -200 * 1e6, 0],
        [-100 * 1e3, 0, 50 * 1e6],
        [-100 * 1e3, -500 * 1e6, 50 * 1e6],
    ]"""
    n_ed = np.random.randint(-2000 * 1e3, 500 * 1e3, 10)
    my_ed = np.random.randint(-350 * 1e6, 400 * 1e6, 10)
    mz_ed = np.random.randint(-100 * 1e6, 100 * 1e6, 10)
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
