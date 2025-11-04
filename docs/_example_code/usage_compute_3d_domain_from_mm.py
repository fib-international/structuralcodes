"""Usage example for creating the full NMM domain using many MM slices."""

import matplotlib.pyplot as plt
import numpy as np

from structuralcodes.codes import set_design_code
from structuralcodes.geometry import (
    RectangularGeometry,
    add_reinforcement_line,
)
from structuralcodes.materials.concrete import create_concrete
from structuralcodes.materials.reinforcement import create_reinforcement
from structuralcodes.sections import GenericSection

# Set parameters
height = 625
width = 400
cover = 50
diameter_rebars = 25
n_bars = 3
fcd = 0.85 * 35
fyd_rebars = 500
ftd_rebars = fyd_rebars  # 520 / 1.15
Es = 200000
eps_ud = 3e-2 / 0.9
d1 = cover + diameter_rebars / 2

# Materials
set_design_code('ec2_2004')
concrete = create_concrete(fcd)
reinforcement = create_reinforcement(
    fyk=fyd_rebars, Es=Es, ftk=ftd_rebars, epsuk=eps_ud
)

# Assemble section
concrete_geometry = RectangularGeometry(width, height, concrete)

z_value = -height / 2 + d1

for z_value in (
    -height / 2 + d1,
    height / 2 - d1,
):
    concrete_geometry = add_reinforcement_line(
        geo=concrete_geometry,
        coords_i=(-width / 2 + d1, z_value),
        coords_j=(width / 2 - d1, z_value),
        diameter=diameter_rebars,
        material=reinforcement,
        n=n_bars,
    )

section = GenericSection(
    concrete_geometry, integrator='fiber', mesh_size=0.001
)

# Perform the analyses and create the plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
n = np.linspace(
    section.section_calculator.n_min, section.section_calculator.n_max, 25
)

for n_i in n:
    res = section.section_calculator.calculate_mm_interaction_domain(n=n_i)
    ax.plot(
        np.ones_like(res.m_y) * res.n * 1e-3, res.m_y * 1e-6, res.m_z * 1e-6
    )
ax.set_xlim(
    xmin=section.section_calculator.n_min * 1e-3,
    xmax=section.section_calculator.n_max * 1e-3,
)
ax.set_ylim(ymin=-600, ymax=600)
ax.set_zlim(zmin=-600, zmax=600)
ax.set_xlabel('N [kN]')
ax.set_ylabel('My [kNm]')
ax.set_zlabel('Mz [kNm]')
