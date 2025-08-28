"""Example code for creating surface geometries."""

import matplotlib.pyplot as plt

from structuralcodes.codes.ec2_2004 import reinforcement_duct_props
from structuralcodes.geometry import (
    RectangularGeometry,
    add_reinforcement_line,
)
from structuralcodes.materials.concrete import ConcreteEC2_2004
from structuralcodes.materials.reinforcement import ReinforcementEC2_2004
from structuralcodes.sections import GenericSection

# Define parameters
fck = 30
fyk = 500
Es = 200000
ductility_class = 'C'
width_web = 250
width_flange = 1000
height_web = 650
height_flange = 150
diameter_bar = 20
cover = 50
n_bars_layer = 3

# Create material
concrete = ConcreteEC2_2004(fck=fck)
reinforcement = ReinforcementEC2_2004(
    fyk=fyk,
    Es=Es,
    **reinforcement_duct_props(fyk=fyk, ductility_class=ductility_class),
)

# Create surface geometries
web_geom = RectangularGeometry(
    width=width_web,
    height=height_web,
    material=concrete,
    origin=(0, height_web / 2),
)
flange_geom = RectangularGeometry(
    width=width_flange,
    height=height_flange,
    material=concrete,
    origin=(0, height_web + height_flange / 2),
)

# Add surface geometries to create a T-shaped geometry
t_geom = web_geom + flange_geom

# Add two layers of reinforcement
for y_coord in (
    cover + diameter_bar / 2,
    2 * cover + 3 * diameter_bar / 2,
):
    t_geom = add_reinforcement_line(
        geo=t_geom,
        coords_i=(
            -width_web / 2 + cover + diameter_bar / 2,
            y_coord,
        ),
        coords_j=(
            width_web / 2 - cover - diameter_bar / 2,
            y_coord,
        ),
        diameter=diameter_bar,
        material=reinforcement,
        n=n_bars_layer,
    )

# Create a section
section_not_translated = GenericSection(geometry=t_geom)

# Use the centroid of the section, given in the gross properties, to re-allign
# the geometry with the origin
t_geom = t_geom.translate(dy=-section_not_translated.gross_properties.cz)
section = GenericSection(geometry=t_geom)

# Calculate the bending strength
bending_strength = section.section_calculator.calculate_bending_strength()

# Calculate the limit axial loads
limit_axial_loads = section.section_calculator.calculate_limit_axial_load()

# Calculate the interaction domain between axial force and bending moment
nm = section.section_calculator.calculate_nm_interaction_domain()

# Calculate the moment-curvature relation
moment_curvature = section.section_calculator.calculate_moment_curvature()

# Visualize nm interaction domain
fig_nm, ax_nm = plt.subplots()
ax_nm.plot(nm.m_y * 1e-6, nm.n * 1e-3, '-k')
ax_nm.grid()
ax_nm.set_xlabel(r'$M_{\mathrm{y}}$ [kNm]')
ax_nm.set_ylabel(r'$N$ [kN]')
ax_nm.set_xlim(xmax=0)
ax_nm.set_ylim(ymax=0)
fig_nm.tight_layout()

# Visualize moment-curvature relation
fig_momcurv, ax_momcurv = plt.subplots()
ax_momcurv.plot(
    -moment_curvature.chi_y * 1e3, -moment_curvature.m_y * 1e-6, '-k'
)
ax_momcurv.grid()
ax_momcurv.set_xlabel(r'$\chi_{\mathrm{y}}$ [1/m]')
ax_momcurv.set_ylabel(r'$M_{\mathrm{y}}$ [kNm]')
ax_momcurv.set_xlim(xmin=0)
ax_momcurv.set_ylim(ymin=0)
fig_momcurv.tight_layout()

fig_momcurv.savefig('moment_curvature.png', dpi=300)
