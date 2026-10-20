"""Example code for creating surface geometries."""

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.animation import FuncAnimation, PillowWriter

from structuralcodes.codes.ec2_2004 import reinforcement_duct_props
from structuralcodes.geometry import (
    RectangularGeometry,
    add_reinforcement_line,
)
from structuralcodes.materials.concrete import ConcreteEC2_2004
from structuralcodes.materials.reinforcement import ReinforcementEC2_2004
from structuralcodes.sections import BeamSection

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
    group_label='web',
)
flange_geom = RectangularGeometry(
    width=width_flange,
    height=height_flange,
    material=concrete,
    origin=(0, height_web + height_flange / 2),
    group_label='flange',
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
        group_label='reinforcement',
    )

# Create a section
section_not_translated = BeamSection(geometry=t_geom)

# Use the centroid of the section, given in the gross properties, to re-allign
# the geometry with the origin
t_geom = t_geom.translate(dy=-section_not_translated.gross_properties.cz)
section = BeamSection(geometry=t_geom)

# Calculate the bending strength
bending_strength = section.section_calculator.calculate_bending_strength()

# Get the stress at a specific point in the section
coord_reinf = (
    section.geometry.point_geometries[0].x,
    section.geometry.point_geometries[0].y,
)
stress = bending_strength.get_point_stress(
    *coord_reinf, group_label='reinforcement'
)
print(f'Stress at point {coord_reinf}: {stress:.2f} MPa')

# Create the detailed_results data structure
bending_strength.create_detailed_result(num_points=10000)

# Visualize bending strength strain field
df = pd.DataFrame(bending_strength.detailed_result.surface_data)
fig_strain, ax_strain = plt.subplots()
norm = mcolors.CenteredNorm(vcenter=0)
abs_max = df['strain'].abs().max()
df.plot.scatter(
    x='y', y='z', s=5, c='strain', ax=ax_strain, cmap='PiYG', norm=norm
)
fig_strain.tight_layout()

fig_strain.savefig('strain_field.png', dpi=300)

# Visualize bending strength stress field
df_p = pd.DataFrame(bending_strength.detailed_result.point_data)
fig_stress, ax_stress = plt.subplots()
# plot as background all concrete fibers in grey
df.plot.scatter(x='y', y='z', s=5, color='lightgray', zorder=1, ax=ax_stress)
# Plot compressed concrete fibers
# Rename just to have a nicer label. Alternative is plotting manually with
# matplotlib and controlling the colorbar
df.rename(columns={'stress': 'stress_c'}).query('stress_c < 0.0').plot.scatter(
    x='y', y='z', s=5, c='stress_c', cmap='Oranges_r', zorder=2, ax=ax_stress
)
# Plot reinforcement fibers
# Rename just to have a nicer label. Alternative is plotting manually with
# matplotlib and controlling the colorbar
df_p.rename(columns={'stress': 'stress_s'}).plot.scatter(
    x='y',
    y='z',
    s=5,
    c='stress_s',
    cmap='Blues',
    zorder=3,
    ax=ax_stress,
    vmin=0,
)

fig_stress.tight_layout()
fig_stress.savefig('stress_field.png', dpi=300)

# Calculate the moment curvature response
moment_curvature = section.section_calculator.calculate_moment_curvature(
    num_pre_yield=10, num_post_yield=50
)

# Plot the moment curvature response
fig_mc, ax_mc = plt.subplots()
ax_mc.plot(-moment_curvature.chi_y, -moment_curvature.m_y * 1e-6, '-k')
ax_mc.set_xlabel(r'$\chi_y$')
ax_mc.set_ylabel(r'$M_y$ [kNm]')
ax_mc.grid(True)
ax_mc.set_xlim(xmin=0)
ax_mc.set_ylim(ymin=0)

fig_mc.tight_layout()
fig_mc.savefig('moment_curvature.png', dpi=300)

# Create the detailed_results data structure
moment_curvature.create_detailed_result(num_points=10000)


# Define a utility function for plotting and saving a figure
def plot_and_save(detailed_result, filename):
    """Simple utility function for plotting and saving stress field."""
    fig, ax = plt.subplots()
    # plot as background all concrete fibers in grey
    df = pd.DataFrame(detailed_result.surface_data)
    df_p = pd.DataFrame(detailed_result.point_data)
    df.plot.scatter(x='y', y='z', s=5, color='lightgray', zorder=1, ax=ax)
    # Plot compressed concrete fibers
    df.rename(columns={'stress': 'stress_c'}).query(
        'stress_c < 0.0'
    ).plot.scatter(
        x='y', y='z', s=5, c='stress_c', cmap='Oranges_r', zorder=2, ax=ax
    )
    # Plot reinforcement fibers
    df_p.rename(columns={'stress': 'stress_s'}).plot.scatter(
        x='y', y='z', s=5, c='stress_s', cmap='Blues', zorder=3, ax=ax, vmin=0
    )

    fig.tight_layout()
    fig.savefig(filename, dpi=300)


# Set the step to the yield moment
moment_curvature.set_step(10)

# Visualize the stress field at yield moment
plot_and_save(
    moment_curvature.detailed_result, 'stress_field_yield_moment.png'
)

# Go the the step of max moment in the detailed results
max_moment_idx = np.argmax(-moment_curvature.m_y)
moment_curvature.set_step(max_moment_idx)

# Visualize the stress field at maximum moment
plot_and_save(moment_curvature.detailed_result, 'stress_field_max_moment.png')

# Plot the stress at rebar for increasing curvature
# Get the stress at a specific point in the section
coord_reinf = (
    section.geometry.point_geometries[0].x,
    section.geometry.point_geometries[0].y,
)
stress = moment_curvature.get_point_stress(
    *coord_reinf, group_label='reinforcement'
)

fig_stress_mc, ax_stress_mc = plt.subplots()
ax_stress_mc.plot(
    -moment_curvature.chi_y,
    stress,
    '-k',
    label=f'Stress at point ({coord_reinf[0]:.0f}, {coord_reinf[1]:.0f}),',
)
ax_stress_mc.set_xlabel(r'$\chi_y$')
ax_stress_mc.set_ylabel('Stress [MPa]')
ax_stress_mc.grid(True)
ax_stress_mc.legend()
ax_stress_mc.set_xlim(xmin=0)
ax_stress_mc.set_ylim(ymin=0)
fig_stress_mc.tight_layout()
fig_stress_mc.savefig('stress_at_rebar.png', dpi=300)

# Create a gif animation of stress field
chi = -moment_curvature.chi_y
mom = -moment_curvature.m_y * 1e-6  # kNm
n_steps = len(chi)

# Precompute axis limits and color limits to avoid jumping during animation
# Section limits from one representative detailed result
moment_curvature.set_step(0)
df0 = pd.DataFrame(moment_curvature.detailed_result.surface_data)
df0_p = pd.DataFrame(moment_curvature.detailed_result.point_data)
ymin = df0['y'].min()
ymax = df0['y'].max()
zmin = df0['z'].min()
zmax = df0['z'].max()
dy = ymax - ymin
dz = zmax - zmin
pad_y = 0.05 * dy if dy > 0 else 1.0
pad_z = 0.05 * dz if dz > 0 else 1.0
# Precompute stress ranges over all steps so the colormap stays fixed
concrete_min = 0.0
steel_max = 0.0

for i in range(n_steps):
    moment_curvature.set_step(i)
    c_min = np.asarray(
        moment_curvature.detailed_result.surface_data['stress']
    ).min()
    concrete_min = min(concrete_min, c_min)
    s_max = np.asarray(
        moment_curvature.detailed_result.point_data['stress']
    ).max()
    steel_max = max(steel_max, s_max)

# Create figure and fixed artists
fig, (ax_mc, ax_sec) = plt.subplots(1, 2, figsize=(14, 6))
# ---- Left: moment-curvature curve
ax_mc.plot(chi, mom, '-k', lw=1.5)
(point_mc,) = ax_mc.plot([], [], 'or', ms=6)

ax_mc.set_xlabel(r'$\chi_y$')
ax_mc.set_ylabel(r'$M_y$ [kNm]')
ax_mc.grid(True)
ax_mc.set_xlim(xmin=0)
ax_mc.set_ylim(ymin=0)

title_mc = ax_mc.set_title('Moment-curvature response')

# ---- Right: section stress field
ax_sec.set_xlabel('y')
ax_sec.set_ylabel('z')
ax_sec.set_aspect('equal')
ax_sec.grid(True)
ax_sec.set_xlim(ymin - pad_y, ymax + pad_y)
ax_sec.set_ylim(zmin - pad_z, zmax + pad_z)

title_sec = ax_sec.set_title('Stress distribution')

# Create empty scatter artists once
ax_sec.scatter(df0['y'], df0['z'], s=4, c='lightgray', zorder=1)
sc_conc = ax_sec.scatter(
    [], [], s=5, c=[], cmap='Oranges_r', zorder=2, vmin=concrete_min, vmax=0.0
)
sc_steel = ax_sec.scatter(
    [], [], s=20, c=[], cmap='Blues', zorder=3, vmin=0.0, vmax=steel_max
)

# Optional colorbars
cbar_conc = fig.colorbar(sc_conc, ax=ax_sec)
cbar_conc.set_label('Concrete stress [MPa]')

cbar_steel = fig.colorbar(sc_steel, ax=ax_sec)
cbar_steel.set_label('Reinforcement stress [MPa]')


# Update function for animation
def update(frame):
    # Update analysis step
    moment_curvature.set_step(frame)
    dr = moment_curvature.detailed_result

    # ---- Left panel: moving point on M-k curve
    point_mc.set_data([chi[frame]], [mom[frame]])
    title_mc.set_text(
        f'Moment-curvature response - step {frame + 1}/{n_steps}'
    )

    # ---- Right panel: section plot
    # Compressed concrete only
    conc_offsets = np.column_stack(
        [np.asarray(dr.surface_data['y']), np.asarray(dr.surface_data['z'])]
    )
    conc_colors = np.asarray(dr.surface_data['stress'])
    conc_offsets = conc_offsets[conc_colors < 0, :]
    conc_colors = conc_colors[conc_colors < 0]

    sc_conc.set_offsets(conc_offsets)
    sc_conc.set_array(conc_colors)

    # Reinforcement fibers
    steel_offsets = np.column_stack(
        [np.asarray(dr.point_data['y']), np.asarray(dr.point_data['z'])]
    )
    steel_colors = np.asarray(dr.point_data['stress'])

    sc_steel.set_offsets(steel_offsets)
    sc_steel.set_array(steel_colors)

    title_sec.set_text(
        f'Stress distribution\n$\\chi_y$={chi[frame]:.4e},  $M_y$={mom[frame]:.2f} kNm'
    )

    return point_mc, sc_conc, sc_steel, title_mc, title_sec


# Build and save animation
fig.tight_layout()

anim = FuncAnimation(
    fig,
    update,
    frames=n_steps,
    interval=150,  # milliseconds between frames
    blit=False,  # safer with scatter + colorbars
    repeat=True,
)

anim.save(
    'moment_curvature_stress_animation.gif',
    writer=PillowWriter(fps=8),
    dpi=150,
)

plt.close(fig)
