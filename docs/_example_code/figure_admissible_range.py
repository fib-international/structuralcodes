"""Plot the light and dark figures of the admissible range of eps_0.

From the repository root, run this with StructuralCodes and Matplotlib
available to write the PNGs to docs/theory::

    py -m docs._example_code.figure_admissible_range

Add --show to preview the figures in Matplotlib (useful for testing in the
local development environment)

This could also be added to the CI/CD pipeline to ensure that the figures are
always up to date with the current implementation of the solver.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure
from matplotlib.patches import Circle, Polygon

from structuralcodes.geometry import RectangularGeometry, add_reinforcement
from structuralcodes.materials.concrete import ConcreteEC2_2004
from structuralcodes.materials.reinforcement import ReinforcementEC2_2004
from structuralcodes.sections import BeamSection

OUTPUT = Path('docs/theory')
# Some parameters for the section
# Section reference system: y horizontal (to the right), z vertical (upwards)
WIDTH, HEIGHT, Z_BARS, Y_BARS, DIAMETER = 250, 500, -192, (-67, 67), 16
# The step of the moment-curvature to illustrate the admissible range of eps_0
STEP = 14


def _data() -> dict:
    """Collect the data for the figure.
    Curvature of a step of the moment-curvature analysis and limits of each
    extreme fibre.
    """
    concrete = ConcreteEC2_2004(fck=45)
    reinforcement = ReinforcementEC2_2004(
        fyk=500, Es=200000, ftk=550, epsuk=0.04
    )
    geo = RectangularGeometry(width=WIDTH, height=HEIGHT, material=concrete)
    for y in Y_BARS:
        geo = add_reinforcement(geo, (y, Z_BARS), DIAMETER, reinforcement)
    section = BeamSection(geo, integrator='marin')
    calculator = section.section_calculator
    mc_res = section.section_calculator.calculate_moment_curvature(n=0)
    curv = mc_res.chi_y[STEP]
    eps_c = concrete.constitutive_law.get_ultimate_strain()
    eps_s = reinforcement.constitutive_law.get_ultimate_strain()
    fibres = [
        ('concrete, top edge', HEIGHT / 2, eps_c),
        ('concrete, bottom edge', -HEIGHT / 2, eps_c),
        ('reinforcement', Z_BARS, eps_s),
    ]
    eps_min, eps_max = calculator._eps0_admissible_bounds(
        section.geometry, curv
    )
    return {
        'polygon': np.array(
            section.geometry.geometries[0].polygon.exterior.coords
        ),
        # coordinates (y, z) of the bars (x and y of the point geometries)
        'bars': [
            (g.x, g.y, g.diameter) for g in section.geometry.point_geometries
        ],
        'curv': curv,
        'fibres': fibres,
        'eps_min': eps_min,
        'eps_max': eps_max,
        'solution': mc_res.eps_a[STEP],
    }


def _draw(fig, data: dict, colors, foreground: str) -> None:
    """Draw the section, the strain profiles and the intervals of eps_0."""
    c_profile, c_steel, c_concrete, c_ok = colors
    curv, fibres = data['curv'], data['fibres']
    eps_min, eps_max = data['eps_min'], data['eps_max']
    color_of = {
        'concrete, top edge': c_concrete,
        'concrete, bottom edge': c_concrete,
        'reinforcement': c_steel,
    }
    # Ratios fine-tuned manually
    ax_sec, ax_eps, ax_rng = fig.subplots(1, 3, width_ratios=[0.6, 1.5, 1.6])
    # Share the vertical axis (z) of the section and of the strain profiles,
    # so that each fibre is at the same height in the two panels.
    ax_sec.sharey(ax_eps)
    z = np.array([-HEIGHT / 2, HEIGHT / 2])

    # (a) Section: polygon of concrete and reinforcement bars
    ax_sec.add_patch(
        Polygon(
            data['polygon'],
            closed=True,
            facecolor=c_concrete,
            alpha=0.3,
            edgecolor=c_concrete,
            lw=1.2,
            clip_on=False,
        )
    )
    for y_bar, z_bar, diameter in data['bars']:
        ax_sec.add_patch(
            Circle((y_bar, z_bar), diameter / 2, color=c_steel, lw=0)
        )
    ax_sec.margins(x=0.15)
    ax_sec.set_aspect('equal', adjustable='datalim')
    ax_sec.set_title('(a) section')
    ax_sec.axis('off')

    # (b) Strain profiles with the same curvature
    x_lim = (-0.008, 0.045)
    for name, z_f, (eps_n, eps_p) in fibres:
        color = color_of[name]
        ax_eps.axhline(z_f, color=color, lw=0.6, ls=':')
        # admissible strains of the fibre (clipped to the plot)
        lo, hi = max(eps_n, x_lim[0]), min(eps_p, x_lim[1])
        ax_eps.plot(
            [lo, hi],
            [z_f, z_f],
            color=color,
            lw=5,
            solid_capstyle='butt',
            alpha=0.35,
        )
        for eps_u, sym in ((eps_n, '-'), (eps_p, '+')):
            if x_lim[0] < eps_u < x_lim[1]:
                ax_eps.plot(eps_u, z_f, '|', color=color, ms=14, mew=2)
                ax_eps.annotate(
                    rf'$\varepsilon_u^{sym}$',
                    (eps_u, z_f),
                    textcoords='offset points',
                    xytext=(4, 6) if z_f < 0 else (-22, -14),
                    fontsize=9,
                    color=color,
                )
    ax_eps.fill_betweenx(
        z,
        eps_min + curv * z,
        eps_max + curv * z,
        color=c_ok,
        alpha=0.2,
        label='admissible profiles',
    )
    for eps_0, ls, label in (
        (eps_min, '-', r'$\varepsilon_0 = \varepsilon_{0,min}$'),
        (eps_max, '--', r'$\varepsilon_0 = \varepsilon_{0,max}$'),
    ):
        ax_eps.plot(
            eps_0 + curv * z, z, ls, color=c_profile, lw=1.8, label=label
        )
    # eps_0 is the strain at z = 0: the admissible range of eps_0 is the
    # segment between the two limit profiles at this height
    ax_eps.axhline(0, color=foreground, lw=0.6, ls=':')
    ax_eps.annotate(
        '',
        xy=(eps_max, 0),
        xytext=(eps_min, 0),
        arrowprops={'arrowstyle': '<->', 'color': c_ok, 'lw': 1.5},
    )
    ax_eps.plot([eps_min, eps_max], [0, 0], 'o', color=c_ok, ms=5, zorder=5)
    for eps_0, label, offset, ha in (
        (eps_min, r'$\varepsilon_{0,min}$', (-4, 6), 'right'),
        (eps_max, r'$\varepsilon_{0,max}$', (4, 6), 'left'),
    ):
        ax_eps.annotate(
            label,
            (eps_0, 0),
            textcoords='offset points',
            xytext=offset,
            ha=ha,
            fontsize=10,
            color=c_ok,
        )
    ax_eps.axvline(0, color=foreground, lw=0.6)
    ax_eps.set(
        xlim=x_lim,
        ylim=(-HEIGHT / 2 - 40, HEIGHT / 2 + 40),
        xlabel=r'$\varepsilon = \varepsilon_0 + \chi \cdot z$',
        ylabel=r'$z$ [mm]',
        title=r'(b) strain profiles with the same curvature $\chi$',
    )
    ax_eps.legend(fontsize=8, loc='upper right')

    # (c) Interval of eps_0 allowed by each fibre and intersection
    rows = [
        (name, eps_n - curv * z_f, eps_p - curv * z_f, color_of[name])
        for name, z_f, (eps_n, eps_p) in fibres
    ]
    lows = [r[1] for r in rows]
    x_lo = min(lows) - 0.004
    x_hi = eps_max + 0.25 * (eps_max - min(lows))
    for k, (name, lo, hi, color) in enumerate(rows):
        row = len(rows) - k
        clipped = hi > x_hi
        ax_rng.plot(
            [lo, min(hi, x_hi)],
            [row, row],
            color=color,
            lw=6,
            solid_capstyle='butt',
            alpha=0.6,
        )
        if clipped:
            ax_rng.annotate(
                '',
                xy=(x_hi, row),
                xytext=(x_hi - 0.03 * (x_hi - x_lo), row),
                arrowprops={'arrowstyle': '->', 'color': color, 'lw': 1.5},
            )
        ax_rng.text(lo, row + 0.22, name, fontsize=8, color=color, va='bottom')
    ax_rng.plot(
        [eps_min, eps_max], [0, 0], color=c_ok, lw=8, solid_capstyle='butt'
    )
    ax_rng.text(
        0.5 * (eps_min + eps_max),
        0.25,
        'admissible range',
        ha='center',
        va='bottom',
        fontsize=8,
        color=c_ok,
    )
    for eps_0 in (eps_min, eps_max):
        ax_rng.axvline(eps_0, color=c_ok, lw=0.8, ls='--')
    ax_rng.plot(
        data['solution'],
        0,
        'D',
        color=foreground,
        mfc='none',
        ms=5,
        label='solution of the current step',
    )
    ax_rng.set(
        xlim=(x_lo, x_hi),
        ylim=(-0.8, len(rows) + 0.8),
        yticks=[],
        xlabel=r'$\varepsilon_0$',
        title=r'(c) $\varepsilon_u^- - \chi \cdot z \leq \varepsilon_0 '
        r'\leq \varepsilon_u^+ - \chi \cdot z$ for each failing material '
        'fibre',
    )
    ax_rng.legend(fontsize=8, loc='lower left')


def main(save: bool = True) -> None:
    """Save the figures to docs/theory or show them in Matplotlib."""
    # Standalone figures can render PNGs without opening GUI windows.
    make_figure = Figure if save else plt.figure
    data = _data()
    for theme, colors in (
        ('light', ('#1f77b4', '#ab4c11', '#64686e', '#267449')),
        ('dark', ('#62b8ff', '#ffad65', '#b6bec7', '#82d5a1')),
    ):
        style = 'dark_background' if theme == 'dark' else 'default'
        background = '#141619' if theme == 'dark' else 'white'
        with (
            plt.style.context(style),
            plt.rc_context({'axes.facecolor': background}),
        ):
            foreground = plt.rcParams['text.color']
            fig = make_figure(figsize=(13, 4.8), layout='constrained')
            _draw(fig, data, colors, foreground)
            if save:
                fig.savefig(
                    OUTPUT / f'FigureAdmissibleRange_{theme}.png',
                    dpi=150,
                    transparent=True,
                )
    if not save:
        plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        '--show',
        action='store_true',
        help='Display the figures in Matplotlib.',
    )
    main(save=not parser.parse_args().show)
