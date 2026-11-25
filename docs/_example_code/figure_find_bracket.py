"""Plot the light and dark figures of the bracket search for a fixed curvature.

From the repository root, run this with StructuralCodes and Matplotlib
available to write the PNGs to docs/theory::

    py -m docs._example_code.figure_find_bracket

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

from structuralcodes.geometry import RectangularGeometry, add_reinforcement
from structuralcodes.materials.concrete import ConcreteEC2_2004
from structuralcodes.materials.reinforcement import ReinforcementEC2_2004
from structuralcodes.sections import BeamSection

OUTPUT = Path('docs/theory')
# The axial force for which the bracket search is illustrated
N_EXT = -1000e3
# The step of the moment-curvature to illustrate the bracket finding
STEP = 3


def _section() -> BeamSection:
    """RC section 250 x 500 with two bars of diameter 16 at the bottom."""
    concrete = ConcreteEC2_2004(fck=45)
    reinforcement = ReinforcementEC2_2004(
        fyk=500, Es=200000, ftk=550, epsuk=0.04
    )
    geo = RectangularGeometry(width=250, height=500, material=concrete)
    for y in (-67, 67):
        geo = add_reinforcement(geo, (y, -192), 16, reinforcement)
    return BeamSection(geo)


def _record_step() -> dict:
    """Inputs of one step of the moment-curvature and points evaluated by
    the bracket search.
    """
    section = _section()
    calculator = section.section_calculator
    mc_res = calculator.calculate_moment_curvature(n=N_EXT)
    chi, eps_a = mc_res.chi_y, mc_res.eps_a
    curv = chi[STEP]
    # Starting point and first step of the search, as computed in
    # calculate_moment_curvature (here eps_0_a is within the admissible range
    # and delta is larger than the lower limit 1e-7 used by the library)
    eps_0_a = eps_a[STEP - 1]
    delta = abs(
        (eps_a[STEP - 1] - eps_a[STEP - 2])
        * (chi[STEP] - chi[STEP - 1])
        / (chi[STEP - 1] - chi[STEP - 2])
    )

    def dn(eps_0):
        """Unbalance of axial force in kN."""
        n_int = calculator.integrator.integrate_strain_response_on_geometry(
            section.geometry, [eps_0, curv, 0]
        )[0]
        return (n_int - N_EXT) * 1e-3

    search = []

    def dn_logged(eps_0):
        """Unbalance recording the points of the bracket search."""
        search.append((eps_0, dn(eps_0)))
        return search[-1][1]

    eps_min, eps_max = calculator._eps0_admissible_bounds(
        section.geometry, curv
    )
    calculator._find_bracket(
        dn_logged,
        eps_0_a,
        dn(eps_0_a),
        eps_min,
        eps_max,
        dn(eps_min),
        dn(eps_max),
        delta,
    )
    return {
        'dn': dn,
        'eps_0_a': eps_0_a,
        'eps_min': eps_min,
        'eps_max': eps_max,
        'search': search,
        'solution': eps_a[STEP],
    }


def _step_label(length: float, delta: float) -> str:
    """Label of a step of the search as a multiple of delta."""
    k = length / delta
    return r'$\delta$' if abs(k - 1) < 1e-6 else rf'${k:.0f}\delta$'


def _draw(fig, step: dict, colors, foreground: str) -> None:
    """Draw the two panels of the figure."""
    c_curve, c_right, c_start, c_left = colors
    dn = step['dn']
    eps_0_a = step['eps_0_a']
    search = step['search']

    # Step of each evaluation from the previous point on the same side
    steps, prev = [], {1: eps_0_a, -1: eps_0_a}
    for eps, _ in search:
        side = 1 if eps > eps_0_a else -1
        steps.append((side, prev[side], eps))
        prev[side] = eps
    # The first step of the search is delta
    delta = abs(steps[0][2] - steps[0][1])

    x_pts = [e for e, _ in search] + [eps_0_a]
    x0, x1 = min(x_pts), max(x_pts)
    pad = 0.15 * (x1 - x0)
    x_zoom = np.linspace(x0 - pad, x1 + pad, 400)
    y_zoom = np.array([dn(e) for e in x_zoom])
    x_all = np.linspace(step['eps_min'], step['eps_max'], 600)
    y_all = np.array([dn(e) for e in x_all])

    ax_all, ax = fig.subplots(1, 2, width_ratios=[1, 1.6])

    # Left: whole admissible range
    ax_all.plot(x_all, y_all, color=c_curve, lw=1.6)
    ax_all.axhline(0, color=foreground, lw=0.6)
    ax_all.axvspan(x0 - pad, x1 + pad, color=c_start, alpha=0.3)
    ax_all.plot(eps_0_a, dn(eps_0_a), 's', color=c_start, ms=7)
    ax_all.set(
        xlim=(step['eps_min'], step['eps_max']),
        xlabel=r'$\varepsilon_0$',
        ylabel=r'$\Delta N(\varepsilon_0)$ [kN]',
    )
    ax_all.set_title(
        r'Admissible range $[\varepsilon_{0,min}, \varepsilon_{0,max}]$',
        fontsize=10,
    )
    # Text in axes coordinates, arrow to the centre of the zoomed region
    ax_all.annotate(
        'zoom on the right',
        xy=(0.5 * (x0 + x1), 0.45),
        xycoords=('data', 'axes fraction'),
        xytext=(0.25, 0.45),
        textcoords='axes fraction',
        va='center',
        color=c_start,
        fontsize=9,
        arrowprops={'arrowstyle': '<-', 'color': c_start, 'lw': 1.2},
    )

    # Right: zoom on the bracket search
    ax.plot(x_zoom, y_zoom, color=c_curve, lw=1.6, zorder=1)
    ax.axhline(0, color=foreground, lw=0.6)
    y_min, y_max = y_zoom.min(), y_zoom.max()
    span = y_max - y_min
    e_l, e_r = sorted(steps[-1][1:])
    ax.axvspan(e_l, e_r, color=c_curve, alpha=0.12, zorder=0)
    ax.text(
        0.5 * (e_l + e_r),
        y_max + 0.05 * span,
        'Found bracket',
        ha='center',
        fontsize=9,
        color=foreground,
    )
    ax.plot(
        eps_0_a,
        dn(eps_0_a),
        's',
        color=c_start,
        ms=9,
        mec=foreground,
        zorder=4,
        label=r'$\varepsilon_{0,A}$ (previous step)',
    )
    for k, ((eps, val), (side, eps_prev, _)) in enumerate(
        zip(search, steps), 1
    ):
        color = c_right if side > 0 else c_left
        ax.plot(eps, val, 'o', color=color, ms=7, zorder=5)
        ax.annotate(
            str(k),
            (eps, val),
            textcoords='offset points',
            xytext=(0, 8),
            ha='center',
            fontsize=9,
            color=color,
        )
        # Arrow of the step below the curve, one row per side
        y_arrow = y_min - (0.18 if side > 0 else 0.34) * span
        ax.annotate(
            '',
            xy=(eps, y_arrow),
            xytext=(eps_prev, y_arrow),
            arrowprops={'arrowstyle': '->', 'color': color, 'lw': 1.2},
        )
        ax.text(
            0.5 * (eps + eps_prev),
            y_arrow + 0.03 * span,
            _step_label(abs(eps - eps_prev), delta),
            ha='center',
            va='bottom',
            fontsize=8,
            color=color,
        )
    ax.plot([], [], 'o', color=c_right, label='search to right')
    ax.plot([], [], 'o', color=c_left, label='search to left')
    ax.plot(
        step['solution'],
        dn(step['solution']),
        'D',
        color=foreground,
        mfc='none',
        ms=10,
        zorder=7,
        label='solution',
    )
    ax.set(
        xlim=(x0 - pad, x1 + pad),
        ylim=(y_min - 0.45 * span, y_max + 0.15 * span),
        xlabel=r'$\varepsilon_0$',
    )
    ax.set_title(
        'Bracket search (alternating sides, growing steps)',
        fontsize=10,
    )
    ax.legend(fontsize=8, loc='upper left')


def main(save: bool = True) -> None:
    """Save the figures to docs/theory or show them in Matplotlib."""
    # Standalone figures can render PNGs without opening GUI windows.
    make_figure = Figure if save else plt.figure
    step = _record_step()
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
            fig = make_figure(figsize=(12, 4.6), layout='constrained')
            _draw(fig, step, colors, foreground)
            if save:
                fig.savefig(
                    OUTPUT / f'FigureFindBracket_{theme}.png',
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
