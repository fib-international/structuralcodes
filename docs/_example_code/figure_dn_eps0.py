"""Plot the light and dark figures of dN(eps_0) at a moment-curvature step.

From the repository root, run this with StructuralCodes and Matplotlib
available to write the PNGs to docs/theory::

    py -m docs._example_code.figure_dn_eps0

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
from scipy.optimize import brentq

from structuralcodes.geometry import RectangularGeometry, add_reinforcement
from structuralcodes.materials.concrete import ConcreteEC2_2004
from structuralcodes.materials.reinforcement import ReinforcementEC2_2004
from structuralcodes.sections import BeamSection

OUTPUT = Path('docs/theory')
# The axial force of the moment-curvature
N_EXT = 0.0
# The step of the moment-curvature to illustrate
STEP = 14


def _section() -> BeamSection:
    """RC section 250 x 500 with two bars of diameter 16 at the bottom."""
    concrete = ConcreteEC2_2004(fck=45)
    reinforcement = ReinforcementEC2_2004(
        fyk=500, Es=200000, ftk=550, epsuk=0.04
    )
    geo = RectangularGeometry(width=250, height=500, material=concrete)
    for y in (-67, 67):
        geo = add_reinforcement(geo, (y, -192), 16, reinforcement)
    return BeamSection(geo, integrator='marin')


def _record_step() -> dict:
    """Inputs of one step of the moment-curvature and points evaluated by
    the bracket search.
    """
    section = _section()
    calc = section.section_calculator
    mc = calc.calculate_moment_curvature(n=N_EXT)
    chi, eps_a = mc.chi_y, mc.eps_a
    # Solution of the previous step and first step of the search, as
    # computed in calculate_moment_curvature
    eps_0_a = eps_a[STEP - 1]
    delta = abs(
        (eps_a[STEP - 1] - eps_a[STEP - 2])
        * (chi[STEP] - chi[STEP - 1])
        / (chi[STEP - 1] - chi[STEP - 2])
    )

    def dn(eps_0):
        """Unbalance of axial force in kN."""
        n_int = calc.integrator.integrate_strain_response_on_geometry(
            section.geometry, [eps_0, chi[STEP], 0]
        )[0]
        return (n_int - N_EXT) * 1e-3

    points = []

    def dn_logged(eps_0):
        """Unbalance recording the evaluated points."""
        points.append((eps_0, dn(eps_0)))
        return points[-1][1]

    eps_min, eps_max = calc._eps0_admissible_bounds(
        section.geometry, chi[STEP]
    )
    # The starting point is clipped to the admissible range
    eps_start = min(max(eps_0_a, eps_min), eps_max)
    eps_l, _, eps_r, _ = calc._find_bracket(
        dn_logged,
        eps_start,
        dn(eps_start),
        eps_min,
        eps_max,
        dn(eps_min),
        dn(eps_max),
        delta,
    )
    # Brent's method with the same tolerances used by the library
    solution = brentq(
        dn,
        *sorted((eps_l, eps_r)),
        xtol=1e-14,
        rtol=4 * np.finfo(float).eps,
    )
    return {
        'dn': dn,
        'eps_0_a': eps_0_a,
        'eps_start': eps_start,
        'bounds': (eps_min, eps_max),
        'bracket': (eps_l, eps_r),
        'search': points,
        'solution': solution,
    }


def _plot_dn(ax, x, dn, bounds, colors) -> None:
    """Plot dN, thicker inside the admissible range."""
    c_curve, _, c_start = colors
    y = np.array([dn(e) for e in x])
    inside = (x >= bounds[0]) & (x <= bounds[1])
    ax.plot(x, y, color=c_start, lw=1)
    ax.plot(x[inside], y[inside], color=c_curve, lw=2)


def _draw(fig, step: dict, colors, foreground: str) -> None:
    """Draw the two panels of the figure."""
    c_curve, c_search, c_start = colors
    dn, eps_0_a, eps_start = step['dn'], step['eps_0_a'], step['eps_start']
    eps_min, eps_max = step['bounds']

    # Style of the markers shared by the two panels: each key is a kind of
    # point and its value holds the keyword arguments of Axes.plot (marker
    # symbol, size, face and edge colors, and the label of the legend, drawn
    # only in the right panel)
    markers = {
        'bound': {'marker': 'o', 'ms': 6, 'mfc': c_curve, 'mec': foreground},
        'start': {
            'marker': 'o',
            'ms': 8,
            'mfc': c_curve,
            'mec': foreground,
            'label': 'start of the search (clipped)',
        },
        'previous': {
            'marker': 's',
            'ms': 7,
            'mfc': c_start,
            'mec': foreground,
            'label': r'$\varepsilon_{0,A}$ (previous step)',
        },
        'search': {
            'marker': 'o',
            'ms': 7,
            'mfc': c_search,
            'mec': c_search,
            'label': 'bracket search to right',
        },
        'solution': {
            'marker': 'D',
            'ms': 5,
            'mfc': 'none',
            'mec': foreground,
            'label': 'solution',
        },
    }

    def mark(axes, eps, kind):
        """Plot the points of dN at eps with the style of kind."""
        axes.plot(
            eps, np.vectorize(dn)(eps), ls='none', zorder=5, **markers[kind]
        )

    ax_all, ax = fig.subplots(1, 2, width_ratios=[1, 1.3])
    for a in (ax_all, ax):
        a.axhline(0, color=foreground, lw=0.6)
        a.set(
            xlabel=r'$\varepsilon_0$', ylabel=r'$\Delta N(\varepsilon_0)$ [kN]'
        )

    # Left: dN in the whole admissible range (and a bit outside)
    margin = 0.1 * (eps_max - eps_min)
    x_all = np.linspace(eps_min - margin, eps_max + margin, 600)
    _plot_dn(ax_all, x_all, dn, step['bounds'], colors)
    ax_all.axvspan(eps_min, eps_max, color=c_curve, alpha=0.08)
    for eps_b, label in zip(
        step['bounds'], (r'$\varepsilon_{0,min}$', r'$\varepsilon_{0,max}$')
    ):
        ax_all.axvline(eps_b, color=c_curve, ls='--', lw=1)
        mark(ax_all, eps_b, 'bound')
        ax_all.annotate(
            label,
            (eps_b, 0),
            textcoords='offset points',
            xytext=(4, 8),
            color=foreground,
        )
    mark(ax_all, eps_0_a, 'previous')
    ax_all.set_title(r'Admissible range of $\varepsilon_0$', fontsize=10)

    # Right: zoom on the start, the bracket search and the solution
    eps_s = [eps for eps, _ in step['search']]
    x0, x1 = min(eps_0_a, *eps_s), max(eps_0_a, *eps_s)
    pad = 0.35 * (x1 - x0)
    _plot_dn(
        ax, np.linspace(x0 - pad, x1 + pad, 300), dn, step['bounds'], colors
    )
    ax.axvspan(*step['bracket'], color=c_curve, alpha=0.12, zorder=0)
    ax.text(
        np.mean(step['bracket']),
        0.95,
        'Found bracket',
        transform=ax.get_xaxis_transform(),
        ha='center',
        va='top',
        fontsize=9,
        color=foreground,
    )
    ax.axvline(
        eps_min, color=c_curve, ls='--', lw=1, label=r'$\varepsilon_{0,min}$'
    )
    mark(ax, eps_start, 'start')
    mark(ax, eps_0_a, 'previous')
    mark(ax, eps_s, 'search')
    mark(ax, step['solution'], 'solution')
    for k, (eps, val) in enumerate(step['search'], 1):
        ax.annotate(
            str(k),
            (eps, val),
            textcoords='offset points',
            xytext=(4, 4),
            color=c_search,
        )
    ax.set_xlim(x0 - pad, x1 + pad)
    ax.set_title('Bracket search and solution', fontsize=10)
    ax.legend(fontsize=8)


def main(save: bool = True) -> None:
    """Save the figures to docs/theory or show them in Matplotlib."""
    # Standalone figures can render PNGs without opening GUI windows.
    make_figure = Figure if save else plt.figure
    step = _record_step()
    # Colors of the curve, of the bracket search and of the points outside
    # the admissible range
    for theme, colors in (
        ('light', ('#1f77b4', '#ab4c11', '#64686e')),
        ('dark', ('#62b8ff', '#ffad65', '#b6bec7')),
    ):
        style = 'dark_background' if theme == 'dark' else 'default'
        background = '#141619' if theme == 'dark' else 'white'
        with (
            plt.style.context(style),
            plt.rc_context({'axes.facecolor': background}),
        ):
            foreground = plt.rcParams['text.color']
            fig = make_figure(figsize=(11, 4.2), layout='constrained')
            _draw(fig, step, colors, foreground)
            if save:
                fig.savefig(
                    OUTPUT / f'FiguredN_eps_0_{theme}.png',
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
