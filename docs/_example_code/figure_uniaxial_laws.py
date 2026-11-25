"""Plot the light and dark figures for the uniaxial constitutive laws.

From the repository root, run this with StructuralCodes and Matplotlib
available to write the PNGs to docs/theory::

    py -m docs._example_code.figure_uniaxial_laws

Add --show to preview the figures in Matplotlib (useful for testing in the
local development environment)

This could also be added to the CI/CD pipeline to ensure that the figures are
always up to date with the current implementation of the constitutive laws.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure

from structuralcodes.materials.constitutive_laws import (
    BilinearCompression,
    Elastic,
    ElasticFragileTension,
    ElasticPlastic,
    InitialStrain,
    ParabolaRectangle,
    Parallel,
    Popovics,
    Sargin,
    UserDefined,
)

OUTPUT = Path('docs/theory')


def _plot_branch(ax, law, start_strain, end_strain, color, **kwargs):
    """Sample one continuous branch from the implemented law."""
    strains = np.linspace(start_strain, end_strain, 400)
    # Some laws snap their input array in place near ultimate strains.
    stresses = law.get_stress(strains.copy())
    ax.plot(strains, stresses, color=color, lw=2, **kwargs)


def _plot_cutoff(ax, law, cutoff_strain, outer_strain, color, *, label=None):
    """Show the closed value at the cutoff and the open limit beyond it."""
    ax.plot(
        [cutoff_strain, outer_strain], [0, 0], color=color, lw=2, label=label
    )
    ax.plot(
        cutoff_strain, law.get_stress(cutoff_strain), 'o', color=color, ms=5
    )
    ax.plot(
        cutoff_strain,
        0,
        'o',
        color=color,
        markerfacecolor=ax.get_facecolor(),
        ms=5,
    )


def _mark_characteristic_strains(ax, ticks, labelsize=None):
    """Align strain symbols and dashed guides with the numeric scale."""
    values, symbols = zip(*ticks)
    for strain in values:
        if strain != 0:
            ax.axvline(
                strain,
                color=plt.rcParams['text.color'],
                alpha=0.35,
                linewidth=0.75,
                linestyle=(0, (3, 3)),
                zorder=0.5,
            )
    strain_axis = ax.secondary_xaxis('top')
    strain_axis.set_xticks(values, ['$' + symbol + '$' for symbol in symbols])
    if labelsize is not None:
        strain_axis.tick_params(labelsize=labelsize)


def _elastic(ax, colors):
    _plot_branch(ax, Elastic(E=200_000), *ax.get_xlim(), colors[0])


def _elastic_plastic(ax, colors):
    modulus, yield_strength, ultimate_strain = 200_000, 450, 0.01
    yield_strain = yield_strength / modulus
    hardening_modulus = (550 - yield_strength) / (
        ultimate_strain - yield_strain
    )
    hardening_law = ElasticPlastic(
        modulus, yield_strength, hardening_modulus, ultimate_strain
    )
    perfect_plastic_law = ElasticPlastic(
        modulus, yield_strength, eps_su=ultimate_strain
    )
    _plot_branch(
        ax,
        perfect_plastic_law,
        -ultimate_strain,
        ultimate_strain,
        colors[1],
        ls='--',
        label='Perfect plasticity',
    )
    _plot_branch(
        ax,
        hardening_law,
        -ultimate_strain,
        ultimate_strain,
        colors[0],
        label='Linear hardening',
    )
    for cutoff_strain, outer_strain in (
        (-ultimate_strain, -0.012),
        (ultimate_strain, 0.012),
    ):
        _plot_cutoff(
            ax,
            hardening_law,
            cutoff_strain,
            outer_strain,
            colors[0],
        )
        ax.plot(
            cutoff_strain,
            perfect_plastic_law.get_stress(cutoff_strain),
            'o',
            color=colors[1],
            ms=5,
        )
    _mark_characteristic_strains(
        ax,
        [
            (-ultimate_strain, r'-\varepsilon_u'),
            (-yield_strain, r'-\varepsilon_y'),
            (yield_strain, r'\varepsilon_y'),
            (ultimate_strain, r'\varepsilon_u'),
        ],
    )
    ax.legend(loc='upper left', frameon=False)


def _fragile_tension(ax, colors):
    law = ElasticFragileTension(fct=3.0, Ec=30_000)
    ax.plot([-0.00005, 0], [0, 0], color=colors[0], lw=2)
    _plot_branch(ax, law, 0, law.eps_ct, colors[0])
    _plot_cutoff(ax, law, law.eps_ct, 0.00016, colors[0])
    _mark_characteristic_strains(ax, [(law.eps_ct, r'\varepsilon_{ct}')])


def _compression_plot(law, ultimate_symbol, peak_symbol):
    ultimate_strain, _ = law.get_ultimate_strain()
    peak_strain, _ = law.get_ultimate_strain(yielding=True)

    def draw(ax, colors):
        xmin, xmax = ax.get_xlim()
        _plot_branch(ax, law, ultimate_strain, 0, colors[0])
        _plot_cutoff(ax, law, ultimate_strain, xmin, colors[0])
        ax.plot([0, xmax], [0, 0], color=colors[0], lw=2)
        _mark_characteristic_strains(
            ax,
            [(ultimate_strain, ultimate_symbol), (peak_strain, peak_symbol)],
        )

    return draw


def _user_defined(ax, colors):
    strains = np.array([-0.004, -0.002, 0, 0.0015, 0.0035])
    stresses = np.array([-36.0, -30.0, 0, 24.0, 32.0])
    xmin, xmax = ax.get_xlim()
    first_strain, last_strain = strains[0], strains[-1]
    law = UserDefined(strains, stresses, flag=0)
    _plot_branch(ax, law, first_strain, last_strain, colors[0])
    _plot_cutoff(ax, law, first_strain, xmin, colors[0])
    _plot_cutoff(
        ax,
        law,
        last_strain,
        xmax,
        colors[0],
        label='0 - Zero stress (default)',
    )
    for flag, label, color, linestyle in (
        (1, 'Constant stress', colors[2], '--'),
        (2, 'Last tangent', colors[1], '-.'),
        (3, 'Last secant', colors[3], ':'),
    ):
        extended_law = UserDefined(strains, stresses, flag=flag)
        _plot_branch(ax, extended_law, xmin, first_strain, color, ls=linestyle)
        _plot_branch(
            ax,
            extended_law,
            last_strain,
            xmax,
            color,
            ls=linestyle,
            label=f'{flag} - {label}',
        )
    ax.plot(strains, stresses, 'o', color=colors[0], ms=5)
    _mark_characteristic_strains(
        ax,
        [
            (strain, rf'\varepsilon_{index}')
            for index, strain in enumerate(strains, 1)
        ],
    )
    ax.legend(loc='upper left', frameon=False, title='Endpoint flag')


def _initial_strain(ax, colors):
    base = Elastic(E=200_000)
    initial_strain = 0.0015
    compatible = InitialStrain(base, initial_strain)
    released = InitialStrain(base, initial_strain, strain_compatibility=False)
    xlim = ax.get_xlim()
    _plot_branch(ax, base, *xlim, colors[2], ls=':', label='Wrapped law')
    _plot_branch(
        ax, compatible, *xlim, colors[0], label='Compatibility enforced'
    )
    _plot_branch(
        ax,
        released,
        *xlim,
        colors[1],
        ls='--',
        label='Compatibility released',
    )
    _mark_characteristic_strains(
        ax,
        [
            (-initial_strain, r'-\varepsilon_{\mathrm{init}}'),
            (initial_strain, r'\varepsilon_{\mathrm{init}}'),
        ],
    )
    ax.legend(loc='upper left', frameon=False)


def _parallel(ax, colors):
    compression = ParabolaRectangle(fc=30, eps_0=-0.002, eps_u=-0.0035)
    tension = ElasticFragileTension(fct=3.0, Ec=30_000)
    law = Parallel([compression, tension])
    ultimate_strain, _ = compression.get_ultimate_strain()
    peak_strain, _ = compression.get_ultimate_strain(yielding=True)
    _plot_branch(ax, law, ultimate_strain, tension.eps_ct, colors[0])
    _plot_cutoff(ax, law, ultimate_strain, -0.004, colors[0])
    _plot_cutoff(ax, law, tension.eps_ct, 0.00045, colors[0])
    _mark_characteristic_strains(
        ax,
        [
            (ultimate_strain, r'\varepsilon_u'),
            (peak_strain, r'\varepsilon_0'),
        ],
    )
    # Magnify the tensile branch so the cracking point is legible.
    inset = ax.inset_axes([0.165, 0.44, 0.25, 0.32])
    inset.axhline(0, color=plt.rcParams['text.color'], lw=1, zorder=0)
    inset.axvline(0, color=plt.rcParams['text.color'], lw=1, zorder=0)
    _plot_branch(inset, law, 0, tension.eps_ct, colors[0])
    _plot_cutoff(inset, law, tension.eps_ct, 0.00015, colors[0])
    inset.set(xlim=(-0.000025, 0.00015), ylim=(-0.6, 3.8))
    inset.text(
        0.18, 0.94, 'Tension', transform=inset.transAxes, va='top', fontsize=8
    )
    inset.ticklabel_format(axis='x', scilimits=(-3, 4), useMathText=True)
    inset.set_yticks([0, 3])
    inset.tick_params(labelsize=8)
    _mark_characteristic_strains(
        inset, [(tension.eps_ct, r'\varepsilon_{ct}')], labelsize=8
    )


def _parallel_springs(ax, colors):
    """Draw two springs sharing a fixed support and a moving rigid bar."""
    foreground = plt.rcParams['text.color']
    initial_bar_y = 0.50
    displacement_arrow_x = 1.67
    ax.plot([-1.4, 1.4], [1.8, 1.8], color=foreground, lw=2.5)
    ax.plot([-1.4, 1.4], [0, 0], color=foreground, lw=3)
    ax.hlines(
        initial_bar_y,
        -1.4,
        displacement_arrow_x,
        color=foreground,
        linestyle='--',
        linewidth=1,
        alpha=0.65,
    )
    ax.hlines(
        0,
        1.4,
        displacement_arrow_x,
        color=foreground,
        linewidth=0.8,
        alpha=0.65,
    )
    ax.text(
        0,
        initial_bar_y + 0.07,
        'Initial position',
        ha='center',
        va='bottom',
        fontsize=9,
    )
    for x in np.linspace(-1.3, 1.3, 12):
        ax.plot([x - 0.13, x], [1.8, 1.93], color=foreground, lw=0.9)

    for center, turns, color, label in (
        (-0.65, 4, colors[0], r'$K_1$'),
        (0.65, 7, colors[1], r'$K_2$'),
    ):
        offsets = np.tile([-0.17, 0.17], turns)
        x = np.r_[center, center, center + offsets, center, center]
        y = np.r_[0, 0.19, np.linspace(0.27, 1.53, 2 * turns), 1.61, 1.8]
        ax.plot(x, y, color=color, lw=2)
        ax.text(center + 0.27, 0.9, label, color=color, fontsize=13)

    for x, y, length, label in (
        (0, 0, 0.47, r'$F$'),
        (displacement_arrow_x, initial_bar_y, initial_bar_y, r'$u$'),
    ):
        ax.annotate(
            '',
            xy=(x, y - length),
            xytext=(x, y),
            arrowprops={
                'arrowstyle': '->',
                'color': foreground,
                'lw': 1.6,
                'shrinkA': 0,
                'shrinkB': 0,
            },
        )
        ax.text(x + 0.1, y - length / 2, label, va='center', fontsize=13)

    ax.set(xlim=(-1.6, 2), ylim=(-0.6, 2.05), aspect='equal')
    ax.set_axis_off()


def _parallel_elastic_plastic(ax, colors):
    first = ElasticPlastic(E=200_000, fy=300, Eh=5_000)
    second = ElasticPlastic(E=100_000, fy=350, Eh=15_000)
    for law, label, linestyle in (
        (first, 'Material 1', '--'),
        (second, 'Material 2', ':'),
    ):
        _plot_branch(ax, law, 0, 0.008, colors[2], ls=linestyle, label=label)
    for weights, color in (
        ((0.50, 0.50), colors[0]),
        ((0.75, 0.25), colors[1]),
    ):
        law = Parallel([first, second], weights)
        label = f'({weights[0]:.2f}, {weights[1]:.2f})'
        _plot_branch(ax, law, 0, 0.008, color, label=label)
    _mark_characteristic_strains(
        ax,
        [
            (
                first.get_ultimate_strain(yielding=True)[1],
                r'\varepsilon_{y,1}',
            ),
            (
                second.get_ultimate_strain(yielding=True)[1],
                r'\varepsilon_{y,2}',
            ),
        ],
    )
    ax.legend(
        loc='lower right',
        bbox_to_anchor=(1, 0.08),
        frameon=False,
        title='Materials and weights',
    )


def _parallel_ultimate_strain(ax, colors):
    first = ElasticPlastic(E=200_000, fy=300, Eh=5_000, eps_su=0.0055)
    second = ElasticPlastic(E=100_000, fy=350, Eh=15_000, eps_su=0.0105)
    parallel = Parallel([first, second], [1.0, 1.0])
    first_ultimate = first.get_ultimate_strain()[1]
    parallel_ultimate = parallel.get_ultimate_strain()[1]

    color = colors[0]
    _plot_branch(
        ax,
        parallel,
        0,
        parallel_ultimate,
        color,
        label='Parallel',
        marker='o',
        markevery=[-1],
        ms=5,
    )
    _plot_cutoff(ax, parallel, parallel_ultimate, ax.get_xlim()[1], color)
    for law, color, label, linestyle in (
        (first, colors[1], 'Material 1', '--'),
        (second, colors[2], 'Material 2', ':'),
    ):
        _plot_branch(
            ax,
            law,
            0,
            law.get_ultimate_strain()[1],
            color,
            label=label,
            ls=linestyle,
        )
    _mark_characteristic_strains(
        ax,
        [
            (first_ultimate, r'\varepsilon_{u,1}'),
            (
                parallel_ultimate,
                r'\varepsilon_{u,2}=\varepsilon_{u,\mathrm{parallel}}',
            ),
        ],
    )
    ax.legend(loc='lower right', bbox_to_anchor=(1, 0.05), frameon=False)


def _parallel_frp(ax, colors):
    fiber_modulus, matrix_modulus = 230_000, 3_000
    fiber, matrix = Elastic(fiber_modulus), Elastic(matrix_modulus)
    _, xmax = ax.get_xlim()
    for phase_law, modulus_mpa, label, color, linestyle in (
        (fiber, fiber_modulus, 'Fiber', colors[2], '--'),
        (matrix, matrix_modulus, 'Matrix', colors[3], ':'),
    ):
        _plot_branch(
            ax,
            phase_law,
            0,
            xmax,
            color,
            ls=linestyle,
            label=f'{label} ({modulus_mpa / 1000:g} GPa)',
        )
    for fiber_fraction, color in ((0.40, colors[0]), (0.50, colors[1])):
        composite = Parallel(
            [fiber, matrix], [fiber_fraction, 1 - fiber_fraction]
        )
        composite_modulus_gpa = composite.get_tangent(0.0) / 1000
        _plot_branch(
            ax,
            composite,
            0,
            xmax,
            color,
            label=(
                f'FRP, {fiber_fraction:.0%} fiber '
                f'({composite_modulus_gpa:.1f} GPa)'
            ),
        )
    ax.legend(loc='upper left', frameon=False)


def main(save: bool = True) -> None:
    """Save the figures to docs/theory or show them in Matplotlib."""
    # Standalone figures can render PNGs without opening GUI windows.
    make_figure = Figure if save else plt.figure
    compression_plot_xlim = (-0.0042, 0.00055)
    figures = [
        ('parallel_springs', _parallel_springs, None),
        ('elastic', _elastic, (-0.003, 0.003)),
        ('elastic_plastic', _elastic_plastic, (-0.012, 0.012)),
        ('elastic_fragile_tension', _fragile_tension, (-0.00005, 0.00016)),
        (
            'bilinear_compression',
            _compression_plot(
                BilinearCompression(30, -0.0018, -0.0035),
                r'\varepsilon_{cu}',
                r'\varepsilon_c',
            ),
            compression_plot_xlim,
        ),
        (
            'parabola_rectangle',
            _compression_plot(
                ParabolaRectangle(30), r'\varepsilon_{cu}', r'\varepsilon_{c0}'
            ),
            compression_plot_xlim,
        ),
        (
            'sargin',
            _compression_plot(
                Sargin(30), r'\varepsilon_{cu1}', r'\varepsilon_{c1}'
            ),
            compression_plot_xlim,
        ),
        (
            'popovics',
            _compression_plot(
                Popovics(30), r'\varepsilon_{cu}', r'\varepsilon_c'
            ),
            compression_plot_xlim,
        ),
        ('user_defined', _user_defined, (-0.007, 0.0065)),
        ('initial_strain', _initial_strain, (-0.003, 0.0025)),
        ('parallel_elastic_plastic', _parallel_elastic_plastic, (0, 0.008)),
        ('parallel_ultimate_strain', _parallel_ultimate_strain, (0, 0.013)),
        ('parallel_frp', _parallel_frp, (0, 0.003)),
        ('parallel', _parallel, (-0.004, 0.00045)),
    ]
    for theme, colors in (
        ('light', ('#1f77b4', '#ab4c11', '#64686e', '#267449')),
        ('dark', ('#62b8ff', '#ffad65', '#b6bec7', '#82d5a1')),
    ):
        style = 'dark_background' if theme == 'dark' else 'default'
        background = '#141619' if theme == 'dark' else 'white'
        with (
            plt.style.context(style),
            plt.rc_context(
                {'axes.facecolor': background, 'figure.max_open_warning': 0}
            ),
        ):
            foreground = plt.rcParams['text.color']
            for name, draw, xlim in figures:
                figure_height = 3 if xlim is None else 4
                fig = make_figure(
                    figsize=(6, figure_height), layout='constrained'
                )
                ax = fig.subplots()
                if xlim is not None:
                    ax.axhline(0, color=foreground, lw=1, zorder=0)
                    ax.axvline(0, color=foreground, lw=1, zorder=0)
                    ax.set(
                        xlim=xlim,
                        xlabel=r'$\varepsilon$',
                        ylabel=r'$\sigma$ [MPa]',
                    )
                    ax.locator_params(axis='x', nbins=6)
                    ax.ticklabel_format(
                        axis='x', scilimits=(-3, 4), useMathText=True
                    )
                draw(ax, colors)
                if save:
                    fig.savefig(
                        OUTPUT / f'Figure_law_{name}_{theme}.png',
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
