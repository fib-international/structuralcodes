"""Plot the effect of a change in load origin on the N–M interaction domain.

From the repository root, run this with StructuralCodes and Matplotlib
available to write the light and dark PNGs to docs/theory::

    py -m docs._example_code.figure_axial_load_limits

Add --show to preview the figures in Matplotlib.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure
from shapely.geometry import LineString

from structuralcodes import set_design_code
from structuralcodes.codes.mc2010 import reinforcement_duct_props
from structuralcodes.geometry import (
    RectangularGeometry,
    add_reinforcement_line,
)
from structuralcodes.materials.concrete import create_concrete
from structuralcodes.materials.reinforcement import create_reinforcement
from structuralcodes.sections import GenericSection

OUTPUT = Path('docs/theory')


def main(save: bool = True) -> None:
    """Compare the same rectangular RC section in two reference systems."""
    set_design_code('mc2010')
    geometry = RectangularGeometry(300, 500, material=create_concrete(25))
    steel = create_reinforcement(
        fyk=450,
        Es=200_000,
        ftk=450,
        epsuk=reinforcement_duct_props(450, 'B')['epsuk'],
    )
    for z in (-210, 210):
        geometry = add_reinforcement_line(
            geometry, (-110, z), (110, z), 20, steel, n=4
        )

    # The new load origin is 500 mm below the centroid.
    eccentricity = -500
    centered = GenericSection(geometry, integrator='fiber')
    shifted = GenericSection(
        geometry.translate(dy=-eccentricity), integrator='fiber'
    )
    n_limits = (
        np.asarray(centered.section_calculator.calculate_limit_axial_load())
        * 1e-3
    )

    cases = []
    for section, e, domain_label, path_label in (
        (centered, eccentricity * 1e-3, 'Centred origin', r'$M_y=Ne$'),
        (shifted, 0, 'Offset origin', r'$M_y=0$'),
    ):
        domain = section.section_calculator.calculate_nm_interaction_domain(
            num=200, complete_domain=True
        )
        n, m_y = domain.n * 1e-3, domain.m_y * 1e-6
        boundary = LineString(np.column_stack((n, m_y)))
        load_path = LineString([(force, e * force) for force in n_limits])
        intersections = boundary.intersection(load_path)
        capacity = np.array(sorted((p.x, p.y) for p in intersections.geoms))
        cases.append((n, m_y, capacity, e, domain_label, path_label))

    centered_capacity, offset_capacity = (case[2] for case in cases)
    make_figure = Figure if save else plt.figure
    for theme, style, colors, capacity_color in (
        ('light', 'default', ('#1f77b4', '#ab4c11'), '#267449'),
        ('dark', 'dark_background', ('#62b8ff', '#ffad65'), '#82d5a1'),
    ):
        with plt.style.context(style):
            fig = make_figure(figsize=(6.4, 4.5), layout='constrained')
            ax = fig.subplots()
            foreground = plt.rcParams['text.color']
            ax.axvline(0, color=foreground, alpha=0.22, lw=0.8)
            ax.axhline(0, color=foreground, alpha=0.22, lw=0.8)
            for limit in n_limits:
                ax.axvline(limit, color=foreground, alpha=0.35, lw=0.9, ls=':')
            limits_axis = ax.secondary_xaxis('top')
            limits_axis.set_xticks(n_limits, [r'$N_{\min}$', r'$N_{\max}$'])

            for case, color in zip(cases, colors):
                n, m_y, capacity, e, domain_label, path_label = case
                ax.plot(n, m_y, color=color, lw=2, label=domain_label)
                path_n = capacity[:, 0] * 1.2
                ax.plot(
                    path_n,
                    e * path_n,
                    color=color,
                    ls='--',
                    label=path_label,
                )

            # Equal N coordinates identify the same capacity in both frames.
            ax.vlines(
                centered_capacity[:, 0],
                centered_capacity[:, 1],
                offset_capacity[:, 1],
                color=capacity_color,
                linestyles=':',
            )
            for capacity in (centered_capacity, offset_capacity):
                ax.plot(*capacity.T, 'o', color=capacity_color, ms=5, zorder=4)
            for point, label, offset, alignment in (
                (centered_capacity[0], r'$N_{\mathrm{min}}(e)$', 9, 'bottom'),
                (centered_capacity[1], r'$N_{\mathrm{max}}(e)$', -9, 'top'),
            ):
                ax.annotate(
                    label,
                    xy=point,
                    xytext=(0, offset),
                    textcoords='offset points',
                    color=capacity_color,
                    ha='center',
                    va=alignment,
                )

            ax.set_xlabel(r'$N$ [kN]')
            ax.set_ylabel(r'$M_y$ [kN·m]')
            ax.legend(loc='lower right', frameon=False)
            if save:
                fig.savefig(
                    OUTPUT / f'Figure_axial_load_limits_{theme}.png',
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
