import marimo

__generated_with = "0.25.1"
app = marimo.App(
    width="medium",
    app_title="Compression + traction par Parallel",
)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Compression + tension, and compression_tension with `Parallel`

    ```python
    compression = ParabolaRectangle(fc=45)      # ou Sargin, Popovics...
    tension = ElasticFragileTension(fct=2.7, Ec=36000)
    compression_tension = Parallel(constitutive_laws=[compression, tension])
    ```
    """)
    return


@app.cell
def _():
    import inspect

    import marimo as mo
    import numpy as np
    from matplotlib import pyplot as plt
    from matplotlib.axes import Axes
    from structuralcodes.core.base import ConstitutiveLaw
    from structuralcodes.materials.concrete import ConcreteEC2_2004
    from structuralcodes.materials.constitutive_laws import (
        BilinearCompression,
        ParabolaRectangle,
        Parallel,
        Popovics,
        Sargin,
        ElasticFragileTension,
    )

    from shapely import Polygon

    from structuralcodes import set_design_code
    from structuralcodes.geometry import SurfaceGeometry, add_reinforcement
    from structuralcodes.materials.concrete import create_concrete
    from structuralcodes.materials.reinforcement import create_reinforcement
    from structuralcodes.sections import BeamSection

    return (
        Axes,
        BeamSection,
        ConstitutiveLaw,
        ElasticFragileTension,
        ParabolaRectangle,
        Parallel,
        Polygon,
        SurfaceGeometry,
        add_reinforcement,
        create_concrete,
        create_reinforcement,
        mo,
        np,
        plt,
        set_design_code,
    )


@app.cell(hide_code=True)
def _(
    Axes,
    compression: "ConstitutiveLaw",
    compression_tension: "ConstitutiveLaw",
    fct,
    np,
    plt,
    tension: "ConstitutiveLaw",
):
    C_COMPRESSION = '#2a78d6'  # bleu
    C_TRACTION = '#eb6834'  # orange
    C_ASSEMBLAGE = '#0b0b0b'  # noir, en pointilles
    C_LIMITE = '#e34948'
    C_SECONDAIRE = '#52514e'
    C_GRILLE = '#d8d7d2'
    C_TEXTE = '#0b0b0b'

    _eps_u_compression = abs(compression.get_ultimate_strain()[0])
    _fig, (_ax, _ax_zoom) = plt.subplots(1, 2, figsize=(11.5, 4.2))


    def _tracer(ax: Axes, eps: np.ndarray) -> None:
        '''Trace les trois courbes sur un jeu de deformations donne.'''
        ax.plot(
            eps * 1000,
            np.asarray(compression.get_stress(eps.copy())),
            color=C_COMPRESSION,
            lw=2,
            label='compression',
        )
        ax.plot(
            eps * 1000,
            np.asarray(tension.get_stress(eps.copy())),
            color=C_TRACTION,
            lw=2,
            label='traction',
        )
        ax.plot(
            eps * 1000,
            np.asarray(compression_tension.get_stress(eps.copy())),
            color=C_ASSEMBLAGE,
            lw=1.6,
            ls=':',
            label='parallel',
        )
        ax.axhline(0, color=C_SECONDAIRE, lw=0.8)
        ax.axvline(0, color=C_SECONDAIRE, lw=0.8)
        ax.set_xlabel('$\\varepsilon$ (‰)', fontsize=9, color=C_SECONDAIRE)
        ax.set_ylabel('$\\sigma$ (MPa)', fontsize=9, color=C_SECONDAIRE)
        ax.grid(True, color=C_GRILLE, lw=0.6)
        ax.set_axisbelow(True)
        for cote in ('top', 'right'):
            ax.spines[cote].set_visible(False)


    _tracer(
        _ax, np.linspace(-1.05 * _eps_u_compression, 4.0 * tension.eps_ct_u, 2000)
    )
    _ax.set_title(
        f'{compression.name} + {tension.name}',
        fontsize=11,
        color=C_TEXTE,
        loc='left',
    )
    _ax.legend(loc='lower right', frameon=False, fontsize=9)

    _tracer(
        _ax_zoom,
        np.linspace(-6.0 * tension.eps_ct_u, 2.5 * tension.eps_ct_u, 2000),
    )
    _ax_zoom.axhline(fct, color=C_LIMITE, lw=1, ls='--')
    _ax_zoom.annotate(
        '$f_{ct}$',
        (_ax_zoom.get_xlim()[0], fct),
        textcoords='offset points',
        xytext=(4, 3),
        fontsize=8,
        color=C_LIMITE,
    )
    _ax_zoom.set_ylim(-3.0 * fct, 1.6 * fct)
    _ax_zoom.set_title('Zoom on traction', fontsize=11, color=C_TEXTE, loc='left')
    _fig.tight_layout()
    _fig
    return


@app.cell(hide_code=True)
def _(
    compression: "ConstitutiveLaw",
    compression_tension: "ConstitutiveLaw",
    mo,
    tension: "ConstitutiveLaw",
):
    mdtext = f"""
    | | compression | tension | compression_tension |
    |---|---|---|---|
    | classe | `{type(compression).__name__}` | `{type(tension).__name__}` | `Parallel` |
    | $\\varepsilon_u$ negative | {compression.get_ultimate_strain()[0] * 1000:.3f} ‰ | {tension.get_ultimate_strain()[0] * 1000:.3f} ‰ | {compression_tension.get_ultimate_strain()[0] * 1000:.3f} ‰ |
    | $\\varepsilon_u$ positive | {compression.get_ultimate_strain()[1] * 1000:.3f} ‰ | {tension.get_ultimate_strain()[1] * 1000:.3f} ‰ | {compression_tension.get_ultimate_strain()[1] * 1000:.3f} ‰ |
    """
    mo.md(mdtext)
    return


@app.cell
def _():
    return


@app.cell
def _(
    ConstitutiveLaw,
    ElasticFragileTension,
    ParabolaRectangle,
    Parallel,
    Polygon,
    SurfaceGeometry,
    add_reinforcement,
    create_concrete,
    create_reinforcement,
    set_design_code,
):
    from operator import ge

    # Concrete
    Ec = 36000
    fct = 2.7
    fck = 45
    set_design_code('ec2_2004')
    compression: ConstitutiveLaw = ParabolaRectangle(fc=fck)
    tension: ConstitutiveLaw = ElasticFragileTension(fct=fct, Ec=Ec)
    compression_tension: ConstitutiveLaw = Parallel(
        constitutive_laws=[compression, tension]
    )

    concrete_c = create_concrete(fck=fck, constitutive_law=compression)
    concrete_t = create_concrete(fck=fck, constitutive_law=tension)

    concrete_ct = create_concrete(fck=fck, constitutive_law=compression_tension)

    # Reinforcement
    fyk = 500
    ftk = 550
    Es = 200000
    epsuk = 0.07

    reinforcement = create_reinforcement(fyk=fyk, Es=Es, ftk=ftk, epsuk=epsuk)

    # Create a rectangular geometry
    width = 1000
    height = 1000
    polygon = Polygon(
        [
            (-width / 2, -height / 2),
            (width / 2, -height / 2),
            (width / 2, height / 2),
            (-width / 2, height / 2),
        ]
    )

    # Concrete geometry sections
    geometry_c = SurfaceGeometry(poly=polygon, material=concrete_c)
    geometry_t = SurfaceGeometry(poly=polygon, material=concrete_t)
    geometry_ct = SurfaceGeometry(poly=polygon, material=concrete_ct)

    # Reinforced concrete geometry sections (copy of concrete geometry sections with added reinforcement)
    diameter_reinf = 40
    cover = 50

    geometry_c_rebar = add_reinforcement(
        geometry_c,
        (
            -width / 2 + cover + diameter_reinf / 2,
            -height / 2 + cover + diameter_reinf / 2,
        ),
        diameter_reinf,
        reinforcement,
    )
    geometry_c_rebar = add_reinforcement(
        geometry_c_rebar,
        (
            width / 2 - cover - diameter_reinf / 2,
            -height / 2 + cover + diameter_reinf / 2,
        ),
        diameter_reinf,
        reinforcement,
    )

    geometry_t_rebar = add_reinforcement(
        geometry_t,
        (
            -width / 2 + cover + diameter_reinf / 2,
            -height / 2 + cover + diameter_reinf / 2,
        ),
        diameter_reinf,
        reinforcement,
    )

    geometry_t_rebar = add_reinforcement(
        geometry_t_rebar,
        (
            width / 2 - cover - diameter_reinf / 2,
            -height / 2 + cover + diameter_reinf / 2,
        ),
        diameter_reinf,
        reinforcement,
    )
    geometry_ct_rebar = add_reinforcement(
        geometry_ct,
        (
            -width / 2 + cover + diameter_reinf / 2,
            -height / 2 + cover + diameter_reinf / 2,
        ),
        diameter_reinf,
        reinforcement,
    )  # The add_reinforcement function returns a CompoundGeometry
    geometry_ct_rebar = add_reinforcement(
        geometry_ct_rebar,
        (
            width / 2 - cover - diameter_reinf / 2,
            -height / 2 + cover + diameter_reinf / 2,
        ),
        diameter_reinf,
        reinforcement,
    )
    return (
        compression,
        compression_tension,
        fct,
        geometry_c,
        geometry_c_rebar,
        geometry_ct,
        geometry_ct_rebar,
        geometry_t,
        geometry_t_rebar,
        tension,
    )


@app.cell
def _(
    BeamSection,
    geometry_c,
    geometry_c_rebar,
    geometry_ct,
    geometry_ct_rebar,
    geometry_t,
    geometry_t_rebar,
):
    # Create section
    section_c = BeamSection(geometry_c)  # , integrator="fiber")
    section_t = BeamSection(geometry_t)  # , integrator="fiber")
    section_ct = BeamSection(geometry_ct)  # , integrator="fiber")
    section_c_rebar = BeamSection(geometry_c_rebar)  # , integrator="fiber")
    section_t_rebar = BeamSection(geometry_t_rebar)  # , integrator="fiber")
    section_ct_rebar = BeamSection(geometry_ct_rebar)  # , integrator="fiber")
    return (
        section_c,
        section_c_rebar,
        section_ct,
        section_ct_rebar,
        section_t,
        section_t_rebar,
    )


@app.cell
def _(
    section_c,
    section_c_rebar,
    section_ct,
    section_ct_rebar,
    section_t,
    section_t_rebar,
):
    myrd_section_t_marin = (
        section_t.section_calculator.calculate_bending_strength().m_y
    )
    myrd_section_c_marin = (
        section_c.section_calculator.calculate_bending_strength().m_y
    )
    myrd_section_ct_marin = (
        section_ct.section_calculator.calculate_bending_strength().m_y
    )
    myrd_section_c_rebar_marin = (
        section_c_rebar.section_calculator.calculate_bending_strength().m_y
    )
    myrd_section_t_rebar_marin = (
        section_t_rebar.section_calculator.calculate_bending_strength().m_y
    )
    myrd_section_ct_rebar_marin = (
        section_ct_rebar.section_calculator.calculate_bending_strength().m_y
    )
    return (
        myrd_section_c_marin,
        myrd_section_c_rebar_marin,
        myrd_section_ct_marin,
        myrd_section_ct_rebar_marin,
        myrd_section_t_marin,
        myrd_section_t_rebar_marin,
    )


@app.cell
def _(
    BeamSection,
    geometry_c,
    geometry_c_rebar,
    geometry_ct,
    geometry_ct_rebar,
    geometry_t,
    geometry_t_rebar,
):
    # Create section
    section_c_fiber = BeamSection(geometry_c, integrator="fiber")
    section_t_fiber = BeamSection(geometry_t, integrator="fiber")
    section_ct_fiber = BeamSection(geometry_ct, integrator="fiber")
    section_c_rebar_fiber = BeamSection(geometry_c_rebar, integrator="fiber")
    section_t_rebar_fiber = BeamSection(geometry_t_rebar, integrator="fiber")
    section_ct_rebar_fiber = BeamSection(geometry_ct_rebar, integrator="fiber")
    return (
        section_c_fiber,
        section_c_rebar_fiber,
        section_ct_fiber,
        section_ct_rebar_fiber,
        section_t_fiber,
        section_t_rebar_fiber,
    )


@app.cell
def _(
    section_c_fiber,
    section_c_rebar_fiber,
    section_ct_fiber,
    section_ct_rebar_fiber,
    section_t_fiber,
    section_t_rebar_fiber,
):
    myrd_section_c_fiber = (
        section_c_fiber.section_calculator.calculate_bending_strength().m_y
    )
    myrd_section_t_fiber = (
        section_t_fiber.section_calculator.calculate_bending_strength().m_y
    )
    myrd_section_ct_fiber = (
        section_ct_fiber.section_calculator.calculate_bending_strength().m_y
    )
    myrd_section_c_rebar_fiber = (
        section_c_rebar_fiber.section_calculator.calculate_bending_strength().m_y
    )
    myrd_section_t_rebar_fiber = (
        section_t_rebar_fiber.section_calculator.calculate_bending_strength().m_y
    )
    myrd_section_ct_rebar_fiber = (
        section_ct_rebar_fiber.section_calculator.calculate_bending_strength().m_y
    )
    return (
        myrd_section_c_fiber,
        myrd_section_c_rebar_fiber,
        myrd_section_ct_fiber,
        myrd_section_ct_rebar_fiber,
        myrd_section_t_fiber,
        myrd_section_t_rebar_fiber,
    )


@app.cell
def _(
    myrd_section_c_marin,
    myrd_section_c_rebar_marin,
    myrd_section_ct_marin,
    myrd_section_ct_rebar_marin,
    myrd_section_t_marin,
    myrd_section_t_rebar_marin,
):
    print(f"Results with marin integrator :")
    print(f"Section C bending strength m_y: {myrd_section_c_marin:.0f}")
    print(f"Section T bending strength m_y: {myrd_section_t_marin:.0f}")
    print(f"Section CT bending strength m_y: {myrd_section_ct_marin:.0f}")
    print(
        f"Section C Rebar bending strength m_y: {myrd_section_c_rebar_marin:.0f}"
    )
    print(
        f"Section T Rebar bending strength m_y: {myrd_section_t_rebar_marin:.0f}"
    )
    print(
        f"Section CT Rebar bending strength m_y: {myrd_section_ct_rebar_marin:.0f}"
    )
    return


@app.cell
def _(
    myrd_section_c_fiber,
    myrd_section_c_rebar_fiber,
    myrd_section_ct_fiber,
    myrd_section_ct_rebar_fiber,
    myrd_section_t_fiber,
    myrd_section_t_rebar_fiber,
):
    print(f"Results with fiber integrator :")
    print(f"Section C bending strength m_y: {myrd_section_c_fiber:.0f}")
    print(f"Section T bending strength m_y: {myrd_section_t_fiber:.0f}")
    print(f"Section CT bending strength m_y: {myrd_section_ct_fiber:.0f}")
    print(
        f"Section C Rebar bending strength m_y: {myrd_section_c_rebar_fiber:.0f}"
    )
    print(
        f"Section T Rebar bending strength m_y: {myrd_section_t_rebar_fiber:.0f}"
    )
    print(
        f"Section CT Rebar bending strength m_y: {myrd_section_ct_rebar_fiber:.0f}"
    )
    return


@app.cell
def _(
    myrd_section_c_fiber,
    myrd_section_c_marin,
    myrd_section_c_rebar_fiber,
    myrd_section_c_rebar_marin,
    myrd_section_ct_fiber,
    myrd_section_ct_marin,
    myrd_section_ct_rebar_fiber,
    myrd_section_ct_rebar_marin,
    myrd_section_t_fiber,
    myrd_section_t_marin,
    myrd_section_t_rebar_fiber,
    myrd_section_t_rebar_marin,
):
    # % of difference between fiber and marin integrator results
    print(
        f"% Difference in Section C bending strength m_y: {(myrd_section_c_fiber - myrd_section_c_marin) / myrd_section_c_marin * 100:.2f}"
    )
    print(
        f"% Difference in Section T bending strength m_y: {(myrd_section_t_fiber - myrd_section_t_marin) / myrd_section_t_marin * 100:.2f}"
    )
    print(
        f"% Difference in Section CT bending strength m_y: {(myrd_section_ct_fiber - myrd_section_ct_marin) / myrd_section_ct_marin * 100:.2f}"
    )
    print(
        f"% Difference in Section C Rebar bending strength m_y: {(myrd_section_c_rebar_fiber - myrd_section_c_rebar_marin) / myrd_section_c_rebar_marin * 100:.2f}"
    )
    print(
        f"% Difference in Section T Rebar bending strength m_y: {(myrd_section_t_rebar_fiber - myrd_section_t_rebar_marin) / myrd_section_t_rebar_marin * 100:.2f}"
    )
    print(
        f"% Difference in Section CT Rebar bending strength m_y: {(myrd_section_ct_rebar_fiber - myrd_section_ct_rebar_marin) / myrd_section_ct_rebar_marin * 100:.2f}"
    )
    return


@app.cell
def _():
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Results FIBER:
    Section C bending strength m_y: 0.0
    Section T bending strength m_y: 0.0
    Section CT bending strength m_y: -469571352.0398918

    Results MARIN:
    Section C bending strength m_y: 0.0
    Section T bending strength m_y: 0.0
    Section CT bending strength m_y: -473490379.7815997
    """)
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
