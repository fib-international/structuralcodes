# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "marimo>=0.25.2",
#     "matplotlib>=3.11.2",
#     "numpy>=2.5.4",
#     "shapely>=2.2.0",
#     "structuralcodes>=0.7.2",
# ]
# ///

import marimo

__generated_with = "0.25.1"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # ElasticFragileTension — les courbes, en interactif

    Toutes les fonctions qui ont servi à tracer les figures du rapport, réunies
    dans un notebook réactif. Bougez les curseurs : tout se retrace.

    - **§1** la loi de comportement — la figure `Figure_const_elasticfragiletension`
    - **§2** ce que renvoient `__marin__` et `__marin_tangent__`
    - **§3** les profils σ(z) et E_t(z) : branches Marin vs loi de référence
    - **§4** l'exactitude de l'intégration : Marin vs fibre vs solution fermée
    - **§5** la combinaison `Parallel(ParabolaRectangle + ElasticFragileTension)`

    Prérequis : `pip install marimo numpy scipy shapely triangle matplotlib`, et
    `structuralcodes` accessible —
    `PYTHONPATH=/chemin/vers/structuralcodes marimo edit notebooks/eft_marin.py`.
    """)
    return


@app.cell
def _():
    import matplotlib

    matplotlib.use("Agg")

    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    from shapely import Polygon

    from structuralcodes.geometry import SurfaceGeometry
    from structuralcodes.materials.basic import GenericMaterial
    from structuralcodes.materials.constitutive_laws import (
        ElasticFragileTension,
        ParabolaRectangle,
        Parallel,
    )
    from structuralcodes.sections import GenericSection

    return (
        ElasticFragileTension,
        GenericMaterial,
        GenericSection,
        ParabolaRectangle,
        Parallel,
        Polygon,
        SurfaceGeometry,
        mo,
        np,
        plt,
    )


@app.cell
def _(plt):
    # Palette et style communs à toutes les figures d'analyse.
    SURFACE, INK, INK2, MUTED, GRID = (
        "#fcfcfb",
        "#0b0b0b",
        "#52514e",
        "#8a8880",
        "#e8e7e3",
    )
    S1, S2, S3, BAND = "#2a78d6", "#eb6834", "#1baf7a", "#cde2fb"

    plt.rcParams.update(
        {
            "figure.facecolor": SURFACE,
            "axes.facecolor": SURFACE,
            "savefig.facecolor": SURFACE,
            "axes.edgecolor": GRID,
            "axes.linewidth": 1.0,
            "axes.labelcolor": INK2,
            "axes.titlecolor": INK,
            "axes.titlesize": 11,
            "axes.titleweight": "bold",
            "axes.labelsize": 9.5,
            "axes.grid": True,
            "grid.color": GRID,
            "grid.linewidth": 0.8,
            "xtick.color": INK2,
            "ytick.color": INK2,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "legend.frameon": False,
            "legend.fontsize": 8.5,
            "legend.labelcolor": INK2,
            "lines.linewidth": 2.0,
            "font.size": 9.5,
        }
    )

    def tidy(ax):
        """Retire les bords haut/droite et passe la grille derrière les tracés."""
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.set_axisbelow(True)
        return ax

    return BAND, INK2, MUTED, S1, S2, S3, SURFACE, tidy


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 1 · La loi de comportement
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    fct_ui = mo.ui.slider(
        0.5, 6.0, 0.1, value=2.7, label="$f_{ct}$ [MPa]", show_value=True
    )
    ec_ui = mo.ui.slider(
        15000, 50000, 1000, value=36000, label="$E_c$ [MPa]", show_value=True
    )
    mo.hstack([fct_ui, ec_ui], justify="start", gap=2)
    return ec_ui, fct_ui


@app.cell
def _(ElasticFragileTension, ec_ui, fct_ui, mo, np):
    law = ElasticFragileTension(fct=fct_ui.value, Ec=ec_ui.value)

    # La loi déclare-t-elle la fissuration comme une ruine ? C'est le correctif
    # discuté dans le rapport : le code patché renvoie (-inf, +inf).
    _ult = law.get_ultimate_strain()
    _patched = not np.isfinite(_ult[1])

    mo.md(
        f"""
    $\\varepsilon_{{ct,u}} = f_{{ct}} / E_c$ = **{law.eps_ct_u:.3e}**

    | | |
    |---|---|
    | `get_ultimate_strain()` | `{_ult}` |
    | `get_ultimate_strain(yielding=True)` | `{law.get_ultimate_strain(yielding=True)}` |
    | état du code | **{"patché" if _patched else "non patché"}** |
    """
    )
    return (law,)


@app.cell
def _(np, plt):
    def plot_constitutive_law(fct, Ec, dark=False):
        """La figure de la doc théorie : σ(ε) avec ses trois domaines.

        Reproduit `Figure_const_elasticfragiletension_{light,dark}.png` dans le
        style des autres figures de `docs/theory/`.
        """
        eps_ct_u = fct / Ec
        ink = "white" if dark else "black"
        x = np.linspace(-2 * eps_ct_u, 3 * eps_ct_u, 2001)
        y = np.where((x > 0) & (x <= eps_ct_u), Ec * x, 0.0)

        fig, ax = plt.subplots(figsize=(8.0, 5.0), dpi=110)
        ax.set_facecolor("none")
        fig.patch.set_alpha(0.0)
        ax.grid(False)
        ax.axvspan(-2 * eps_ct_u, 0, color="#e8e8e8")
        ax.axvspan(0, eps_ct_u, color="#d6eaf1")
        ax.axvspan(eps_ct_u, 3 * eps_ct_u, color="#e8e8e8")
        ax.plot(x, y, color="tab:blue", lw=2.5, solid_capstyle="butt")
        ax.axhline(0, color=ink, lw=0.9)
        ax.axvline(0, color=ink, lw=0.9)

        ax.text(-eps_ct_u, 0.87 * fct, r"$\varepsilon < 0$", ha="center", color=ink)
        ax.text(
            eps_ct_u / 2,
            0.87 * fct,
            r"$0 < \varepsilon \leq \varepsilon_{ct,u}$",
            ha="center",
            color=ink,
        )
        ax.text(
            2 * eps_ct_u,
            0.87 * fct,
            r"$\varepsilon > \varepsilon_{ct,u}$",
            ha="center",
            color=ink,
        )

        ax.set_xlim(-2 * eps_ct_u, 3 * eps_ct_u)
        ax.set_ylim(-0.13 * fct, 1.11 * fct)
        ax.set_xticks([0, eps_ct_u])
        ax.set_xticklabels(["0", r"$\varepsilon_{ct,u}$"])
        ax.set_yticks([0, fct])
        ax.set_yticklabels(["0", r"$f_{ct}$"])
        ax.set_ylabel(r"$\sigma$ [MPa]", color=ink)
        ax.tick_params(colors=ink)
        for spine in ax.spines.values():
            spine.set_color(ink)
        fig.tight_layout()
        return fig

    return (plot_constitutive_law,)


@app.cell(hide_code=True)
def _(mo):
    dark_ui = mo.ui.checkbox(value=False, label="variante sombre (texte blanc)")
    dark_ui
    return (dark_ui,)


@app.cell
def _(dark_ui, ec_ui, fct_ui, plot_constitutive_law):
    plot_constitutive_law(fct_ui.value, ec_ui.value, dark=dark_ui.value)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 2 · Ce que renvoient `__marin__` et `__marin_tangent__`

    Le plan de déformation est $\varepsilon(z) = \varepsilon_a + \chi_y \cdot z$.
    Mettez $\chi_y$ à zéro pour emprunter l'autre chemin de code — celui de la
    déformation uniforme, qui renvoie `strains = None`.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    eps_a_ui = mo.ui.slider(
        -30, 30, 1, value=-4, label=r"$\varepsilon_a$ [$\times 10^{-5}$]", show_value=True
    )
    chi_ui = mo.ui.slider(
        -30, 30, 1, value=6, label=r"$\chi_y$ [$\times 10^{-7}$ mm$^{-1}$]", show_value=True
    )
    h_ui = mo.ui.slider(200, 1000, 50, value=500, label="h [mm]", show_value=True)
    b_ui = mo.ui.slider(150, 600, 50, value=300, label="b [mm]", show_value=True)
    mo.hstack([eps_a_ui, chi_ui, b_ui, h_ui], justify="start", gap=1.5)
    return b_ui, chi_ui, eps_a_ui, h_ui


@app.cell
def _(chi_ui, eps_a_ui):
    eps_a = eps_a_ui.value * 1e-5
    chi_y = chi_ui.value * 1e-7
    return chi_y, eps_a


@app.cell
def _(mo):
    def marin_summary(law, eps_a, chi_y):
        """Tableau de ce que les deux méthodes renvoient pour ce plan.

        `__marin__` mute `strain[0]` quand la courbure est nulle : on lui passe
        donc une liste neuve à chaque appel.
        """
        rows = []
        for name, method in (
            ("`__marin__`", law.__marin__),
            ("`__marin_tangent__`", law.__marin_tangent__),
        ):
            strains, coeffs = method(strain=[eps_a, chi_y])
            if strains is None:
                rows.append(f"| {name} | `None` (toute la géométrie) | `{coeffs[0]}` |")
            else:
                for (e0, e1), coeff in zip(strains, coeffs):
                    pretty = tuple(round(c, 6) for c in coeff)
                    rows.append(f"| {name} | `({e0:.3g}, {e1:.3g})` | `{pretty}` |")
        return mo.md(
            "| méthode | plage de déformation | coefficients $(a_0, a_1, \\ldots)$ |\n"
            "|---|---|---|\n" + "\n".join(rows)
        )

    return (marin_summary,)


@app.cell
def _(chi_y, eps_a, law, marin_summary):
    marin_summary(law, eps_a, chi_y)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 3 · Les branches Marin reproduisent-elles la loi ?

    On reconstruit $\sigma(z)$ et $E_t(z)$ à partir des seules branches
    renvoyées — exactement comme le fait l'intégrateur — et on compare à
    `get_stress` / `get_tangent`. Ce qui n'est pas renvoyé n'est pas intégré :
    la zone comprimée et la zone fissurée valent zéro.
    """)
    return


@app.cell
def _(np):
    def marin_profile(law, eps_a, chi_y, z, tangent=False):
        """Évalue les branches Marin sur la hauteur de section.

        Pour chaque plage de déformation renvoyée, on convertit les bornes en
        cotes z, puis on somme le polynôme $\\sum a_n z^n$ sur cette tranche.
        Hors des plages renvoyées : zéro.
        """
        method = law.__marin_tangent__ if tangent else law.__marin__
        strains, coeffs = method(strain=[eps_a, chi_y])
        z = np.atleast_1d(np.asarray(z, dtype=float))
        out = np.zeros_like(z)

        if strains is None:
            for order, a_n in enumerate(coeffs[0]):
                out += a_n * z**order
            return out

        for (eps_start, eps_end), coeff in zip(strains, coeffs):
            z_0, z_1 = sorted(
                ((eps_start - eps_a) / chi_y, (eps_end - eps_a) / chi_y)
            )
            inside = (z >= z_0) & (z <= z_1)
            for order, a_n in enumerate(coeff):
                out[inside] += a_n * z[inside] ** order
        return out

    return (marin_profile,)


@app.cell
def _(BAND, INK2, MUTED, S1, S2, marin_profile, np, plt, tidy):
    def plot_profiles(law, eps_a, chi_y, h):
        """σ(z) et E_t(z) : loi de référence vs branches Marin."""
        z = np.linspace(-h / 2, h / 2, 4001)
        eps = eps_a + chi_y * z

        fig, axs = plt.subplots(1, 2, figsize=(10.5, 4.6))
        panels = (
            (axs[0], law.get_stress, False, r"$\sigma$  [MPa]", "(a)  Contraintes", 1.0),
            (axs[1], law.get_tangent, True, r"$E_t$  [GPa]", "(b)  Module tangent", 1e-3),
        )

        for ax, reference, tangent, xlabel, title, scale in panels:
            tidy(ax)
            if chi_y != 0:
                z_0, z_1 = sorted((-eps_a / chi_y, (law.eps_ct_u - eps_a) / chi_y))
                ax.axhspan(
                    max(z_0, -h / 2), min(z_1, h / 2), color=BAND, alpha=0.45, lw=0
                )
            # get_stress recale le tableau en place : on lui passe une copie
            ax.plot(
                reference(eps.copy()) * scale,
                z,
                color=S1,
                lw=3.2,
                label="loi de référence",
            )
            ax.plot(
                marin_profile(law, eps_a, chi_y, z, tangent=tangent) * scale,
                z,
                color=S2,
                lw=1.6,
                ls=(0, (5, 3)),
                label="branches Marin",
            )
            ax.set_xlabel(xlabel)
            ax.set_ylabel("z  [mm]")
            ax.set_title(title)
            ax.legend(loc="lower right")

        gap = np.max(
            np.abs(
                law.get_stress(eps.copy())
                - marin_profile(law, eps_a, chi_y, z, tangent=False)
            )
        )
        note = (
            "accord exact"
            if gap < 1e-9
            else f"écart max {gap:.3g} MPa — recalage $\\pm$1e-6 de get_stress"
        )
        fig.text(0.5, 0.005, note, ha="center", color=INK2 if gap < 1e-9 else MUTED)
        fig.tight_layout(rect=[0, 0.035, 1, 1])
        return fig

    return (plot_profiles,)


@app.cell
def _(chi_y, eps_a, h_ui, law, plot_profiles):
    plot_profiles(law, eps_a, chi_y, h_ui.value)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 4 · Exactitude de l'intégration

    Pour une section rectangulaire et cette loi, les résultantes ont une forme
    fermée : la bande élastique est un trapèze de contrainte sur
    $[z_0, z_1]$, donc $N$ et $M_y$ s'intègrent à la main. On compare Marin et
    fibre à cette solution.
    """)
    return


@app.cell
def _(GenericMaterial, GenericSection, Polygon, SurfaceGeometry, np):
    def make_section(law, integrator, b, h, mesh_size=0.0002):
        """Section rectangulaire b×h centrée sur (0,0), faite de cette seule loi."""
        material = GenericMaterial(density=2500, constitutive_law=law)
        polygon = Polygon(
            (
                (-b / 2, -h / 2),
                (b / 2, -h / 2),
                (b / 2, h / 2),
                (-b / 2, h / 2),
            )
        )
        geometry = SurfaceGeometry(polygon, material)
        return GenericSection(geometry, integrator=integrator, mesh_size=mesh_size)

    def integrate(section, strain, quantity="stress"):
        """Appelle l'intégrateur de la section pour un plan (eps_a, chi_y, chi_z)."""
        result = (
            section.section_calculator.integrator.integrate_strain_response_on_geometry(
                section.geometry,
                list(strain),
                integrate=quantity,
                mesh_size=section.section_calculator.mesh_size,
            )
        )
        return np.array(result[:3]) if quantity == "stress" else result[0]

    def exact_resultants(law, eps_a, chi_y, b, h):
        """Solution fermée de N et M_y pour la bande élastique 0 < ε ≤ ε_ct,u."""
        eps_ct_u = law.eps_ct_u
        if chi_y == 0:
            stress = law.Ec * eps_a if 0 < eps_a <= eps_ct_u else 0.0
            return np.array([stress * b * h, 0.0, 0.0])

        z_0, z_1 = sorted((-eps_a / chi_y, (eps_ct_u - eps_a) / chi_y))
        z_0, z_1 = max(z_0, -h / 2), min(z_1, h / 2)
        if z_1 <= z_0:
            return np.zeros(3)

        return np.array(
            [
                law.Ec * b * (eps_a * (z_1 - z_0) + chi_y * (z_1**2 - z_0**2) / 2),
                law.Ec
                * b
                * (eps_a * (z_1**2 - z_0**2) / 2 + chi_y * (z_1**3 - z_0**3) / 3),
                0.0,
            ]
        )

    def exact_stiffness(law, eps_a, chi_y, b, h):
        """Solution fermée de EA, ES, EI pour la même bande."""
        eps_ct_u = law.eps_ct_u
        if chi_y == 0:
            modulus = law.Ec if 0 < eps_a <= eps_ct_u else 0.0
            return np.array([modulus * b * h, 0.0, modulus * b * h**3 / 12])

        z_0, z_1 = sorted((-eps_a / chi_y, (eps_ct_u - eps_a) / chi_y))
        z_0, z_1 = max(z_0, -h / 2), min(z_1, h / 2)
        if z_1 <= z_0:
            return np.zeros(3)

        return np.array(
            [
                law.Ec * b * (z_1 - z_0),
                law.Ec * b * (z_1**2 - z_0**2) / 2,
                law.Ec * b * (z_1**3 - z_0**3) / 3,
            ]
        )

    return exact_resultants, integrate, make_section


@app.cell
def _(
    INK2,
    S1,
    S2,
    S3,
    SURFACE,
    exact_resultants,
    integrate,
    make_section,
    np,
    plt,
    tidy,
):
    def plot_accuracy(law, b, h, n_points=14):
        """N(χ), M_y(χ) et l'erreur relative : Marin, fibre, solution fermée."""
        section_marin = make_section(law, "marin", b, h)
        section_fiber = make_section(law, "fiber", b, h)
        curvatures = np.logspace(-7.4, -4.5, n_points)

        marin = np.array([integrate(section_marin, [0.0, k, 0]) for k in curvatures])
        fiber = np.array([integrate(section_fiber, [0.0, k, 0]) for k in curvatures])
        exact = np.array([exact_resultants(law, 0.0, k, b, h) for k in curvatures])

        fig, axs = plt.subplots(1, 3, figsize=(14.0, 4.3))

        for ax, column, scale, ylabel, title in (
            (axs[0], 0, 1e-3, "N  [kN]", "(a)  Effort normal"),
            (axs[1], 1, 1e-6, r"$M_y$  [kNm]", "(b)  Moment"),
        ):
            tidy(ax)
            ax.plot(
                curvatures,
                exact[:, column] * scale,
                color=INK2,
                lw=4.0,
                alpha=0.35,
                label="solution exacte",
            )
            ax.plot(curvatures, marin[:, column] * scale, color=S1, label="Marin")
            ax.plot(
                curvatures,
                fiber[:, column] * scale,
                color=S2,
                ls=(0, (5, 3)),
                label="Fibre",
            )
            ax.set_xscale("log")
            ax.set_xlabel(r"courbure $\chi_y$  [mm$^{-1}$]")
            ax.set_ylabel(ylabel)
            ax.set_title(title)
            ax.legend(loc="upper right")

        ax = tidy(axs[2])
        safe = np.where(exact[:, 0] == 0, np.nan, exact[:, 0])
        for values, color, marker, label in (
            (marin[:, 0], S1, "o", "Marin"),
            (fiber[:, 0], S2, "s", "Fibre"),
        ):
            error = np.abs(100 * (values - exact[:, 0]) / safe)
            ax.plot(
                curvatures,
                np.maximum(error, 1e-15),
                color=color,
                marker=marker,
                ms=4.5,
                mfc=color,
                mec=SURFACE,
                mew=1.2,
                ls="-" if color == S1 else (0, (5, 3)),
                label=label,
            )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylim(3e-16, 3e2)
        ax.set_xlabel(r"courbure $\chi_y$  [mm$^{-1}$]")
        ax.set_ylabel("|erreur|  [%]")
        ax.set_title("(c)  Erreur relative sur N")
        ax.legend(loc="center left")
        ax.axhline(1.0, color=S3, lw=1, ls=(0, (2, 3)))

        fig.tight_layout()
        return fig

    return (plot_accuracy,)


@app.cell(hide_code=True)
def _(mo):
    sweep_ui = mo.ui.run_button(
        label="Lancer le balayage (l'intégrateur fibre prend ~2 s)"
    )
    sweep_ui
    return (sweep_ui,)


@app.cell
def _(b_ui, h_ui, law, mo, plot_accuracy, sweep_ui):
    mo.stop(
        not sweep_ui.value,
        mo.md("*Cliquez sur le bouton pour lancer le balayage en courbure.*"),
    )
    plot_accuracy(law, b_ui.value, h_ui.value)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 5 · La combinaison `Parallel`

    La loi est faite pour être associée à une loi de compression.
    `Parallel._split_marin` fusionne les plages des deux lois : la bande de
    traction ne porte que la loi de traction, les branches de compression
    restent intactes.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    fc_ui = mo.ui.slider(
        20, 80, 5, value=30, label="$f_c$ [MPa] (ParabolaRectangle)", show_value=True
    )
    fc_ui
    return (fc_ui,)


@app.cell
def _(BAND, S1, S2, marin_profile, np, plt, tidy):
    def plot_parallel(combined, eps_ct_u, eps_a, chi_y, h):
        """σ(z) de la loi combinée : branches Marin vs get_stress."""
        z = np.linspace(-h / 2, h / 2, 4001)
        eps = eps_a + chi_y * z

        fig, ax = plt.subplots(figsize=(7.0, 5.0))
        tidy(ax)
        if chi_y != 0:
            z_0, z_1 = sorted((-eps_a / chi_y, (eps_ct_u - eps_a) / chi_y))
            ax.axhspan(
                max(z_0, -h / 2), min(z_1, h / 2), color=BAND, alpha=0.45, lw=0
            )
        ax.plot(combined.get_stress(eps.copy()), z, color=S1, lw=3.2, label="get_stress")
        ax.plot(
            marin_profile(combined, eps_a, chi_y, z),
            z,
            color=S2,
            lw=1.6,
            ls=(0, (5, 3)),
            label="branches Marin",
        )
        ax.set_xlabel(r"$\sigma$  [MPa]")
        ax.set_ylabel("z  [mm]")
        ax.set_title("Parallel(ParabolaRectangle + ElasticFragileTension)")
        ax.legend(loc="lower right")
        fig.tight_layout()
        return fig

    return (plot_parallel,)


@app.cell
def _(
    ParabolaRectangle,
    Parallel,
    chi_y,
    eps_a,
    fc_ui,
    h_ui,
    law,
    plot_parallel,
):
    combined_law = Parallel([ParabolaRectangle(fc=fc_ui.value), law])
    plot_parallel(combined_law, law.eps_ct_u, eps_a, chi_y, h_ui.value)
    return (combined_law,)


@app.cell
def _(chi_y, combined_law, eps_a, marin_summary):
    marin_summary(combined_law, eps_a, chi_y)
    return


if __name__ == "__main__":
    app.run()
