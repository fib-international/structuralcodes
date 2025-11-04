(theory-crs-sign-convention)=
# Coordinate system and sign convention

(theory-grs)=
## Coordinate system
StructuralCodes adopts a global reference system (GRS) with the following conventions:

**$x$-axis**
: Points out of the screen (toward the observer).

**$y$-axis**
: Points to the right on the screen.

**$z$-axis**
: Points upward on the screen.

Note that this coordinate system is consistent with graphical representations, where Shapely operates in screen coordinates $xy$, mapped to **GRS** coordinates as $yz$.

The figure [below](#theory-fig-coordinate-system) illustrates the coordinate system in a rectangular geometry.

(theory-fig-coordinate-system)=
:::::{grid}
::::{grid-item}
:class: caption-text, sd-text-center

:::{image} FigureGRS_light.png
:align: center
:class: only-light
:::

:::{image} FigureGRS_dark.png
:align: center
:class: only-dark
:::

The coordinate system adopted in StructuralCodes shown on a rectangular geometry.

::::
:::::


(theory-sign-convention)=
## Sign conventions

The following sign conventions apply.

**Forces**
: Are negative when in compression.

**Moments**
: Follow the right-hand rule as illustrated in the figure [below](#theory-fig-moment-signs):
   - $M_{\textrm{y}}$, bending about the $y$-axis, is positive when top fibers are stretched, and bottom fibers are compressed.
   - $M_{\textrm{z}}$, bending about the $z$-axis, is positive when left fibers are stretched, and right fibers are compressed.

**Stresses** and **strains**
: Are positive in tension and negative in compression.

**Loads**
: Act in the origin of the **GRS**.

(theory-fig-moment-signs)=
:::::{grid}
::::{grid-item}
:class: caption-text, sd-text-center

:::{image} FigureSigns_light.png
:align: center
:class: only-light
:::

:::{image} FigureSigns_dark.png
:align: center
:class: only-dark
:::

The definition of positive moments.

::::
:::::

::::::{admonition} Loads act in the origin!
:class: caution

Pay particular attention to **loads**. When the section is subjected to axial load, this is considered acting in the origin, i.e. `(0, 0)`. If the center of gravity of the section is not aligned with the origin, offset moments are generated, as illustrated in the figure [below](#theory-fig-geometry-not-in-origin). If you want the load to act on the center of the geometry, translate the geometry in order to have the center in `(0, 0)`.

(theory-fig-geometry-not-in-origin)=
:::::{grid}
::::{grid-item}
:class: caption-text, sd-text-center

:::{image} FigureLoadOrigin_light.png
:align: center
:class: only-light
:width: 60%
:::

:::{image} FigureLoadOrigin_dark.png
:align: center
:class: only-dark
:width: 60%
:::

Illustration of the offset moments that are generated when a geometry which is not aligned with the origin is subjected to an axial force.

::::
:::::

::::::