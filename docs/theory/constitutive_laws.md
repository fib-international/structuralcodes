(theory-constitutive-laws)=
# Constitutive laws

## General concept

A uniaxial constitutive law in StructuralCodes gives the stress $\sigma$ for a given strain $\varepsilon$. Following the [sign convention](theory-sign-convention), stress and strain are positive in tension and negative in compression.

:::{math}
:label: eq:constitutive-stress
\sigma(\varepsilon) = f(\varepsilon)
:::

The tangent modulus $E_{\mathrm{t}}$ is the slope of the stress–strain curve on a differentiable branch.

:::{math}
:label: eq:constitutive-tangent
E_{\mathrm{t}}(\varepsilon) = \frac{\mathrm{d}\sigma}{\mathrm{d}\varepsilon}.
:::

The secant modulus $E_{\mathrm{s}}$ is obtained by dividing the stress by the strain. At zero strain, the implementation returns the tangent as the secant modulus.

:::{math}
:label: eq:constitutive-moduli
E_{\mathrm{s}}(\varepsilon) =
\begin{cases}
\sigma(\varepsilon)/\varepsilon, & \varepsilon \ne 0, \\
E_{\mathrm{t}}(0), & \varepsilon = 0.
\end{cases}
:::

At a branch junction, the value returned by the tangent method follows the convention of the individual law. Some laws also snap strains close to a reported ultimate strain to that limit before evaluating the response.

The laws below implement the common uniaxial interface, which can also be used to introduce new laws. See the [constitutive laws API](api-constitutive-laws) for the class interfaces.

(theory-constitutive-elastic)=
## Linear elastic

The {class}`Elastic <structuralcodes.materials.constitutive_laws.Elastic>` law has a constant elastic modulus $E$ in tension and compression.

:::{math}
:label: eq:constitutive-elastic
\sigma(\varepsilon) = E\varepsilon.
:::

:::{math}
:label: eq:constitutive-elastic-tangent
E_{\mathrm{t}}(\varepsilon) = E.
:::

Ultimate strains can be specified for section calculations. The stress–strain curve extends linearly through these limits.

(figure-theory-law-elastic)=
:::::{grid}
::::{grid-item}
:class: caption-text, sd-text-center

:::{image} Figure_law_elastic_light.png
:align: center
:alt: Linear elastic response in compression and tension.
:class: only-light
:::

:::{image} Figure_law_elastic_dark.png
:align: center
:alt: Linear elastic response in compression and tension.
:class: only-dark
:::

Linear elastic response in compression and tension with $E=200\,\mathrm{GPa}$.

::::
:::::

(theory-constitutive-elastic-plastic)=
## Elastic-plastic

The {class}`ElasticPlastic <structuralcodes.materials.constitutive_laws.ElasticPlastic>` law is symmetric in tension and compression. The initial slope is the elastic modulus $E$. Yielding starts at the strain magnitude $\varepsilon_y=f_y/E$, where $f_y$ is the yield strength. The post-yield slope is the hardening modulus $E_h$. In the stress function below, $\varepsilon_u$ is the specified ultimate strain magnitude, or $\infty$ if `eps_su` is omitted.

:::{math}
:label: eq:constitutive-elastic-plastic
\sigma(\varepsilon) =
\begin{cases}
E\varepsilon, & |\varepsilon| \le \varepsilon_y, \\
\operatorname{sgn}(\varepsilon)
\bigl[f_y + E_h(|\varepsilon|-\varepsilon_y)\bigr],
& \varepsilon_y < |\varepsilon| \le \varepsilon_u, \\
0, & |\varepsilon| > \varepsilon_u.
\end{cases}
:::

The tangent is $E$ up to and including yield, $E_h$ after yield up to and including the ultimate strain, and zero beyond it. With $E_h=0$, the law is elastic-perfectly plastic. When `eps_su` is omitted, the law reports $\pm 2\varepsilon_y$ as the ultimate strains for section calculations.

(figure-theory-law-elastic-plastic)=
:::::{grid}
::::{grid-item}
:class: caption-text, sd-text-center

:::{image} Figure_law_elastic_plastic_light.png
:align: center
:alt: Elastic-plastic response with linear hardening and with perfect plasticity.
:class: only-light
:::

:::{image} Figure_law_elastic_plastic_dark.png
:align: center
:alt: Elastic-plastic response with linear hardening and with perfect plasticity.
:class: only-dark
:::

Elastic-plastic response with $E=200\,\mathrm{GPa}$, $f_y=450\,\mathrm{MPa}$, and $\varepsilon_u=0.01$ for both curves. The hardening modulus is $E_h\approx12.90\,\mathrm{GPa}$ for the linear hardening curve, giving $\sigma(\varepsilon_u)=550\,\mathrm{MPa}$, and $E_h=0$ for the perfectly-plastic curve.

::::
:::::

(theory-constitutive-bilinear-compression)=
## Bilinear compression

The {class}`BilinearCompression <structuralcodes.materials.constitutive_laws.BilinearCompression>` law is elastic-perfectly plastic in compression. With the signed compressive strength $f_c<0$, the strain at the start of the plateau $\varepsilon_c<0$, and the ultimate strain $\varepsilon_{cu}\le\varepsilon_c$, the elastic modulus is $E_c=f_c/\varepsilon_c$.

:::{math}
:label: eq:constitutive-bilinear-compression
\sigma(\varepsilon) =
\begin{cases}
E_c\varepsilon, & \varepsilon_c < \varepsilon \le 0, \\
f_c, & \varepsilon_{cu} \le \varepsilon \le \varepsilon_c, \\
0, & \text{otherwise}.
\end{cases}
:::

Within the elastic branch the tangent is $E_c$. It is zero on the plateau, in tension, and beyond the ultimate strain. At the branch endpoints, a scalar call returns the elastic tangent at zero and zero at $\varepsilon_c$. An array call returns zero at zero strain and the elastic tangent at $\varepsilon_c$.

(figure-theory-law-bilinear-compression)=
:::::{grid}
::::{grid-item}
:class: caption-text, sd-text-center

:::{image} Figure_law_bilinear_compression_light.png
:align: center
:alt: Linear compression branch followed by a constant-stress plateau.
:class: only-light
:::

:::{image} Figure_law_bilinear_compression_dark.png
:align: center
:alt: Linear compression branch followed by a constant-stress plateau.
:class: only-dark
:::

Bilinear compression with $f_c=-30\,\mathrm{MPa}$, $\varepsilon_c=-0.0018$, and $\varepsilon_{cu}=-0.0035$.

::::
:::::

(theory-constitutive-parabola-rectangle)=
## Parabola-rectangle

The {class}`ParabolaRectangle <structuralcodes.materials.constitutive_laws.ParabolaRectangle>` law uses a power function from zero strain to the peak strain $\varepsilon_0<0$, followed by a plateau at the signed compressive strength $f_c<0$ down to the ultimate strain $\varepsilon_u\le\varepsilon_0$. The exponent $n$ defaults to $2$, which gives a parabola.

:::{math}
:label: eq:constitutive-parabola-rectangle
\sigma(\varepsilon) =
\begin{cases}
f_c\bigl[1-(1-\varepsilon/\varepsilon_0)^n\bigr],
& \varepsilon_0 \le \varepsilon \le 0, \\
f_c, & \varepsilon_u \le \varepsilon < \varepsilon_0, \\
0, & \text{otherwise}.
\end{cases}
:::

On the first branch, the tangent is

:::{math}
:label: eq:constitutive-parabola-rectangle-tangent
E_{\mathrm{t}}(\varepsilon)
= \frac{n f_c}{\varepsilon_0}(1-\varepsilon/\varepsilon_0)^{n-1},
\qquad \varepsilon_0 \le \varepsilon \le 0.
:::

The tangent is zero elsewhere.

(figure-theory-law-parabola-rectangle)=
:::::{grid}
::::{grid-item}
:class: caption-text, sd-text-center

:::{image} Figure_law_parabola_rectangle_light.png
:align: center
:alt: Parabola-rectangle response with exponent 2.
:class: only-light
:::

:::{image} Figure_law_parabola_rectangle_dark.png
:align: center
:alt: Parabola-rectangle response with exponent 2.
:class: only-dark
:::

Parabola-rectangle response with $f_c=-30\,\mathrm{MPa}$, $\varepsilon_0=-0.002$, $\varepsilon_u=-0.0035$, and $n=2$.

::::
:::::

(theory-constitutive-sargin)=
## Sargin

The {class}`Sargin <structuralcodes.materials.constitutive_laws.Sargin>` law uses a rational function for the ascending and descending branches of the concrete compression curve.[^sargin-1971] Let $f_c<0$ be the signed compressive strength, $\varepsilon_{c1}<0$ the strain at peak stress, $\varepsilon_{cu1}\le\varepsilon_{c1}$ the ultimate strain, and $k$ the parameter controlling the shape of the curve. With $\eta=\varepsilon/\varepsilon_{c1}$, the stress and tangent below apply over $\varepsilon_{cu1} \le \varepsilon \le 0$.

:::{math}
:label: eq:constitutive-sargin
\sigma(\varepsilon)
= f_c\frac{k\eta-\eta^2}{1+(k-2)\eta}.
:::

:::{math}
:label: eq:constitutive-sargin-tangent
E_{\mathrm{t}}(\varepsilon)
= \frac{f_c}{\varepsilon_{c1}}
  \frac{k-2\eta+(2-k)\eta^2}{[1+(k-2)\eta]^2}.
:::

The initial tangent is $k f_c/\varepsilon_{c1}$, and $\eta=1$ gives the peak compressive stress $f_c$. The stress and tangent are zero in tension and beyond $\varepsilon_{cu1}$ in compression.

(figure-theory-law-sargin)=
:::::{grid}
::::{grid-item}
:class: caption-text, sd-text-center

:::{image} Figure_law_sargin_light.png
:align: center
:alt: Sargin compression curve with ascending and descending branches.
:class: only-light
:::

:::{image} Figure_law_sargin_dark.png
:align: center
:alt: Sargin compression curve with ascending and descending branches.
:class: only-dark
:::

Sargin compression curve with $f_c=-30\,\mathrm{MPa}$, $\varepsilon_{c1}=-0.0023$, $\varepsilon_{cu1}=-0.0035$, and $k=2.04$.

::::
:::::

(theory-constitutive-popovics)=
## Popovics

The {class}`Popovics <structuralcodes.materials.constitutive_laws.Popovics>` law describes the concrete compression curve in terms of the signed compressive strength $f_c<0$, the peak strain $\varepsilon_c<0$, and the initial modulus $E_c$.[^popovics-1973] The secant modulus at the peak, $E_{\mathrm{p}}=f_c/\varepsilon_c$, defines the shape parameter $n$.

:::{math}
:label: eq:constitutive-popovics-parameters
\eta = \frac{\varepsilon}{\varepsilon_c},
\qquad
n = \frac{E_c}{E_c-E_{\mathrm{p}}}.
:::

The stress and tangent below apply on the compression branch $\varepsilon_{cu}\le\varepsilon\le 0$, where $\varepsilon_{cu}\le\varepsilon_c$ is the ultimate strain.

:::{math}
:label: eq:constitutive-popovics
\sigma(\varepsilon)
= f_c\frac{n\eta}{n-1+\eta^n}.
:::

:::{math}
:label: eq:constitutive-popovics-tangent
E_{\mathrm{t}}(\varepsilon)
= \frac{f_c}{\varepsilon_c}
  \frac{n(n-1)(1-\eta^n)}{(n-1+\eta^n)^2}.
:::

The stress and tangent are zero in tension and beyond $\varepsilon_{cu}$ in compression. If $E_c$ is omitted, the implementation uses the relation $E_c=5000\sqrt{|f_c|}$ from Mander et al. [^mander-1988], with $f_c$ and $E_c$ expressed in MPa.

(figure-theory-law-popovics)=
:::::{grid}
::::{grid-item}
:class: caption-text, sd-text-center

:::{image} Figure_law_popovics_light.png
:align: center
:alt: Popovics compression curve with ascending and descending branches.
:class: only-light
:::

:::{image} Figure_law_popovics_dark.png
:align: center
:alt: Popovics compression curve with ascending and descending branches.
:class: only-dark
:::

Popovics compression curve with $f_c=-30\,\mathrm{MPa}$, $\varepsilon_c=-0.002$, and $\varepsilon_{cu}=-0.0035$. The default relation gives $E_c\approx27.39\,\mathrm{GPa}$.

::::
:::::

(theory-constitutive-user-defined)=
## User defined

The {class}`UserDefined <structuralcodes.materials.constitutive_laws.UserDefined>` law interpolates a series of stress–strain points $(\varepsilon_i,\sigma_i)$ ordered by increasing strain. For a curve supplied from the origin into tension, the implementation reflects those points through the origin to construct the compression branch.

:::{math}
:label: eq:constitutive-user-defined
\begin{aligned}
m_i &= \frac{\sigma_{i+1}-\sigma_i}{\varepsilon_{i+1}-\varepsilon_i}, \\
\sigma(\varepsilon) &= \sigma_i + m_i(\varepsilon-\varepsilon_i),
\qquad \varepsilon_i \le \varepsilon \le \varepsilon_{i+1}.
\end{aligned}
:::

The tangent is the segment slope $m_i$. At an interior point it is the slope of the segment on the lower-strain side. At the first and last interpolation points, the tangent is taken from the adjacent segment.

The `flag` parameter controls the response beyond the original endpoints. For an end segment, let $(\varepsilon_b,\sigma_b)$ be the corresponding endpoint and $m_b$ the slope of the adjacent segment.

| `flag` | Stress beyond the endpoint |
| --- | --- |
| `0` (default) | Zero |
| `1` | $\sigma_b$ |
| `2` | $\sigma_b + m_b(\varepsilon-\varepsilon_b)$ |
| `3` | $\sigma_b + (\sigma_b/\varepsilon_b)(\varepsilon-\varepsilon_b)$ |

For flags `1`–`3`, the added segment at each end runs to $100\varepsilon_b$. Beyond the resulting interpolation range, both stress and tangent are zero. The reported ultimate strains are the original endpoints by default and can be set independently using `eps_u` initializing the object.

The example below uses a smaller positive slope beyond the first change in stiffness on each side of the origin. The stress magnitude therefore continues to increase with strain magnitude in both compression and tension. All four flags give the same curve between the supplied endpoints. Their differences appear in the extensions.

(figure-theory-law-user-defined)=
:::::{grid}
::::{grid-item}
:class: caption-text, sd-text-center

:::{image} Figure_law_user_defined_light.png
:align: center
:alt: A piecewise linear curve with hardening in compression and tension. The four flags give zero stress, constant stress, last tangent, or last secant beyond each endpoint.
:class: only-light
:::

:::{image} Figure_law_user_defined_dark.png
:align: center
:alt: A piecewise linear curve with hardening in compression and tension. The four flags give zero stress, constant stress, last tangent, or last secant beyond each endpoint.
:class: only-dark
:::

Hardening in compression and tension with strain–stress points $(-0.004,-36)$, $(-0.002,-30)$, $(0,0)$, $(0.0015,24)$, and $(0.0035,32)$, with stress in MPa. The four `flag` settings control the extensions beyond the first and last points.

::::
:::::

(theory-constitutive-initial-strain)=
## Initial strain

The {class}`InitialStrain <structuralcodes.materials.constitutive_laws.InitialStrain>` law wraps an existing law and assigns an initial strain $\varepsilon_{\mathrm{init}}$. When strain compatibility is enforced, the wrapped law is evaluated at the sum of the imposed strain and the initial strain.

:::{math}
:label: eq:constitutive-initial-strain
\sigma(\varepsilon)
= \sigma_{\mathrm{w}}(\varepsilon+\varepsilon_{\mathrm{init}}).
:::

:::{math}
:label: eq:constitutive-initial-strain-tangent
E_{\mathrm{t}}(\varepsilon)
= E_{\mathrm{t,w}}(\varepsilon+\varepsilon_{\mathrm{init}}).
:::

With strain compatibility released, as for an unbonded tendon, the pointwise stress remains at its initial value $\sigma_{\mathrm{w}}(\varepsilon_{\mathrm{init}})$. The tangent method returns $10^{-6}E_{\mathrm{t,w}}(0)$ in this case. Both reported ultimate strains of the wrapped law are shifted by $-\varepsilon_{\mathrm{init}}$.

(figure-theory-law-initial-strain)=
:::::{grid}
::::{grid-item}
:class: caption-text, sd-text-center

:::{image} Figure_law_initial_strain_light.png
:align: center
:alt: An elastic law with a positive initial strain, with strain compatibility enforced and released.
:class: only-light
:::

:::{image} Figure_law_initial_strain_dark.png
:align: center
:alt: An elastic law with a positive initial strain, with strain compatibility enforced and released.
:class: only-dark
:::

Linear elastic law with $E=200\,\mathrm{GPa}$ and initial strain $\varepsilon_{\mathrm{init}}=0.0015$. The two initial-strain curves show the response with strain compatibility enforced and released.

::::
:::::

(theory-constitutive-elastic-fragile-tension)=
## Elastic-fragile tension

The {class}`ElasticFragileTension <structuralcodes.materials.constitutive_laws.ElasticFragileTension>` law represents the tensile response of concrete up to cracking. The tensile strength $f_{ct}$ and the elastic modulus $E_c$ define the cracking strain $\varepsilon_{ct}=f_{ct}/E_c$.

:::{math}
:label: eq:constitutive-elastic-fragile-tension
\sigma(\varepsilon) =
\begin{cases}
E_c\varepsilon, & 0 < \varepsilon \le \varepsilon_{ct}, \\
0, & \text{otherwise}.
\end{cases}
:::

:::{math}
:label: eq:constitutive-elastic-fragile-tension-tangent
E_{\mathrm{t}}(\varepsilon) =
\begin{cases}
E_c, & 0 < \varepsilon \le \varepsilon_{ct}, \\
0, & \text{otherwise}.
\end{cases}
:::

The stress reaches $f_{ct}$ at $\varepsilon_{ct}$ and drops to zero immediately beyond it. Both stress and tangent are zero at the origin and in compression. A concrete compression law can be added using the [parallel law](theory-constitutive-parallel).

(figure-theory-law-elastic-fragile-tension)=
:::::{grid}
::::{grid-item}
:class: caption-text, sd-text-center

:::{image} Figure_law_elastic_fragile_tension_light.png
:align: center
:alt: Elastic response in tension up to cracking. The open marker shows the zero-stress limit immediately after cracking.
:class: only-light
:::

:::{image} Figure_law_elastic_fragile_tension_dark.png
:align: center
:alt: Elastic response in tension up to cracking. The open marker shows the zero-stress limit immediately after cracking.
:class: only-dark
:::

Elastic response in tension with $E_c=30\,\mathrm{GPa}$ and $f_{ct}=3\,\mathrm{MPa}$, giving $\varepsilon_{ct}=0.0001$. The open marker shows the zero-stress limit immediately after cracking.

::::
:::::

(theory-constitutive-parallel)=
## Parallel

### General concepts

Two springs attached to the same rigid bar undergo the same displacement $u$. Each spring develops a force proportional to its own stiffness, and the applied force $F$ is the sum of the two spring forces.

(figure-theory-law-parallel-springs)=
:::::{grid}
::::{grid-item}
:class: caption-text, sd-text-center

:::{image} Figure_law_parallel_springs_light.png
:align: center
:alt: Two springs share a fixed support and a rigid bar. The displacement arrow runs from the dashed initial bar position to the displaced bar. The force arrow starts at the displaced bar.
:class: only-light
:::

:::{image} Figure_law_parallel_springs_dark.png
:align: center
:alt: Two springs share a fixed support and a rigid bar. The displacement arrow runs from the dashed initial bar position to the displaced bar. The force arrow starts at the displaced bar.
:class: only-dark
:::

Two springs in parallel carry the applied force $F$ together. The displacement $u$ is measured from the dashed initial position of the common bar to its displaced position.

::::
:::::

:::{math}
:label: eq:constitutive-parallel-springs
F = F_1 + F_2 = K_1 u + K_2 u = (K_1+K_2)u,
\qquad K_{\mathrm{eq}} = K_1+K_2.
:::

The {class}`Parallel <structuralcodes.materials.constitutive_laws.Parallel>` law applies the same principle to uniaxial constitutive laws. It evaluates each component at the same strain $\varepsilon$ and adds the stress and tangent contributions with weights $w_i$.

:::{math}
:label: eq:constitutive-parallel
\sigma(\varepsilon) = \sum_{i=1}^{N} w_i \sigma_i(\varepsilon).
:::

:::{math}
:label: eq:constitutive-parallel-tangent
E_{\mathrm{t}}(\varepsilon) = \sum_{i=1}^{N} w_i E_{\mathrm{t},i}(\varepsilon).
:::

The weights default to $1$, which gives a direct sum of the component stresses and tangents. The implementation uses supplied weights as given. For a composite of perfectly bonded components with cross-sectional areas $A_i$, the total axial force is $\sum_i A_i \sigma_i$. To express that force as stress referred to an area $A_{\mathrm{ref}}$, set $w_i=A_i/A_{\mathrm{ref}}$. The weights sum to $1$ when the reference area is the sum of the component areas.

The combined law reports the smallest lower ultimate strain and the largest upper ultimate strain among its components.

### Two elastic-plastic materials

Consider two [elastic-plastic laws](theory-constitutive-elastic-plastic) at the same tensile strain, with the following properties.

| Material | Elastic modulus $E$ [GPa] | Yield strength $f_y$ [MPa] | Hardening modulus $E_h$ [GPa] | Yield strain $\varepsilon_y$ |
| --- | ---: | ---: | ---: | ---: |
| 1 | 200 | 300 | 5 | 0.0015 |
| 2 | 100 | 350 | 15 | 0.0035 |

Material 1 yields first. Between the two yield strains, the combined tangent is $w_1 E_{h,1}+w_2 E_2$ because material 1 is hardening while material 2 is still elastic. After both materials yield, the tangent is $w_1 E_{h,1}+w_2 E_{h,2}$. The weights change the combined slopes while the strain at which each component yields is set by its own law.

For equal weights $(0.50,0.50)$, the three successive tangent moduli are $150$, $52.5$, and $10\,\mathrm{GPa}$. With weights $(0.75,0.25)$, they are $175$, $28.75$, and $7.5\,\mathrm{GPa}$. Giving more weight to material 1 raises the initial stiffness and reduces the stiffness after yield because this material has the higher elastic modulus and the lower hardening modulus.

(figure-theory-law-parallel-elastic-plastic)=
:::::{grid}
::::{grid-item}
:class: caption-text, sd-text-center

:::{image} Figure_law_parallel_elastic_plastic_light.png
:align: center
:alt: Tensile stress–strain curves for two elastic-plastic materials and their parallel combinations at weights 0.50 and 0.50, and 0.75 and 0.25.
:class: only-light
:::

:::{image} Figure_law_parallel_elastic_plastic_dark.png
:align: center
:alt: Tensile stress–strain curves for two elastic-plastic materials and their parallel combinations at weights 0.50 and 0.50, and 0.75 and 0.25.
:class: only-dark
:::

Two elastic-plastic materials in parallel with weights $(w_1,w_2)=(0.50,0.50)$ and $(0.75,0.25)$. Material 1 uses $E_1=200\,\mathrm{GPa}$, $f_{y,1}=300\,\mathrm{MPa}$, and $E_{h,1}=5\,\mathrm{GPa}$. Material 2 uses $E_2=100\,\mathrm{GPa}$, $f_{y,2}=350\,\mathrm{MPa}$, and $E_{h,2}=15\,\mathrm{GPa}$. Neither material has a specified ultimate strain.

::::
:::::

### Separate tension and compression branches

A specific application of parallel material is to combine a compression law for concrete (e.g. [ParabolaRectangle](theory-constitutive-parabola-rectangle)) with a tensile law for concrete (e.g. [ElasticFragileTension](theory-constitutive-elastic-fragile-tension)). With unit weights, each component supplies the stress on its own side of the origin, giving a single curve with distinct tension and compression branches.

(figure-theory-law-parallel)=
:::::{grid}
::::{grid-item}
:class: caption-text, sd-text-center

:::{image} Figure_law_parallel_light.png
:align: center
:alt: Parabola-rectangle compression and elastic-fragile tension acting in parallel. The inset enlarges the tensile branch.
:class: only-light
:::

:::{image} Figure_law_parallel_dark.png
:align: center
:alt: Parabola-rectangle compression and elastic-fragile tension acting in parallel. The inset enlarges the tensile branch.
:class: only-dark
:::

Parabola-rectangle compression with $f_c=-30\,\mathrm{MPa}$, $\varepsilon_0=-0.002$, $\varepsilon_u=-0.0035$, and $n=2$, combined with elastic-fragile tension with $E_c=30\,\mathrm{GPa}$ and $f_{ct}=3\,\mathrm{MPa}$. Both weights are $1$. The inset enlarges the tensile branch.

::::
:::::

### Longitudinal stiffness of a fiber-reinforced polymer

Consider continuous, aligned fibers bonded to a polymer matrix and loaded along the fiber direction. In the equal-strain model, each phase is represented by a [linear elastic law](theory-constitutive-elastic). For fibers and matrix extending over the same length, their volume fractions equal their area fractions. Using these fractions as weights gives the longitudinal rule of mixtures.[^parallel-rule-of-mixtures]

:::{math}
:label: eq:constitutive-parallel-rule-of-mixtures
\begin{aligned}
w_f &= V_f,\qquad w_m = V_m = 1-V_f, \\
E_L &= V_f E_f + (1-V_f)E_m.
\end{aligned}
:::

:::{math}
:label: eq:constitutive-parallel-longitudinal-stress
\sigma_L(\varepsilon) = E_L\varepsilon.
:::

For illustrative moduli $E_f=230\,\mathrm{GPa}$ and $E_m=3\,\mathrm{GPa}$, the slope of each composite stress–strain curve returned by `Parallel` is its longitudinal modulus $E_L$. The same values follow from the rule of mixtures.

| Fiber volume fraction $V_f$ | Weights $(w_f,w_m)$ | $E_L$ from `Parallel` [GPa] | Rule of mixtures [GPa] |
| ---: | ---: | ---: | ---: |
| 40% | $(0.40,0.60)$ | 93.8 | $0.40(230)+0.60(3)=93.8$ |
| 50% | $(0.50,0.50)$ | 116.5 | $0.50(230)+0.50(3)=116.5$ |

(figure-theory-law-parallel-frp)=
:::::{grid}
::::{grid-item}
:class: caption-text, sd-text-center

:::{image} Figure_law_parallel_frp_light.png
:align: center
:alt: Longitudinal stress–strain curves for fiber, matrix, and composites containing 40 and 50 percent fiber by volume.
:class: only-light
:::

:::{image} Figure_law_parallel_frp_dark.png
:align: center
:alt: Longitudinal stress–strain curves for fiber, matrix, and composites containing 40 and 50 percent fiber by volume.
:class: only-dark
:::

Longitudinal stress–strain curves for the fiber with $E_f=230\,\mathrm{GPa}$, the matrix with $E_m=3\,\mathrm{GPa}$, and the two composites. Weights $(w_f,w_m)=(0.40,0.60)$ and $(0.50,0.50)$ give $E_L=93.8\,\mathrm{GPa}$ and $116.5\,\mathrm{GPa}$, respectively.

::::
:::::

[^sargin-1971]: Sargin, M. (1971). [*Stress-strain relationships for concrete and the analysis of structural concrete sections*](https://books.google.ca/books/about/Stress_strain_relationships_for_concrete.html?id=FN7ZAAAAIAAJ). Study No. 4, Solid Mechanics Division, University of Waterloo, Ontario, Canada.
[^popovics-1973]: Popovics, S. (1973). [A numerical approach to the complete stress-strain curve of concrete](https://doi.org/10.1016/0008-8846%2873%2990096-3). *Cement and Concrete Research*, 3(5), 583–599.
[^mander-1988]: Mander, J. B., Priestley, M. J. N., and Park, R. (1988). [Theoretical stress-strain model for confined concrete](https://doi.org/10.1061/%28ASCE%290733-9445%281988%29114%3A8%281804%29). *Journal of Structural Engineering*, 114(8), 1804–1826.
[^parallel-rule-of-mixtures]: [Experiments to determine elastic moduli for filament-wound cylinders](https://ntrs.nasa.gov/citations/19660001965), NASA Technical Reports Server.

