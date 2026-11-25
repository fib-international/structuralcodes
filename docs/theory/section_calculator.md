(theory-section-calculator)=
# Section calculator

## General concept
Section calculators evaluate the structural response based on section properties and applied loads. They integrate results from the geometry and material models to provide insights into section performance by using a specific [*section integrator*](theory-section-integrators).

The section calculator contains the methods that permits to handle a structual analysis for a cross section. For instance it is possible to compute the bending strength given the axial force, or to compute the moment curvature, or to compute the strength domain. The outputs of those methods are represented as *ad-hoc* defined results objects which contains the data and the methods to further process them.

:::{Attention}
Note that in the current implementation both material constitutive laws and section are *unitless*. It means that user should choose a consistent set of units, like N for forces and mm for lengths; in this way the stresses will be in MPa and the moments in Nmm.
Pay attention that default-defined constitutive laws by material classes work with the units defined by the considered code, e.g. MPa for EC2 or MC2010.
:::

## Neutral axis orientation (`theta`)

As better detailed below, most of the methods require the angle $\theta$ (`theta` input argument) for defining the neutral axis inclination $\theta$.

The parameter $\theta$ represents the angle between the original axis $y$ of the section reference system $(y,z)$ and the rotated axis $y^*$ of the local reference system $(y^*,z^*)$ used internally to calculate the section response (see figure [below](theory-fig-bending-calc-rotated-system)).

The angle `theta` is assumed positive if counter-clockwise (CCW), negative if clockwise (CW).

Therefore:
- `theta` = 0 → axis $y^*$ coincides with axis $y$.
- `theta` > 0 → the rotated system $(y^*, z^*)$ is obtained by a CCW rotation.
- `theta` < 0 → the rotated system $(y^*, z^*)$ is obtained by a CW rotation.

The neutral axis is orthogonal to the $z^*$ direction. Hence, theta directly controls the bending direction and allows the treatment of:
- uniaxial bending about principal or non-principal axes,
- biaxial bending,

Internally, once `theta` is defined, coordinates are transformed according to the standard rotation:

$$
\begin{aligned}
y^* &= y \cos\theta + z \sin\theta \\
z^* &= -y \sin\theta + z \cos\theta
\end{aligned}
$$

In the rotated system $(y^*, z^*)$ the bending is unixial, therefore the strain field is expressed in the rotated system as:

$$
\varepsilon(y^*) = \varepsilon_0 + \kappa^* \, z^*  
$$

where:

- $\varepsilon_0$ is the axial strain,
- $\kappa^*$ is the curvature,
- $z^*$ is the coordinate measured in the direction perpendicular to the neutral axis.

This formulation provides a fully general description of section behaviour under combined axial force and biaxial bending.

(theory-compute-bending-strength)=
## Compute bending strength
With this algorithm, StructuralCodes computes the bending strength of the section given the axial load (positive in tension and negative in compression) and an angle of the neutral axes respect to the y axis. 

(theory-fig-bending-calc-rotated-system)=
:::::{grid}
::::{grid-item}
:class: caption-text, sd-text-center

:::{image} FigureBendingRotated_light.png
:align: center
:class: only-light
:::

:::{image} FigureBendingRotated_dark.png
:align: center
:class: only-dark
:::

The reference system used for computing bending strength; C indicates the compressed portion, T indicates the stretched portion.

::::
:::::

In the rotated reference system **$y^*z^*$**, the bending strength in terms of ***negative*** moment $M_{y^*}$ is computed. This means that the moment calculated is the one stretching bottom fiber in the rotated system **$y^*z^*$**.

After the bending strength is computed in the rotated reference system **$y^*z^*$**, the result is rotated back to the original coordinate system **$yz$** with moments signs following our sign convention described [here](theory-crs-sign-convention).

::::::{Note}
According to such definition, to compute the bending strength for a section with top fibers in tension and bottom fibers in compression, the angle theta should be equal to $\pi$.

(theory-fig-system-uniaxial-bending)=
:::::{grid}
::::{grid-item}
:class: caption-text, sd-text-center

:::{image} FigureBendingTopCompression_light.png
:align: center
:class: only-light
:::

:::{image} FigureBendingTopCompression_dark.png
:align: center
:class: only-dark
:::

Rotated Coordinate system to be used for computing uniaxial bending with bottom fibers compressed and top fibers stretched.

::::
:::::

::::::

According to classic RC theory[^leonhardt1973], the deformations domains could be represented as following (with domains from 1 to 6 moving respectivelly from pure tension to pure compression).

(theory-fig-deformation-domain)=
:::::{grid}
::::{grid-item}
:class: caption-text, sd-text-center

:::{image} FigureDeformationDomain_light.png
:align: center
:class: only-light
:::

:::{image} FigureDeformationDomain_dark.png
:align: center
:class: only-dark
:::

Ultimate deformation domains for uniaxial bending with or without axial force [^leonhardt1973].

::::
:::::

The algorithm developed can be summarized as follows:
1. **Rotate the section**: the section is rotated by the given angle theta. In this new CRS *$y^*z^*$*, the problem becomes uniaxial bending about the *$y^*$* axis
2. **Axial load check**: verify if axial load is within the admissible range of axial loads (in tension and compression). This ensures the section can withstand the applied aixal load without failure.
3. **Ultimate strain profile**: find a strain profile that reaches the utimate strain for at least one of the materials. The found strain profile must guarantee equilibrium with external axial load. This is computed with an iterative algorithm based on bisection method.

    a. The internal axial load, defined by the balanced failure condition (i.e., the simultaneous reaching of ultimate strain in both the stretched and compressed materials), is evaluated by integrating the strain profile.

    b. If the internal axial load is greater than the external axial load, the neutral axis needs to be lowered, indicating excessive tension in the section. If the internal axial load is lower, the neutral axis should be raised to reduce compression.

    c. The strain profile is then adjusted by changing the curvature, pivoting on either the top or bottom chord. The goal is to balance the axial load and reach a solution that satisfies equilibrium. This is done solving the equation $\Delta N(\chi_{y^*})=0$ where $\Delta N(\chi_{y^*})=N_{ext}-N_{int}(\chi_{y^*})$. For instance the function could be something like depicted in the following picture.

    (theory-fig-dN_chi)=
    :::::{grid}
    ::::{grid-item}
    :class: caption-text, sd-text-center

    :::{image} FiguredN_chi_light.png
    :align: center
    :class: only-light
    :::

    :::{image} FiguredN_chi_dark.png
    :align: center
    :class: only-dark
    :::

    The function $\Delta N(\chi_y^*)$ the searched curvature is the intersection with the horizontal axis (i.e. when $\Delta N(\chi_y^*) = 0$).

    ::::
    :::::

    Bisection algorithm permits to find the zero of the function within some iterations.

(theory-fig-bisection-dem)=
:::::{grid}
::::{grid-item}
:class: caption-text, sd-text-center

:::{image} bisection_demo_light.gif
:align: center
:class: only-light
:::

:::{image} bisection_demo_dark.gif
:align: center
:class: only-dark
:::

Bisection algorithm for finding the strain plain that is in equilibrium with extenal axial force.

::::
:::::


4. **Final Computation of Bending Strength**: once the equilibrium strain profile is found, the bending strength is calculated by integrating the obtained strain profile. This results in the final bending strength in the rotated coordinate system described in step 1, which is then transformed back to the original coordinate system.

## Compute Moment-curvature relation

With this algorithm, StructuralCodes computes the moment-curvature of the section given the axial load (positive in tension and negative in compression) and an angle of the neutral axes respect to the y axis (see figure [above](theory-fig-bending-calc-rotated-system)). 

The algorithm works with the following steps:

1. **Rotate the section**: the section is rotated by the given angle theta. In this new CRS *$y^*z^*$*, the problem becomes uniaxial bending about the *$y^*$* axis.
2. **Axial load check**: verify if axial load is within the admissible range of axial loads (in tension and compression). This ensures the section can withstand the applied aixal load without failure.
3. **If not provided, compute array of curvatures**: if curvatures are not provided, the array of curvatures is computed. First the yield and ultimate curvatures are computed and the array of curvatures is computed with some points before yielding, and other points from yielding to ultimate condition. Ultimate curvature represents the last value computed automatically. If the user wants to go beyond that point, he should provide a custom curvatures array. For more information, refert to the [api](#api-section-calculator).
4. **For each curvature value compute moment**: for each value of the curvature, the algorithm finds the strain profile that is in equilibrium with the external axial force and computes the corresponding value of moment.
To do so, for a fixed value of curvature $\chi$ the strain is $\varepsilon(z^*) = \varepsilon_0 + \chi\, z^*$, where $\varepsilon_0$ is the axial strain at coordinates $(0, 0)$ of the section, and the equation $\Delta N(\varepsilon_0)=0$ is solved, where $\Delta N(\varepsilon_0)=N_{int}(\varepsilon_0)-N_{ext}$. The solution is searched only in the range of $\varepsilon_0$ where all materials are within their ultimate strains: first a bracket $[\varepsilon_{0,l}, \varepsilon_{0,r}]$ containing the zero closest to the solution of the previous step is found, then the zero is computed with Brent's method (see [below](theory-find-equilibrium-fixed-curvature)).
If for a value of curvature no equilibrium can be found (e.g. the curvature is beyond the ultimate one or, for constitutive laws with softening and high compression, beyond a limit point), a `NoConvergenceWarning` is raised and the analysis is stopped at the last converged step.

::::::{note}
(theory-find-equilibrium-fixed-curvature)=
### Find the equilibrium for a fixed curvature

**Admissible range of $\varepsilon_0$.** For a fixed curvature the strain is linear over the section, so for each geometry $g$ (surface or point) with ultimate strains $\varepsilon_{u,g}^-$ and $\varepsilon_{u,g}^+$ it is sufficient to check its extreme coordinates $z^*_{g,min}$ and $z^*_{g,max}$. The admissible range is the intersection of the ranges of all the geometries:

$$\varepsilon_{0,min} = \max_{g} \max_{z^* \in \{z^*_{g,min},\, z^*_{g,max}\}} \left(\varepsilon_{u,g}^- - \chi\, z^*\right) \qquad \varepsilon_{0,max} = \min_{g} \min_{z^* \in \{z^*_{g,min},\, z^*_{g,max}\}} \left(\varepsilon_{u,g}^+ - \chi\, z^*\right)$$

Outside this range at least one material has failed. Restricting the search to $[\varepsilon_{0,min}, \varepsilon_{0,max}]$ excludes, for instance, the spurious solutions found when the reinforcement is broken and the section is completely in tension ($\Delta N = 0$ with zero moment), and the discontinuities of $\Delta N(\varepsilon_0)$ due to the rupture of a material, which are located at the bounds of the range. If $\varepsilon_{0,min} > \varepsilon_{0,max}$ the curvature exceeds the ultimate one and no equilibrium exists.

The figure [below](#theory-fig-admissible-range) shows the admissible range for a step of the plastic branch of the moment-curvature of a rectangular section ($N = 0$). For a fixed curvature all the strain profiles are parallel lines (b): moving the profile, i.e. changing $\varepsilon_0$, the first limit reached on the compressed side is the ultimate strain of concrete at the top edge ($\varepsilon_{0,min}$), and on the tension side the ultimate strain of the reinforcement ($\varepsilon_{0,max}$). Each extreme fibre allows an interval of $\varepsilon_0$ (c), and the admissible range is the intersection of these intervals.

(theory-fig-admissible-range)=
:::::{grid}
::::{grid-item}
:class: caption-text, sd-text-center

:::{image} FigureAdmissibleRange_light.png
:align: center
:class: only-light
:::

:::{image} FigureAdmissibleRange_dark.png
:align: center
:class: only-dark
:::

Admissible range of $\varepsilon_0$ for a fixed curvature. (a) Section. (b) Strain profiles with the same curvature: the bands are the strains allowed in each extreme fibre (the concrete has no limit in tension), the shaded area contains all the admissible profiles. (c) Interval of $\varepsilon_0$ allowed by each fibre and their intersection, which contains the solution of the step.

::::
:::::

**Starting point and bounds.** The solution of the previous step $\varepsilon_{0,A}$, clipped to the admissible range, is the starting point. $\Delta N$ is computed at $\varepsilon_{0,A}$ and at the two bounds: if one of these values is already in equilibrium (within the tolerance on $|\Delta N|$) it is the solution. This is the case of the last step of a moment-curvature computed automatically: at the ultimate curvature one material is exactly at its ultimate strain, so the solution lies on a bound of the admissible range.

**Bracket search.** Starting from $\varepsilon_{0,A}$ the range is explored in both directions, alternating the sides, until a change of sign of $\Delta N$ is found. The steps start from $\delta$ and grow exponentially:

$$\varepsilon_{0,i} = \varepsilon_{0,i-1} \pm \min\left(\delta \cdot 2^{i-1},\ \Delta_{max}\right) \qquad \Delta_{max} = \max\left(\frac{\varepsilon_{0,max}-\varepsilon_{0,min}}{100},\ \delta\right)$$

where $\delta$ is the estimate of the change of $\varepsilon_0$ obtained from the two previous steps, scaled with the increment of curvature:

$$\delta = \left| \left(\varepsilon_{0,j-1} - \varepsilon_{0,j-2}\right) \frac{\chi_j - \chi_{j-1}}{\chi_{j-1} - \chi_{j-2}} \right|$$

($\delta = 10^{-4}$ for the first two steps). In this way:

- the bracket is found with few evaluations, since the solution is generally close to the one of the previous step;
- the zero closest to $\varepsilon_{0,A}$ is found, on whichever side it is. This can be relevant for constitutive laws with softening (e.g. Sargin and Popovics) and high compression, where $\Delta N(\varepsilon_0)$ could be not monotic having more than one zero in the admissible range: following the zero closest to the previous solution keeps the analysis on the same branch;

If no change of sign is found in the whole admissible range, no equilibrium exists for the given curvature and a `NoConvergenceWarning` is raised.

The figure [below](#theory-fig-find-bracket) shows the bracket search at the third step of the moment-curvature of a rectangular section 250 x 500 mm ($f_{ck} = 45$ MPa, two bars of diameter 16 mm at the bottom, $\varepsilon_{uk} = 4\%$) with $N = -1000$ kN. Here the estimate $\delta$ is much smaller than the actual move of the solution, so the evaluations alternate on the two sides of $\varepsilon_{0,A}$ (1, 2, 3, ...) with steps $\delta$, $2\delta$, $4\delta$, ... until the first change of sign is found, between the points 7 and 9 on the right side. Brent's method then finds the zero inside this bracket with a few evaluations.

(theory-fig-find-bracket)=
:::::{grid}
::::{grid-item}
:class: caption-text, sd-text-center

:::{image} FigureFindBracket_light.png
:align: center
:class: only-light
:::

:::{image} FigureFindBracket_dark.png
:align: center
:class: only-dark
:::

Bracket search for a fixed curvature. Left: $\Delta N(\varepsilon_0)$ in the whole admissible range. Right: zoom on the search; the points are numbered in the order of evaluation (orange: increasing $\varepsilon_0$, green: decreasing $\varepsilon_0$) and the arrows show the step from the previous point on the same side.

::::
:::::

**Brent's method.** Once the bracket is found, the zero is computed with Brent's method (`scipy.optimize.brentq`). At each iteration it keeps a bracket of the solution and tries a step with inverse quadratic interpolation (or with the secant method), falling back to a bisection step when the interpolation is not reliable. For the piecewise smooth function $\Delta N(\varepsilon_0)$ it converges in a few iterations, never losing the bracket. If the bracket collapses on a point where $\Delta N$ jumps across zero without an equilibrium (e.g. a brittle material that suddenly stops carrying stress), the residual $|\Delta N|$ remains larger than the tolerance and a `NoConvergenceWarning` is raised.

Graphically this can be seen in the figure [below](#theory-fig-dN_eps0), where the function $\Delta N(\varepsilon_0)$ is plotted for a step of the plastic branch of a moment-curvature. In this case $\varepsilon_{0,A}$ is outside the admissible range (with the new curvature the concrete would exceed its ultimate strain) and the search starts from $\varepsilon_{0,min}$.

(theory-fig-dN_eps0)=
:::::{grid}
::::{grid-item}
:class: caption-text, sd-text-center

:::{image} FiguredN_eps_0_light.png
:align: center
:class: only-light
:::

:::{image} FiguredN_eps_0_dark.png
:align: center
:class: only-dark
:::

The function $\Delta N(\varepsilon_0)$ during the $i$-th step of moment-curvature computation. Left: admissible range of $\varepsilon_0$ and values at the bounds. Right: the previous solution $\varepsilon_{0,A}$ is outside the admissible range, so the search starts from $\varepsilon_{0,min}$; the first point of the bracket search already gives a change of sign, and Brent's method finds the solution inside the bracket.

::::
:::::

::::::

(theory-nm-domain)=
## Compute NM interaction domain

With this algorithm, StructuralCodes computes the two dimensional $N$-$M_y^*$ representing uniaxial bending respect to the rotated reference system *$y^*z^*$*.

The algorithm works with the following steps:

1. **Rotate the section**: the section is rotated by the given angle theta. In this new CRS *$y^*z^*$*, the problem becomes uniaxial bending about the *$y^*$* axis.
2. **Compute the list of ultimate strain profiles**: according to classic RC theory, the failure domain can be grouped into 6 differen domain, from 1 to 6, moving respectivelly from pure tension to pure compression. (see figure [here](theory-fig-deformation-domain)). The algorithm created a list of all possible failure strain profiles within the six domains. User can customatize the discretization used (see [api](#api-section-calculator)).
3. **Integrate stresses for all computed ultimate strain profiles**: for all the strain profiles determined in step 2, the stresses are integrated computing the pair $(N, M_y^*)$ that represent a point in the ultimate domain profile.
This is represented by the following image that shows for a rectangular RC section the ultimate strain profiles computed and the corresponding MN interaction domain. The different colors correspond to different domains. The same is visible in the animation [here](theory-anim-domain_nm)


(theory-fig-domain_nm)=
:::::{grid}
::::{grid-item}
:class: caption-text, sd-text-center

:::{image} Figure_domain_mn_light.png
:align: center
:class: only-light
:::

:::{image} Figure_domain_mn_dark.png
:align: center
:class: only-dark
:::

Ultimate strain profiles for uniaxial bending and corresponding interaction domain.

::::
:::::


(theory-anim-domain_nm)=
:::::{grid}
::::{grid-item}
:class: caption-text, sd-text-center

:::{image} Figure_domain_mn_light.gif
:align: center
:class: only-light
:::

:::{image} Figure_domain_mn_dark.gif
:align: center
:class: only-dark
:::

Animation of ultimate strain profiles for uniaxial bending and corresponding interaction domain.

::::
:::::

4. **Rotate back to original $yz$ coordinate system**: the pair of values $(N, M_y^*)$ is rotated back obtaining the triplet $(N, M_y, M_z)$

(theory-nmm-domain)=
## Compute NMM interaction domain

With this algorithm, StructuralCodes computes the full three dimensional $N$-$M_y$-$M_z$ interaction domain.

To do so, the algorithm loops different values of angle $\theta_i$ of neutral axes respect axis $y$, going from $0$ to $2\pi$ (the number of discretizations can be customized with the parameter `num_theta` whose default value is equal to 32).
For each value of the angle $\theta_i$, the algorithm proceeds as the case of unixial bending ([here](theory-nm-domain)).
The final result is represented by the data structure collecting for all points the triplet $(N, M_y, M_z)$ and the corresponding ultimate strain profile.

(theory-figure-nmm)=
:::::{grid}
::::{grid-item}
:class: caption-text, sd-text-center

:::{image} Figure_domain_nmm_light.png
:align: center
:class: only-light
:::

:::{image} Figure_domain_nmm_dark.png
:align: center
:class: only-dark
:::

Left: ultimate domain NMM. Right: corresponding points in terms of N, $\chi_y$, $\chi_z$.

::::
:::::

## Compute MM interaction domain

With this algorithm, StructuralCodes computes the two dimensional $M_y$, $M_z$ interaction domain for a given value of external axial load  $N$.

The algorithm works with the following steps:
1. **Create an array of angles**: the algorithm creates an array of linearly spaced values of angle $\theta$ of neutral axis respect to axis $y$. 
2. **Compute bending strength**: for each value of the angle $\theta_i$, the bending strength (with the iterative algorithm described [here](theory-compute-bending-strength)) for the given value of external axial load $N$ is computed obtaining a pair $(M_{y,i}, M_{z,i})$ that represent the point in the failure domain. Doing this for all values, the entire $M_yM_z$ failure domain is created.

:::{note}
In the previous cases (i.e. [NM domain](theory-nm-domain) and [NMM domain](theory-nmm-domain)) we simply integrate a bunch of strain profiles obtaining a set of points in the interaction domain. This is a very fast and no iterative solutions are needed. A side effect of this choice is that the obtained domain points are not controlled in terms of a constant external axial forces, see e.g. [this figure](theory-figure-nmm). 

For this reason the computation of MM interaction domain for a given value of axial load $N$ does not work in the same way. The side effect of this algorithm is that it is significantly slower since each point of the failure domain is determined with an iterative process.
::: 

::::::{hint}
Using this method it is possible to loop from minimum to maximum axial load that the section is able to carry and construct a nice structured version of the failure domain. Of course this operation is pretty slow compared to the direct computation of [NMM domain](theory-nmm-domain).

::::{dropdown-syntax}
:::{literalinclude} ../_example_code/usage_compute_3d_domain_from_mm.py
:lines: 3-
:caption: Create many slices of MM interaction domains to have a full NMM domain.
:::
::::

The obtain figure is shown in the figure [below](#theory-figure-many-mm).

(theory-figure-many-mm)=
:::::{grid}
::::{grid-item}
:class: caption-text, sd-text-center

:::{image} Figure-many-mm-3d_light.png
:align: center
:class: only-light
:::

:::{image} Figure-many-mm-3d_dark.png
:align: center
:class: only-dark
:::

Ultimate domain NMM obtained as a set of slices computed as MM interaction domains.

::::
:::::



::::::

## Compute strain profile for given external forces

With this algorithm, StructuralCodes  computes the strain profile (represented by the triplet of values $\varepsilon_0$, $\chi_y$, $\chi_z$ that are respectively the axial strain corresponding to point $(0,0)$, and the curvatures with respect to axes $y$ and $z$) for a given set of external forces $N_{\mathrm{ext}}$, $M_{\mathrm{y,ext}}$, $M_{\mathrm{z,ext}}$, with an iterative solution algorithm based on Newton-Rhapson method.

The vector containing the internal force $N$, $M_{\mathrm{y}}$, $M_{\mathrm{z}}$ is obtained by integrating the stresses over the cross-section area. These quantities can be expressed as follows:

:::{math}
\begin{aligned}
N  &= \int_A \sigma(y,z) \, dA \\
M_y &= \int_A \sigma(y,z) \, z \, dA \\
M_z &= -\int_A \sigma(y,z) \, y \, dA
\end{aligned}
:::

or, in vector form, as:

:::{math}
:label: eq:section_forces
\begin{bmatrix}
N \\
M_y \\
M_z
\end{bmatrix}
=
\int_A
\begin{bmatrix}
1 \\
z \\
-y
\end{bmatrix}
\sigma(y,z)\, dA
:::

The strain $\varepsilon(y,z)$ can be expressed as a function of the axial deformation $\varepsilon_0$ and the curvatures $\chi_y$ and $\chi_z$ as:

:::{math}
:label: eq:strains_determination
\varepsilon(y,z) = \varepsilon_0 + \chi_y \cdot z - \chi_z \cdot y
:::

or, in vector form, as:

:::{math}
:label: eq:strains_determination_vector
\varepsilon(y,z) =
\begin{bmatrix} 1 & z & -y \end{bmatrix} \cdot 
\begin{bmatrix}
\varepsilon_0 \\
\chi_y \\
\chi_z
\end{bmatrix}
:::

The row vector 
$\mathbf{B}(y,z) = \begin{bmatrix} 1 & z & -y \end{bmatrix}$
is the **strain-displacement vector** that maps the generalized strain vector 
$\boldsymbol{\varepsilon}_g = \begin{bmatrix} \varepsilon_0 & \chi_y & \chi_z \end{bmatrix}^T$
to the pointwise axial strain according to {eq}`eq:strains_determination_vector`.

:::{hint}
The sign appearing in equation {eq}`eq:strains_determination` are explained considering the sign convention shown [here](theory-crs-sign-convention). 

According to such convention a positive $\chi_y$ stretches top fibers and compresses bottom fibers (i.e. hogging), so when multiplying positive $z$ coordinates (points towars top edge) the final strain is positive. A negative $\chi_y$ stretches instead bottom fibers (i.e. sagging), so when multiplying positive $z$ coordinates (points towards top edge) the final strain is negative. 

A similar consideration can be made for the other curvature $\chi_z$
:::

The stress $\sigma(y,z)$ can be related to the strain $\varepsilon(y,z)$ through the consitutive law (in general non-linear), according to:

:::{math}
:label: eq:stress-strain
\sigma(y,z) = f[\varepsilon(y,z)] = f[\varepsilon_0 + \chi_y \cdot z - \chi_z \cdot y]
:::

or, in vector form, as:

:::{math}
:label: eq:stress-strain_vector
\sigma(y,z)
=
f\left(
\begin{bmatrix} 1 & z & -y \end{bmatrix} \cdot 
\begin{bmatrix}
\varepsilon_0 \\
\chi_y \\
\chi_z
\end{bmatrix} 
\right)
:::

Substituting {eq}`eq:stress-strain` into {eq}`eq:section_forces`, the internal forces vector can be written as:

:::{math}
:label: eq:forces-strains
\begin{bmatrix}
N \\
M_y \\
M_z
\end{bmatrix}
=
\int_A
\begin{bmatrix}
1 \\
z \\
-y
\end{bmatrix}
f[\varepsilon_0 + \chi_y \cdot z - \chi_z \cdot y]
\, dA
:::

or, in vector form, as:

:::{math}
:label: eq:forces-strains-vector
\begin{bmatrix}
N \\
M_y \\
M_z
\end{bmatrix}
=
\int_A
\begin{bmatrix}
1 \\
z \\
-y
\end{bmatrix}
f\left(
\begin{bmatrix} 1 & z & -y \end{bmatrix} \cdot 
\begin{bmatrix}
\varepsilon_0 \\
\chi_y \\
\chi_z
\end{bmatrix} 
\right)
\, dA
:::

In compact matrix notation, {eq}`eq:forces-strains-vector` can be written as:

:::{math}
:label: eq:forces-strains-matrix
\mathbf{f}^{int}
=
\int_A
\mathbf{B}^T
\cdot
f\left(
\mathbf{B} \cdot \boldsymbol{\varepsilon}_g
\right)
:::

The vector $\boldsymbol{\varepsilon}_g$ represent the unknown strain profile that we want to find for the given external forces $\mathbf{f}^{ext}$. To do so, a numeric method like Newton-Raphson can be used.

To this aim, the residual vector can be written as:

:::{math}
:label: eq:residuals
\begin{aligned}
\mathbf{R}(\boldsymbol{\varepsilon}_g) &= \mathbf{f}^{ext} - \mathbf{f}^{int} \\
  &= \mathbf{f}^{ext} - \int_A
\mathbf{B}^T
\cdot
f\left(
\mathbf{B} \cdot \boldsymbol{\varepsilon}_g
\right)
\end{aligned} 
:::

and linearizing it:

:::{math}
:label: eq:residuals_linearized
\mathbf{R}(\boldsymbol{\varepsilon}_g^{k+1}) = \mathbf{R}(\boldsymbol{\varepsilon}_g^k) + \mathbf{J}(\boldsymbol{\varepsilon}_g^k) \Delta \boldsymbol{\varepsilon}_g
:::

It is possible then to solve the current step:

:::{math}
:label: eq:nr_solve_current_step
\Delta \boldsymbol{\varepsilon}_g
= - \mathbf{J}^{-1}(\boldsymbol{\varepsilon}_g^k) 
\cdot
\mathbf{R}(\boldsymbol{\varepsilon}_g^k)
:::

and therefore obtaining the new iteration of $\boldsymbol{\varepsilon}_g$:

:::{math}
:label: eq:nr_new_step
\boldsymbol{\varepsilon}_g^{k+1}
=
\boldsymbol{\varepsilon}_g^k +
\Delta \boldsymbol{\varepsilon}_g
:::

The iterative process can continue up to a certain convergence criterion is met, like for instance up to when the norm of the residual vector $\|\mathbf{R}(\boldsymbol{\varepsilon}_g^{k+1})\|$ is smaller then a selected tolerance.

The Jacobian matrix $\mathbf{J}$ is basically the tanget stiffness of the section and is a 3x3 matrix given by:

:::{math}
:label: eq:jacobian
\mathbf{J} = \int_A
\dfrac{d \sigma}{d \boldsymbol{\varepsilon}_g}
\cdot
\begin{bmatrix}
1 \\
z \\
-y
\end{bmatrix}
\, dA
=
\int_A
\dfrac{d \sigma}{d \boldsymbol{\varepsilon}_g}
\cdot
\mathbf{B}^T
\, dA
:::

where:

:::{math}
:label: eq:derivative_sigma_epsilon_g
\dfrac{d \sigma}{d \boldsymbol{\varepsilon}_g} =
\dfrac{d \sigma}{d \varepsilon} \cdot
\dfrac{d \varepsilon}{d \boldsymbol{\varepsilon}_g}
:::

Therefore the Jacobian matrix can be written as:

:::{math}
:label: eq:jacobian_full
\begin{aligned}
\mathbf{J} 
&=
\int_A
\dfrac{d \sigma}{d \varepsilon} 
\begin{bmatrix}
1 \\
z \\
-y
\end{bmatrix} \cdot
\begin{bmatrix} 1 && z && -y \end{bmatrix}
\, dA =
\\
&=
\int_A
\dfrac{d \sigma}{d \varepsilon}
\begin{bmatrix}
1 & z & -y \\
z & z^2 & -y z \\
-y & -y z & y^2
\end{bmatrix}
\, dA
\end{aligned}
:::

where $\dfrac{d \sigma}{d \varepsilon}$ represents the tangent modulus of the constitutive law.

::::{hint}
For a linear elastic material, $\dfrac{d \sigma}{d \varepsilon}$ is constant and equal to $E$, therefore the tangent stiffness matrix of the section can be written as the usual stiffness matrix:

:::{math}
\mathbf{K} = E
\begin{bmatrix}
A & S_z & S_y \\
S_z & I_{yy} & I_{yz} \\
S_y & I_{yz} & I_{zz}
\end{bmatrix}
:::

where:

| Symbol | Definition |
|--------|------------|
| $A$      | Cross-sectional area: $A = \int_A dA $ |
| $S_y$    | First moment about z-axis: $S_y = \int_A y \, dA$ |
| $S_z$    | First moment about y-axis: $S_z = \int_A z \, dA$ |
| $I_{yy}$ | Moment of inertia about y-axis: $I_{yy} = \int_A z^2 \, dA$ |
| $I_{zz}$ | Moment of inertia about z-axis: $I_{zz} = \int_A y^2 \, dA$ |
| $I_{yz}$ | Product of inertia: $I_{yz} = \int_A y z \, dA$ |
| $E$      | Young's modulus (constant for linear elastic material) |

::::

The algorithm implemented in StructuralCodes is basically as above described, with the possibility of using also the modified Newton-Rhapson with initial tangent (setting `initial = True`). Further the user can set the maximum number of iterations (set to 10 by default) and the tolerance (set to 1e-6 by default).

[^leonhardt1973]: Leonhardt, F., Mönning, E. Vorlesungen über Massivbau, Springer, Verlag, Berlin/Heidelberg, 1973.