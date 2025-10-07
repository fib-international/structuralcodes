# Section calculator


## General concept
Section calculators evaluate the structural response based on section properties and applied loads. They integrate results from the geometry and material models to provide insights into section performance by using a specific [*section integrator*](theory-section-integrators).

The section calculator contains the methods that permits to handle a structual analysis for a cross section. For instance it is possible to compute the bending strength given the axial force, or to compute the moment curvature, or to compute the strength domain. The outputs of those methods are represented as *ad-hoc* defined results objects which contains the data and the methods to further process them.

:::{Attention}
Note that in the current implementation both material constitutive laws and section are *unitless*. It means that user should choose a consistent set of units, like N for forces and mm for lengths; in this way the stresses will be in MPa and the moments in Nmm.
Pay attention that default-defined constitutive laws by material classes work with the units defined by the considered code, e.g. MPa for EC2 or MC2010.
:::

## Compute bending strength
With this algorithm, *structuralcodes* computes the bending strength of the section given the axial load (positive in tension and negative in compression) and an angle of the neutral axes respect to the y axis. 

(theory-fig-bending-calc-rotated-system)=
:::{figure} FigureBendingRotated.png

The reference system used for computing bending strength; C indicates the compressed portion, T indicates the stretched portion.
:::


In the rotated reference system **$y^*z^*$**, the bending strength in terms of positive $M_{y^*}$ is computed.

:::{Note}
According to such definition, to compute the bending strength for a section with top fibers in compression and bottom fibers in tension, the angle theta should be equal to $\pi$.

(theory-fig-system-uniaxial-bending)=
:::{figure} FigureBendingTopCompression.png

Rotated Coordinate system to be used for computing uniaxial bending with bottom fibers stretched and top fibers compressed.
:::

:::

According to classic RC theory, the deformations domains could be represented as following (with domains from 1 to 6 moving respectivelly from pure tension to pure compression).

(theory-fig-deformation-domain)=
:::{figure} Figure_DeformationDomain.png

Ultimate deformation domains for uniaxial bending with or without axial force.
:::

The algorithm developed can be summarized as follows:
1. **Rotate the section**: the section is rotated by the given angle theta. In this new CRS *$y^*z^*$*, the problem becomes uniaxial bending about the *$y^*$* axis
2. **Axial load check**: verify if axial load is within the admissible range of axial loads (in tension and compression). This ensures the section can withstand the applied aixal load without failure.
3. **Ultimate strain profile**: find a strain profile that reaches the utimate strain for at least one of the materials. The found strain profile must guarantee equilibrium with external axial load. This is computed with an iterative algorithm based on bisection method.
    a. The internal axial load, defined by the balanced failure condition (i.e., the simultaneous reaching of ultimate strain in both the stretched and compressed materials), is evaluated by integrating the strain profile.
    b. If the internal axial load is greater than the external axial load, the neutral axis needs to be lowered, indicating excessive tension in the section. If the internal axial load is lower, the neutral axis should be raised to reduce compression.
    c. The strain profile is then adjusted by changing the curvature, pivoting on either the top or bottom chord. The goal is to balance the axial load and reach a solution that satisfies equilibrium. This is done solving the equation $\Delta N(k_{y^*})=0$ where $\Delta N(k_{y^*})=N_{ext}-N_{int}(k_{y^*})$. For instance the function could be something like depicted in the following picture:

    Bisection algorithm permits to find the zero of the function within some iterations.

(theory-fig-bisection-dem)=
:::{figure} bisection_demo.gif

Bisection algorithm for finding the strain plain that is in equilibrium with extenal axial force.
:::

4. **Final Computation of Bending Strength**: once the equilibrium strain profile is found, the bending strength is calculated by integrating the obtained strain profile. This results in the final bending strength in the rotated coordinate system, which is then transformed back to the original coordinate system.

## Compute Moment-curvature relation

With this algorithm, *Structuralcodes* computes the moment-curvature of the section given the axial load (positive in tension and negative in compression) and an angle of the neutral axes respect to the y axis. 

(theory-fig-bending-calc-rotated-system-mcurv)=
:::{figure} FigureBendingRotated.png

The reference system used for computing moment-curvature.
:::

The algorithm works with the following steps:

1. **Rotate the section**: the section is rotated by the given angle theta. In this new CRS *$y^*z^*$*, the problem becomes uniaxial bending about the *$y^*$* axis.
2. **Axial load check**: verify if axial load is within the admissible range of axial loads (in tension and compression). This ensures the section can withstand the applied aixal load without failure.
3. **If not provided, compute array of curvatures**: if curvatures are not provided, the array of curvatures is computed. First the yield and ultimate curvatures are computed and the array of curvatures is computed with some points before yielding, and other points from yielding to ultimate condition. Ultimate curvature represents the last value computed automatically. If the user wants to go beyond that poit, he should provide a custom curvatures array. For more information, refert to the [api](api-section-calculator).
4. **For each curvature value compute moment**: for each value of the curvature, the algorithms finds the strain profile that is in equlibrium with external axial force and computes the corresponding value of moment. 
To do so, the algorithm proceeds iterativelly, for a fixed value of curvature, solving the equation $\Delta N(\varepsilon_0)=0$ where $\Delta N(\varepsilon_0)=N_{ext}-N_{int}(\varepsilon_0)$ and $\varepsilon_0$ is the axial strain at coordinates $(0, 0)$ of the section. The iterative solution is performed with a bisection algorithm applied to the range $[\varepsilon_{0,A}, \varepsilon_{0,B}]$ where $\varepsilon_{0,A}$ is the axial strain from the last step and $\varepsilon_{0,B}$ is found by a quick [pre-computation algorithm](theory-pre-find-range) whose aim is finding the range where there is one solution. The bisection algorithm is then applied and the variable $\varepsilon_0$ is detemined within a fixed tolerance.

(theory-pre-find-range)=
### Find the range to which apply bisection
Starting with $\varepsilon_{0,A}$ and the corresponding value of $Delta N(\varepsilon_{0,A})$, the value $\varepsilon_{0,B}'$ is found using a fast pre-find algorithm. This algorithm uses the following exponential law for evaluating the tentative value of $\varepsilon_{0,B,i}$:

$$\varepsilon_{0,B,i} = \varepsilon_{0,B,i-1} + \delta \cdot r^i$$

where $\delta = 1e-3$ and $r = 2$. 

In this way, $\varepsilon_{0,B}$ increases in fast way within very few iterations (generally 2-3 iterations are enough).

The value of $\Delta N(\varepsilon_{0,B,i})$ is computed and whenever $\Delta N(\varepsilon_{0,A}) \cdot \Delta N(\varepsilon_{0,B,i}) < 0$ the algorithm stops having found the range $[\varepsilon_{0,A}, \varepsilon_{0,B}]$ where the bisection algorithm is applied. The use of this algorithm permits to find a good range where bisection algorithm is applied in an efficient way without too many iterations.

Graphically this can be seen in the figure [below](#theory-fig-dN_eps0), where the function $\Delta N(\varepsilon_0)$ is plotted indicating the value $\varepsilon_{0,A}$ from the previous step and the two iterations needed to find $\varepsilon_{0,B}$. 

(theory-fig-dN_eps0)=
:::{figure} FiguredN_eps_0.png

The function $\Delta N(\varepsilon_0)$ during $i$-th step of moment-curvatura computation: application of the pre-find algorithm to select the range where bisection is applied.
:::

Then, the usual bisection algorithm is adopted to find the value $\varepsilon_0$ for which $\Delta N(\varepsilon_0) = N_{ext} - N_{int}(\varepsilon_0)$ is sufficiently near to 0.

