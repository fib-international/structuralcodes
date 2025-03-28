# Section calculator


## General concept
Section calculators evaluate the structural response based on section properties and applied loads. They integrate results from the geometry and material models to provide insights into section performance by using a specific *section integrator*.

The section calculator contains the methods that permits to handle a structual analysis for a cross section. For instance it is possible to compute the bending strength given the axial force, or to compute the moment curvatura, or to compute the strength domain. The outputs of those methods are represented as *ad-hoc* defined results objects which contains the data and the methods to further process them.

:::{Attention}
Note that in the current implementation both material constitutive laws and section are *unitless*. It means that user should choose a consistent set of units, like N for forces and mm for lengths; in this way the stresses will be in MPa and the moments in Nmm.
Pay attention that default-defined constitutive laws by material classes work with the units defined by the considered code, e.g. MPa for EC2 or MC2010.
:::

## Compute bending strength
With this algorithm, *structuralcodes* computes the bending strength of the section given the axial load (positive in tension and negative in compression) and an angle of the neutral axes respect to the y axis. 

![Bending Strength Rotated](FigureBendingRotated.png)

In the rotated reference system **y\*z\***, the bending strength in terms of positive {math}`M_{\textrm{y}^*}` is computed.

:::{Note}
According to such definition, to compute the bending strength for a section with top fibers in compression and bottom fibers in tension, the angle theta should be equal to {math}`\pi`.

![Bending Strength Bottom Tension](FigureBendingTopCompression.png)
:::

According to classic RC theory, the deformations domains could be represented as following (with domains from 1 to 6 moving respectivelly from pure tension to pure compression):

![Deformation Domain](Figure_DeformationDomain.jpg)

The algorithm developed can be summarized as follows:
1. **Rotate the section**: the section is rotated by the given angle theta. In this new CRS (**y\*z\***), the problem becomes uniaxial bending about the **y\*** axis
2. **Axial load check**: verify if axial load is within the admissible range of axial loads (in tension and compression). This ensures the section can withstand the applied aixal load without failure.
3. **Ultimate strain profile**: find a strain profile that reaches the utimate strain for at least one of the materials. The found strain profile must guarantee equilibrium with external axial load. This is computed with an iterative algorithm based on bisection method.
    a. The internal axial load, defined by the balanced failure condition (i.e., the simultaneous reaching of ultimate strain in both the stretched and compressed materials), is evaluated by integrating the strain profile.
    b. If the internal axial load is greater than the external axial load, the neutral axis needs to be lowered, indicating excessive tension in the section. If the internal axial load is lower, the neutral axis should be raised to reduce compression.
    c. The strain profile is then adjusted by changing the curvature, pivoting on either the top or bottom chord. The goal is to balance the axial load and reach a solution that satisfies equilibrium. This is done solving the equation {math}`\Delta N(k_{y^*})=0` where {math}`\Delta N(k_{y^*})=N_{ext}-N_{int}(k_{y^*})`. For instance the function could be something like depicted in the following picture:

    ![dN vs chi](FiguredN_chi.png)

    Bisection algorithm permits to find the zero of the function within some iterations.

    ![Bisection demo](bisection_demo.gif)

4. **Final Computation of Bending Strength**: once the equilibrium strain profile is found, the bending strength is calculated by integrating the obtained strain profile. This results in the final bending strength in the rotated coordinate system, which is then transformed back to the original coordinate system.

## Compute Moment-curvature relation

With this algorithm, *structuralcodes* computes the moment-curvature of the section given the axial load (positive in tension and negative in compression) and an angle of the neutral axes respect to the y axis. 

![Bending Strength Rotated](FigureBendingRotated.png)

The algorithm works with the following steps:

1. **If not provided, compute array of curvatures**: if curvatures are not provided, the array of curvatures corresponding to which computing the moments is computed. First the yield and 