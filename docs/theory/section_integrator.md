(theory-section-integrators)=
# Section integrators
The algorithms described in [Section calculators](theory-section-calculator) require integration of stress or stiffness over the cross-section. This operation is performed using proper `SectionIntegrator` objects. 

In *structuralcodes* we have currently two distinct integrators: *fiber integrator* that discretize the cross section in a bunch of fibers (whose behavior is uniaxial), and *marin integrator* that is based on the computation of moments of area.

## General concept
Section integrators compute properties such as stiffness, stress distribution, or strain compatibility over a section. They operate on the geometry and material definitions to derive these results.

(theory-marin-integrator)=
## Marin integrator
Marin integrator implements the idea of prof. Joaquín Marín from University of Venezuala that he published during the '80s [^marin1984].

The method performs an exact numerical solution of double integrals for polynomial functions acting over any plane polygonal surface (with or without holes) and is based on the computation of moments of area.

Within this context the moment of area of order $m, n$ is described by the basic double integral:

:::{math}
:label: eq:marin_moment_area
\mathcal{M}_{m,n} = \int_A y^m z^n \, dA
:::

where it is implicitly understood that $\int_A ... \, dA$ represents a double integral.

If we consider a polynomial function with maximum exponents $M$ and $N$ respectivelly for $y$ and $z$, expressed in the following form, having indicated $a_{m,n}$ any real coefficient for the term $y^m z^n$:

:::{math}
:label: eq:marin_polynomial
P(y,z) = 
\sum_{m=0}^M \sum_{n=0}^N a_{m,n} y^m z^n
:::

the integral of the polynomial function over the cross-section, can be written as:

:::{math}
:label: eq:marin_integral_polyniamo
\int_A P(y,z) \, dA
=
\sum_{m=0}^M \sum_{n=0}^N a_{m,n} \mathcal{M}_{m,n}
:::

When the cross-section (i.e. the integration domain) is represented by a closed simply connected polygonal contour with $P$ vertices, the moment of area of order $m, n$ defined by {eq}`eq:marin_moment_area` can be shown to be computed by the following simple relation:

:::{math}
:label: eq:marin_moment_area_polygon
\mathcal{M}_{m,n} = \dfrac{m! n!}{\left(m+n+2\right)!}
\sum_{i=1}^P w_i \sum_{j=0}^m \sum_{k=0}^n 
\binom{j+k}{k} \cdot \binom{m+n-j-k}{n-k}
y_i^{m-j} y_{i+1}^j z_i^{n-k} y_{i+1}^k
:::

where $\binom{n}{k}$ indicates the usual binomial coefficients:

:::{math}
\binom{n}{k} = 
\dfrac{n!}{k!(n-k)!}
:::

and (note that here we assume that $P+1$ restarts from $1$):

:::{math}
w_i = y_i z_{i+1} - y_{i+1} z_i
\qquad
\text{for } i = 1, 2, \ldots, P
:::

:::{hint}
As mentioned, the integral is exact for polygonal surfaces and polynomial functions to be integrated over the surface.

If the cross-section is not polynomial, it must be discretized into a polygon (for instance a circular shape must be discretized as a polygon with $P$ vertices) and therefore the integral is not exact anymore.

If the function to be integrated over the surface is not polynomial, the function can be discretized in linear branches and the cross section must be discretized in small polygons for each polynomial branch; also in this case the integral is not exact anymore.
:::

The main application of this algorithm in *structuralcodes* is for integrating the stresses over the cross-section.

In this context, let's consider a case of unixial bending, for which the stress function can be considered variable only with respect to coordinate $z$:

:::{math}
:label: eq:marin_polynom_stress_unixial_bend
P(y,z) = \sigma(z) = \sum_{n=0}^N a_n z^n = a_0 + a_1 z + a_2 z^2 + \ldots + a_N z^N
:::

In this case the internal forces $N$, $M_y$, $M_z$ can be written as:

<span style="background-color:yellow; color:red">**Note**: check signs! </span>
:::{math}
:label: eq:marin_stress_integration
\begin{aligned}
N &= \int_A \sigma(z) \, dA = \sum_{n=0}^N a_n \mathcal{M}_{0,n}\\
M_y &= \int_A z \sigma(z) \, dA = \sum_{n=0}^N a_n \mathcal{M}_{0,n+1}\\
M_z &= - \int_A y \sigma(z) \, dA = - \sum_{n=0}^N a_n \mathcal{M}_{1,n}
\end{aligned}
:::

::::{hint}
When the bending is not unixial anymore, one can consider the bidimensional normal stress polynomial:

:::{math}
P(y,z) = \sigma(y,z) = \sum_{m=0}^M \sum_{n=0}^N a_{m,n} y^m z^n 
:::

An alternative approach is to rotate the section in order to have uniaxial bending in the new reference system **$y^*z^*$**, see figure [below](theory-fig-marin-rotation). In *structuralcodes* we are adopting the latter approach.

(theory-fig-marin-rotation)=
:::{figure} FigureBendingRotated.png

Rotation of reference system for having uniaxial bending.
:::

::::

The coeficcients $a_n$ are dependent on the stress function $\sigma(z)$, therefore they depend on the constitutive law $\sigma(\varepsilon)$. For this reason, each constitutive law that is used with marin integration, must implement a special `__marin__` method that returns the coefficients $a_n$.

The determination of the coefficients for the constitutive laws implemented in *structuralcodes* is reported in the following subsections.

### Linear elastic material

The constitutive law of a linear elastic material can be written as:

:::{math}
:label: eq:marin-linear-constitutive
\sigma(\varepsilon) = E \cdot \varepsilon
:::

For uniaxial bending the strain $varepsilon$ is equal to:

:::{math}
:label: eq:marin-linear-strain
\varepsilon(z) = \varepsilon_0 + \chi_y \cdot z
:::

The stress function is easily obtained substituting {eq}`eq:marin-linear-strain` in {eq}`eq:marin-linear-constitutive`:

:::{math}
:label: eq:marin-linear-stress-function
\sigma(z) = E \left( \varepsilon_0 + \chi_y \cdot z \right)
:::

Therefore the marin coefficients $a_n$ for linear elastic material are:

:::{math}
:label: eq:marin-linear-coefficients
\begin{aligned}
a_0 &= E \cdot \varepsilon_0 \\
a_1 &= E \cdot \chi_z
\end{aligned}
:::

### Elastic - plastic material

The constitutive law of an elastic plastic material with hardening can be written as:

:::{math}
:label: eq:marin-ep-constitutive
:::

In this case, the constitutive law can be expressed as two polynomial branches (each of degree 1): one for the elastic part and the other for the plastic part.




[^marin1984]: Marín, J. Computing columns, footings and gates through moments of area, Computers & Structures, 18(2), 343-349, 1984.

(theory-fiber-integrator)=
## Fiber integrator
Uses a discretized approach where the section is divided into fibers, each representing a material point. This discretization is performed by a triangulation of the geometry.
