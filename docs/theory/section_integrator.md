(theory-section-integrators)=
# Section integrators
The algorithms described in [Section calculators](theory-section-calculator) require integration of stress or stiffness over the cross-section. This operation is performed using proper `SectionIntegrator` objects. 

In *structuralcodes* we have currently two distinct integrators: *fiber integrator* that discretize the cross section in a bunch of fibers (whose behavior is uniaxial), and *marin integrator* that is based on the computation of moments of area.

## General concept
Section integrators compute properties such as stiffness, stress distribution, or strain compatibility over a section. They operate on the geometry and material definitions to derive these results.

(theory-fiber-integrator)=
## Fiber integrator
This integrator uses a discretized approach where the section is divided into fibers, each representing a material point. This discretization is performed by a triangulation of the geometry.
Each fiber is defined as the centre point of each triangle and is characterized by its competence area.

According to this integration, the integrals over the cross section are simply sums for all fibers. For instance, internal forces are determined as:


:::{math}
:label: eq:fiber_stress_integration
\begin{aligned}
N &= \int_A \sigma(z) \, dA  \approx \sum_{i=1}^{N_{fibers}} A_i \sigma (\varepsilon_i) \\
M_y &= \int_A z \sigma(z) \, dA  \approx \sum_{i=1}^{N_{fibers}} z A_i \sigma (\varepsilon_i)\\
M_z &= - \int_A y \sigma(z) \, dA  \approx - \sum_{i=1}^{N_{fibers}} y A_i \sigma (\varepsilon_i)
\end{aligned}
:::

where $\varepsilon_i$ is determined using equation {eq}`eq:marin-linear-strain` for each fiber:

:::{math}
:label: eq:fiber_strain_fibers

\varepsilon_i (y_i, z_i) = \varepsilon_0 + \chi_y \cdot z_i - \chi_z \cdot y_i
:::

:::{note}
The triangulation is optimized in order to be executed only the first time a calculation on the section is performed. All fibers information (position $y_i$, $z_i$ and area $A_i$) are stored in numpy arrays. 
Therefore the application of {eq}`eq:fiber_strain_fibers` and of {eq}`eq:fiber_stress_integration` is extremely fast making Fiber integrator the fastest integrator in *structuralcodes*.
:::

(theory-marin-integrator)=
## Marin integrator
Marin integrator implements the idea of prof. Joaquín Marín from University of Venezuela that he published during the '80s [^marin1984].

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

When the cross-section (i.e. the integration domain) is represented by a closed simply connected polygonal contour with $P$ vertices, the moment of area of order $m, n$ defined by {eq}`eq:marin_moment_area` can be shown to be computed by the following relation:

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

:::{note}
As mentioned, the integral is exact for polygonal surfaces and polynomial functions to be integrated over the surface.

If the cross-section is not polynomial, it must be discretized into a polygon (for instance a circular shape must be discretized as a polygon with $P$ vertices) and therefore the integral is not exact anymore.

If the function to be integrated over the surface is not polynomial, the function can be discretized in linear branches and the cross section must be discretized in small polygons for each polynomial branch; also in this case the integral is not exact anymore.
:::

The main application of this algorithm in *structuralcodes* is for integrating the stress or modulus over the cross-section.

In this context, let's consider a case of unixial bending, for which the stress function can be considered variable only with respect to coordinate $z$:

:::{math}
:label: eq:marin_polynom_stress_unixial_bend
P(y,z) = \sigma(z) = \sum_{n=0}^N a_n z^n = a_0 + a_1 z + a_2 z^2 + \ldots + a_N z^N
:::

In this case the internal forces $N$, $M_y$, $M_z$ can be written as:

:::{math}
:label: eq:marin_stress_integration
\begin{aligned}
N &= \int_A \sigma(z) \, dA = \sum_{n=0}^N a_n \mathcal{M}_{0,n}\\
M_y &= \int_A z \sigma(z) \, dA = \sum_{n=0}^N a_n \mathcal{M}_{0,n+1}\\
M_z &= - \int_A y \sigma(z) \, dA = - \sum_{n=0}^N a_n \mathcal{M}_{1,n}
\end{aligned}
:::

::::{hint}
When the bending is not uniaxial anymore, one can consider the bidimensional normal stress polynomial:

:::{math}
P(y,z) = \sigma(y,z) = \sum_{m=0}^M \sum_{n=0}^N a_{m,n} y^m z^n 
:::

An alternative approach is to rotate the section in order to have uniaxial bending in the new reference system **$y^*z^*$**, see figure [below](theory-fig-marin-rotation). In *structuralcodes* we are adopting the latter approach.

(theory-fig-marin-rotation)=
:::{figure} FigureBendingRotated.png

Rotation of reference system for having uniaxial bending.
:::

::::

The coeficients $a_n$ are dependent on the stress function $\sigma(z)$, therefore they depend on the constitutive law $\sigma(\varepsilon)$. For this reason, each constitutive law that is used with Marin integration, must implement a special `__marin__` method that returns the coefficients $a_n$ for each polynomial branch that defines the constitutive law.

In determination of stiffness matrix, we need to integrate the tangent modulus over the cross-section. In this case we are integrating the function $\sigma'$, where $\cdot'$ indicates the derivative (i.e. the tangent modulus).

The determination of the coefficients for the constitutive laws implemented in *structuralcodes* is reported in the following subsections.

:::{note}
This means that if you write your custom constitutive law you must calculate the Marin coefficients for that law?

Don't worry, we have done some work for you! Every material has a base version of `__marin__` method, so you don't have to implement it mandatorily. Keep in mind that the free `__marin__` method works under a very simple idea: we linearize in several branches your law, obtaining a discretize piecewise constitutive law. For this discretized law we have the coefficients computed for the piece-wise linear law described [below](theory-marin-userdefined-law).

If your law can be represented (at least in branches) as polynomial functions it is way more efficient writing your custom `__marin__` method since the integrator will split the polygon into as few as possible sub-polygons to perform the intergration.

So it is up to you! You can avoid bothering calculating and implementing Marin coefficients for stress and stiffness, but you can (and probably should) do it if you can optimize the calculation speed and if your constitutive law permits it.
:::

### Linear elastic material

The constitutive law of a linear elastic material can be written as:

:::{math}
:label: eq:marin-linear-constitutive
\sigma(\varepsilon) = E \cdot \varepsilon
:::

For uniaxial bending the strain $\varepsilon$ is equal to:

:::{math}
:label: eq:marin-linear-strain
\varepsilon(z) = \varepsilon_0 + \chi_y \cdot z
:::

The stress function is easily obtained substituting {eq}`eq:marin-linear-strain` in {eq}`eq:marin-linear-constitutive`:

:::{math}
:label: eq:marin-linear-stress-function
\sigma(z) = E \left( \varepsilon_0 + \chi_y \cdot z \right)
:::

Therefore the Marin coefficients $a_n$ for linear elastic material are:

:::{math}
:label: eq:marin-linear-coefficients
\begin{aligned}
a_0 &= E \cdot \varepsilon_0 \\
a_1 &= E \cdot \chi_y
\end{aligned}
:::

The stiffness is simply $E$. Therefore:

:::{math}
:label: eq:marin-linear-stiffness-function
\sigma'(z) = E
:::

And the only Marin coefficient is:

:::{math}
:label: eq:marin-linear-coefficients-stiffness
a_0 = E 
:::

### Elastic-plastic material

The constitutive law of an elastic plastic material with linear hardening is represented as in the figure [below](figure-theory-marin-ep-constitutive). 

(figure-theory-marin-ep-constitutive)=
:::{figure} Figure_const_ep.png

Elastic plastic with linear hardening law
:::

In this case, the constitutive law cannot be expressed as a whole with a single polynomial function, but looking at the picture [above](figure-theory-marin-ep-constitutive) we can recognize 5 different portions, each one representing a different polynomial function. The constitutive law can then be written for the five branches as:

:::{math}
:label: eq:const-law-ep

\sigma(\varepsilon) =
\begin{cases}
0 & \text{for } \varepsilon < -\varepsilon_u \\
E_h \cdot \varepsilon - \Delta \sigma & \text{for } -\varepsilon_u \le \varepsilon < \varepsilon_y \\
E \cdot \varepsilon & \text{for } -\varepsilon_y \le \varepsilon < \varepsilon_y \\
E_h \cdot \varepsilon + \Delta \sigma & \text{for } \varepsilon_y \le \varepsilon < \varepsilon_u \\
0 & \text{for } \varepsilon > \varepsilon_u
\end{cases} 
:::

where $E_h$ is the hardening modulus (i.e. the slope of the hardning branch) and $\Delta \sigma = f_y \cdot \left(1 - \dfrac{E_h}{E}\right)$.

Substituting {eq}`eq:marin-linear-strain` in {eq}`eq:const-law-ep`, one can compute the Marin coefficients for each branch. The coefficients are reported in the following table.

:::{list-table} Marin coefficients for stress integration for the different branches of elastic-plastic constitutive law
:header-rows: 1


* - Branch name
  - Strain range
  - $a_0$
  - $a_1$
* - Elastic portion
  - $-\varepsilon_y < \varepsilon \le \varepsilon_y$
  - $E \cdot \varepsilon_0$
  - $E \cdot \chi_y$
* - Hardening (positive)
  - $\varepsilon_y < \varepsilon \le \varepsilon_u$
  - $E_h \cdot \varepsilon_0 + \Delta \sigma$
  - $E_h \cdot \chi_y$
* - Hardening (negative)
  - $-\varepsilon_u < \varepsilon \le -\varepsilon_y$
  - $E_h \cdot \varepsilon_0 - \Delta \sigma$
  - $E_h \cdot \chi_y$
* - Zero
  - otherwise
  - 0.0
  - n.a.
:::

The `__marin__` method returns the coefficients reported in the table above and the strain limits for each portion.

Then the integrator discretizes the polygon cutting it into the portions needed (depending on the strain distribution in the polygon) and for each subpolygon the exact Marin integration is performed.

The tangent stiffness can be represented by the following equation:

:::{math}
:label: eq:const-law-ep-stiffness

\sigma'(\varepsilon) =
\begin{cases}
0 & \text{for } \varepsilon < -\varepsilon_u \\
E_h  & \text{for } -\varepsilon_u \le \varepsilon < \varepsilon_y \\
E  & \text{for } -\varepsilon_y \le \varepsilon < \varepsilon_y \\
E_h  & \text{for } \varepsilon_y \le \varepsilon < \varepsilon_u \\
0 & \text{for } \varepsilon > \varepsilon_u
\end{cases} 
:::

Therefore the marin coefficients for stiffness function are:

:::{list-table} Marin coefficients for stiffness integration for the different branches of elastic-plastic constitutive law
:header-rows: 1

* - Branch name
  - Strain range
  - $a_0$
* - Elastic portion
  - $-\varepsilon_y < \varepsilon \le \varepsilon_y$
  - $E$
* - Hardening (positive)
  - $\varepsilon_y < \varepsilon \le \varepsilon_u$
  - $E_h$
* - Hardening (negative)
  - $-\varepsilon_u < \varepsilon \le -\varepsilon_y$
  - $E_h$
* - Zero
  - otherwise
  - 0.0
:::

### Bilinear compression

The bilinear compression is characterized by no strength in tension and a linear behavior in compression up to $\varepsilon_{c0}$ and then a constant stress equal to $f_c$ (negative). After $\varepsilon_{cu}$, the stress drops to 0. The constituive law is depicted the figure [below](figure-theory-marin-bilin-compr-constitutive). 

(figure-theory-marin-bilin-compr-constitutive)=
:::{figure} Figure_const_bilinearcompression.png

Elastic plastic with linear hardening law
:::

The constitutive law can be written in the following branches:

:::{math}
:label: eq:const-law-bilin-compr

\sigma(\varepsilon) =
\begin{cases}
0 & \text{for } \varepsilon \ge 0 \\
E \cdot \varepsilon & \text{for } \varepsilon_{c0} \le \varepsilon < 0 \\
f_c & \text{for } \varepsilon_{cu} \le \varepsilon < \varepsilon_{c0} \\
0 & \text{for } \varepsilon < \varepsilon_{cu}
\end{cases} 
:::

Substituting {eq}`eq:marin-linear-strain` in {eq}`eq:const-law-bilin-compr`, one can compute the Marin coefficients for each branch. The coefficients are reported in the following table.

:::{list-table} Marin coefficients for stress integration for the different branches of bilinear compression constitutive law
:header-rows: 1


* - Branch name
  - Strain range
  - $a_0$
  - $a_1$
* - Elastic portion
  - $\varepsilon_{c0} \le \varepsilon < 0$
  - $E \cdot \varepsilon_0$
  - $E \cdot \chi_y$
* - Constant portion
  - $\varepsilon_{cu} \le \varepsilon < \varepsilon_{c0}$
  - $f_c$
  - n.a.
* - Zero
  - otherwise
  - 0.0
  - n.a.
:::

The Marin coefficients for each branch for stiffness function are reported in the following table.

:::{list-table} Marin coefficients for stiffness integration for the different branches of bilinear compression constitutive law
:header-rows: 1


* - Branch name
  - Strain range
  - $a_0$
* - Elastic portion
  - $\varepsilon_{c0} \le \varepsilon < 0$
  - $E$
* - Constant portion
  - $\varepsilon_{cu} \le \varepsilon < \varepsilon_{c0}$
  - 0.0
* - Zero
  - otherwise
  - 0.0
:::

### Parabola rectangle

The parabola rectangle law, typically used for concrete, cannot be written as a whole with a polynomial law. Despite that, looking at the constitutive law plotted in the figure [below](figure-theory-marin-paraborect-constitutive), we can recognize a portion for the parabolic branch, a portion for the constant part and finally a portion for zero stress (both in tension and for excessive compression)


(figure-theory-marin-paraborect-constitutive)=
:::{figure} Figure_const_parabolarectangle.png

Parabola rectangle law
:::

The constitutive law can therefore be written as:

:::{math}
:label: eq:const-law-paraborect

\sigma(\varepsilon) =
\begin{cases}
0 & \text{for } \varepsilon < \varepsilon_{cu} \\
f_c & \text{for } \varepsilon_{cu} \le \varepsilon < \varepsilon_{c0} \\
f_c \left[ 1 - \left( 1 - \dfrac{\varepsilon}{\varepsilon_{c0}} \right) ^2 \right] & \text{for } \varepsilon_{c0} \le \varepsilon < 0 \\
0, & \text{for } \varepsilon > 0 \\
\end{cases} 
:::

Substituting  {eq}`eq:marin-linear-strain` in {eq}`eq:const-law-paraborect` the Marin coefficients for each branch can be determined. They are reported in the following table.

:::{list-table} Marin coefficients for stress integration for the different branches of parabola rectangle constitutive law
:header-rows: 1


* - Branch name
  - Strain range
  - $a_0$
  - $a_1$
  - $a_2$
* - Parabolic portion
  - $\varepsilon_{c0} \le \varepsilon < 0$
  - $\dfrac{f_c \varepsilon_{0}}{\varepsilon_{c0}}  \left( 2 - \dfrac{\varepsilon_0}{\varepsilon_{c0}} \right)$
  - $\dfrac{2 f_c \chi_y}{\varepsilon_{c0}} \left( 1 - \dfrac{\varepsilon_0}{\varepsilon_{c0}} \right)$
  - $- \dfrac{f_c \chi_y^2}{\varepsilon_{c0}}$
* - Constant portion
  - $\varepsilon_{cu} \le \varepsilon < \varepsilon_{c0}$
  - $f_c$
  - n.a.
  - n.a.
* - Zero
  - otherwise
  - 0.0
  - n.a.
  - n.a.
:::

The `__marin__` method return the coefficients reported in the table above and the strain limits for each portion.

Then the integrator discretizes the polygon cutting it into the portions needed (depending on the strain distribution in the polygon) and for each subpolygon the exact Marin integration is performed.

If the function to be integrated over the cross-section is the stiffness, represented by the following function (only corresponding to the parabolic portion, being 0 elsewhere):

:::{math}
:label: eq:const-law-paraborect-stiffness

\sigma'(\varepsilon) = \dfrac{2 f_c}{\varepsilon_{c0}} \left( 1 - \dfrac{\varepsilon}{\varepsilon_{c0}}\right)
:::

one can write the stiffness as a function of $z$:

:::{math}
:label: eq:const-law-paraborect-stiffness-fun

\sigma'(z) = \dfrac{2 f_c}{\varepsilon_{c0}} \left( 1 - \dfrac{\varepsilon_0}{\varepsilon_{c0}}\right) - \dfrac{2 f_c}{\varepsilon_{c0}} \dfrac{\chi_y}{\varepsilon_{c0}} \cdot z
:::

The marin coefficients for stiffness function integration are therefore:

:::{list-table} Marin coefficients for stiffness integration for the different branches of parabola rectangle constitutive law
:header-rows: 1


* - Branch name
  - Strain range
  - $a_0$
  - $a_1$
* - Parabolic portion
  - $\varepsilon_{c0} \le \varepsilon < 0$
  - $\dfrac{2 f_c}{\varepsilon_{c0}} \left( 1 - \dfrac{\varepsilon_0}{\varepsilon_{c0}}\right)$
  - $- \dfrac{2 f_c}{\varepsilon_{c0}} \dfrac{\chi_y}{\varepsilon_{c0}}$
* - Zero
  - otherwise
  - 0.0
  - n.a.
:::

(theory-marin-userdefined-law)=
### UserDefined constitutive law

A `UserDefined` constitutive law (see [api documentation]()) permits the user to input any custom law discretizing it in a piecewise linear function.

In this case, for each branch, the `__marin__` function return the coefficients $a_0$, $a_1$ and the strain limits for the branch.

Then the integrator discretizes the polygon cutting it into how many portions needed for the several branches (depending on the strain distribution in the polygon) and for each subpolygon the exact Marin integration is performed

For each $i$-th branch, the constitutive law can be written as:

:::{math}
:label: eq:marin-piecewise-linear-constitutive
\sigma(\varepsilon) = \sigma_{i-1} + \dfrac{\sigma_i - \sigma_{i-1}}{\varepsilon_i - \varepsilon_{i-1}} {\varepsilon - \varepsilon_{i-1}}
:::

And substituting {eq}`eq:marin-linear-strain` in {eq}`eq:marin-piecewise-linear-constitutive`, the stress function can be written as:

:::{math}
:label: eq:marin-piecewise-linear-stress-function

\begin{aligned}
\sigma(z) &= \sigma_{i-1} + \dfrac{\sigma_i - \sigma_{i-1}}{\varepsilon_i - \varepsilon_{i-1}} {\varepsilon_0 + \chi_y \cdot z - \varepsilon_{i-1}} \\
&= \sigma_{i-1} + K_i \left( \varepsilon_0 + \chi_y \cdot z - \varepsilon_{i-1} \right)
\end{aligned}
:::

Having introduced the symbol $K_i$ for the stiffness of the $i$-th branch.

Therefore the Marin coefficients $a_n$ for each branch of a piecewise-linear elastic material are:

:::{math}
:label: eq:marin-piecewise-linear-coefficients
\begin{aligned}
a_0 &= \sigma_{i-1} + K_i \left( \varepsilon_0 - \varepsilon_i-1 \right) \\
a_1 &= K_i \cdot \chi_y
\end{aligned}
:::

The Marin coefficient for stiffness function integration for each branch of a piecewise-linear elastic material is:

:::{math}
:label: eq:marin-piecewise-linear-coefficients-stiffness
a_0 = K_i
:::

[^marin1984]: Marín, J. Computing columns, footings and gates through moments of area, Computers & Structures, 18(2), 343-349, 1984.


