(theory-shell-incremental)=
# Incremental-iterative formulation

The constitutive laws for shells are in general non-linear. To calculate the generalized strains for a given set of stress resultants, we therefore need an incremental-iterative formulation.

The stress resultants of the shell section are calculated for a given set of generalized strains in Eq. [](eq:shell-stress-resultants). The error in this prediction, the residual $\mathbf{R}_{\mathrm{res}}$ is defined as

:::{math}

\mathbf{R}_{\mathrm{res}}
=
\mathbf{R}_{\mathrm{ext}}
-
\mathbf{R}_{\mathrm{int}}

:::

where $\mathbf{R}_{\mathrm{ext}}$ are the applied stress resultants on the section. The residual can furthermore be linearized using a Taylor series expansion around the last incremental value of the generalized strains $\boldsymbol{\varepsilon}_0^k$

:::{math}

\mathbf{R}_{\mathrm{res}}(\boldsymbol{\varepsilon}_0^{k+1})
=
\mathbf{R}_{\mathrm{res}}(\boldsymbol{\varepsilon}_0^{k})
+
\mathbf{J}(\boldsymbol{\varepsilon}_0^{k}) \Delta(\boldsymbol{\varepsilon}_0)

:::

where $\mathbf{J}$ is the Jacobian of the residual. Since we are seeking a solution where the residual is zero, the generalized strain increment can be expressed as

:::{math}
\Delta\boldsymbol{\varepsilon}_0
=
-\mathbf{J}(\boldsymbol{\varepsilon}_0^{k})^{-1}
\mathbf{R}_{\mathrm{res}}(\boldsymbol{\varepsilon}_0^{k})
:::

and the updated generalized strains are calculated as following.

:::{math}

\boldsymbol{\varepsilon}_0^{k+1}
=
\boldsymbol{\varepsilon}_0^k
+
\Delta\boldsymbol{\varepsilon}_0

:::

The Jacobian is expressed as the partial derivative of the residual with respect to the generalized strains.

:::{math}
:label: eq:shell-jacobian-intermediate

\mathbf{J}
=
\dfrac{\partial \mathbf{R}_{\mathrm{res}}}{\partial \boldsymbol{\varepsilon}_0}
=
-\dfrac{\partial \mathbf{R}_{\mathrm{int}}}{\partial \boldsymbol{\varepsilon}_0}
=
-\dfrac{\partial}{\partial \boldsymbol{\varepsilon}_0}
\int_t \mathbf{B}^T \boldsymbol{\sigma}(z)dz

:::

The partial derivative of the stresses with respect to the generalized strains can be derived by differentiating by parts:

:::{math}
:label: eq:shell-partial-derivative-stress

\dfrac{\partial\boldsymbol{\sigma}}{\partial \boldsymbol{\varepsilon}_0}
=
\dfrac{\partial\boldsymbol{\sigma}}{\partial \boldsymbol{\varepsilon}}
\dfrac{\partial\boldsymbol{\varepsilon}}{\partial \boldsymbol{\varepsilon}_0}
=
\dfrac{\partial\boldsymbol{\sigma}}{\partial \boldsymbol{\varepsilon}}
\mathbf{B}

:::

Inserting [](eq:shell-partial-derivative-stress) into [](eq:shell-jacobian-intermediate) we get the following general expression for the Jacobian expressed as the tangent stiffness matrix $\mathbf{K}_T$ of the section, where $\mathbf{D}_T$ is the tangent constitutive matrix.

:::{math}

\mathbf{J}
=
-\int_t \mathbf{B}^T \dfrac{\partial\boldsymbol{\sigma}}{\partial\boldsymbol{\varepsilon}} \mathbf{B}dz
=
-\int_t \mathbf{B}^T \mathbf{D}_T \mathbf{B}dz
=
-\int_t
\begin{bmatrix}
\mathbf{D}_T & -z\mathbf{D}_T\\
-z\mathbf{D}_T & z^2\mathbf{D}_T\\
\end{bmatrix}
dz
=
-\mathbf{K}_T

:::

The incremental-iterative formulation can alternatively be derived using a secant approach. Assume that the section has a size $a$ in both in-plane directions. A set of generalized displacements can be obtained as:

:::{math}
\mathbf{r} = a \boldsymbol{\varepsilon}_0
:::

The stresses can furthermore be expressed as

:::{math}

\boldsymbol{\sigma} = \mathbf{D}_S \boldsymbol{\varepsilon}

:::

where $\mathbf{D}_S$ is the secant constitutive matrix.

We can impose a set of virtual generalized strains $\delta\boldsymbol{\varepsilon}_0$ and obtain the external and internal virtual work, $\delta W_{\mathrm{ext}}$ and $\delta W_{\mathrm{int}}$, respectively.

:::{math}

\delta W_{\mathrm{ext}}
=
\delta \mathbf{r}^T \mathbf{R}_{\mathrm{ext}}a
=
\delta \boldsymbol{\varepsilon}_0^T \mathbf{R}_{\mathrm{ext}}a^2

:::

:::{math}

\delta W_{\mathrm{int}}
=
\int_V
\delta\boldsymbol{\varepsilon}^T
\boldsymbol{\sigma}
dV
=
\delta \boldsymbol{\varepsilon}_0^T
\int_t
\mathbf{B}^T
\mathbf{D}_S
\mathbf{B}
dt
\boldsymbol{\varepsilon}_0
a^2

:::

Equating the external and internal virtual work, we obtain the following relation between the applied stress resultants and the generalized strains, where $\mathbf{K}_S$ is the secant stiffness matrix of the section.

:::{math}
:label: eq:shell-secant-formulation

\mathbf{R}_{\mathrm{ext}}
=
\int_t
\begin{bmatrix}
\mathbf{D}_S & -z\mathbf{D}_S\\
-z\mathbf{D}_S & z^2\mathbf{D}_S\\
\end{bmatrix}
dz
\boldsymbol{\varepsilon}_0
=
\mathbf{K}_S
\boldsymbol{\varepsilon}_0

:::

Finally, the generalized strains are obtained by inverting Eq. [](eq:shell-secant-formulation):

:::{math}

\boldsymbol{\varepsilon}_0
=
\mathbf{K}_S^{-1}
\mathbf{R}_{\mathrm{ext}}

:::

The formulation given by [](eq:shell-secant-formulation) is the default solution method implemented for the shell section in StructuralCodes.