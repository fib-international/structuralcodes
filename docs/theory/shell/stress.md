(theory-shell-stress)=
# Stress resultants

The stresses at a level $z$ through the thickness are calculated using suitable two-dimensional constitutive laws:

:::{math}

\begin{bmatrix}
\sigma_x(z) \\
\sigma_y(z) \\
\tau_{xy}(z)
\end{bmatrix}
=
\boldsymbol{\sigma}(z)
=
\boldsymbol{\sigma}(\boldsymbol{\varepsilon(z)})

:::

The stress resultants of the shell section are obtained through the following integrals.

:::{math}

n_x = \int_t \sigma_x dz

:::

:::{math}

n_y = \int_t \sigma_y dz

:::

:::{math}

n_{xy} = \int_t \tau_{xy} dz

:::

:::{math}

m_x = -\int_t \sigma_x z dz

:::

:::{math}

m_y = -\int_t \sigma_y z dz

:::

:::{math}

m_{xy} = -\int_t \tau_{xy} z dz

:::

The stress resultants can alternatively be written in the following more compact form.

:::{math}
:label: eq:shell-stress-resultants

\mathbf{R}_{\mathrm{int}}
=
\int_t
\begin{bmatrix}
\begin{array}{c:c}
\mathbf{1} & -z \mathbf{1}
\end{array}
\end{bmatrix}^T
\boldsymbol{\sigma}(z)
dz
=
\int_t \mathbf{B}^T \boldsymbol{\sigma}(z)dz

:::