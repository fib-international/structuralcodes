(theory-shell-displacement)=
# Displacement formulation

Assuming Kirchhoff-Love plate behaviour which is applicable to thin plates with insignificant transverse shear effects, where it is reasonable to assume that plane sections through the thickness remain plane under deformation, the deformed position at mid-section is properly defined by the three displacement components $\{u, v, w\}^T$.

The displacements at any point through the thickness can be expressed as

:::{math}
\begin{bmatrix}
u \\
v \\
w
\end{bmatrix}
=
\begin{bmatrix}
u_0 \\
v_0 \\
w_0
\end{bmatrix}
-z\begin{bmatrix}
\dfrac{\partial w}{\partial x} \\
\dfrac{\partial w}{\partial y} \\
0
\end{bmatrix}
:::

where $\{u_0, v_0, w_0\}$ are the displacements at mid-section. The in-plane strains at an arbitrary level are given by the following relation.

:::{math}

\boldsymbol{\varepsilon}
=
\begin{bmatrix}
\varepsilon_x \\
\varepsilon_y \\
\gamma_{xy}
\end{bmatrix}
=
\begin{bmatrix}
\dfrac{\partial u}{\partial x} \\
\dfrac{\partial v}{\partial y} \\
\dfrac{\partial u}{\partial y} + \dfrac{\partial v}{\partial x}
\end{bmatrix}
-
z
\begin{bmatrix}
\dfrac{\partial^2 w}{\partial x^2} \\
\dfrac{\partial^2 w}{\partial y^2} \\
2\dfrac{\partial^2 w}{\partial x \partial y}
\end{bmatrix}
=
\begin{bmatrix}
\varepsilon_{0x} \\
\varepsilon_{0y} \\
\gamma_{0xy}
\end{bmatrix}
-
z
\begin{bmatrix}
\chi_x \\
\chi_y \\
\chi_{xy}
\end{bmatrix}

:::

The relation between the generalized strains and the strains at an arbitrary level can furthermore be written on the following compact form

:::{math}

\boldsymbol{\varepsilon}
=
\begin{bmatrix}
\begin{array}{c:c}
\mathbf{1} & -z\mathbf{1}
\end{array}
\end{bmatrix}
\begin{bmatrix}
\varepsilon_{0x} \\
\varepsilon_{0y} \\
\gamma_{0xy}\\
\chi_x \\
\chi_y \\
\chi_{xy}
\end{bmatrix}
=
\mathbf{B} \boldsymbol{\varepsilon_0}

:::

where $\mathbf{1}$ is the 3x3 unit matrix, $\mathbf{B}$ is the strain-displacement array, and $\boldsymbol{\varepsilon_0}$ is the array of generalized strains.