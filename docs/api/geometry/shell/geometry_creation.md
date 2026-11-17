(api-shell-geometry-creation)=
# Geometry creation

## Shell geometry

```{eval-rst}
.. autoclass:: structuralcodes.geometry.ShellGeometry

    .. automethod:: __init__

    .. autoproperty:: thickness
    .. autoproperty:: material
    .. autoproperty:: reinforcement

    .. automethod:: add_reinforcement

    .. automethod:: from_geometry

```

## Shell reinforcement

```{eval-rst}
.. autoclass:: structuralcodes.geometry.ShellReinforcement

    .. automethod:: __init__

    .. autoproperty:: z
    .. autoproperty:: n_bars
    .. autoproperty:: cc_bars
    .. autoproperty:: diameter_bar
    .. autoproperty:: phi
    .. autoproperty:: area
    .. autoproperty:: T
    .. autoproperty:: local_stress_to_global_force
    .. autoproperty:: local_modulus_to_global_stiffness
    

    .. automethod:: generalized_strain_to_longitudinal

```