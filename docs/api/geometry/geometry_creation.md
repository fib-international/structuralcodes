(api-geometry-creation)=
# Geometry creation

```{eval-rst}
.. autoclass:: structuralcodes.geometry.Geometry

    .. automethod:: __init__

    .. autoproperty:: name
    .. autoproperty:: group_label

```

```{eval-rst}
.. autoclass:: structuralcodes.geometry.PointGeometry

    .. automethod:: __init__

    .. autoproperty:: diameter
    .. autoproperty:: area
    .. autoproperty:: material
    .. autoproperty:: density
    .. autoproperty:: x
    .. autoproperty:: y
    .. autoproperty:: point

    .. automethod:: translate
    .. automethod:: rotate
    .. automethod:: from_geometry

```

```{eval-rst}
.. autoclass:: structuralcodes.geometry.SurfaceGeometry

    .. automethod:: __init__

    .. autoproperty:: area
    .. autoproperty:: centroid
    .. autoproperty:: density

    .. automethod:: calculate_extents
    .. automethod:: split
    .. automethod:: split_two_lines
    .. automethod:: translate
    .. automethod:: rotate

    .. automethod:: from_geometry

```

```{eval-rst}
.. autoclass:: structuralcodes.geometry.CompoundGeometry

    .. automethod:: __init__

    .. autoproperty:: reinforced_concrete
    .. autoproperty:: area

    .. automethod:: calculate_extents
    .. automethod:: translate
    .. automethod:: rotate
    .. automethod:: from_geometry


```

```{eval-rst}
.. autofunction:: structuralcodes.geometry.create_line_point_angle

```

## Functions for adding reinforcement

```{eval-rst}
.. autofunction:: structuralcodes.geometry.add_reinforcement

```

```{eval-rst}
.. autofunction:: structuralcodes.geometry.add_reinforcement_line

```
