(api-geometry-creation)=
# Geometry creation

## Point geometry

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

## SurfaceGeometry

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

## Compound geometry

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

## Line object

```{eval-rst}
.. autofunction:: structuralcodes.geometry.create_line_point_angle

```

:::{note}

This function is useful for creating a line which can be used with the {func}`split() <structuralcodes.geometry.SurfaceGeometry.split>` and {func}`split_two_lines() <structuralcodes.geometry.SurfaceGeometry.split_two_lines>` methods.

:::

## Common geometries

In this section the classes and methods for creating special and common geometries are described. Generally these are simply wrappers of base geometries.

```{eval-rst}
.. autoclass:: structuralcodes.geometry.RectangularGeometry

    .. automethod:: __init__

    .. autoproperty:: height
    .. autoproperty:: width

```

```{eval-rst}
.. autoclass:: structuralcodes.geometry.CircularGeometry

    .. automethod:: __init__

    .. autoproperty:: radius
    .. autoproperty:: diameter

```

## Functions for adding reinforcement

```{eval-rst}
.. autofunction:: structuralcodes.geometry.add_reinforcement

```

```{eval-rst}
.. autofunction:: structuralcodes.geometry.add_reinforcement_line

```

```{eval-rst}
.. autofunction:: structuralcodes.geometry.add_reinforcement_circle

```

## Base geometry class

```{eval-rst}
.. autoclass:: structuralcodes.geometry.Geometry

    .. automethod:: __init__

    .. autoproperty:: name
    .. autoproperty:: group_label

```
