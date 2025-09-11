(usage-geometries)=
# Geometries

## General

Geometry objects are responsible for storing a polygon, a point, or a _compound_ of two or more of the former, and their respective material objects. Geometry objects are therefore the physical representation of a [section](#usage-sections).

## Surface geometries

A {class}`SurfaceGeometry <structuralcodes.geometry.SurfaceGeometry>` is initialized by passing a {class}`shapely.Polygon` with an arbitrary shape, and a [material object](#api-materials). See the example below where we create a T-shaped geometry.

Note that adding two or more {class}`SurfaceGeometry <structuralcodes.geometry.SurfaceGeometry>` objects results in a {class}`CompoundGeometry <structuralcodes.geometry.CompoundGeometry>` object.

(code-usage-surface-geometries)=
::::{dropdown-syntax}
:::{literalinclude} ../../_example_code/usage_create_surface_geometries.py
:lines: 3-
:caption: Create surface geometries.
:::
::::

The final geometry is shown in the figure [below](#fig-usage-t-shaped-geometry).

(fig-usage-t-shaped-geometry)=
:::{figure} t_shaped_geometry.svg

The t-shaped geometry created with the code [above](#code-usage-surface-geometries).
:::

:::{tip}

If you work in a Jupyter notebook, and return the geometry object in the final line of a code cell, the shape is visualized in the output below the code cell ✨
:::

:::{tip}

In the code [above](#code-usage-surface-geometries) we used [Shapely](https://shapely.readthedocs.io) to model the shape of the geometry. We created the t-shaped geometry by taking the union of two polygons, but we could just as well have created one t-shaped polygon directly. You can read more about geometry creation with Shapely [here](https://shapely.readthedocs.io/en/stable/geometry.html).

If the t-shaped geometry was made with two different materials, say one material for the flange and one for the web, we could have created two geometry objects, e.g. one {class}`SurfaceGeometry <structuralcodes.geometries.SurfaceGeometry>` for the flange and one {class}`SurfaceGeometry <structuralcodes.geometries.SurfaceGeometry>` for the web, and added these together (`+`) to create a {class}`CompoundGeometry <structuralcodes.geometries.CompoundGeometry>`.

:::

:::{admonition} Creating holes
:class: tip

If you ever need to create a hole in a geometry, simply create a separate geometry object for the hole, and subtract it from the other. This operation returns a new {class}`CompoundGeometry <structuralcodes.geometry.CompoundGeometry>` with a hole ✅

If you prefer to do the modelling with Shapely, this is also possible. Simply subtract the polygon of the hole from the polygon of the base geometry, and use the result as input to a {class}`SurfaceGeometry <structuralcodes.geometries.SurfaceGeometry>`.
:::

## Rectangular and circular geometries

To simplify creating common geometrical shapes, we can use the {class}`RectangularGeometry <structuralcodes.geometry.RectangularGeometry>` or {class}`CircularGeometry <structuralcodes.geometry.CircularGeometry>` classes. These classes are wrappers for creating surface geometries with either a rectangular or a circular shape.

Using {class}`RectangularGeometry <structuralcodes.geometry.RectangularGeometry>`, the above example can be simplified to the example below, where the changed lines are highlighted.

(code-usage-surface-geometries-simple)=
::::{dropdown-syntax}
:::{literalinclude} ../../_example_code/usage_create_surface_geometries_simple.py
:lines: 3-
:emphasize-lines: 1, 14-26
:caption: Create surface geometries using wrapper class for rectangular shape.
:::
::::

## Point geometries and ways to add reinforcement

A {class}`PointGeometry <structuralcodes.geometry.PointGeometry>` is a geometrical object defined by a coordinate, a diameter, and a material. Point geometries are used to represent reinforcement, and can be included in arbitrary locations, and in arbitrary numbers.

There are several ways to create point geometries and add these to other geometry objects as reinforcement. Each of the following options results in a {class}`CompoundGeometry <structuralcodes.geometry.CompoundGeometry>` object holding the geometry object(s) that are reinforced, and the point geometries that act as reinforcement.

- Create {class}`PointGeometry <structuralcodes.geometry.PointGeometry>` objects and add (`+`) them together with the geometry they reinforce.
- Use the {func}`add_reinforcement() <structuralcodes.geometry.add_reinforcement>` function to add one {class}`PointGeometry <structuralcodes.geometry.PointGeometry>` to the geometry.
- Use the {func}`add_reinforcement_line() <structuralcodes.geometry.add_reinforcement_line>` function to add point geometries along a line. This is useful for reinforcing geometries with straight edges.
- Use the {func}`add_reinforcement_circle() <structuralcodes.geometry.add_reinforcement_circle>` function to add point geometries along a circle. This is useful for reinforcing geometries with curved edges.

The example below demonstrates how to add two layers of longitudinal reinforcement to the T-shaped geometry created above. The changed lines are highlighted.

Notice how {func}`add_reinforcement_line() <structuralcodes.geometry.add_reinforcement_line>` returns a new {class}`CompoundGeometry <structuralcodes.geometry.CompoundGeometry>` representing the T-shaped surface geometry _and_ the point geometries for the longitudinal reinforcement. The final geometry with reinforcement is shown in the figure [below](#fig-usage-t-shaped-geometry-with-reinforcement).

(code-usage-surface-geometries-with-reinforcement)=
::::{dropdown-syntax}
:::{literalinclude} ../../_example_code/usage_create_surface_geometries_with_reinforcement.py
:lines: 5-11, 13-70
:emphasize-lines: 1-5, 7, 12-15, 18-20, 24-28, 47-65
:caption: Create surface geometries with reinforcement.
:::
::::

(fig-usage-t-shaped-geometry-with-reinforcement)=
:::{figure} t_shaped_geometry_with_reinforcement.svg

The t-shaped geometry with reinforcement created with the code [above](#code-usage-surface-geometries-with-reinforcement).
:::

## Profiles

StructuralCodes comes with a set of predefined common profiles. A profile class is a wrapper around a {class}`shapely.Polygon` and exposes the most common elastic and plastic section properties.

The following families of profiles are available:

- {class}`IPE <structuralcodes.geometry.profiles.IPE>`
- {class}`HE <structuralcodes.geometry.profiles.HE>`
- {class}`UB <structuralcodes.geometry.profiles.UB>`
- {class}`UC <structuralcodes.geometry.profiles.UC>`
- {class}`UBP <structuralcodes.geometry.profiles.UBP>`
- {class}`IPN <structuralcodes.geometry.profiles.IPN>`
- {class}`UPN <structuralcodes.geometry.profiles.UPN>`

The available profiles within each family can be listed by calling the `.profiles` method on the respective class. The code below shows how to list the available profiles in the {class}`HE <structuralcodes.geometry.profiles.HE>` family.

(code-usage-list-steel-profiles)=
::::{dropdown-syntax}
:::{literalinclude} ../../_example_code/usage_steel_profiles.py
   :lines: 3-5
   :caption: List the available profiles in the HE family.
:::
::::

To create an HEA100 profile, simply pass the name of the profile to the class as shown below.

(code-usage-create-steel-profile)=
::::{dropdown-syntax}
:::{literalinclude} ../../_example_code/usage_steel_profiles.py
   :lines: 3, 6-7
   :caption: List the available profiles in the HE family.
:::
::::

The figure [below](#fig-usage-hea100) shows the shape of the profile.

(fig-usage-hea100)=
:::{figure} hea100.svg

The shape of the profile created with the code [above](#code-usage-create-steel-profile).
:::

Notice how all the predefined profiles expose thicknesses, widths, heights, radii, etc. as properties along with elastic and plastic section properties like {py:attr}`Iy <structuralcodes.geometry.HE.Iy>`, {py:attr}`Wely <structuralcodes.geometry.HE.Wely>`, {py:attr}`Wply <structuralcodes.geometry.HE.Wply>`, and {py:attr}`iy <structuralcodes.geometry.HE.iy>`.

:::::{tip}

The profiles expose the underlying {class}`shapely.Polygon` as the `.polygon` property. This can be passed to a {class}`SurfaceGeometry <structuralcodes.geometry.SurfaceGeometry>` along with a material for use in further calculations. See for example the following code where calculate the bending strength of an IPE100 profile.

(code-usage-steel-profile-in-geometry)=
::::{dropdown-syntax}
:::{literalinclude} ../../_example_code/usage_steel_profile_in_geometry.py
   :lines: 3-
   :caption: Use a profile in a geometry and calculate the bending strength.
:::
::::

:::::
