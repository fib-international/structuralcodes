(quickstart)=
# Quickstart

This example shows how to use `structuralcodes` to calculate the response of a rectangular reinforced concrete section. Follow the example step-by-step, or {ref}`skip to the end <quickstart-full-example>` if you are in a hurry.

Import relevant functions and classes:

```{eval-rst}
.. literalinclude:: quickstart_example.py
   :lines: 3-9
```

:::{seealso}
{ref}`Library structure <api-structure>`.
:::

Set the active design code to Eurocode 2, 2004 (`ec2_2004`):

```{eval-rst}
.. literalinclude:: quickstart_example.py
   :lines: 11-12
```

:::{seealso}
{func}`set_design_code() <structuralcodes.codes.set_design_code>`
:::

Create a concrete and a reinforcement material:

```{eval-rst}
.. literalinclude:: quickstart_example.py
   :lines: 14-24
```

:::{seealso}
{ref}`Material reference <api-materials>`.

The concrete factory {func}`create_concrete() <structuralcodes.materials.concrete.create_concrete>`.

The reinforcement {func}`create_reinforcement() <structuralcodes.materials.reinforcement.create_reinforcement>`.
:::

Create a {class}`SurfaceGeometry <structuralcodes.geometry.SurfaceGeometry>` based on a {class}`shapely.Polygon` and the {class}`Concrete <structuralcodes.materials.concrete.Concrete>` created with {func}`create_concrete() <structuralcodes.materials.concrete.create_concrete>`:

```{eval-rst}
.. literalinclude:: quickstart_example.py
   :lines: 26-39
```

:::{seealso}
{class}`shapely.Polygon`

{ref}`Geometry reference <api-geometry>`.
:::

Add reinforcement to the geometry:

```{eval-rst}
.. literalinclude:: quickstart_example.py
   :lines: 41-62
```

Create a {class}`GenericSection <structuralcodes.sections.GenericSection>` based on the geometry:

```{eval-rst}
.. literalinclude:: quickstart_example.py
   :lines: 64-65
```

:::{seealso}
{ref}`Section reference <api-sections>`
:::

Call the {func}`.calculate_moment_curvature() <structuralcodes.sections.GenericSectionCalculator.calculate_moment_curvature>` method on the {class}`GenericSectionCalculator <structuralcodes.sections.GenericSectionCalculator>` to calculate the moment-curvature relation:

```{eval-rst}
.. literalinclude:: quickstart_example.py
   :lines: 67-68
```

:::{seealso}
{ref}`Section calculator reference <api-section-calculator>`

{ref}`Section integrator reference <api-section-integrator>`

{ref}`Section results reference <api-section-results>`
:::

Full example:

(quickstart-full-example)=
```{eval-rst}
.. literalinclude:: quickstart_example.py
   :linenos:
```
