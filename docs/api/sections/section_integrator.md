(api-section-integrator)=
# Section integrator

## Base class for section integrators

```{eval-rst}
.. autoclass:: structuralcodes.sections.SectionIntegrator

    .. automethod:: prepare_input
    .. automethod:: integrate_stress
    .. automethod:: integrate_modulus
    .. automethod:: integrate_strain_response_on_geometry

```

## Integrator factory

```{eval-rst}
.. autofunction:: structuralcodes.sections.integrator_factory

```

## The Marin integrator

```{eval-rst}
.. autoclass:: structuralcodes.sections.MarinIntegrator

    .. automethod:: prepare_input
    .. automethod:: integrate_stress
    .. automethod:: integrate_modulus
    .. automethod:: integrate_strain_response_on_geometry

```

```{eval-rst}
.. autofunction:: structuralcodes.sections.marin_integration

```

## The fiber integrator

```{eval-rst}
.. autoclass:: structuralcodes.sections.FiberIntegrator

    .. automethod:: prepare_input
    .. automethod:: integrate_stress
    .. automethod:: integrate_modulus
    .. automethod:: integrate_strain_response_on_geometry
    .. automethod:: prepare_triangulation

```
