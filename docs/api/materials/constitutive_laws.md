(api-constitutive-laws)=
# Constitutive laws

(api-constitutive-law-factory)=
## Constitutive law factory

```{eval-rst}
.. autofunction:: structuralcodes.materials.constitutive_laws.get_constitutive_laws_list

```

```{eval-rst}
.. autofunction:: structuralcodes.materials.constitutive_laws.create_constitutive_law

```

## Linear elastic

This section describes the api of {class}`Elastic`. For the theory reference see [this page](theory-constitutive-elastic).

```{eval-rst}
.. autoclass:: structuralcodes.materials.constitutive_laws.Elastic

   .. automethod:: __init__
   .. automethod:: get_stress
   .. automethod:: get_tangent
   .. automethod:: get_ultimate_strain

```

## Elastic plastic

This section describes the api of {class}`ElasticPlastic`. For the theory reference see [this page](theory-constitutive-elastic-plastic).

```{eval-rst}
.. autoclass:: structuralcodes.materials.constitutive_laws.ElasticPlastic

   .. automethod:: __init__
   .. automethod:: get_stress
   .. automethod:: get_tangent
   .. automethod:: get_ultimate_strain

```

## Parabola rectangle

This section describes the api of {class}`ParabolaRectangle`. For the theory reference see [this page](theory-constitutive-parabola-rectangle).

```{eval-rst}
.. autoclass:: structuralcodes.materials.constitutive_laws.ParabolaRectangle

   .. automethod:: __init__
   .. automethod:: get_stress
   .. automethod:: get_tangent
   .. automethod:: get_ultimate_strain

```

## Bilinear compression

This section describes the api of {class}`BilinearCompression`. For the theory reference see [this page](theory-constitutive-bilinear-compression).

```{eval-rst}
.. autoclass:: structuralcodes.materials.constitutive_laws.BilinearCompression

   .. automethod:: __init__
   .. automethod:: get_stress
   .. automethod:: get_tangent
   .. automethod:: get_ultimate_strain

```

## Sargin

This section describes the api of {class}`Sargin`. For the theory reference see [this page](theory-constitutive-sargin).

```{eval-rst}
.. autoclass:: structuralcodes.materials.constitutive_laws.Sargin

   .. automethod:: __init__
   .. automethod:: get_stress
   .. automethod:: get_tangent
   .. automethod:: get_ultimate_strain

```

## Popovics

This section describes the api of {class}`Popovics`. For the theory reference see [this page](theory-constitutive-popovics).

```{eval-rst}
.. autoclass:: structuralcodes.materials.constitutive_laws.Popovics

   .. automethod:: __init__
   .. automethod:: get_stress
   .. automethod:: get_tangent
   .. automethod:: get_ultimate_strain

```

(api-elasticfragiletension-law)=
## Linear elastic fragile tension law

This section describes the api of {class}`ElasticFragileTension`. For the theory reference see [this page](theory-constitutive-elastic-fragile-tension).

```{eval-rst}
.. autoclass:: structuralcodes.materials.constitutive_laws.ElasticFragileTension

   .. automethod:: __init__
   .. automethod:: get_stress
   .. automethod:: get_tangent
   .. automethod:: get_ultimate_strain

 ```

(api-userdefined-law)=
## Userdefined

This section describes the api of {class}`UserDefined`. For the theory reference see [this page](theory-constitutive-user-defined).

```{eval-rst}
.. autoclass:: structuralcodes.materials.constitutive_laws.UserDefined

   .. automethod:: __init__
   .. automethod:: get_stress
   .. automethod:: get_tangent
   .. automethod:: get_ultimate_strain

```

## Initial strain

This section describes the api of {class}`InitialStrain`. For the theory reference see [this page](theory-constitutive-initial-strain).

```{eval-rst}
.. autoclass:: structuralcodes.materials.constitutive_laws.InitialStrain

   .. automethod:: __init__
   .. automethod:: get_stress
   .. automethod:: get_tangent
   .. automethod:: get_ultimate_strain

   .. autoproperty:: strain_compatibility
   .. autoproperty:: wrapped_law

```

## Parallel

This section describes the api of {class}`Parallel`. For the theory reference see [this page](theory-constitutive-parallel).

```{eval-rst}
.. autoclass:: structuralcodes.materials.constitutive_laws.Parallel

   .. automethod:: __init__
   .. automethod:: get_stress
   .. automethod:: get_tangent
   .. automethod:: get_ultimate_strain

   .. autoproperty:: weights
   .. autoproperty:: wrapped_laws
```
