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

```{eval-rst}
.. autoclass:: structuralcodes.materials.constitutive_laws.Elastic

   .. automethod:: __init__
   .. automethod:: get_stress
   .. automethod:: get_tangent
   .. automethod:: get_ultimate_strain

```

## Elastic plastic

```{eval-rst}
.. autoclass:: structuralcodes.materials.constitutive_laws.ElasticPlastic

   .. automethod:: __init__
   .. automethod:: get_stress
   .. automethod:: get_tangent
   .. automethod:: get_ultimate_strain

```

## Parabola rectangle

```{eval-rst}
.. autoclass:: structuralcodes.materials.constitutive_laws.ParabolaRectangle

   .. automethod:: __init__
   .. automethod:: get_stress
   .. automethod:: get_tangent
   .. automethod:: get_ultimate_strain

```

## Bilinear elastic perfectly plastic

```{eval-rst}
.. autoclass:: structuralcodes.materials.constitutive_laws.BilinearCompression

   .. automethod:: __init__
   .. automethod:: get_stress
   .. automethod:: get_tangent
   .. automethod:: get_ultimate_strain

```

## Sargin

```{eval-rst}
.. autoclass:: structuralcodes.materials.constitutive_laws.Sargin

   .. automethod:: __init__
   .. automethod:: get_stress
   .. automethod:: get_tangent
   .. automethod:: get_ultimate_strain

```

## Popovics

```{eval-rst}
.. autoclass:: structuralcodes.materials.constitutive_laws.Popovics

   .. automethod:: __init__
   .. automethod:: get_stress
   .. automethod:: get_tangent
   .. automethod:: get_ultimate_strain

```

## Userdefined

```{eval-rst}
.. autoclass:: structuralcodes.materials.constitutive_laws.UserDefined

   .. automethod:: __init__
   .. automethod:: get_stress
   .. automethod:: get_tangent
   .. automethod:: get_ultimate_strain

```

## Initial strain

```{eval-rst}
.. autoclass:: structuralcodes.materials.constitutive_laws.InitialStrain

   .. automethod:: __init__
   .. automethod:: get_stress
   .. automethod:: get_tangent
   .. automethod:: get_ultimate_strain

   .. autoproperty:: strain_compatibility
   .. autoproperty:: wrapped_law

```