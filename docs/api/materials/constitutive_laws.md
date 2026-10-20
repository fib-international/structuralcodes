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

(api-userdefined-law)=
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

## Parallel

```{eval-rst}
.. autoclass:: structuralcodes.materials.constitutive_laws.Parallel

   .. automethod:: __init__
   .. automethod:: get_stress
   .. automethod:: get_tangent
   .. automethod:: get_ultimate_strain

   .. autoproperty:: weights
   .. autoproperty:: wrapped_laws
```

## Concrete smeared cracking

```{eval-rst}
.. autoclass:: structuralcodes.materials.constitutive_laws.ConcreteSmearedCracking

   .. automethod:: __init__
   .. automethod:: get_stress
   .. automethod:: get_tangent
   .. automethod:: get_secant
   .. automethod:: calculate_effective_principal_strains

   .. autoproperty:: uniaxial_law
   .. autoproperty:: strength_reduction_lateral_cracking
   .. autoproperty:: poisson_reduction
   .. autoproperty:: cracking_criterion
   .. autoproperty:: initial_modulus_compression

```

```{eval-rst}
.. autoclass:: structuralcodes.materials.constitutive_laws.GeneralVecchioCollins

   .. automethod:: __init__
   .. automethod:: reduction

   .. autoproperty:: c_1
   .. autoproperty:: c_2
   .. autoproperty:: minimum_value

```

```{eval-rst}
.. autoclass:: structuralcodes.materials.constitutive_laws.ConstantPoissonReduction

   .. automethod:: __init__
   .. automethod:: nu
   .. automethod:: poisson_matrix

```

```{eval-rst}
.. autoclass:: structuralcodes.materials.constitutive_laws.NoTension

   .. automethod:: cracked

```

## Elastic 2D

```{eval-rst}
.. autoclass:: structuralcodes.materials.constitutive_laws.Elastic2D

   .. automethod:: __init__
   .. automethod:: get_stress
   .. automethod:: get_secant

   .. autoproperty:: E
   .. autoproperty:: nu
   .. autoproperty:: stiffness_matrix

```

## Utility functions

```{eval-rst}

.. autofunction:: structuralcodes.materials.constitutive_laws.calculate_principal_strains

```

```{eval-rst}

.. autofunction:: structuralcodes.materials.constitutive_laws.establish_strain_transformation_matrix

```