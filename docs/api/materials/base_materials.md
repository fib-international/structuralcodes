(api-base-materials)=
# Base material classes

## Base material class

```{eval-rst}
.. autoclass:: structuralcodes.core.base.Material

   .. autoproperty:: constitutive_law
   .. autoproperty:: name
   .. autoproperty:: density
```

## Base concrete class

```{eval-rst}
.. autoclass:: structuralcodes.materials.concrete.Concrete

   .. autoproperty:: fck
   .. autoproperty:: constitutive_law
   .. autoproperty:: gamma_c
```

## Base reinforcement steel class

```{eval-rst}
.. autoclass:: structuralcodes.materials.reinforcement.Reinforcement

   .. autoproperty:: fyk
   .. autoproperty:: Es
   .. autoproperty:: ftk
   .. autoproperty:: epsuk
   .. autoproperty:: epsyk
   .. autoproperty:: epsyd
   .. autoproperty:: constitutive_law
```

## Base constitutive law class

```{eval-rst}
.. autoclass:: structuralcodes.core.base.ConstitutiveLaw

   .. autoproperty:: name
   .. automethod:: preprocess_strains_with_limits
   .. automethod:: get_secant
```