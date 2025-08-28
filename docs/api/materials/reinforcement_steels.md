(api-reinforcement-steel-materials)=
# Reinforcement steel materials

## Reinforcement steel factory

```{eval-rst}
.. autofunction:: structuralcodes.materials.reinforcement.create_reinforcement
```

(api-reinforcement-ec2-2004)=
## Eurocode 2 (2004)

```{eval-rst}
.. autoclass:: structuralcodes.materials.reinforcement.ReinforcementEC2_2004

   .. automethod:: __init__

   .. automethod:: fyd
   .. automethod:: ftd
   .. automethod:: epsud
   .. autoproperty:: fyk
   .. autoproperty:: Es
   .. autoproperty:: ftk
   .. autoproperty:: epsuk
   .. autoproperty:: epsyk
   .. autoproperty:: epsyd
   .. autoproperty:: constitutive_law
   .. autoproperty:: gamma_s
   .. autoproperty:: gamma_eps
   .. autoproperty:: name
   .. autoproperty:: density
```

(api-reinforcement-mc2010)=
## _fib_ Model Code 2010

```{eval-rst}
.. autoclass:: structuralcodes.materials.reinforcement.ReinforcementMC2010

   .. automethod:: __init__
   
   .. automethod:: fyd
   .. automethod:: ftd
   .. automethod:: epsud
   .. autoproperty:: fyk
   .. autoproperty:: Es
   .. autoproperty:: ftk
   .. autoproperty:: epsuk
   .. autoproperty:: epsyk
   .. autoproperty:: epsyd
   .. autoproperty:: constitutive_law
   .. autoproperty:: gamma_s
   .. autoproperty:: gamma_eps
   .. autoproperty:: name
   .. autoproperty:: density
```

(api-reinforcement-ec2-2023)=
## Eurocode 2 (2023)

```{eval-rst}
.. autoclass:: structuralcodes.materials.reinforcement.ReinforcementEC2_2023

   .. automethod:: __init__
   
   .. automethod:: fyd
   .. automethod:: ftd
   .. automethod:: epsud
   .. autoproperty:: fyk
   .. autoproperty:: Es
   .. autoproperty:: ftk
   .. autoproperty:: epsuk
   .. autoproperty:: epsyk
   .. autoproperty:: epsyd
   .. autoproperty:: constitutive_law
   .. autoproperty:: gamma_s
   .. autoproperty:: name
   .. autoproperty:: density
```