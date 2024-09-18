(api-reinforcement-steel-materials)=
# Reinforcement steel materials

## Reinforcement steel factory

```{eval-rst}
.. autofunction:: structuralcodes.materials.reinforcement.create_reinforcement
```

## Eurocode 2 (2004)

```{eval-rst}
.. autoclass:: structuralcodes.materials.reinforcement.ReinforcementEC2_2004

   .. automethod:: fyd
   .. automethod:: ftd
   .. automethod:: epsud
   .. autoproperty:: gamma_s
   .. autoproperty:: gamma_eps
```

## _fib_ Model Code 2010

```{eval-rst}
.. autoclass:: structuralcodes.materials.reinforcement.ReinforcementMC2010

   .. automethod:: fyd
   .. automethod:: ftd
   .. automethod:: epsud
   .. autoproperty:: gamma_s
   .. autoproperty:: gamma_eps
```

## Eurocode 2 (2023)

```{eval-rst}
.. autoclass:: structuralcodes.materials.reinforcement.ReinforcementEC2_2023

   .. automethod:: fyd
   .. automethod:: ftd
   .. automethod:: epsud
   .. autoproperty:: gamma_s
```