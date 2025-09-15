(api-concrete-materials)=
# Concrete materials

## Concrete factory

```{eval-rst}
.. autofunction:: structuralcodes.materials.concrete.create_concrete
```

(api-concrete-ec2-2004)=
## Eurocode 2 (2004)

```{eval-rst}
.. autoclass:: structuralcodes.materials.concrete.ConcreteEC2_2004

   .. automethod:: __init__

   .. automethod:: fcd
   .. autoproperty:: fck
   .. autoproperty:: fcm
   .. autoproperty:: fctm
   .. autoproperty:: fctk_5
   .. autoproperty:: fctk_95
   .. autoproperty:: Ecm
   .. autoproperty:: gamma_c
   .. autoproperty:: alpha_cc
   .. autoproperty:: eps_c1
   .. autoproperty:: eps_cu1
   .. autoproperty:: k_sargin
   .. autoproperty:: eps_c2
   .. autoproperty:: eps_cu2
   .. autoproperty:: n_parabolic_rectangular
   .. autoproperty:: eps_c3
   .. autoproperty:: eps_cu3
   .. autoproperty:: constitutive_law
   .. autoproperty:: name
   .. autoproperty:: density
```

(api-concrete-mc2010)=
## _fib_ Model Code 2010

```{eval-rst}
.. autoclass:: structuralcodes.materials.concrete.ConcreteMC2010

   .. automethod:: __init__
   
   .. automethod:: fcd
   .. autoproperty:: fck
   .. autoproperty:: fcm
   .. autoproperty:: Eci
   .. autoproperty:: fctm
   .. autoproperty:: fctkmin
   .. autoproperty:: fctkmax
   .. autoproperty:: Gf
   .. autoproperty:: gamma_c
   .. autoproperty:: alpha_cc
   .. autoproperty:: eps_c1
   .. autoproperty:: eps_cu1
   .. autoproperty:: k_sargin
   .. autoproperty:: eps_c2
   .. autoproperty:: eps_cu2
   .. autoproperty:: n_parabolic_rectangular
   .. autoproperty:: eps_c3
   .. autoproperty:: eps_cu3
   .. autoproperty:: constitutive_law
   .. autoproperty:: name
   .. autoproperty:: density
```

(api-concrete-ec2-2023)=
## Eurocode 2 (2023)

```{eval-rst}
.. autoclass:: structuralcodes.materials.concrete.ConcreteEC2_2023

   .. automethod:: __init__

   .. automethod:: fcd
   .. automethod:: fctd
   .. autoproperty:: fck
   .. autoproperty:: fcm
   .. autoproperty:: fctm
   .. autoproperty:: fctk_5
   .. autoproperty:: fctk_95
   .. autoproperty:: Ecm
   .. autoproperty:: gamma_c
   .. autoproperty:: eps_c1
   .. autoproperty:: eps_cu1
   .. autoproperty:: k_sargin
   .. autoproperty:: eps_c2
   .. autoproperty:: eps_cu2
   .. autoproperty:: n_parabolic_rectangular
   .. autoproperty:: constitutive_law
   .. autoproperty:: name
   .. autoproperty:: density
```