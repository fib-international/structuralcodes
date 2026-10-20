(api-basic-materials)=
# Basic materials

(api-elastic-material)=
## Elastic material

```{eval-rst}
.. autoclass:: structuralcodes.materials.basic.ElasticMaterial

  .. automethod:: __init__
  
  .. autoproperty:: E
  .. autoproperty:: constitutive_law
  .. autoproperty:: name
  .. autoproperty:: density
  .. automethod:: from_material

```

(api-elastic-plastic-material)=
## Elastic-plastic material

```{eval-rst}
.. autoclass:: structuralcodes.materials.basic.ElasticPlasticMaterial

  .. automethod:: __init__
  
  .. autoproperty:: E
  .. autoproperty:: fy
  .. autoproperty:: Eh
  .. autoproperty:: eps_su
  .. autoproperty:: constitutive_law
  .. autoproperty:: name
  .. autoproperty:: density

```

(api-generic-material)=
## Generic material

```{eval-rst}
.. autoclass:: structuralcodes.materials.basic.GenericMaterial

  .. automethod:: __init__
  
  .. autoproperty:: constitutive_law
  .. autoproperty:: name
  .. autoproperty:: density

```
