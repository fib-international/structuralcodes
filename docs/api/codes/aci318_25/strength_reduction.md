# Strength reduction factors

Strength reduction factors according to ACI 318-25, Chapter 21.

ACI 318 uses LRFD -- strength reduction factors (phi) are applied at the
member capacity level, not at the material level. Design functions in
this library return nominal strengths; phi is applied externally by the caller.

## Moment, axial force, or combined

```{eval-rst}
.. autofunction:: structuralcodes.codes.aci318_25.phi_flexure
```

```{eval-rst}
.. autofunction:: structuralcodes.codes.aci318_25.section_classification
```

## Shear

```{eval-rst}
.. autofunction:: structuralcodes.codes.aci318_25.phi_shear
```

## Torsion

```{eval-rst}
.. autofunction:: structuralcodes.codes.aci318_25.phi_torsion
```

## Bearing

```{eval-rst}
.. autofunction:: structuralcodes.codes.aci318_25.phi_bearing
```
