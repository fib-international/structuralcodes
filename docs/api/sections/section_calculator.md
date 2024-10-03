(api-section-calculator)=
# Section calculator

## Generic section calculator

```{eval-rst}
.. autoclass:: structuralcodes.sections.GenericSectionCalculator

    .. autoproperty:: n_min
    .. autoproperty:: n_max


    .. automethod:: calculate_limit_axial_load
    .. automethod:: check_axial_load
    .. automethod:: integrate_strain_profile
    .. automethod:: calculate_bending_strength
    .. automethod:: calculate_moment_curvature
    .. automethod:: calculate_nm_interaction_domain
    .. automethod:: calculate_nmm_interaction_domain
    .. automethod:: calculate_mm_interaction_domain

    .. automethod:: get_balanced_failure_strain
    .. automethod:: find_equilibrium_fixed_pivot
    .. automethod:: find_equilibrium_fixed_curvature

```

(api-section-results)=
## Section results

```{eval-rst}
.. autoclass:: structuralcodes.core._section_results.GrossProperties
```

```{eval-rst}
.. autoclass:: structuralcodes.core._section_results.CrackedProperties
```

```{eval-rst}
.. autoclass:: structuralcodes.core._section_results.MomentCurvatureResults
```

```{eval-rst}
.. autoclass:: structuralcodes.core._section_results.UltimateBendingMomentResults
```

```{eval-rst}
.. autoclass:: structuralcodes.core._section_results.NMInteractionDomain

    .. autoproperty:: n
    .. autoproperty:: m_y
    .. autoproperty:: e_a
    .. autoproperty:: k_y
```

```{eval-rst}
.. autoclass:: structuralcodes.core._section_results.NMMInteractionDomain

    .. autoproperty:: n
    .. autoproperty:: m_y
    .. autoproperty:: m_z
    .. autoproperty:: e_a
    .. autoproperty:: k_y
    .. autoproperty:: k_z
```

```{eval-rst}
.. autoclass:: structuralcodes.core._section_results.MMInteractionDomain

    .. autoproperty:: n
    .. autoproperty:: m_y
    .. autoproperty:: m_z
    .. autoproperty:: e_a
    .. autoproperty:: k_y
    .. autoproperty:: k_z
```
