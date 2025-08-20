(usage-sections)=
# Sections

## General

In StructuralCodes, a _Section_ is responsible for connecting a [geometry](#usage-geometries) object with an API for performing calculations on the geometry.

The section contains a _SectionCalculator_ which contains the API for doing calculations on the geometry, including for example lower level methods for integrating a strain profile to stress resultants, or calculating a strain profile based on stress resultants, or higher level methods like calculating a moment-curvature relation, or an interaction diagram between moments and axial force.

The SectionCalculator uses a _SectionIntegrator_ to integrate the strain response on a geometry.

## The generic beam section

The {class}`GenericSection <structuralcodes.sections.GenericSection>` takes a {class}`SurfaceGeometry <structuralcodes.geometry.SurfaceGeometry>` or a {class}`CompoundGeometry <structuralcodes.geometry.CompoundGeometry>` as input, and is capable of calculating the response of an arbitrarily shaped geometry with arbitrary reinforcement layout, subject to stresses in the direction of the beam axis.

In the example [below](#code-usage-generic-section), we continue the example with the T-shaped geometry.

(usage-sections-gross-properties-tip)=
:::{tip}
If you are looking for the gross properties of the section, these are available at the `.gross_properties` property on the `Section` object ✅
:::

Since the axial force and the moments are all relative to the origin, we start by translating the geometry such that the centroid is alligned with the origin.

Notice how we can use {py:meth}`.calculate_bending_strength() <structuralcodes.sections.GenericSectionCalculator.calculate_bending_strength>` and {py:meth}`.calculate_limit_axial_load() <structuralcodes.sections.GenericSectionCalculator.calculate_limit_axial_load>` to calculate the bending strength and the limit axial loads in tension and compression.

Furthermore, we have the following methods:

{py:meth}`.integrate_strain_profile() <structuralcodes.sections.GenericSectionCalculator.integrate_strain_profile>`
: Calculate the stress resultants for a given strain profile.

{py:meth}`.calculate_strain_profile() <structuralcodes.sections.GenericSectionCalculator.calculate_strain_profile>`
: Calculate the strain profile for given stress resultants.

{py:meth}`.calculate_moment_curvature() <structuralcodes.sections.GenericSectionCalculator.calculate_moment_curvature>`
: Calculate the moment-curvature relation.

{py:meth}`.calculate_nm_interaction_domain() <structuralcodes.sections.GenericSectionCalculator.calculate_nm_interaction_domain>`
: Calculate the interaction domain between axial load and bending moment.

{py:meth}`.calculate_nmm_interaction_domain() <structuralcodes.sections.GenericSectionCalculator.calculate_nmm_interaction_domain>`
: Calculate the interaction domain between axial load and biaxial bending.

{py:meth}`.calculate_mm_interaction_domain() <structuralcodes.sections.GenericSectionCalculator.calculate_mm_interaction_domain>`
: Calculate the interaction domain between biaxial bending for a given axial load.

See the {class}`GenericSectionCalculator <structuralcodes.sections.GenericSectionCalculator>` for a complete list.

(code-usage-generic-section)=
::::{dropdown-syntax}
:::{literalinclude} ../../_example_code/usage_create_surface_geometries_with_reinforcement.py
:lines: 12-13, 72-90
:caption: Create a section and perform calculations.
:::
::::

Notice how for example {py:meth}`.calculate_moment_curvature() <structuralcodes.sections.GenericSectionCalculator.calculate_moment_curvature>` returns a custom dataclass of type {class}`MomentCurvatureResults <structuralcodes.core._section_results.MomentCurvatureResults>`. If we inspect this class further, we find among its attributes `chi_y` and `m_y`. These are the curvature and the moment about the selected neutral axis. We also find the curvature and moment about the axis orthogonal to the neutral axis `chi_z` and `m_z`, and the axial strain at the level of the neutral axis `eps_axial`. All these attributes are stored as arrays, ready for visualization or further processing.

:::::{tip}

Use your favourite plotting library to visualize the results from the {class}`GenericSectionCalculator <structuralcodes.sections.GenericSectionCalculator>`. The code below shows how to plot the moment-curvature relation in the figure [below](#fig-usage-moment-curvature) with [Matplotlib](https://matplotlib.org/). Notice how we are plotting the negative values of the curvatures and moments due to the sign convention.

(code-usage-visualize-moment-curvature)=
::::{dropdown-syntax}
:::{literalinclude} ../../_example_code/usage_create_surface_geometries_with_reinforcement.py
:lines: 3, 101-111
:caption: Visualize the moment-curvature relation.
:::
::::

(fig-usage-moment-curvature)=
:::{figure} moment_curvature.png
:width: 75%

The moment-curvature relation computed with the code [above](#code-usage-generic-section).
:::

:::::