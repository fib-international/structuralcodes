(usage-sections)=
# Sections

## General

In StructuralCodes, a _Section_ is responsible for connecting a [geometry](#usage-geometries) object with an API for performing calculations on the geometry.

The section contains a _SectionCalculator_ which contains the API for doing calculations on the geometry, including for example lower level methods for integrating a strain profile to stress resultants, or calculating a strain profile based on stress resultants, or higher level methods like calculating a moment-curvature relation, or an interaction diagram between moments and axial force.

The SectionCalculator uses a _SectionIntegrator_ to integrate the strain response on a geometry.

:::{tip}
See the theory reference for a guide on the [sign convention](#theory-sign-convention) used in StructuralCodes.
:::

:::{tip}
The theory reference provides an overview of the theory behind the [section calculator](#theory-section-calculator) and the [section integrators](#theory-section-integrators).
:::

## The beam section

The {class}`BeamSection <structuralcodes.sections.BeamSection>` takes a {class}`SurfaceGeometry <structuralcodes.geometry.SurfaceGeometry>` or a {class}`CompoundGeometry <structuralcodes.geometry.CompoundGeometry>` as input, and is capable of calculating the response of an arbitrarily shaped geometry with arbitrary reinforcement layout, subject to stresses in the direction of the beam axis.

In the example [below](#code-usage-beam-section), we continue the example with the T-shaped geometry.

(usage-sections-gross-properties-tip)=
:::{tip}
If you are looking for the gross properties of the section, these are available at the `.gross_properties` property on the `Section` object ✅
:::

Since the axial force and the moments are all relative to the origin, we start by translating the geometry such that the centroid is alligned with the origin.

Notice how we can use {py:meth}`.calculate_bending_strength() <structuralcodes.sections.BeamSectionCalculator.calculate_bending_strength>` and {py:meth}`.calculate_limit_axial_load() <structuralcodes.sections.BeamSectionCalculator.calculate_limit_axial_load>` to calculate the bending strength and the limit axial loads in tension and compression.

Furthermore, we have the following methods:

{py:meth}`.integrate_strain_profile() <structuralcodes.sections.BeamSectionCalculator.integrate_strain_profile>`
: Calculate the stress resultants for a given strain profile.

{py:meth}`.calculate_strain_profile() <structuralcodes.sections.BeamSectionCalculator.calculate_strain_profile>`
: Calculate the strain profile for given stress resultants.

{py:meth}`.calculate_moment_curvature() <structuralcodes.sections.BeamSectionCalculator.calculate_moment_curvature>`
: Calculate the moment-curvature relation.

{py:meth}`.calculate_nm_interaction_domain() <structuralcodes.sections.BeamSectionCalculator.calculate_nm_interaction_domain>`
: Calculate the interaction domain between axial load and bending moment.

{py:meth}`.calculate_nmm_interaction_domain() <structuralcodes.sections.BeamSectionCalculator.calculate_nmm_interaction_domain>`
: Calculate the interaction domain between axial load and biaxial bending.

{py:meth}`.calculate_mm_interaction_domain() <structuralcodes.sections.BeamSectionCalculator.calculate_mm_interaction_domain>`
: Calculate the interaction domain between biaxial bending for a given axial load.

See the {class}`BeamSectionCalculator <structuralcodes.sections.BeamSectionCalculator>` for a complete list.

(usage-sections-complete-nm-note)=
:::{note}
Notice that a call to {py:meth}`.calculate_nm_interaction_domain() <structuralcodes.sections.BeamSectionCalculator.calculate_nm_interaction_domain>` returns the interaction domain for a negative bending moment, i.e. a bending moment that gives compression at the top of the cross section according to the [sign convention](#theory-sign-convention). To obtain the interaction domain for the positive bending moment, the neutral axis for the calculation should be rotated an angle {math}`\theta = \pi`, i.e. the method should be called with the keyword argument `theta = np.pi` or `theta = math.pi`. Alternatively, to return the complete domain, the method could be called with the keyword argument `complete_domain = True`.
:::

(code-usage-beam-section)=
::::{dropdown-syntax}
:::{literalinclude} ../../_example_code/usage_create_surface_geometries_with_reinforcement.py
:lines: 12-13, 72-90
:caption: Create a section and perform calculations.
:::
::::

Notice how for example {py:meth}`.calculate_moment_curvature() <structuralcodes.sections.BeamSectionCalculator.calculate_moment_curvature>` returns a custom dataclass of type {class}`MomentCurvatureResults <structuralcodes.core._section_results.MomentCurvatureResults>`. If we inspect this class further, we find among its attributes `chi_y` and `m_y`. These are the curvature and the moment about the selected global axis of the section. We also find the curvature and moment about the axis orthogonal to the global axis `chi_z` and `m_z`, and the axial strain at the level of the global axis `eps_axial`. All these attributes are stored as arrays, ready for visualization or further processing.

:::::{tip}

Use your favourite plotting library to visualize the results from the {class}`BeamSectionCalculator <structuralcodes.sections.BeamSectionCalculator>`. The code below shows how to plot the moment-curvature relation in the figure [below](#fig-usage-moment-curvature) with [Matplotlib](https://matplotlib.org/). Notice how we are plotting the negative values of the curvatures and moments due to the sign convention.

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

The moment-curvature relation computed with the code [above](#code-usage-beam-section).
:::

:::::

## Inspecting further with results objects

The results objects contain methods for inspecting and processing further the detailed results on the section after a calculation is performed.

For instance, calling the method {py:meth}`.calculate_bending_strength() <structuralcodes.sections.BeamSectionCalculator.calculate_bending_strength>` returns an object of class {class}`UltimateBendingMomentResults <structuralcodes.core._section_results.UltimateBendingMomentResults>`.

The object contains the following methods:

{py:meth}`.create_detailed_result() <structuralcodes.core._section_results.UltimateBendingMomentResults.create_detailed_result>`
: Create an object of class {class}`SectionDetailedResultState <structuralcodes.core._section_results.SectionDetailedResultState>` with the datastructure storing strain and stress fields across the section.

{py:meth}`.get_point_strain() <structuralcodes.core._section_results.UltimateBendingMomentResults.get_point_strain>`
: Get the strain at a point.

{py:meth}`.get_point_stress() <structuralcodes.core._section_results.UltimateBendingMomentResults.get_point_stress>`
: Get the stress at a point.

In the example [below](#code-usage-beam-section-result-strength), we continue the example with the T-shaped geometry and we create the detailed_result.

(code-usage-beam-section-result-strength)=
::::{dropdown-syntax}
:::{literalinclude} ../../_example_code/usage_bending_strength_results.py
:lines: 87-88, 100-101
:caption: Create the detailed result datastructure.
:::
::::

If we inspect this result object, we find among its attributes `detailed_result` that has two useful datastructures: one for surface geometries and one for point geometries. The datastructure for surface geometries contains all points samples from the geometries. These data structures can be easily input for creating a DataFrame object permitting further processing. For instance in the example [below](#code-usage-beam-section-result-strength-plot), we use pandas plotting abilities to create a visualization of strain and stress fields.

(code-usage-beam-section-result-strength-plot)=
::::{dropdown-syntax}
:::{literalinclude} ../../_example_code/usage_bending_strength_results.py
:lines: 103-141
:caption: Create a plot of strain and stress fields for bending strength.
:::
::::

(fig-usage-strain-field)=
:::{figure} strain_field.png
:width: 75%

The strain field plotted with the code [above](#code-usage-beam-section-result-strength-plot).
:::

(fig-usage-stress-field)=
:::{figure} stress_field.png
:width: 75%

The stress field plotted with the code [above](#code-usage-beam-section-result-strength-plot).
:::

If we only want to know the strain and/or stress at a specific point we can directly use methods {py:meth}`.get_point_strain() <structuralcodes.core._section_results.UltimateBendingMomentResults.get_point_strain>` or {py:meth}`.get_point_stress() <structuralcodes.core._section_results.UltimateBendingMomentResults.get_point_stress>` without the need of creating the detailed result. For instance with the code [below](#code-usage-beam-section-result-strength-point-stress) we query the stress at a point for any geometry whose `group_label` matches the pattern `"reinforcement"`. Note that a pattern can be given with wildcards (like `'*'` or `'?'`); for more information refer to {py:meth}`.get_point_stress() <structuralcodes.core._section_results.UltimateBendingMomentResults.get_point_stress>` and {py:meth}`.group_filter()<structuralcodes.geometry.CompoundGeometry.group_filter>`.

(code-usage-beam-section-result-strength-point-stress)=
::::{dropdown-syntax}
:::{literalinclude} ../../_example_code/usage_bending_strength_results.py
:lines: 90-98
:caption: Get the stress at the coordinates of one reinforcement
:::
::::

When dealing with analysis that do not involve a single section state, like for instance a Moment-Curvature analysis, the `detailed_result` object is created for the `current_step`. The `current_step` is initialized at 0. Then methods are available for navigating through steps:

{py:meth}`.next_step() <structuralcodes.core._section_results.MomentCurvatureResults.next_step>`
: Navigate to the next step.

{py:meth}`.previous_step() <structuralcodes.core._section_results.MomentCurvatureResults.previous_step>`
: Navigate to the previous step.

{py:meth}`.set_step() <structuralcodes.core._section_results.MomentCurvatureResults.set_step>`
: Navigate to a specific step.

If we are not interested in getting the detailed results but we want to get the stresses and/or strains at specific points the same methods {py:meth}`.get_point_stress() <structuralcodes.core._section_results.MomentCurvatureResults.get_point_stress>` and {py:meth}`.get_point_strain() <structuralcodes.core._section_results.MomentCurvatureResults.get_point_strain>` can be used. In this case the return will be an array of stresses (or strains) for the point through the moment-curvature analysis. For instance with the code [below](#code-usage-beam-section-result-mcurv-point-stress) it is possible to plot the stress for one rebar at a given position for the increasing curvature.

(code-usage-beam-section-result-mcurv-point-stress)=
::::{dropdown-syntax}
:::{literalinclude} ../../_example_code/usage_bending_strength_results.py
:lines: 202-226
:caption: Plot the stress at reinforcement position for increasing curvature.
:::
::::

(fig-usage-stress-curvature)=
:::{figure} stress_at_rebar.png
:width: 75%

The stress for increasing curvature plotted with the code [above](#code-usage-beam-section-result-mcurv-point-stress).
:::

Using in a more advanced way matplotlib, it is possible to create an animation (like for instance a *gif* file) showing the stress field for the steps of the moment-curvature analysis. The code [below](code-usage-beam-section-result-mcurv-anim) create an animation with the moment curvature diagram on the left and the cross section stress field on the right.

(code-usage-beam-section-result-mcurv-anim)=
::::{dropdown-syntax}
:::{literalinclude} ../../_example_code/usage_bending_strength_results.py
:lines: 228-359
:caption: Generate an animation for the moment-curvature section response.
:::
::::

(fig-usage-stress-animation)=
:::{figure} moment_curvature_stress_animation.gif
:width: 75%

The animation of the stress field with increasing curvature created with the code [above](#code-usage-beam-section-result-mcurv-anim).
:::