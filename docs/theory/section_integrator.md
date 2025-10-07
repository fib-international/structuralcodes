(theory-section-integrators)=
# Section integrators
## General concept
Section integrators compute properties such as stiffness, stress distribution, or strain compatibility over a section. They operate on the geometry and material definitions to derive these results.

## Marin integrator
Handles sections with predefined shapes and material properties, ideal for standard geometries. It applies assumptions about material distribution to simplify calculations.

## Fiber integrator
Uses a discretized approach where the section is divided into fibers, each representing a material point. This discretization is performed by a triangulation of the geometry.

**Figure Placeholder**: Schematic comparison of Marin integrator and Fiber integrator approaches.
