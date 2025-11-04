"""A collection of functions related to reinforced concrete sections."""

import typing as t

import structuralcodes.core._section_results as s_res
from structuralcodes.geometry import CompoundGeometry, SurfaceGeometry
from structuralcodes.materials.basic import ElasticMaterial, GenericMaterial
from structuralcodes.materials.constitutive_laws import UserDefined
from structuralcodes.sections import GenericSection
from structuralcodes.sections.section_integrators import FiberIntegrator


def _create_surface_geometries(polygons_list, material):
    """Process shapely polygons to SurfaceGeometries."""
    # Create an empty list to store SurfaceGeometry objects
    surface_geometries = []

    # Iterate over the list of polygons and create SurfaceGeometry for each
    for polygon in polygons_list:
        # Create a new SurfaceGeometry for the current polygon
        surface_geometry = SurfaceGeometry(polygon, material)
        # Add the new SurfaceGeometry to the list
        surface_geometries.append(surface_geometry)

    return CompoundGeometry(surface_geometries)


def calculate_elastic_cracked_properties(
    section: GenericSection,
    theta: float = 0,
    return_cracked_section: bool = False,
) -> t.Union[
    s_res.SectionProperties, t.Tuple[s_res.SectionProperties, CompoundGeometry]
]:
    """Calculates the cracked section properties of a reinforced concrete
    section. (GenericSection). Materials in surface geometries and point
    geometries are elastic-linear in order to make the cracking properties
    independent of the stress state. Tension in all surface geometries is
    neglected.

    Args:
        section (GenericSection): The section to use as basis for the
            calculation.
        theta (float): Angle of the neutral axis to the horizontal in radians.
            theta=0 implies upper compression block.
        return_cracked_section (bool): If true, returns also the cracked
            section.

    Returns:
        t.Union[SectionProperties, t.Tuple[SectionProperties, GenericSection]]:
        SectionProperties data of a cracked section (i.e cracked section
        properties) and (if return_cracked_section is True) the cracked
        section.
    """
    if not section.geometry.reinforced_concrete:
        return None

    # Get the integrator properties
    if isinstance(section.section_calculator.integrator, FiberIntegrator):
        mesh_size = getattr(
            section.section_calculator.integrator, 'mesh_size', 0.01
        )
        integrator = 'fiber'
    else:
        mesh_size = None
        integrator = 'marin'

    rotated_geometry = section.geometry.rotate(-theta)

    processed_geometries = []
    for geo in rotated_geometry.geometries:
        E = geo.material.constitutive_law.get_tangent(eps=0)
        density = geo.material.density
        if geo.concrete:
            # It is concrete we only give compression behavior
            elastic_law = UserDefined([-100, 0], [-100 * E, 0])
            elastic_law = GenericMaterial(
                density=density,
                constitutive_law=elastic_law,
                name='elastic concrete',
            )
        else:
            # It is not concrete we give tension and compression behavior
            elastic_law = ElasticMaterial(
                E=E, density=density, name='elastic material'
            )
        processed_geometries.append(
            geo.from_geometry(geo, new_material=elastic_law)
        )

    for pg in rotated_geometry.point_geometries:
        Es = pg.material.constitutive_law.get_tangent(eps=0)
        density = pg.material.density
        elastic_steel = ElasticMaterial(
            E=Es, density=density, name='elastic steel'
        )
        processed_geometries.append(
            pg.from_geometry(geo=pg, new_material=elastic_steel)
        )

    # Create a new temporary section
    geo = CompoundGeometry(geometries=processed_geometries)
    elastic_section = GenericSection(
        geo, integrator=integrator, mesh_size=mesh_size
    )
    rotated_geometry = elastic_section.geometry

    curv = -1e-5  # Any curvature should return the same mechanical properties.
    # Find the equilibrium with fixed curvature
    eps = elastic_section.section_calculator.find_equilibrium_fixed_curvature(
        rotated_geometry, 0, curv, 0
    )[0]
    z_na = -eps / curv  # distance to neutral fibre

    # Cutting concrete geometries and retaining the compressed block
    cut_geom = None
    for part in rotated_geometry.geometries:
        upper_div, lower_div = part.split(((0, z_na), 0))
        # Convert to SurfaceGeometry
        subpart = _create_surface_geometries(upper_div, part.material)
        if cut_geom is None:
            cut_geom = subpart
        else:
            cut_geom += subpart

    # Add reinforcement geometries
    for reinf in rotated_geometry.point_geometries:
        cut_geom += reinf

    # return geometry to original rotation
    cracked_geom = cut_geom.rotate(theta)

    # Define the cracked section
    cracked_sec = GenericSection(
        cracked_geom, integrator=integrator, mesh_size=mesh_size
    )

    cracked_prop = cracked_sec.gross_properties

    return (
        (cracked_prop, cracked_geom)
        if return_cracked_section
        else cracked_prop
    )
