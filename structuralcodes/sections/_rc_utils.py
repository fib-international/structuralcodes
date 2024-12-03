from structuralcodes.geometry import CompoundGeometry, SurfaceGeometry
from structuralcodes.materials.constitutive_laws import Elastic, UserDefined
from structuralcodes.sections import GenericSection


def calculate_elastic_cracked_properties(section: GenericSection, theta=0):
    """Calculates the cracked section properties of a reinforced concrete
    section.  (GenericSection). Materials in surface geometries and point
    geometries are elastic-linear  in order to make the cracking properties
    independent of the stress state.  Tension in all surface geometries is
    neglected.

    Args:
        section: GenericSection
        theta: Angle of the neutral axis to the horizontal. theta=0 implies
               upper compression block.

    Returns:
        cracked_prop : GrossProperties data of cracked_sec (i.e cracked
                       section properties)
        cracked_sec: GenericSection removing tension areas in geometries.
    """

    def create_surface_geometries(polygons_list, material):
        """Process shapely polygons to SurfaceGeometries."""
        # Create an empty list to store SurfaceGeometry objects
        surface_geometries = []

        # Iterate over the list of polygons and create SurfaceGeometry for each
        for polygon in polygons_list:
            # Create a new SurfaceGeometry for the current polygon
            surface_geometry = SurfaceGeometry(polygon, material)
            # Add the new SurfaceGeometry to the list
            surface_geometries.append(surface_geometry)

        return surface_geometries

    if not section.geometry.reinforced_concrete:
        return None

    rotated_geometry = section.geometry.rotate(theta)

    for geo in rotated_geometry.geometries:
        Ec = geo.material.get_tangent(eps=0)
        elastic_concrete = UserDefined([-100, 0], [-100 * Ec, 0])
        geo._material = elastic_concrete

    for pg in rotated_geometry.point_geometries:
        Es = pg.material.get_tangent(eps=0)
        elastic_steel = Elastic(Es, 'elastic steel')
        pg._material = elastic_steel

    curv = -1e-5  # Any curvature should return the same mechanical properties.
    # Find the equilibrium with fixed curvature
    eps = section.section_calculator.find_equilibrium_fixed_curvature(
        rotated_geometry, 0, curv, 0
    )[0]
    z_na = -eps / curv  # distance to neutral fibre

    # Cutting concrete geometries and retaining the compressed block
    cut_geom = CompoundGeometry(None, None)
    for i, part in enumerate(rotated_geometry.geometries):
        min_x, max_x, min_y, max_y = part.calculate_extents()
        upper_div, lower_div = part.split(((min_x, z_na), 0))
        # Convert to SurfaceGeometry
        subpart_sg = create_surface_geometries(upper_div, part.material)
        if i == 0:
            cut_geom = CompoundGeometry(subpart_sg)
        else:
            cut_geom += subpart_sg

    # Add reinforcement geometries
    for reinf in rotated_geometry.point_geometries:
        cut_geom += reinf

    # return geoemtry to original rotation
    cracked_geom = cut_geom.rotate(-theta)

    # Define the cracked section
    cracked_sec = GenericSection(cracked_geom)
    cracked_prop = cracked_sec.gross_properties

    print(f'{cracked_prop:cracked}')
    return cracked_prop, cracked_sec
