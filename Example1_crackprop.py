import math

from shapely import Polygon
from shapely.geometry import Polygon

from structuralcodes import codes, materials
from structuralcodes.geometry import CompoundGeometry, SurfaceGeometry
from structuralcodes.plots.section_plots import draw_section
from structuralcodes.sections._generic import GenericSection
from structuralcodes.sections._reinforcement import add_reinforcement_line

# from structuralcodes.plots.section_plots import draw_section_response,draw_section,get_stress_point

### CREATE THE SECTION ########
codes.set_design_code(design_code='ec2_2004')
concrete = materials.concrete.create_concrete(fck=25)
reinforcemnet = materials.reinforcement.create_reinforcement(
    fyk=500,
    Es=200000,
    density=7850,
    ftk=500,
    epsuk=0.07,
)
poly = Polygon(((0, 0), (350, 0), (350, 500), (0, 500)))
geo = SurfaceGeometry(poly, concrete)
geo = add_reinforcement_line(
    geo, (50, 450), (300, 450), 12, reinforcemnet, n=3
)
geo = add_reinforcement_line(geo, (50, 50), (300, 50), 20, reinforcemnet, n=6)
sec = GenericSection(
    geo,
)
# sec.geometry = sec.geometry.translate(-sec.gross_properties.cy, -sec.gross_properties.cz)
print(sec.geometry.reinforced_concrete)
sec.geometry = sec.geometry.translate(-175, -250)
print(sec.geometry.reinforced_concrete)
sec.geometry._reinforced_concrete = (
    True  # TODO ISSUE. Traslation modifies reinforced_concrete attribure
)
### CREATE THE SECTION ########

print(sec.cracked_properties)
# print(sec.cracked_properties.i_yy_2 / sec.gross_properties.i11)


def cracked_properties(section: GenericSection, n, neg_bending=False):
    """TODO."""
    if not section.geometry.reinforced_concrete:
        return None

    if neg_bending:
        r = section.section_calculator.calculate_bending_strength(theta=0, n=n)
        curv = -r.chi_y / 10  # small curvature to get n.a
    else:
        r = section.section_calculator.calculate_bending_strength(
            theta=math.pi, n=n
        )
        curv = r.chi_y / 10  # small curvature to get n.a

    eps = section.section_calculator.find_equilibrium_fixed_curvature(
        section.geometry, n, curv, 0
    )[0]
    y = -eps / curv  # distance to neutral fibre

    def create_surface_geometries(polygons_list, material):
        """Process shapely polygons to SurfaceGeometries."""
        # Create an empty list to store SurfaceGeometry objects
        surface_geometries = []

        # Iterate over the list of polygons and create a SurfaceGeometry for each one
        for polygon in polygons_list:
            # Create a new SurfaceGeometry for the current polygon
            surface_geometry = SurfaceGeometry(polygon, material)
            # Add the new SurfaceGeometry to the list
            surface_geometries.append(surface_geometry)

        return surface_geometries

    cut_geom = CompoundGeometry(None, None)
    for i, part in enumerate(section.geometry.geometries):
        min_x, max_x, min_y, max_y = part.calculate_extents()
        above_div, below_div = part.split(((min_x, y), 0))
        if neg_bending:
            subpart_poly = below_div
        else:
            subpart_poly = above_div
        # convert to SurfaceGeometry
        subpart_sg = create_surface_geometries(subpart_poly, part.material)
        if True:
            if i == 0:
                cut_geom = CompoundGeometry(subpart_sg)
            else:
                cut_geom += subpart_sg
    for reinf in section.geometry.point_geometries:
        cut_geom += reinf
    cracked_sec = GenericSection(cut_geom)
    draw_section(cracked_sec)
    # draw_section_response(cracked_sec, eps, curv)
    # cracked_sec = cut_generic_section_upper_part(section, x)
    return cracked_sec.gross_properties


"""print('----------')
res = cracked_properties(sec, 0, True)
print(res)
print('----------')
res = cracked_properties(sec, 0, False)
print(res)
print('----------')"""
