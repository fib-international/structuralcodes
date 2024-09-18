from shapely import Polygon
from shapely.geometry import Polygon

from structuralcodes import codes, materials
from structuralcodes.geometry import SurfaceGeometry
from structuralcodes.sections._generic import GenericSection
from structuralcodes.sections._reinforcement import add_reinforcement_line

# from structuralcodes.plots.section_plots import draw_section_response,draw_section

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

# print(sec.cracked_properties)
print(sec.cracked_properties.i_yy_1 / sec.gross_properties.i11)


"""print('----------')
res = cracked_properties(sec, 0, True)
print(res)
print('----------')
res = cracked_properties(sec, 0, False)
print(res)
print('----------')"""
