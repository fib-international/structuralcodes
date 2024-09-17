from shapely import Polygon

from structuralcodes import codes, materials
from structuralcodes.geometry import SurfaceGeometry
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
sec.geometry = sec.geometry.translate(-175, -250)
### CREATE THE SECTION ########


### TEST calculate_strain_profile_biaxial  ########
n = -50 * 1e3
my = -210 * 1e6  # -150
mz = 47 * 1e6  # -100


res = sec.section_calculator.calculate_strain_profile(n, my, mz)
print(
    'eps',
    round(res[0] * 1e3, 2),
    '  chiy',
    round(res[1] * 1e6, 2),
    '  chiz',
    round(res[2] * 1e6, 2),
)
forces = sec.section_calculator.integrate_strain_profile(res)
print(
    'N',
    round(forces[0] / 1e3, 2),
    '  My',
    round(forces[1] / 1e6, 2),
    '  Mz',
    round(forces[2] / 1e6, 2),
)
