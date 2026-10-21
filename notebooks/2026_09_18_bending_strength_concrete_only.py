from shapely import Polygon

from structuralcodes import set_design_code
from structuralcodes.core.base import ConstitutiveLaw
from structuralcodes.geometry import SurfaceGeometry
from structuralcodes.materials.concrete import create_concrete
from structuralcodes.materials.constitutive_laws import (
    ElasticFragileTension,
    ParabolaRectangle,
    Parallel,
)
from structuralcodes.sections import BeamSection

Ec = 36000
fct = 2.7
fck = 45
set_design_code('ec2_2004')
compression: ConstitutiveLaw = ParabolaRectangle(fc=fck)
tension: ConstitutiveLaw = ElasticFragileTension(fct=fct, Ec=Ec)
compression_tension: ConstitutiveLaw = Parallel(
    constitutive_laws=[compression, tension]
)

concrete_c = create_concrete(fck=fck, constitutive_law=compression)
concrete_t = create_concrete(fck=fck, constitutive_law=tension)

concrete_ct = create_concrete(fck=fck, constitutive_law=compression_tension)


# Create a rectangular geometry
width = 10
height = 10
polygon = Polygon(
    [
        (-width / 2, -height / 2),
        (width / 2, -height / 2),
        (width / 2, height / 2),
        (-width / 2, height / 2),
    ]
)

geometry_c = SurfaceGeometry(poly=polygon, material=concrete_c)
geometry_t = SurfaceGeometry(poly=polygon, material=concrete_t)
geometry_ct = SurfaceGeometry(poly=polygon, material=concrete_ct)

# Create section
section_c = BeamSection(geometry_c, integrator='fiber')
section_t = BeamSection(geometry_t, integrator='fiber')
section_ct = BeamSection(geometry_ct, integrator='fiber')
print('Results :')
section_c.section_calculator.calculate_limit_axial_load()
print(
    f'Section C bending strength m_y: {section_c.section_calculator.calculate_bending_strength().m_y}'
)
print(
    f'Section T bending strength m_y: {section_t.section_calculator.calculate_bending_strength().m_y}'
)
print(
    f'Section CT bending strength m_y: {section_ct.section_calculator.calculate_bending_strength().m_y}'
)
print(section_c.section_calculator.calculate_bending_strength().m_y)