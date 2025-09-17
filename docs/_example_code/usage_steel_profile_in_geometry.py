"""Example code for using a steel profile in a geometry."""

from structuralcodes.geometry import SurfaceGeometry
from structuralcodes.geometry.profiles import IPE
from structuralcodes.materials.basic import ElasticPlasticMaterial
from structuralcodes.sections import GenericSection

# Create a profile
ipe100 = IPE('IPE100')

# Create an elastic perfect plastic material
steel = ElasticPlasticMaterial(E=210000, fy=275 / 1.05, density=7850)

# Create a geometry and a section
geom = SurfaceGeometry(poly=ipe100.polygon, material=steel)
section = GenericSection(geometry=geom)

# Use the section calculator to calculate the bending strength
bending_strength = section.section_calculator.calculate_bending_strength()
