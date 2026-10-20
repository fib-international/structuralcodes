"""Example code for creating surface geometries."""

from structuralcodes.geometry import RectangularGeometry
from structuralcodes.materials.concrete import ConcreteEC2_2004

# Define parameters
fck = 45
width_web = 250
width_flange = 1000
height_web = 650
height_flange = 150

# Create material
concrete = ConcreteEC2_2004(fck=fck)

# Create surface geometries
web_geom = RectangularGeometry(
    width=width_web,
    height=height_web,
    material=concrete,
    origin=(0, height_web / 2),
)
flange_geom = RectangularGeometry(
    width=width_flange,
    height=height_flange,
    material=concrete,
    origin=(0, height_web + height_flange / 2),
)

# Add surface geometries to create a T-shaped geometry
t_geom = web_geom + flange_geom
