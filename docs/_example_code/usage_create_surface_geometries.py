"""Example code for creating surface geometries."""

from shapely import Polygon

from structuralcodes.geometry import SurfaceGeometry
from structuralcodes.materials.concrete import ConcreteEC2_2004

# Define parameters
fck = 45
width_web = 250
width_flange = 1000
height_web = 650
height_flange = 150

# Create material
concrete = ConcreteEC2_2004(fck=fck)

# Create polygons
polygon_web = Polygon(
    (
        (-width_web / 2, 0),
        (width_web / 2, 0),
        (width_web / 2, height_web),
        (-width_web / 2, height_web),
    )
)
polygon_flange = Polygon(
    (
        (-width_flange / 2, height_web),
        (width_flange / 2, height_web),
        (width_flange / 2, height_web + height_flange),
        (-width_flange / 2, height_web + height_flange),
    ),
)

# Create surface geometries
web_geom = SurfaceGeometry(poly=polygon_web, material=concrete)
flange_geom = SurfaceGeometry(poly=polygon_flange, material=concrete)

# Add surface geometries to create a T-shaped geometry
t_geom = web_geom + flange_geom
