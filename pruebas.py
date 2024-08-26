from shapely import Polygon

from structuralcodes import codes, materials
from structuralcodes.geometry import SurfaceGeometry
from structuralcodes.plots.section_plots import draw_section_response3D
from structuralcodes.sections._generic import GenericSection
from structuralcodes.sections._reinforcement import add_reinforcement_line

# from structuralcodes.plots.section_plots import draw_section_response,draw_section,get_stress_point


codes.set_design_code(design_code='ec2_2004')
# Create materials
concrete = materials.concrete.create_concrete(fck=25)
reinforcemnet = materials.reinforcement.create_reinforcement(
    fyk=500,
    Es=200000,
    density=7850,
    ftk=500,
    epsuk=0.07,
)
# Create section
poly = Polygon(((0, 0), (350, 0), (350, 500), (0, 500)))
geo = SurfaceGeometry(poly, concrete)
geo = add_reinforcement_line(geo, (50, 50), (300, 50), 12, reinforcemnet, n=3)
geo = add_reinforcement_line(
    geo, (50, 440), (300, 440), 20, reinforcemnet, n=6
)
sec = GenericSection(geo)


n = 0 * 1e3
my = 200 * 1e6
mz = 100 * 1e6
res = sec.section_calculator.calculate_strain_profile(n, my, mz)
print(res)

draw_section_response3D(sec, res[0], res[1], 0)

"""
#######################33
from shapely.geometry import Polygon

# Define un polígono simple
polygon1 = Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
polygon2 = Polygon([(1, 1), (2, 1), (2, 2), (1, 2)])
def get_polygon_coordinates(polygon):
    x, y = polygon.exterior.xy
    return x, y

x1, y1 = get_polygon_coordinates(polygon1)
x2, y2 = get_polygon_coordinates(polygon2)
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# Crear una figura y un eje 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Función para graficar un polígono en 3D
def plot_polygon_3d(ax, x, y, z=0):
    verts = [list(zip(x, y, [z]*len(x)))]
    poly = Poly3DCollection(verts, alpha=0.5, edgecolor='r')
    ax.add_collection3d(poly)

# Graficar los polígonos
plot_polygon_3d(ax, x1, y1, z=0)  # Graficar el primer polígono en z=0
plot_polygon_3d(ax, x2, y2, z=1)  # Graficar el segundo polígono en z=1

# Configuración de los ejes
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_xlim([0, 3])
ax.set_ylim([0, 3])
ax.set_zlim([0, 2])

plt.show()
"""
######################### s
######################### s
######################### s
######################### superficie desde puntos


######################### s
######################### s
######################### s
########################### s####################### s
######################### s
#################### sacar puntos de un polygon shapely

#########################################3
