import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata
from scipy.spatial import Delaunay
from shapely.geometry import Point, Polygon
from shapely.geometry.polygon import orient


def create_surface_from_points(x, y, z):
    # Crear una cuadrícula de valores para la superficie
    grid_x, grid_y = np.mgrid[min(x) : max(x) : 100j, min(y) : max(y) : 100j]

    # Interpolar los puntos z en la cuadrícula
    grid_z = griddata((x, y), z, (grid_x, grid_y), method='cubic')

    return grid_x, grid_y, grid_z


def get_border_points(polygon, num_points=100):
    # Asegurarse de que el polígono esté orientado correctamente
    polygon = orient(polygon)
    exterior_coords = np.array(polygon.exterior.coords)

    # Calcular distancias acumuladas a lo largo del borde del polígono
    distances = np.sqrt(np.sum(np.diff(exterior_coords, axis=0) ** 2, axis=1))
    cumulative_distances = np.insert(np.cumsum(distances), 0, 0)

    # Crear una interpolación lineal de los puntos del borde
    total_length = cumulative_distances[-1]
    interpolated_distances = np.linspace(0, total_length, num_points)

    border_points = np.empty((num_points, 2))
    for i, d in enumerate(interpolated_distances):
        index = np.searchsorted(cumulative_distances, d) - 1
        t = (d - cumulative_distances[index]) / distances[index]
        border_points[i] = (1 - t) * exterior_coords[
            index
        ] + t * exterior_coords[index + 1]

    return border_points


def get_interior_points(polygon, num_points=100):
    min_x, min_y, max_x, max_y = polygon.bounds
    points = []
    while len(points) < num_points:
        random_point = Point(
            np.random.uniform(min_x, max_x), np.random.uniform(min_y, max_y)
        )
        if polygon.contains(random_point):
            points.append(random_point)
    return np.array([point.coords[0] for point in points])


# Ejemplo de uso
poly = Polygon([(0, 0), (2, 0), (1.5, 1), (2, 2), (0, 2)])

border_points = get_border_points(poly, num_points=50)
interior_points = get_interior_points(poly, num_points=500)

# Combinar puntos de borde e interior
all_points = np.vstack((border_points, interior_points))


x = all_points[:, 0]
y = all_points[:, 1]
z = np.sin(x**2 + y**2)
# Crear la superficie a partir de los puntos
grid_x, grid_y, grid_z = create_surface_from_points(x, y, z)
print((grid_x))
print((x))
# Visualizar la superficie
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(grid_x, grid_y, grid_z, cmap='viridis')

# También se pueden graficar los puntos originales
ax.scatter(x, y, z, color='r')

plt.show()


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Graficar los puntos en 3D
ax.scatter(x, y, z, color='blue', s=10)

# Configuración de los ejes
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

# Ajustar los límites de los ejes
ax.set_xlim([x.min(), x.max()])
ax.set_ylim([y.min(), y.max()])
ax.set_zlim([z.min(), z.max()])
plt.show()


import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from shapely.geometry import Polygon

# Crear la malla Delaunay en 3D
points = np.vstack((x, y)).T
tri = Delaunay(points)

x = np.array([0, 2, 1.5, 2, 0])
y = np.array([0, 0, 1, 2, 2])
z = np.array([0, 1, 2, 1, 0])

# Graficar los puntos y la malla filtrada
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, y, z, color='r')

# Recorrer cada triángulo de la malla Delaunay
for simplex in tri.simplices:
    # Calcular el punto central del triángulo en 2D
    centroid_2d = np.mean(points[simplex], axis=0)

    # Verificar si el punto central está dentro del polígono
    if poly.contains(Point(centroid_2d)):
        vertices_3d = np.array([x[simplex], y[simplex], z[simplex]]).T
        poly = Poly3DCollection(
            [vertices_3d], color='b', alpha=0.5, edgecolor='k'
        )
        ax.add_collection3d(poly)

# Ajustar los límites de la gráfica
ax.set_xlim([x.min(), x.max()])
ax.set_ylim([y.min(), y.max()])
ax.set_zlim([z.min(), z.max()])

# Mostrar la gráfica
plt.show()
