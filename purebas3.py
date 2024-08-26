import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata
from shapely.geometry import Point, Polygon
from shapely.geometry.polygon import orient


def create_surface_from_points(x, y, z):
    # Crear una cuadrícula de valores para la superficie
    grid_x, grid_y = np.mgrid[min(x) : max(x) : 100j, min(y) : max(y) : 100j]

    # Interpolar los puntos z en la cuadrícula
    grid_z = griddata((x, y), z, (grid_x, grid_y), method='cubic')

    return grid_x, grid_y, grid_z


# Ejemplo de puntos desordenados
x = np.random.random(100)
y = np.random.random(100)
z = np.sin(x**2 + y**2)

# Crear la superficie a partir de los puntos
grid_x, grid_y, grid_z = create_surface_from_points(x, y, z)

# Visualizar la superficie
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(grid_x, grid_y, grid_z, cmap='viridis')

# También se pueden graficar los puntos originales
ax.scatter(x, y, z, color='r')

plt.show()


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

# Visualizar los puntos
import matplotlib.pyplot as plt

plt.figure(figsize=(8, 8))
plt.plot(*poly.exterior.xy, 'k-')  # Contorno del polígono
plt.scatter(
    border_points[:, 0],
    border_points[:, 1],
    color='blue',
    label='Border Points',
)
plt.scatter(
    interior_points[:, 0],
    interior_points[:, 1],
    color='red',
    label='Interior Points',
)
plt.legend()
plt.show()

x = all_points[:, 0]
y = all_points[:, 1]
z = np.sin(x**2 + y**2)
# Crear la superficie a partir de los puntos
grid_x, grid_y, grid_z = create_surface_from_points(x, y, z)

"""# Visualizar la superficie
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(grid_x, grid_y, grid_z, cmap='viridis')

# También se pueden graficar los puntos originales
ax.scatter(x, y, z, color='r')

plt.show()"""


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
