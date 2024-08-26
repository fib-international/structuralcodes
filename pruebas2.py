import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import Point, Polygon
from shapely.geometry.polygon import orient

# Definir el polígono
polygon = Polygon([(0, 0), (2, 0), (2, 2), (0, 2)])

# Asegurarse de que el polígono esté orientado correctamente
polygon = orient(polygon)


# Función que transforma Y y Z para obtener X
def func(y, z):
    return np.sin(y) + np.cos(z)


# Generar una cuadrícula de puntos dentro del polígono
num_points = 500
min_x, min_y, max_x, max_y = polygon.bounds
points = []

while len(points) < num_points:
    random_point = Point(
        np.random.uniform(min_x, max_x), np.random.uniform(min_y, max_y)
    )
    if polygon.contains(random_point):
        points.append(random_point)

# Extraer las coordenadas Y y Z de los puntos
coords = np.array([point.coords[0] for point in points])
y_points = coords[:, 0]
z_points = coords[:, 1]

# Aplicar la función a los puntos para calcular X
x_points = func(y_points, z_points)

# Crear la figura y el eje 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Graficar los puntos en 3D
ax.scatter(x_points, y_points, z_points, color='blue', s=10)

# Configuración de los ejes
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

# Ajustar los límites de los ejes
ax.set_xlim([x_points.min(), x_points.max()])
ax.set_ylim([y_points.min(), y_points.max()])
ax.set_zlim([z_points.min(), z_points.max()])

plt.show()
