"""Classes for circular geometries.

The class `CircularGeometry` represents a circular SurfaceGeometry with
homogenous material.
This class is simply a wrapper of `SurfaceGeometry` class and permits an easy
input by the user.
"""

import typing as t

import numpy as np
from shapely import Polygon

from structuralcodes.core.base import ConstitutiveLaw, Material

from ._geometry import SurfaceGeometry


def _create_circle(radius, npoints=20):
    """Create a circle with a given radius."""
    phi = np.linspace(0, 2 * np.pi, npoints + 1)
    x = radius * np.cos(phi)
    y = radius * np.sin(phi)
    points = np.transpose(np.array([x, y]))
    return Polygon(points)


class CircularGeometry(SurfaceGeometry):
    """This is a wrapper class for defining a `SurfaceGeometry` of circular
    shape with a homogeneous material.
    """

    _radius: float

    def __init__(
        self,
        diameter: float,
        material: t.Union[Material, ConstitutiveLaw],
        n_points: int = 20,
        density: t.Optional[float] = None,
        concrete: bool = False,
    ) -> None:
        """Initialize a CircularGeometry.

        Arguments:
            diameter (float): The diameter of the geometry.
            material (Union(Material, ConstitutiveLaw)): A Material or
                ConsitutiveLaw class applied to the geometry.
            n_points (int): The number of points used to discretize the
                circle as a shapely `Polygon` (default = 20).
            density (Optional(float)): When a ConstitutiveLaw is passed as
                material, the density can be provided by this argument. When
                material is a Material object the density is taken from the
                material.
            concrete (bool): Flag to indicate if the geometry is concrete.

        Note:
            The CircularGeometry is simply a wrapper for a SurfaceGeometry
            object.
        """
        # Check that size is strictly positive
        if diameter <= 0:
            raise ValueError('Diameter must be a positive number.')
        # Manage size as radius or diameter (default)
        self._radius = diameter / 2.0
        # Create the shapely polygon
        polygon = _create_circle(radius=self._radius, npoints=n_points)
        # Pass everything to the base class
        super().__init__(
            poly=polygon, material=material, density=density, concrete=concrete
        )

    @property
    def radius(self):
        """Returns the radius of the circle."""
        return self._radius

    @property
    def diameter(self):
        """Return diameter of the circle."""
        return 2 * self._radius
