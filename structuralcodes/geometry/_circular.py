"""Classes for circular geometries.

The class `CircularGeometry` represents a circular SurfaceGeometry with
homogenous material.
This class is simply a wrapper of `SurfaceGeometry` class and permits an easy
input by the user.
"""

import typing as t

import numpy as np
from numpy.typing import ArrayLike
from shapely import Polygon

from structuralcodes.core.base import ConstitutiveLaw, Material

from ._geometry import SurfaceGeometry


def _create_circle(radius, npoints=20, origin: t.Optional[ArrayLike] = None):
    """Create a circle with a given radius."""
    origin = origin if origin is not None else (0.0, 0.0)
    phi = np.linspace(0, 2 * np.pi, npoints + 1)
    x = radius * np.cos(phi) + origin[0]
    y = radius * np.sin(phi) + origin[1]
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
        origin: t.Optional[ArrayLike] = None,
        name: t.Optional[str] = None,
        group_label: t.Optional[str] = None,
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
            origin (Optional(ArrayLike)): The center point of the circle.
                (0.0, 0.0) is used as default.
            name (Optional(str)): The name to be given to the object.
            group_label (Optional(str)): A label for grouping several objects.

        Note:
            The CircularGeometry is simply a wrapper for a SurfaceGeometry
            object.
        """
        # Check that size is strictly positive
        if diameter <= 0:
            raise ValueError('Diameter must be a positive number.')
        # Manage size as radius or diameter (default)
        self._radius = diameter / 2.0
        # Parse origin
        if origin is not None and len(origin) != 2:
            raise ValueError('origin must be an ArrayLike with len == 2')
        origin = origin if origin is not None else (0.0, 0.0)
        # Create the shapely polygon
        polygon = _create_circle(
            radius=self._radius, npoints=n_points, origin=origin
        )
        # Pass everything to the base class
        super().__init__(
            poly=polygon,
            material=material,
            density=density,
            concrete=concrete,
            name=name,
            group_label=group_label,
        )

    @property
    def radius(self):
        """Returns the radius of the circle."""
        return self._radius

    @property
    def diameter(self):
        """Return diameter of the circle."""
        return 2 * self._radius
