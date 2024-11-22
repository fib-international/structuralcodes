"""Classes for circular sections.

The class `Circular` represents a circular section with homogenous material.
The class `CircularRC` represents a reinforced concrete circular section
with `Concrete` material for the `SurfaceGeometry` and `Reinforcement` material
for the `PointGeometries`.

These classes are simply wrappers of `GenericSection` class and permits an easy
input by the user, but all calculations are performed in the under-the-hood
`GenericSection` class.

In the future this class will contain utility methods for automatically define
confined concrete according to some confinement models. For now a not
implemented error is raised if this feature is used.
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
    points = [(this_x, this_y) for this_x, this_y in zip(x, y)]
    return Polygon(points)


class CircularGeometry(SurfaceGeometry):
    """This is a wrapper class for defining a SurfaceGeometry of circular shape
    with a homogeneous material.
    """

    _radius: float

    def __init__(
        self,
        size: float,
        material: t.Union[Material, ConstitutiveLaw],
        is_radius: bool = False,
        n_points: int = 20,
        density: t.Optional[float] = None,
        concrete: bool = False,
    ) -> None:
        """Initialize a CircularSection.

        Arguments:
            size (float): The size of the geometry. By default this number
                is assumed as the diameter of the section. To interpret this
                as the radius `is_radius` should be set as `True`.
            material (Union(Material, ConstitutiveLaw)): A Material or
                ConsitutiveLaw class applied to the geometry.
                is_radius (bool): Indicates if size is interpreted as radius
                    (default = False).
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
        if size <= 0:
            raise ValueError('Size must be a positive number.')
        # Manage size as radius or diameter (default)
        self._radius = size if is_radius else size / 2.0
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
