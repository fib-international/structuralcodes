"""Classes for circular geometries.

The class `CircularGeometry` represents a circular SurfaceGeometry with
homogenous material.
This class is simply a wrapper of `SurfaceGeometry` class and permits an easy
input by the user.
"""

import typing as t
from itertools import chain

import numpy as np
from shapely import Polygon

from structuralcodes.core.base import ConstitutiveLaw, Material

from ._geometry import CompoundGeometry, SurfaceGeometry
from ._reinforcement import add_reinforcement_circle


def _create_circle(radius, npoints=20):
    """Create a circle with a given radius."""
    phi = np.linspace(0, 2 * np.pi, npoints + 1)
    x = radius * np.cos(phi)
    y = radius * np.sin(phi)
    points = np.transpose(np.array([x, y]))
    return Polygon(points)


class CircularGeometry(SurfaceGeometry):
    """This is a wrapper class for defining a SurfaceGeometry of circular shape
    with a homogeneous material.
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


class CircularRCGeometry(CompoundGeometry):
    """This is a wrapper class for defining a CompoundGeometry made of a
    surface with circular shape and with a set of bars around the circle.
    """

    _radius: float

    def __init__(
        self,
        size: float,
        surface_material: t.Union[Material, ConstitutiveLaw],
        diameter: float,
        number: int,
        reinforcement_material: t.Union[Material, ConstitutiveLaw],
        cover: float,
        is_radius: bool = False,
        n_points: int = 20,
        surface_density: t.Optional[float] = None,
        concrete: bool = False,
    ) -> None:
        """Initialize a CircularRCGeometry.

        Arguments:
            size (float): The size of the geometry. By default this number
                is assumed as the diameter of the section. To interpret this
                as the radius `is_radius` should be set as `True`.
            suface_material (Union(Material, ConstitutiveLaw)): A Material or
                ConsitutiveLaw class applied to the surface geometry.
            diameter (float): The diameter of the bars.
            number (int): The number of bars to distribute.
            reinforcement_material (Union(Material, ConstitutiveLaw)): A
                Material or ConsitutiveLaw class applied to the point geometry.
            cover (float): The cover of concrete intended as the distance
                between outer line and barycenter of longitudinal bars.
            is_radius (bool): Indicates if size is interpreted as radius
                    (default = False).
            n_points (int): The number of points used to discretize the
                circle as a shapely `Polygon` (default = 20).
            surface_density (Optional(float)): When a ConstitutiveLaw is passed
                as surface_material, the density can be provided by this
                argument. When surface_material is a Material object the
                density is taken from the material.
            concrete (bool): Flag to indicate if the surface geometry is
                concrete.

        Note:
            The CircularRCGeometry is simply a wrapper for a CompoundGeometry
            object.
        """
        # Check that size is strictly positive
        if size <= 0:
            raise ValueError('Size must be a positive number.')
        # Manage size as radius or diameter (default)
        self._radius = size if is_radius else size / 2.0
        # Create the CircularGeometry
        geometry = CircularGeometry(
            size,
            surface_material,
            is_radius,
            n_points,
            surface_density,
            concrete,
        )
        # Add the reinforcement
        r = geometry.radius - cover
        rc_geometry = add_reinforcement_circle(
            geometry, (0, 0), r, diameter, reinforcement_material, number
        )

        # Pass everything to the base class CompoundGeometry
        super().__init__(
            list(chain(rc_geometry.geometries, rc_geometry.point_geometries))
        )

    @property
    def radius(self):
        """Returns the radius of the circle."""
        return self._radius

    @property
    def diameter(self):
        """Return diameter of the circle."""
        return 2 * self._radius
