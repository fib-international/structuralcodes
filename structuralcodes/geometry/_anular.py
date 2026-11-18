"""Classes for circular geometries.

The class `CircularGeometry` represents a circular SurfaceGeometry with
homogenous material.
This class is simply a wrapper of `SurfaceGeometry` class and permits an easy
input by the user.
"""

import typing as t

from numpy.typing import ArrayLike

from ..core.base import Material
from ._circular import _create_circle
from ._geometry import SurfaceGeometry


class AnularGeometry(SurfaceGeometry):
    """This is a wrapper class for defining a `SurfaceGeometry` of circular
    shape with a homogeneous material.
    """

    _outer_radius: float
    _inner_radius: float

    def __init__(
        self,
        outer_diameter: float,
        inner_diameter: float,
        material: Material,
        n_points: int = 20,
        concrete: bool = False,
        origin: t.Optional[ArrayLike] = None,
        name: t.Optional[str] = None,
        group_label: t.Optional[str] = None,
    ) -> None:
        """Initialize a CircularGeometry.

        Arguments:
            outer_diameter (float): The external diameter of the geometry.
            inner_diameter (float): The internal diameter of the geometry.
            material (Material): A Material class applied to the geometry.
            n_points (int): The number of points used to discretize the
                circle as a shapely `Polygon` (default = 20). Note that the
                number of points is rounded to the nearest multiple of 4.
            concrete (bool): Flag to indicate if the geometry is concrete.
            origin (Optional(ArrayLike)): The center point of the circle.
                (0.0, 0.0) is used as default.
            name (Optional(str)): The name to be given to the object.
            group_label (Optional(str)): A label for grouping several objects.

        Note:
            The AnularGeometry is simply a wrapper for a SurfaceGeometry
            object.
        """
        # Check that size is strictly positive
        if outer_diameter <= 0 or inner_diameter <= 0:
            raise ValueError('Diameters must be positive numbers.')
        if outer_diameter <= inner_diameter:
            raise ValueError(
                'Outer diameter should be greater than inner diameter.'
            )
        # Manage size as radius or diameter
        self._outer_radius = outer_diameter / 2.0
        self._inner_radius = inner_diameter / 2.0
        # Parse origin
        if origin is not None and len(origin) != 2:
            raise ValueError('origin must be an ArrayLike with len == 2')
        origin = origin if origin is not None else (0.0, 0.0)
        # Create the shapely polygon
        # Outer polygon
        polygon_out = _create_circle(
            radius=self._outer_radius, npoints=n_points, origin=origin
        )
        # Inner polygon
        polygon_in = _create_circle(
            radius=self._inner_radius, npoints=n_points, origin=origin
        )

        # Pass everything to the base class
        super().__init__(
            poly=polygon_out - polygon_in,
            material=material,
            concrete=concrete,
            name=name,
            group_label=group_label,
        )

    @property
    def outer_radius(self):
        """Returns the outer radius of the anular geometry."""
        return self._outer_radius

    @property
    def outer_diameter(self):
        """Return outer diameter of the anular geometry."""
        return 2 * self._outer_radius

    @property
    def inner_radius(self):
        """Returns the inner radius of the anular geometry."""
        return self._inner_radius

    @property
    def inner_diameter(self):
        """Return inner diameter of the anular geometry."""
        return 2 * self._inner_radius

    @property
    def thickness(self):
        """Return the thickness of the anular geometry."""
        return self._outer_radius - self._inner_radius
