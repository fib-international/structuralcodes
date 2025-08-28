"""Classes for rectangular geometries.

The class `RectangularGeometry` represents a rectangular SurfaceGeometry with
homogenous material.
This class is simply a wrapper of `SurfaceGeometry` class and permits an easy
input by the user.
"""

import typing as t

from numpy.typing import ArrayLike
from shapely import Polygon

from structuralcodes.core.base import Material

from ._geometry import SurfaceGeometry


class RectangularGeometry(SurfaceGeometry):
    """This is a wrapper class for defining a `SurfaceGeometry` of rectangular
    shape with a homogeneous material.
    """

    _width: float
    _height: float

    def __init__(
        self,
        width: float,
        height: float,
        material: Material,
        concrete: bool = False,
        origin: t.Optional[ArrayLike] = None,
        name: t.Optional[str] = None,
        group_label: t.Optional[str] = None,
    ) -> None:
        """Initialize a RectangularGeometry.

        Arguments:
            width (float): The width of the geometry.
            height (float): The height of the geometry.
            material (Material): A Material class applied to the geometry.
            concrete (bool): Flag to indicate if the geometry is concrete. When
                passing a Material as material, this is automatically inferred.
            origin (Optional(ArrayLike)): The center point of the rectangle.
                (0.0, 0.0) is used as default.
            name (Optional(str)): The name to be given to the object.
            group_label (Optional(str)): A label for grouping several objects.

        Note:
            The RectangularGeometry is simply a wrapper for a SurfaceGeometry
            object.
        """
        # Check that size is strictly positive
        if width <= 0:
            raise ValueError('Width must be a positive number.')
        if height <= 0:
            raise ValueError('Height must be a positive number.')

        self._width = width
        self._height = height

        # Parse origin
        if origin is not None and len(origin) != 2:
            raise ValueError('origin must be an ArrayLike with len == 2')
        origin = origin if origin is not None else (0.0, 0.0)

        # Create the shapely polygon
        polygon = Polygon(
            (
                (-width / 2 + origin[0], -height / 2 + origin[1]),
                (width / 2 + origin[0], -height / 2 + origin[1]),
                (width / 2 + origin[0], height / 2 + origin[1]),
                (-width / 2 + origin[0], height / 2 + origin[1]),
            )
        )
        # Pass everything to the base class
        super().__init__(
            poly=polygon,
            material=material,
            concrete=concrete,
            name=name,
            group_label=group_label,
        )

    @property
    def width(self):
        """Returns the width of the rectangle."""
        return self._width

    @property
    def height(self):
        """Return the height of the rectangle."""
        return self._height
