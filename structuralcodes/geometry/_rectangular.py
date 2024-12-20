"""Classes for rectangular geometries.

The class `RectangularGeometry` represents a rectangular SurfaceGeometry with
homogenous material.
This class is simply a wrapper of `SurfaceGeometry` class and permits an easy
input by the user.
"""

import typing as t

from shapely import Polygon

from structuralcodes.core.base import ConstitutiveLaw, Material

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
        material: t.Union[Material, ConstitutiveLaw],
        density: t.Optional[float] = None,
        concrete: bool = False,
    ) -> None:
        """Initialize a RectangularGeometry.

        Arguments:
            width (float): The width of the geometry.
            height (float): The height of the geometry.
            material (Union(Material, ConstitutiveLaw)): A Material or
                ConsitutiveLaw class applied to the geometry.
            density (Optional(float)): When a ConstitutiveLaw is passed as
                material, the density can be provided by this argument. When
                material is a Material object the density is taken from the
                material.
            concrete (bool): Flag to indicate if the geometry is concrete. When
                passing a Material as material, this is automatically inferred.

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

        # Create the shapely polygon
        polygon = Polygon(
            (
                (-width / 2, -height / 2),
                (width / 2, -height / 2),
                (width / 2, height / 2),
                (-width / 2, height / 2),
            )
        )
        # Pass everything to the base class
        super().__init__(
            poly=polygon, material=material, density=density, concrete=concrete
        )

    @property
    def width(self):
        """Returns the width of the rectangle."""
        return self._width

    @property
    def height(self):
        """Return the height of the rectangle."""
        return self._height
