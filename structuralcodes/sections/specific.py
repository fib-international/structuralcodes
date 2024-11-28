"""Specific class section implementation."""

from __future__ import annotations

import typing as t

from shapely import Polygon

from structuralcodes.core.base import Material
from structuralcodes.geometry import SurfaceGeometry
from structuralcodes.sections._generic import GenericSection


class RectangularSection(GenericSection):
    """This is the specific implementation of a rectangular class section.

    The section is a 2D geometry where Y axis is horizontal while Z axis is
    vertical. The geometry has its center in origin.

    The moments and curvatures around Y and Z axes are assumed positive
    according to RHR.

    Attributes: TODO
        geometry (Union(SurfaceGeometry, CompoundGeometry)): The geometry of
            the section.
        name (str): The name of the section.
        section_calculator (GenericSectionCalculator): The object responsible
            for performing different calculations on the section (e.g. bending
            strength, moment curvature, etc.).
    """

    def __init__(
        self,
        height: float,
        width: float,
        material: Material,
        name: t.Optional[str] = None,
        integrator: t.Literal['marin', 'fiber'] = 'marin',
    ) -> None:
        """Initialize a rectangular section.

        Arguments:
            height (float): height of the rectangular section.
            width (float): width of the rectagnular section.
            material (Material): material of the section.
            name (str): The name of the section. Default value is set to None.
            integrator (str): Name of the integrator. Default value: marin.
        """
        poly: Polygon = Polygon(
            (
                (-width / 2, -height / 2),
                (width / 2, -height / 2),
                (width / 2, height / 2),
                (-width / 2, height / 2),
            )
        )
        geometry: SurfaceGeometry = SurfaceGeometry(poly, material)

        if name is None:
            name: str = 'RectangularSection'
        super().__init__(geometry, name, integrator)

        self.height: float = height
        self.width: float = width
