""" Specific class section implementation"""

from __future__ import annotations

import typing as t
import structuralcodes.core._section_results as s_res
from shapely import Polygon
from structuralcodes.sections._generic import GenericSection
from structuralcodes.geometry import SurfaceGeometry
from structuralcodes.materials.concrete import Concrete
from structuralcodes.core.base import Material

class RectangularSection(GenericSection): 
    """This is the specific implementation of a rectangular class section.

    The section is a 2D geometry where Y axis is horizontal while Z axis is
    vertical. The geometry has its center in origo

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

    def __init__(self, height: float, width: float, material: Material, name: t.Optional[str] = None, integrator: t.Literal['marin', 'fiber'] = 'marin') -> None:
        """Initialize a GenericSection.

        Arguments:
            geometry (Union(SurfaceGeometry, CompoundGeometry)): The geometry
                of the section.
            name (str): The name of the section. Default value is set to None
            integrator (str): The name of the SectionIntegrator to use. Default value is set to marin
        """
        poly: Polygon = Polygon(((-width/2, -width/2), (width/2, -width/2), (width/2, width/2), (-width/2, width/2)))
        geometry: SurfaceGeometry = SurfaceGeometry(poly, material)

        if name is None:
            name: str = 'RectangularSection'
        super().__init__(geometry, name, integrator)

        self.height: float = height
        self.width: float = width
        
    def calculate_gross_section_properties(self):
        # Calculating gross section properties
        # Gross area
        area: float = self.height * self.width
        perimeter: float = 2 * (self.height + self.width)
        ea: float = 0
        mass: float = 0
        area_reinforcement: float = 0
        for geo in self.geometry.geometries:
            ea += area * geo.material.get_tangent(eps=0)[0]
            if geo.density is not None:
                mass += area * geo.density * 1e-9

        for geo in self.geometry.point_geometries:
            ea += area * geo.material.get_tangent(eps=0)[0]
            area_reinforcement += geo.area
            if geo.density is not None:
                mass += area * geo.density * 1e-9


        gp = s_res.GrossProperties()
        gp.area = area
        gp.perimeter = perimeter
        gp.ea = ea
        gp.mass = mass
        self.gross_properties = gp
        # TODO: more gross sectional properties, and store them in self._gross_properties

    @property
    def area(self) -> float:
        """Returns Area of the cross section.
        Only the gross area of concrete is computed.

        Returns:
            float: The area of concrete section.
        """
        return self.width * self.height


class RectangularRC(Section):
    """Class for RC rectangular cross section.
    The class assumes width and height respectively
    along local x and y axes.

    Input parameters:
    width: width of the section
    height: height of the section
    cover: net cover (distance from outer edget of the section to
    edge of the stirrup)
    reinforcement_data: ReinforcementData
    concrete: Concrete
    reinforcement: Reinforcement (TODO!!!!!! for now const law)
    confined_concrete: Concrete (default = None) When set to None the
    same concrete is adopted for concrete core. When using the string
    "auto" the confinement will be automatically computed
    stirrups_data: StirrupsData (default = None) to use when need
    automatic computation of concrete confinement
    name: optional string for identidying the section

    Axes definition:

                     z ▲
                       │
                       │
            ┌──────────┼───────────┐
            │          │           │
            │          │           │
            │          │           │
            │          │           │
            │          │           │
            │          │           │
            │          │           │
            │          └───────────┼─────►
            │                      │     y
            │                      │
            │                      │
            │                      │
            │                      │
            │                      │
            │                      │
            └──────────────────────┘

            RectangularRCSection.fromGeometry(Geometry(itPolygons,itMaterials))


    """

    def __init__(
        self,
        width: float,
        height: float,
        cover: float,
        reinforcement_data: ReinforcementData,
        concrete: Concrete,
        reinforcement: ElasticPlastic,
        confined_concrete: t.Union[Concrete, str, None] = None,
        stirrups_data: t.Optional[StirrupsData] = None,
        name: t.Optional[str] = None,
    ) -> None:
        if name is None:
            name = 'RectangularRCSection'
        super().__init__(name)

        # section geometry
        self.height = height
        self.width = width
        self.cover = cover

        # materials
        self.concrete = concrete
        self.reinforcement = reinforcement
        if confined_concrete is not None:
            if isinstance(confined_concrete, str):
                # provided a string, check if it auto
                if (
                    confined_concrete.lower() == 'automatic'
                    or confined_concrete.lower() == 'auto'
                ):
                    # TODO: compute automatic confinement
                    if stirrups_data is None:
                        raise ValueError(
                            'Asked for automatic confinement \
                            computation, but stirrups_data is None \
                            \nPlease provide a proper stirrups_data'
                        )
                    self.confined_concrete = concrete
                else:
                    raise ValueError(f'Value "{confined_concrete}" unknown')
            elif isinstance(confined_concrete, Concrete):
                # provided a Concrete material, used that one
                self.confined_concrete = confined_concrete
        else:
            self.confined_concrete = concrete

        # stirrups
        if stirrups_data is None:
            # If no stirrups are provided by default it is assumed the
            # diameter is 8 mm - The spacing is not significant
            stirrups_data = StirrupsData(8, 150)
        self.stirrups_data = stirrups_data

        # longitudinal reinforcement
        self.reinforcement_data = reinforcement_data
        # arrange in the section the reinforcement
        self._arrange_reinforcement_in_section()
