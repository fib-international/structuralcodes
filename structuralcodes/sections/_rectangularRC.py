"""Rectangular RC cross section"""
from dataclasses import dataclass, field
from math import cos, sin, pi
import warnings
import typing as t
# import numpy as np
# import numpy.typing as npt
from structuralcodes.core.base import Section
from structuralcodes.materials.concrete import Concrete
# For now this is just a trick! refactor later TODO!
from structuralcodes.materials.constitutive_laws import ElasticPlastic


@dataclass(frozen=True)
class StirrupsData:
    """dataclass representing stirrups data for RectanularRC
    section.
    diameter: diameter of transverse reinforcement
    spacing: spacing between consecutive stirrups
    number_legs_x: number of legs parallel to x local axis
    number_legs_y: number of legs parallel to y local axis
    alfa: angle between stirrups and longitudinal axis of beam (rad)"""

    diameter: float
    spacing: float
    number_legs_x: float = 2
    number_legs_y: float = 2
    alfa: float = pi / 2
    sin_alfa: float = field(init=False, repr=False)
    ctg_alfa: float = field(init=False, repr=False)

    def __post_init__(self):
        if self.alfa > 2 * pi:
            # Warning: I assume that this are given in degrees
            radians = (self.alfa * pi) / 180
            warnings.warn(
                "StirrupsData: alfa should be given in radians. "
                + f"Since the value alfa={self.alfa} has been "
                + "provided it is assumed that is given in degrees."
                + f"Converted to {radians:.2f} radians."
            )
            object.__setattr__(self, 'alfa', radians)
        sin_alfa = sin(self.alfa)
        if sin_alfa == 0:
            raise ValueError("StirrupsData: the sin of alfa cannot be zero.")
        object.__setattr__(self, 'sin_alfa', sin_alfa)
        object.__setattr__(self, 'ctg_alfa', cos(self.alfa) / self.sin_alfa)


@dataclass(frozen=True)
class ReinforcementData:
    """dataclass representing stirrups data for RectangularRC
    section.
    diameter_angle: diameter of longitudinal bar in each angle of the section
    diameter_bottom: diameter of bars in bottom side (default = None)
    number_bottom: number of bars in bottom side (default = 0)
    diameter_top: diameter of bars in top side (default = None)
    number_top: number of bars in top side (default = 0)
    diameter_left: diameter of bars in left side (default = None)
    number_left: number of bars in left side (default = 0)
    diameter_right: diameter of bars in right side (default = None)
    number_right: number of bars in right side (default = 0)

    if no diameter is provided for sides, it is assumed that same diameter
    of angle bars is mantained for all other sides."""
    diameter_corner: float
    diameter_bottom: t.Optional[float] = None
    number_bottom: int = 0
    diameter_top: t.Optional[float] = None
    number_top: int = 0
    diameter_left: t.Optional[float] = None
    number_left: int = 0
    diameter_right: t.Optional[float] = None
    number_right: int = 0

    def __post_init__(self):
        if self.diameter_bottom is None:
            object.__setattr__(self, 'diameter_bottom', self.diameter_corner)
        if self.diameter_top is None:
            object.__setattr__(self, 'diameter_top', self.diameter_corner)
        if self.diameter_left is None:
            object.__setattr__(self, 'diameter_left', self.diameter_corner)
        if self.diameter_right is None:
            object.__setattr__(self, 'diameter_right', self.diameter_corner)


class RectangularRC(Section):
    """Class for RC rectangular cross section.
    The class assumes width and height respectively
    along local x and y axes

    Input parameters:
    width: width of the section
    height: height of the section
    cover: net cover (distance from outer edget of the section to
    edge of the stirrup)
    reinforcement_data: ReinforcementData
    concrete: Concrete
    reinforcement: Reinforcement (TODO!!!!!! )
    confined_concrete: Concrete (default = None) When set to None the
    same concrete is adopted for concrete core. When using the string
    "auto" the confinement will be automatically computed
    stirrups_data: StirrupsData (default = None) to use when need
    automatic computation of concrete confinement
    name: optional string for identidying the section"""

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
            name: t.Optional[str] = None) -> None:
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
                if confined_concrete.lower() == 'automatic':
                    # TODO: compute automatic confinement
                    if stirrups_data is None:
                        raise ValueError(
                            'Asked for automatic confinement \
                            computation, but stirrups_data is None \
                            \nPlease provide a proper stirrups_data')
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

    @property
    def area(self) -> float:
        """Returns Area of the cross section.
        Only the gross area of concrete is computed

        Returns:
            float: The area of concrete section.
        """
        return self.width * self.height

    def bending_strength_xp(self, N: float = 0) -> float:
        """Returns the beding strength in x+ direction for a given
        value of axial force (+: tension, -: compression)"""
        return 0.0

    def bending_strength_xn(self, N: float = 0) -> float:
        """Returns the beding strength in x- direction for a given
        value of axial force (+: tension, -: compression)"""
        return 0.0

    def bending_strength_yp(self, N: float = 0) -> float:
        """Returns the beding strength in y+ direction for a given
        value of axial force (+: tension, -: compression)"""
        return 0.0

    def bending_strength_yn(self, N: float = 0) -> float:
        """Returns the beding strength in y- direction for a given
        value of axial force (+: tension, -: compression)"""
        return 0.0

    def _arrange_reinforcement_in_section(self):
        """the method creates the arrays for reinforcements
        starting from reinforcement_data"""
        delta = 2 * (self.cover + self.stirrups_data.diameter / 2.0)
        hc = self.height - delta
        wc = self.width - delta
        delta = self.reinforcement_data.diameter_corner + \
            self.stirrups_data.diameter
        hcc = hc - delta
        wcc = wc - delta
        xmin = -hcc / 2.0
        xmax = hcc / 2.0
        ymin = -wcc / 2.0
        ymax = wcc / 2.0
        self.xs = [xmin, xmax, xmax, xmin]
        self.ys = [ymin, ymin, ymax, ymax]
        self.d = [self.reinforcement_data.diameter_corner] * 4
        # bottom bars
        for i in range(self.reinforcement_data.number_bottom):
            self.xs.append(xmin + (xmax - xmin) / (self.reinforcement_data.number_bottom + 1) * (i + 1))
            self.ys.append(ymin)
            self.d.append(self.reinforcement_data.diameter_bottom)

