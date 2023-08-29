""""Generic class section implemenetation"""
import typing as t
from structuralcodes.core.base import Section

# Generic cross section implementation
# Characteristics:
# outer polygon (define a polygon class as a list of points)
# inner polygons (for holes)
# list of reinforcement positions, area and material (point reinforcement)
# list of reinforcement poistions, thickness and material (line reinforcement)?


class GenericSection(Section):
    """This is the implementation of the generic class section"""

    def __init__(self, name: t.Optional[str] = None) -> None:
        if name is None:
            name = 'GenericSection'
        super().__init__(name)
        
    @property
    def area(self) -> float:
        """Returns the area of concrete"""
        return NotImplemented
    
    def bending_strength_xp(self, N: float = 0) -> float:
        """Returns the beding strength in x+ direction for a given
        value of axial force (+: tension, -: compression)"""
        return NotImplemented

    def bending_strength_xn(self, N: float = 0) -> float:
        """Returns the beding strength in x- direction for a given
        value of axial force (+: tension, -: compression)"""
        return NotImplemented

    def bending_strength_yp(self, N: float = 0) -> float:
        """Returns the beding strength in y+ direction for a given
        value of axial force (+: tension, -: compression)"""
        return NotImplemented

    def bending_strength_yn(self, N: float = 0) -> float:
        """Returns the beding strength in y- direction for a given
        value of axial force (+: tension, -: compression)"""
        return NotImplemented
