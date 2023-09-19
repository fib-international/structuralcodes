""""Generic class section implemenetation"""
import typing as t
from structuralcodes.core.base import Section, SectionCalculator

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
        sectionanlyzer: GenericSectionCalculator


class GenericSectionCalculator(SectionCalculator):
    '''Calculator class implementing analysis algortims for '''

    def __init__(self) -> None:
        super().__init__()
