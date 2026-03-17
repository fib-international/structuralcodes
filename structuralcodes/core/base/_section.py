"""Abstract base classes for sections."""

from __future__ import annotations  # To have clean hints of ArrayLike in docs

import abc
import typing as t

from ._geometry import Geometry


class Section(abc.ABC):
    """Abstract base class for a cross secion.
    The section is defined by local axes y and z (mapped to x and y cartesian
    plane respectively).
    """

    geometry: Geometry
    section_counter: t.ClassVar[int] = 0
    section_calculator: SectionCalculator

    def __init__(self, name: t.Optional[str] = None) -> None:
        """Initialize a Section object."""
        self.id = self.section_counter
        self._name = name if name is not None else f'Section_{self.id}'
        self._increase_global_counter()

    @property
    def name(self):
        """Returns the name of the section."""
        return self._name

    @classmethod
    def _increase_global_counter(cls):
        cls.section_counter += 1


class SectionCalculator(abc.ABC):
    """Abstract class for SectionCalculators
    defining the interface.
    """

    def __init__(self, section: Section) -> None:
        """Initialization of SectionCalculator object, providing
        a Section object.
        """
        self.section = section
