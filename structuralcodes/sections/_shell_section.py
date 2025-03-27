"""A concrete implementation of a shell section."""

import typing as t

from numpy.typing import ArrayLike

from ..core.base import Section, SectionCalculator
from ..geometry._shell_geometry import ShellGeometry
from .section_integrators._shell_integrator import ShellFiberIntegrator


class ShellSection(Section):
    """A shell section."""

    def __init__(self, geometry: ShellGeometry):
        """Initialize a shell section."""
        super().__init__()
        self._geometry = geometry
        self.section_calculator = ShellSectionCalculator(section=self)

    @property
    def geometry(self):
        return self._geometry


class ShellSectionCalculator(SectionCalculator):
    """A calculator for shell sections."""

    section: ShellSection

    def __init__(self, section: ShellSection):
        super().__init__(section=section)
        self.integrator = ShellFiberIntegrator()

    def integrate_strain_profile(
        self,
        strain: ArrayLike,
        integrate: t.Literal['stress', 'modulus'] = 'stress',
    ):
        """Integrate a strain profile returning stress resultants or tangent
        section stiffness matrix.
        """
        return self.integrator.integrate_strain_response_on_geometry(
            geo=self.section.geometry, strain=strain, integrate=integrate
        )

    def calculate_strain_profile(
        self,
        nx: float,
        ny: float,
        nxy: float,
        mx: float,
        my: float,
        mxy: float,
    ):
        """Get the strain plane for a given set of stress resultants."""
        raise NotImplementedError
