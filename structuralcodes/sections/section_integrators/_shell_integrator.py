"""A concrete implementation of a shell integrator."""

import typing as t

from _section_integrator import SectionIntegrator
from numpy.typing import ArrayLike

from ...geometry._shell_geometry import ShellGeometry


class ShellFiberIntegrator(SectionIntegrator):
    """A concrete implementation of a fiber integrator for shell sections."""

    def prepare_input(
        self,
        geo: ShellGeometry,
        strain: ArrayLike,
        integrate: t.Literal['stress', 'modulus'] = 'stress',
    ):
        raise NotImplementedError

    def integrate_stress(self, *prepared_input):
        raise NotImplementedError

    def integrate_modulus(self, *prepared_input):
        raise NotImplementedError

    def integrate_strain_response_on_geometry(
        self,
        geo: ShellGeometry,
        strain: ArrayLike,
        integrate: t.Literal['stress', 'modulus'] = 'stress',
    ):
        prepared_input = self.prepare_input(
            geo=geo, strain=strain, integrate=integrate
        )

        # Return the calculated response
        if integrate == 'stress':
            return self.integrate_stress(prepared_input)
        if integrate == 'modulus':
            return self.integrate_modulus(prepared_input)
        raise ValueError(f'Unknown integrate type: {integrate}')
