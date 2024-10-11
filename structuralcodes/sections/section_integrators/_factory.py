"""Factory for integrators."""

import typing as t

from ._fiber_integrator import FiberIntegrator
from ._marin_integrator import MarinIntegrator
from ._section_integrator import SectionIntegrator

integrator_registry = {'marin': MarinIntegrator, 'fiber': FiberIntegrator}


class IntegratorFactory:
    """A factory for integrators."""

    instances: t.Dict  # A dict of integrators that have been created.
    registry: t.Dict[str, SectionIntegrator]

    def __init__(self, registry: t.Dict[str, SectionIntegrator]):
        self.instances = {}
        self.registry = registry

    def __call__(
        self, method: t.Literal['marin', 'fiber']
    ) -> SectionIntegrator:
        """Create an integrator based on its name.

        Arguments:
            method (str): The name of the integrator to use.

        Returns:
            SectionIntegrator: The section integrator.
        """
        self.instances.setdefault(
            method.lower(), self.registry.get(method.lower(), MarinIntegrator)
        )
        # Here we should throw a warning if the user input an integrator not
        # given in the registry.
        return self.instances.get(method.lower())


integrator_factory = IntegratorFactory(integrator_registry)
