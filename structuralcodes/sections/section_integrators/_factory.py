"""Factory for integrators."""

import typing as t

from ._fiber_integrator import FiberIntegrator
from ._marin_integrator import MarinIntegrator
from ._section_integrator import SectionIntegrator


integrator_registry = {'Marin': MarinIntegrator, 'Fiber': FiberIntegrator}


class IntegratorFactory:
    """A factory for integrators."""

    instances: t.Dict  # A dict of integrators that have been created.
    registry: t.Dict[str, SectionIntegrator]

    def __init__(self, registry: t.Dict[str, SectionIntegrator]):
        self.instances = {}
        self.registry = registry

    def __call__(self, method: str) -> SectionIntegrator:
        """Create an integrator based on its name."""
        self.instances.setdefault(
            method, self.registry.get(method, MarinIntegrator)
        )
        # Here we should throw a warning if the user input an integrator not
        # given in the registry.
        return self.instances.get(method)


integrator_factory = IntegratorFactory(integrator_registry)
