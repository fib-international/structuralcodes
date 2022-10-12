"""Main class and builder for fib Model Code 2010"""
import typing as t

from structuralcodes.core.base import DesignCode

from ._concrete_material_properties import fcm, fctm

FUNCTIONS = (fcm, fctm)


class MC2010(DesignCode):
    """Main class collecting models from fib Model Code 2010."""

    def __init__(self) -> None:
        super().__init__()
        for fun in FUNCTIONS:
            setattr(self, fun.__name__, fun)


class MC2010Builder:
    """Build the instance of MC2010"""

    def __init__(self) -> None:
        """Initialize the builder."""
        self._instance: t.Optional[MC2010] = None

    def __call__(self, **_ignored) -> MC2010:
        """Build the code"""
        del _ignored
        if self._instance is None:
            self._instance = MC2010()
        return self._instance
