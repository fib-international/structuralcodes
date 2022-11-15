"""Main class and builder for fib Model Code 2010"""
import typing as t

from structuralcodes.core.base import DesignCode

from ._concrete_material_properties import fcm, fctm, fctkmin, fctkmax, Gf, fcd

FUNCTIONS = (fcm, fctm, fctkmin, fctkmax, Gf, fcd)


class MC2020(DesignCode):
    """Main class collecting models from fib Model Code 2010."""

    def __init__(self) -> None:
        super().__init__(
            name='MC 2020', release_year=0000, materials=('concrete')
        )
        for fun in FUNCTIONS:
            setattr(self, fun.__name__, fun)


class MC2020Builder:
    """Build the instance of MC2020"""

    def __init__(self) -> None:
        """Initialize the builder."""
        self._instance: t.Optional[MC2020] = None

    def __call__(self, **_ignored) -> MC2020:
        """Build the code"""
        del _ignored
        if self._instance is None:
            self._instance = MC2020()
        return self._instance
