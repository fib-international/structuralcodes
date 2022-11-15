"""Main class and builder for Eurocode 2 EN1992:2004"""
import typing as t

from structuralcodes.core.base import DesignCode

from ._concrete_material_properties import fcm, fctm, fctkmin, fctkmax, Gf, fcd

FUNCTIONS = (fcm, fctm, fctkmin, fctkmax, Gf, fcd)


class EN1992_2004(DesignCode):
    """Main class collecting models from Eurocode 2 EN1992:2004"""

    def __init__(self) -> None:
        super().__init__(
            name='EN1992-1-1', release_year=2004, materials=('concrete')
        )
        for fun in FUNCTIONS:
            setattr(self, fun.__name__, fun)


class EN1992_2004Builder:
    """Build the instance of EN1992-1-1:2004"""

    def __init__(self) -> None:
        """Initialize the builder."""
        self._instance: t.Optional[EN1992_2004] = None

    def __call__(self, **_ignored) -> EN1992_2004:
        """Build the code"""
        del _ignored
        if self._instance is None:
            self._instance = EN1992_2004()
        return self._instance
