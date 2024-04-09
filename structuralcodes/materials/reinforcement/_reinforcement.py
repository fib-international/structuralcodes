"""Core implementation of the reinforcement material."""

import typing as t

from structuralcodes.core.base import Material


class Reinforcement(Material):
    """The abstract reinforcement material."""

    _fyk: float

    def __init__(
        self, fyk: float, density: float, name: t.Optional[str] = None
    ) -> None:
        """Initializes an abstract reinforcement material."""
        name = name if name is not None else 'Reinforcement'
        super().__init__(density, name)

        self._fyk = abs(fyk)

    @property
    def fyk(self) -> float:
        """Returns fyk in MPa."""
        return self._fyk

    @fyk.setter
    def fyk(self, fyk: float) -> None:
        """Setter for fyk (in MPa)."""
        self._fck = abs(fyk)
