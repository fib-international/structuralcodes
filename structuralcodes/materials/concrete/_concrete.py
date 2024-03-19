"""Core implementation of the concrete material."""

import abc
import typing as t

from structuralcodes.core.base import Material


class Concrete(Material):
    """The abstract concrete material."""

    _fck: float
    _existing: bool

    def __init__(
        self,
        fck: float,
        name: t.Optional[str] = None,
        density: float = 2400,
        existing: t.Optional[bool] = False,
    ) -> None:
        """Initializes an abstract concrete material."""
        name = name if name is not None else 'Concrete'
        super().__init__(density=density, name=name)

        self._fck = abs(fck)
        if existing:
            raise NotImplementedError(
                'Existing concrete feature not implemented yet'
            )
        self._existing = existing

    @property
    def fck(self) -> float:
        """Returns fck in MPa."""
        return self._fck

    @fck.setter
    def fck(self, fck: float) -> None:
        """Setter for fck (in MPa)."""
        self._fck = abs(fck)
        self._reset_attributes()

    @abc.abstractmethod
    def _reset_attributes(self):
        """Each concrete should define its own _reset_attributes method
        This is because fck setting, reset the object arguments.
        """
