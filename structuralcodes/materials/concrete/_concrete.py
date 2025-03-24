"""Core implementation of the concrete material."""

import abc
import typing as t

from structuralcodes.core.base import ConstitutiveLaw, Material


class Concrete(Material):
    """The abstract concrete material."""

    _fck: float
    _gamma_c: t.Optional[float] = None
    _existing: bool
    _constitutive_law: t.Optional[ConstitutiveLaw]

    def __init__(
        self,
        fck: float,
        name: t.Optional[str] = None,
        density: float = 2400,
        gamma_c: t.Optional[float] = None,
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
        self._gamma_c = gamma_c

    @property
    def fck(self) -> float:
        """Returns fck in MPa."""
        return self._fck

    @property
    @abc.abstractmethod
    def gamma_c(self) -> float:
        """Each concrete should implement its own getter for the partial factor
        in order to interact with the globally set national annex.
        """

    @abc.abstractmethod
    def fcd(self) -> float:
        """Each concrete should implement its own method for calculating the
        design compressive strength.
        """
