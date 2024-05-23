"""Core implementation of the concrete material."""

import abc
import typing as t

from structuralcodes.core.base import ConstitutiveLaw, Material
from structuralcodes.materials.constitutive_laws import ParabolaRectangle


class Concrete(Material):
    """The abstract concrete material."""

    _fck: float
    _gamma_c: t.Optional[float] = None
    _existing: bool

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
        self._constitutive_law = ParabolaRectangle(
            self._fck, name=name + '_ConstLaw'
        )
        self._gamma_c = gamma_c

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

    @property
    def constitutive_law(self) -> ConstitutiveLaw:
        """Returns the constitutive law object."""
        return self._constitutive_law

    @constitutive_law.setter
    def constitutive_law(self, constitutive_law: ConstitutiveLaw) -> None:
        """Setter for constitutive law."""
        if 'concrete' in constitutive_law.__materials__:
            self._constitutive_law = constitutive_law
        else:
            raise ValueError(
                'The constitutive law selected is not suitable '
                'for being used with a concrete material.'
            )

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
