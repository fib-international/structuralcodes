"""Core implementation of the concrete material."""

import abc
import typing as t

from structuralcodes.core.base import ConstitutiveLaw, Material
from structuralcodes.materials.constitutive_laws import (
    ParabolaRectangle,
    create_constitutive_law,
)


class Concrete(Material):
    """The abstract concrete material."""

    _fck: float
    _gamma_c: t.Optional[float] = None
    _existing: bool
    _constitutive_law: t.Optional[ConstitutiveLaw] = None

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
        self._constitutive_law = None

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
        if self._constitutive_law is None:
            self._constitutive_law = ParabolaRectangle(
                self.fcd(), name=self.name + '_ConstLaw'
            )
        return self._constitutive_law

    @constitutive_law.setter
    def constitutive_law(
        self,
        constitutive_law: t.Optional[ConstitutiveLaw] = None,
        constitutive_law_name: t.Optional[str] = None,
    ) -> None:
        """Setter for constitutive law.

        Args:
        consitutive_law (ConstitutiveLaw): a valid ConstitutiveLaw object
                for concrete (optional)
        constitutive_law_name (str): a string defining a valid constitutive law
                type for concrete. If a name is provided a new constitutive law
                of that class will be created ignoring if the user provided
                also the constitutive_law argument
                (valid options: 'elastic', 'parabolarectangle', 'sargin',
                'popovics').
        """
        if constitutive_law is None and constitutive_law_name is None:
            raise ValueError(
                'At least a constitutive law or a string defining the '
                'constitutive law must be provided.'
            )
        if constitutive_law_name is not None:
            constitutive_law = create_constitutive_law(
                constitutive_law_name=constitutive_law_name, material=self
            )

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
