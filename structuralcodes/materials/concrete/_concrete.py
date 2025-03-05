"""Core implementation of the concrete material."""

import abc
import typing as t

from structuralcodes.core.base import ConstitutiveLaw, Material
from structuralcodes.materials.constitutive_laws import create_constitutive_law


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
        """Each concrete should define its own _reset_attributes method. This
        is because fck setting, reset the object arguments.
        """

    @property
    def constitutive_law(self) -> ConstitutiveLaw:
        """Returns the constitutive law object."""
        if self._constitutive_law is None:
            self.constitutive_law = 'parabolarectangle'
        return self._constitutive_law

    @constitutive_law.setter
    def constitutive_law(
        self,
        constitutive_law: t.Union[
            ConstitutiveLaw,
            t.Literal['elastic', 'parabolarectangle', 'sargin', 'popovics'],
        ],
    ) -> None:
        """Setter for constitutive law.

        Arguments:
            consitutive_law (ConstitutiveLaw | str): A valid ConstitutiveLaw
                object for concrete or a string defining a valid constitutive
                law type for concrete. (valid options for string: 'elastic',
                'parabolarectangle', 'bilinearcompression', 'sargin',
                'popovics').
        """
        if constitutive_law is None:
            raise ValueError(
                'At least a constitutive law or a string defining the '
                'constitutive law must be provided.'
            )
        if isinstance(constitutive_law, str):
            constitutive_law = create_constitutive_law(
                constitutive_law_name=constitutive_law, material=self
            )

        if isinstance(constitutive_law, ConstitutiveLaw):
            if 'concrete' in constitutive_law.__materials__:
                self._constitutive_law = constitutive_law
            else:
                raise ValueError(
                    'The constitutive law selected is not suitable '
                    'for being used with a concrete material.'
                )
        else:
            raise ValueError(
                f'The constitutive law {constitutive_law} could not be created'
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
