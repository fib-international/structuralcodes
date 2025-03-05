"""Core implementation of the reinforcement material."""

import abc
import typing as t

from structuralcodes.core.base import ConstitutiveLaw, Material
from structuralcodes.materials.constitutive_laws import (
    ElasticPlastic,
    create_constitutive_law,
)


class Reinforcement(Material):
    """The abstract reinforcement material."""

    _fyk: float
    _Es: float
    _ftk: float
    _epsuk: float
    _gamma_s: t.Optional[float] = None

    def __init__(
        self,
        fyk: float,
        Es: float,
        density: float,
        ftk: float,
        epsuk: float,
        gamma_s: t.Optional[float] = None,
        name: t.Optional[str] = None,
    ) -> None:
        """Initializes an abstract reinforcement material."""
        name = name if name is not None else 'Reinforcement'
        super().__init__(density, name)

        self._fyk = abs(fyk)
        self._Es = abs(Es)
        self._ftk = abs(ftk)
        self._epsuk = abs(epsuk)
        self._gamma_s = gamma_s

        # Calculate the plastic hardening modulus
        if self.epsuk - self.fyk / self.Es > 0.0 and self.ftk - self.fyk > 0:
            Eh = (self.ftk - self.fyk) / (self.epsuk - self.fyk / self.Es)
        else:
            Eh = 0.0
        self._constitutive_law = ElasticPlastic(
            E=self.Es, fy=self.fyd(), Eh=Eh, eps_su=self.epsud()
        )

    @property
    def fyk(self) -> float:
        """Returns fyk in MPa."""
        return self._fyk

    @fyk.setter
    def fyk(self, fyk: float) -> None:
        """Setter for fyk (in MPa)."""
        self._fyk = abs(fyk)

    @property
    def Es(self) -> float:
        """Returns Es in MPa."""
        return self._Es

    @Es.setter
    def Es(self, Es: float) -> None:
        """Setter for Es (in MPa)."""
        self._Es = abs(Es)

    @property
    def ftk(self) -> float:
        """Returns ftk in MPa."""
        return self._ftk

    @ftk.setter
    def ftk(self, ftk: float) -> None:
        """Setter for ftk (in MPa)."""
        self._ftk = abs(ftk)

    @property
    def epsuk(self) -> float:
        """Returns epsuk."""
        return self._epsuk

    @epsuk.setter
    def epsuk(self, epsuk: float) -> None:
        """Setter for epsuk."""
        self._epsuk = abs(epsuk)

    @property
    def epsyk(self) -> float:
        """Returns characteristic yield strain epsyk."""
        return self.fyk / self.Es

    @property
    def constitutive_law(self) -> ConstitutiveLaw:
        """Returns the constitutive law object."""
        return self._constitutive_law

    @constitutive_law.setter
    def constitutive_law(
        self,
        constitutive_law: t.Union[
            ConstitutiveLaw,
            t.Literal[
                'elastic',
                'elasticperfecltyplastic',
                'elasticplastic',
            ],
        ],
    ) -> None:
        """Setter for constitutive law.

        Arguments:
            consitutive_law (ConstitutiveLaw | str): a valid ConstitutiveLaw
                object for reinforcement or  a string defining a valid
                constitutive law type for reinforcement. (valid options:
                'elastic', 'elasticperfectlyplastic', 'elasticplastic').
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
            if 'rebars' in constitutive_law.__materials__:
                self._constitutive_law = constitutive_law
            else:
                raise ValueError(
                    'The constitutive law selected is not suitable '
                    'for being used with a Reinforcement material.'
                )
        else:
            raise ValueError(
                f'The constitutive law {constitutive_law} could not be created'
            )

    @property
    @abc.abstractmethod
    def gamma_s(self) -> float:
        """Each reinforcement should implement its own getter for the partial
        factor in order to interact with the globally set national annex.
        """

    @abc.abstractmethod
    def fyd(self) -> float:
        """Each reinforcement should implement its own method for calculating
        the design yield strength.
        """

    @abc.abstractmethod
    def ftd(self) -> float:
        """Each reinforcement should implement its own method for calculating
        the design ultimate strength.
        """

    @abc.abstractmethod
    def epsud(self) -> float:
        """Each reinforcement should implement its own method for calculating
        the design ultimate strain.
        """

    @property
    def epsyd(self) -> float:
        """Return the design yield strain."""
        return self.fyd() / self.Es
