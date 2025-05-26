"""Core implementation of the reinforcement material."""

import abc
import typing as t

from structuralcodes.core.base import ConstitutiveLaw, Material


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

    @property
    def fyk(self) -> float:
        """Returns fyk in MPa."""
        return self._fyk

    @property
    def Es(self) -> float:
        """Returns Es in MPa."""
        return self._Es

    @property
    def ftk(self) -> float:
        """Returns ftk in MPa."""
        return self._ftk

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
