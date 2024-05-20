"""Core implementation of the reinforcement material."""

import typing as t

from structuralcodes.core.base import ConstitutiveLaw, Material
from structuralcodes.materials.constitutive_laws import ElasticPlastic


class Reinforcement(Material):
    """The abstract reinforcement material."""

    _fyk: float
    _Es: float
    _ftk: float
    _epsuk: float

    def __init__(
        self,
        fyk: float,
        Es: float,
        density: float,
        ftk: float,
        epsuk: float,
        name: t.Optional[str] = None,
    ) -> None:
        """Initializes an abstract reinforcement material."""
        name = name if name is not None else 'Reinforcement'
        super().__init__(density, name)

        self._fyk = abs(fyk)
        self._Es = abs(Es)
        self._ftk = abs(ftk)
        self._epsuk = abs(epsuk)

        # Calculate the plastic hardening modulus
        if self.epsuk - self.fyk / self.Es > 0.0 and self.ftk - self.fyk > 0:
            Eh = (self.ftk - self.fyk) / (self.epsuk - self.fyk / self.Es)
        else:
            Eh = 0.0
        self._constitutive_law = ElasticPlastic(
            E=self.Es, fy=self.fyk, Eh=Eh, eps_su=self.epsuk
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
    def constitutive_law(self) -> ConstitutiveLaw:
        """Returns the constitutive law object."""
        return self._constitutive_law

    @constitutive_law.setter
    def constitutive_law(self, constitutive_law: ConstitutiveLaw) -> None:
        """Setter for constitutive law."""
        if 'rebars' in constitutive_law.__materials__:
            self._constitutive_law = constitutive_law
        else:
            raise ValueError(
                'The constitutive law selected is not suitable '
                'for being used with a reinforcement material.'
            )
