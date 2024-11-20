"""The concrete class for Model Code 2010 Reinforcement Material."""

import typing as t

from structuralcodes.codes import mc2010
from structuralcodes.core._units import UnitSet

from ._reinforcement import Reinforcement


class ReinforcementMC2010(Reinforcement):
    """Reinforcement implementation for MC 2010."""

    _default_units = UnitSet(length='mm', force='N')

    def __init__(
        self,
        fyk: float,
        Es: float,
        ftk: float,
        epsuk: float,
        gamma_s: t.Optional[float] = None,
        gamma_eps: t.Optional[float] = None,
        name: t.Optional[str] = None,
        density: float = 7850.0,
        units: t.Optional[UnitSet] = None,
    ):
        """Initializes a new instance of Reinforcement for MC2010.

        Args:
            fyk (float): Characteristic yield strength.
            Es (float): The Young's modulus.
            ftk (float): Characteristic ultimate strength.
            epsuk (float): The characteristik strain at the ultimate stress
                level.
            gamma_s (Optional(float)): The partial factor for reinforcement.
                Default value is 1.15.

        Keyword Args:
            name (str): A descriptive name for the reinforcement.
            density (float): Density of material in kg/m3 (default: 7850).
            units (Optional[UnitSet]): The selected set of units to work in.
                The default is length=m and force=N.

        Note:
            The arguments should be provided compatible with the selected set
            of units.
        """
        if name is None:
            name = f'Reinforcement{round(fyk):d}'
        self._gamma_eps = gamma_eps
        super().__init__(
            fyk=fyk,
            Es=Es,
            name=name,
            density=density,
            ftk=ftk,
            epsuk=epsuk,
            gamma_s=gamma_s,
            units=units,
        )

    def fyd(self) -> float:
        """The design yield strength."""
        return self.unit_converter.convert_stress_forwards(
            mc2010.fyd(
                self.unit_converter.convert_stress_backwards(self.fyk),
                self.gamma_s,
            )
        )

    @property
    def gamma_s(self) -> float:
        """The partial factor for reinforcement."""
        return self._gamma_s or 1.15

    def ftd(self) -> float:
        """The design ultimate strength."""
        return self.unit_converter.convert_stress_forwards(
            mc2010.fyd(
                self.unit_converter.convert_stress_backwards(self.ftk),
                self.gamma_s,
            )
        )

    def epsud(self) -> float:
        """The design ultimate strain."""
        return mc2010.epsud(self.epsuk, self.gamma_eps)

    @property
    def gamma_eps(self) -> float:
        """The partial factor for ultimate strain."""
        return self._gamma_eps or 0.9

    def __elastic__(self) -> dict:
        """Returns kwargs for creating an elastic constitutive law."""
        return {'E': self.Es}

    def __elasticperfectlyplastic__(self) -> dict:
        """Returns kwargs for ElasticPlastic constitutive law with no strain
        hardening.
        """
        return {
            'E': self.Es,
            'fy': self.fyd(),
            'eps_su': self.epsud(),
        }

    def __elasticplastic__(self) -> dict:
        """Returns kwargs for ElasticPlastic constitutive law with strain
        hardening.
        """
        Eh = (self.ftd() - self.fyd()) / (self.epsuk - self.epsyd)
        return {
            'E': self.Es,
            'fy': self.fyd(),
            'Eh': Eh,
            'eps_su': self.epsud(),
        }
