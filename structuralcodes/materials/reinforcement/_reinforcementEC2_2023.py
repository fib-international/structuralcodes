"""The concrete class for Eurocode 2 2023 Reinforcement Material."""

import typing as t

from structuralcodes.codes import ec2_2023
from structuralcodes.core._units import UnitSet

from ._reinforcement import Reinforcement


class ReinforcementEC2_2023(Reinforcement):  # noqa: N801
    """Reinforcement implementation for EC2 2023."""

    _default_units = UnitSet(length='mm', force='N')

    def __init__(
        self,
        fyk: float,
        Es: float,
        ftk: float,
        epsuk: float,
        gamma_s: t.Optional[float] = None,
        name: t.Optional[str] = None,
        density: float = 7850.0,
        units: t.Optional[UnitSet] = None,
    ):
        """Initializes a new instance of Reinforcement for EC2 2023.

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
            ec2_2023.fyd(
                self.unit_converter.convert_stress_backwards(self.fyk),
                self.gamma_s,
            )
        )

    def ftd(self) -> float:
        """The design ultimate strength."""
        return self.unit_converter.convert_stress_forwards(
            ec2_2023.fyd(
                self.unit_converter.convert_stress_backwards(self.ftk),
                self.gamma_s,
            )
        )

    def epsud(self) -> float:
        """The design ultimate strain."""
        return ec2_2023.eps_ud(self.epsuk, self.gamma_s)

    @property
    def gamma_s(self) -> float:
        """The partial factor for reinforcement."""
        # Here we should implement the interaction with the globally set
        # national annex. For now, we simply return the default value.
        return self._gamma_s or 1.15

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
