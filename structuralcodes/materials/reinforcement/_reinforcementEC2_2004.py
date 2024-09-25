"""The concrete class for Eurocode 2 2004 Reinforcement Material."""

import typing as t

from structuralcodes.codes import ec2_2004

from ._reinforcement import Reinforcement


class ReinforcementEC2_2004(Reinforcement):  # noqa: N801
    """Reinforcement implementation for EC2 2004."""

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
    ):
        """Initializes a new instance of Reinforcement for EC2 2004.

        Arguments:
            fyk (float): Characteristic yield strength in MPa.
            Es (float): The Young's modulus in MPa.
            ftk (float): Characteristic ultimate strength in MPa.
            epsuk (float): The characteristik strain at the ultimate stress
                level.
            gamma_s (Optional(float)): The partial factor for reinforcement.
                Default value is 1.15.

        Keyword Arguments:
            name (str): A descriptive name for the reinforcement.
            density (float): Density of material in kg/m3 (default: 7850).
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
        )

    def fyd(self) -> float:
        """The design yield strength."""
        return ec2_2004.fyd(self.fyk, self.gamma_s)

    @property
    def gamma_s(self) -> float:
        """The partial factor for reinforcement."""
        # Here we should implement the interaction with the globally set
        # national annex. For now, we simply return the default value.
        return self._gamma_s or 1.15

    def ftd(self) -> float:
        """The design ultimate strength."""
        return ec2_2004.fyd(self.ftk, self.gamma_s)

    def epsud(self) -> float:
        """The design ultimate strain."""
        return ec2_2004.epsud(self.epsuk, self.gamma_eps)

    @property
    def gamma_eps(self) -> float:
        """The partial factor for ultimate strain."""
        # Here we should implement the interaction with the globally set
        # national annex. For now, we simply return the default value.
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
