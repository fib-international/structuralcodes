"""The concrete class for Eurocode 2 2004 Reinforcement Material."""

import typing as t

from ._reinforcement import Reinforcement


class ReinforcementEC2_2004(Reinforcement):  # noqa: N801
    """Reinforcement implementation for EC2 2004."""

    def __init__(
        self,
        fyk: float,
        Es: float,
        name: t.Optional[str] = None,
        density: float = 7850.0,
    ):
        """Initializes a new instance of Reinforcement for EC2 2004.

        Args:
            fyk (float): Characteristic yield strength in MPa.
            Es (float): The Young's modulus in MPa.

        Keyword Args:
            name (str): A descriptive name for the reinforcement.
            desnsity (float): Density of material in kg/m3 (default: 7850).
        """
        if name is None:
            name = f'Reinforcement{round(fyk):d}'
        super().__init__(fyk=fyk, Es=Es, name=name, density=density)
