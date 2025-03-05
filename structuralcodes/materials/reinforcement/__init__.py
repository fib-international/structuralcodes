"""Reinforcement material."""

import typing as t

from structuralcodes.codes import _use_design_code

from ._reinforcement import Reinforcement
from ._reinforcementEC2_2004 import ReinforcementEC2_2004
from ._reinforcementEC2_2023 import ReinforcementEC2_2023
from ._reinforcementMC2010 import ReinforcementMC2010

__all__ = [
    'create_reinforcement',
    'Reinforcement',
    'ReinforcementMC2010',
    'ReinforcementEC2_2004',
    'ReinforcementEC2_2023',
]

REINFORCEMENTS: t.Dict[str, Reinforcement] = {
    'fib Model Code 2010': ReinforcementMC2010,
    'EUROCODE 2 1992-1-1:2004': ReinforcementEC2_2004,
    'EUROCODE 2 1992-1-1:2023': ReinforcementEC2_2023,
}


def create_reinforcement(
    fyk: float,
    Es: float,
    ftk: float,
    epsuk: float,
    gamma_s: t.Optional[float] = None,
    name: t.Optional[str] = None,
    density: float = 7850,
    design_code: t.Optional[str] = None,
) -> t.Optional[Reinforcement]:
    """A factory function to create the correct type of reinforcement based on
    the desired design code.

    Arguments:
        fyk (float): Characteristic yield strength in MPa.
        Es (float): The Young's modulus in MPa.
        ftk (float): Characteristic ultimate strength in MPa.
        epsuk (float): The characteristik strain at the ultimate stress level.

    Keyword Arguments:
        gamma_s (Optional(float)): The partial factor for reinforcement.
        density (float): Density of the material in kg/m3 (default: 7850)
        design_code (str): Optional string (default: None) indicating the
            desired standard. If None (default) the globally used design
            standard will be adopted. Otherwise the design standard specified
            will be used for the instance of the material.

    Raises:
        ValueError: If the design code is not valid or does not cover
            reinforcement as a material.
    """
    # Get the code from the global variable
    _code = _use_design_code(design_code)

    # Check if the code is a proper code for reinforcement
    code = None
    if _code is not None and 'reinforcement' in _code.__materials__:
        code = _code
    if code is None:
        raise ValueError(
            'The design code is not set, either use '
            'structuralcodes.code.set_designcode, or provide a valid '
            'string in the function.'
        )

    # Create the proper reinforcement object
    current_reinforcement = REINFORCEMENTS.get(code.__title__, None)
    if current_reinforcement is not None:
        return current_reinforcement(
            fyk=fyk,
            Es=Es,
            name=name,
            density=density,
            ftk=ftk,
            epsuk=epsuk,
            gamma_s=gamma_s,
        )
    return None
