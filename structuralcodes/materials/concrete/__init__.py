"""Concrete material."""

import typing as t

from structuralcodes.codes import _use_design_code

from ._concrete import Concrete
from ._concreteEC2_2004 import ConcreteEC2_2004
from ._concreteEC2_2023 import ConcreteEC2_2023
from ._concreteMC2010 import ConcreteMC2010

__all__ = [
    'create_concrete',
    'Concrete',
    'ConcreteMC2010',
    'ConcreteEC2_2023',
    'ConcreteEC2_2004',
]

CONCRETES: t.Dict[str, Concrete] = {
    'fib Model Code 2010': ConcreteMC2010,
    'EUROCODE 2 1992-1-1:2004': ConcreteEC2_2004,
    'EUROCODE 2 1992-1-1:2023': ConcreteEC2_2023,
}


def create_concrete(
    fck: float,
    name: t.Optional[str] = None,
    density: float = 2400.0,
    gamma_c: t.Optional[float] = None,
    existing: bool = False,
    design_code: t.Optional[str] = None,
    **kwargs,
) -> t.Optional[Concrete]:
    """A factory function to create the correct type of concrete based on the
    desired design code.

    Arguments:
        fck (float): Characteristic strength of concrete in MPa. (if existing
            it is intended as the mean strength).

    Keyword Arguments:
        density (float): Density of Concrete in kg/m3 (default: 2400).
        gamma_c (Optional(float)): The partial factor for concrete.
        existing (bool): Boolean indicating if the concrete is of an existing
            structure (default: False).
        design_code (str): Optional string (default: None) indicating the
            desired standard. If None (default) the globally used design
            standard will be adopted. Otherwise the design standard specified
            will be used for the instance of the material.

    Raises:
        ValueError: if the design code is not valid or does not cover concrete
            as a material.
    """
    # Get the code from the global variable
    _code = _use_design_code(design_code)

    # Check if the code is a proper concrete code
    code = None
    if _code is not None and 'concrete' in _code.__materials__:
        code = _code
    if code is None:
        raise ValueError(
            'The design code is not set, either use '
            'structuralcodes.code.set_designcode, or provide a valid '
            'string in the function.'
        )

    # Create the proper concrete object
    current_concrete = CONCRETES.get(code.__title__, None)
    if current_concrete is not None:
        return current_concrete(
            fck=fck,
            name=name,
            density=density,
            gamma_c=gamma_c,
            existing=existing,
            **kwargs,
        )
    return None
