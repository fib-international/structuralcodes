"""Concrete material"""
import typing as t
from structuralcodes.code import _use_design_code
from ._concrete import Concrete
from ._concreteMC2010 import ConcreteMC2010


def create_concrete(
    fck: float,
    name: t.Optional[str] = None,
    density: float = 2400.0,
    existing: bool = False,
    design_code: t.Optional[str] = None,
) -> t.Optional[Concrete]:
    """
    A factory function to create the correct type of concrete based on the
    desired design code.

    Args:
        fck (float): Characteristic strength of concrete in MPa.
            (if existing it is intended as the mean strength)

    Kwargs:
        density (float): Density of Concrete in kg/m3 (default: 2400)
        existing (bool): Boolean indicating if the concrete is of an
            existing structure (default: False)
        deisgn_code (str): Optional string (default: None) indicating the
            desired standard. If None (default) the globally used design
            standard will be adopted. Otherwise the design standard specified
            will be used for the instance of the material.
            Currently available codes: 'mc2010'

    Raises:
        ValueError: if the design code is not valid or does not cover
            concrete as a material.
    """
    # Get the code from the global variable
    _code = _use_design_code(design_code)
    # Check if the code is a proper concrete code
    code = _code if 'concrete' in _code.__materials__ else None
    if code is None:
        raise ValueError(
            'The design code is not set, either use '
            'structuralcodes.code.set_designcode, or provide a valid '
            'string in the function.'
        )
    # Create the proper concrete object
    if code.__title__ == 'fib Model Code 2010':
        return ConcreteMC2010(fck, name, density, existing)
    return None


__all__ = [
    'create_concrete',
    'Concrete',
    'ConcreteMC2010',
]
