"""Collection of functions related to design codes."""

import types
import typing as t

from . import ec2_2004, ec2_2023, mc2010, mc2020

__all__ = [
    'mc2010',
    'mc2020',
    'ec2_2023',
    'ec2_2004',
    'set_design_code',
    'get_design_codes',
    'set_national_annex',
]

# Global code object used by material classes
_CODE: t.Optional[types.ModuleType] = None

# Global national annex object
_NATIONAL_ANNEX: t.Optional[str] = None

# Design code registry
_DESIGN_CODES = {
    'mc2010': mc2010,
    'mc2020': mc2020,
    'ec2_2004': ec2_2004,
    'ec2_2023': ec2_2023,
}


def set_design_code(
    design_code: t.Optional[t.Union[str, types.ModuleType]] = None,
) -> None:
    """Set the current design code globally.

    Args:
        design_code (Union[str, Moduletype]): The abbreviation of the code
            (str), or a module that represents the code (ModuleType).

    Note:
        Call get_design_codes() to get a list of the available codes.
    """
    global _CODE  # pylint: disable=W0603
    if design_code is None:
        # Reset to None
        _CODE = None
    elif isinstance(design_code, str) and design_code.lower() in _DESIGN_CODES:
        # The design code abbreviation is valid
        _CODE = _DESIGN_CODES.get(design_code.lower())
    elif isinstance(design_code, types.ModuleType) and all(
        name in dir(design_code)
        for name in ('__title__', '__year__', '__materials__')
    ):
        # The module is a valid design code
        _CODE = design_code
    else:
        raise ValueError(
            f'{design_code} is not a valid abbreviation for a design code, or'
            ' a valid module representing a design code.\nType '
            'get_design_codes() to list the available design codes.'
        )


def get_design_codes() -> t.List[str]:
    """Get a list of the available design codes."""
    return list(_DESIGN_CODES.keys())


def set_national_annex(national_annex: str) -> None:
    """Set the current national annex globally.

    Args:
        national_annex (str): The abbreviation of the national annex.

    Note:
        Call get_national_annexes() on the relevant design code to see a list
        of available national annexes.
    """
    global _NATIONAL_ANNEX  # pylint: disable=W0603
    _NATIONAL_ANNEX = national_annex.lower()


def _use_design_code(
    design_code: t.Optional[str] = None,
) -> t.Optional[types.ModuleType]:
    """Use a design code in a class.

    Kwargs:
        design_code (str): The abbreviation of the code.

    Note:
        Call get_design_codes() to get a list of the available codes.
    """
    if design_code is None:
        # Returned the globally set design code
        return _CODE

    # Set design code before returning
    set_design_code(design_code)
    return _CODE
