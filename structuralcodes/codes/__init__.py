"""Collection of functions related to design codes."""

import types
import typing as t

from . import ec2_2004, ec2_2023, mc2010

__all__ = [
    'mc2010',
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
    'ec2_2004': ec2_2004,
    'ec2_2023': ec2_2023,
}


def set_design_code(design_code: t.Optional[str] = None) -> None:
    """Set the current design code globally.

    Args:
        design_code (str): The abbreviation of the code.

    Note:
        Call get_design_codes() to get a list of the available codes.
    """
    global _CODE  # pylint: disable=W0603
    if design_code is not None:
        _CODE = _DESIGN_CODES.get(design_code.lower())
    else:
        _CODE = None


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
        return _CODE
    return _DESIGN_CODES.get(design_code.lower())
