"""Collection of functions related to design codes"""
import typing as t

from structuralcodes.core.base import DesignCode

from ._code_factory import code_factory

# Global code object used by material classes
_CODE: t.Optional[DesignCode] = None


def set_design_code(design_code: str) -> None:
    """Set the current design code globally.

    Args:
        design_code (str): The abbreviation of the code.

    Note:
        Call get_design_codes() to get a list of the available codes.
    """
    global _CODE  # pylint: disable=W0603
    _CODE = code_factory.create(design_code)


def get_design_codes() -> t.List[str]:
    """Get a list of the available design codes."""
    return code_factory.get_registered_codes()


def set_national_annex(national_annex: str) -> None:
    """Set the current national annex globally."""
    raise NotImplementedError


def _use_design_code(
    design_code: t.Optional[str] = None,
) -> t.Optional[DesignCode]:
    """Use a design code in a class"""
    if design_code is None:
        return _CODE
    return code_factory.create(design_code)
