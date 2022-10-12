"""Core implementation of the concrete material"""
import typing as t

# from structuralcodes.code import _CODE
from structuralcodes.code import _use_design_code
from structuralcodes.core.base import Material, DesignCode


REQUIRED_FUNCTIONS = (
    'fcm',
    'fctm',
)


class Concrete(Material):
    """The concrete material."""

    def __init__(self, design_code: t.Optional[str] = None) -> None:
        super().__init__()

        _code = _use_design_code(design_code)
        code = _code if isinstance(_code, DesignCode) else None

        if code is None:
            raise Exception(
                'The design code is not set, either use '
                'structuralcodes.code.set_designcode, or provide a valid '
                'string in the material constructor.'
            )

        # Set attributes from the design code
        for fun in REQUIRED_FUNCTIONS:
            setattr(self, fun, code.__getattribute__(fun))
