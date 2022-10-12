"""Factory for design codes"""
import typing as t

from structuralcodes.core.base import DesignCode

from .mc2010.mc2010 import MC2010Builder


class DesignCodeFactory:
    """Factory for design codes."""

    def __init__(self) -> None:
        """Initialize factory."""
        self._builders: t.Dict = {}

    def register_code(self, code: str, builder: t.Any) -> None:
        """Register the code builder of a new code.

        Args:
            code (str): The abbreviation of the code.
            builder (-): The builder for the code.
        """
        self._builders[code] = builder

    def create(self, code: str, **kwargs) -> DesignCode:
        """Create the code using the code builder.

        Args:
            code (str): The abbreviation of the code.

        Kwargs:
            **kwargs: The keyword arguments of the code builder for the code.
        """
        builder = self._builders.get(code)
        if builder is None:
            raise ValueError(code)
        return builder(**kwargs)

    def get_registered_codes(self) -> t.List[str]:
        """Get a list of the registered design code builders."""
        return list(self._builders.keys())


# Instantiate the code factory and register codes
code_factory = DesignCodeFactory()
code_factory.register_code('mc2010', MC2010Builder())
