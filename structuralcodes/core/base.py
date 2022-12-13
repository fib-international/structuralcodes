"""Abstract base classes"""
import abc
import typing as t


class Material(abc.ABC):
    """Abstract base class for materials."""

    def __init__(self, density: float, name: t.Optional[str] = None) -> None:
        """Initializes an instance of a new material

        Args:
            density (float): density of the material in kg/m3

        Keyword Args:
            name (Optional[str]): descriptive name of the material
        """
        self._density = abs(density)
        self._name = name if name is not None else "Material"

    @property
    def name(self):
        """Returns the name of the material"""
        return self._name

    @property
    def density(self):
        """Returns the density of the material in kg/m3"""
        return self._density
