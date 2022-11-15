"""Abstract base classes"""
import abc


class Material(abc.ABC):
    """Abstract base class for materials."""

    def __init__(self, name: str, density: float):
        '''Initialize a material'''
        self._name = name
        self._density = density


class DesignCode(abc.ABC):
    """Abstract base class for design codes."""

    def __init__(self, name: str, release_year: int, materials: tuple):
        self._name = name
        self._release_year = release_year
        self._materials = materials

    @property
    def name(self) -> str:
        '''Getting the name of the standard'''
        return self._name

    @property
    def release_year(self) -> int:
        '''Getting the release year of the standard'''
        return self._release_year

    @property
    def materials(self) -> tuple:
        '''Getting the materials to which stadard refer'''
        return self._materials
