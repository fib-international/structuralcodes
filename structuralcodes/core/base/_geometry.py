"""Abstract base class for geometries."""

from __future__ import annotations  # To have clean hints of ArrayLike in docs

import typing as t
from fnmatch import fnmatchcase

from ._material import Material


class _GroupMixin:
    """A mixin class for handling functionality related to group labels."""

    _group_label: t.Optional[str] = None

    @property
    def group_label(self):
        """Returns the group_label."""
        return self._group_label

    def _group_matches(
        self, pattern: str, *, case_sensitive: bool = True
    ) -> bool:
        """Checks if the group_label matches a pattern.

        Arguments:
            pattern (str): the string pattern to be checked

        Keyword Arguments:
            case_sensitive (bool, optional): if True (default) the check is
                case sensitive.

        Returns:
            (bool): Returns True if the group_label matches the pattern.

        Note:
            The matching permits to use:
                - "*" any chars
                - "?" single char
                - "[abc]" character set

        Examples:
            >>> geo.group_matches("nametos*")
            >>> geo.group_matches("*pier*")
            >>> geo.group_matches("Abutment??", case_senstive=False)
        """
        if self.group_label is None:
            return False
        if not case_sensitive:
            return fnmatchcase(self.group_label.casefold(), pattern.casefold())
        return fnmatchcase(self.group_label, pattern)


class Geometry:
    """Base class for a geometry object."""

    _geometry_counter: t.ClassVar[int] = 0
    id: int

    def __init__(
        self, name: t.Optional[str] = None, base_name: str = 'Geometry'
    ) -> None:
        """Initializes a geometry object.

        The name and grouplabel serve for filtering in a compound object. By
        default it creates a new name each time.

        Arguments:
            name (Optional(str)): The name to be given to the object.
            base_name (str): If name is not given, use this argument together
                with a global counter to create a unique name for the geometry.
        """
        self.id = Geometry.return_global_counter_and_increase()
        self._name = name if name is not None else f'{base_name}_{self.id}'

    @property
    def name(self):
        """Returns the name of the Geometry."""
        return self._name

    @classmethod
    def _increase_global_counter(cls):
        """Increases the global counter by one."""
        cls._geometry_counter += 1

    @classmethod
    def return_global_counter_and_increase(cls):
        """Returns the current counter and increases it by one."""
        counter = cls._geometry_counter
        cls._increase_global_counter()
        return counter

    @staticmethod
    def from_geometry(
        geo: Geometry,
        new_material: t.Optional[Material] = None,
    ) -> Geometry:
        """Create a new geometry with a different material."""
        raise NotImplementedError(
            'This method should be implemented by subclasses'
        )

    def _name_matches(
        self, pattern: str, *, case_sensitive: bool = True
    ) -> bool:
        """Checks if the name matches a pattern.

        Arguments:
            pattern (str): the string pattern to be checked

        Keyword Arguments:
            case_sensitive (bool, optional): if True (default) the check is
                case sensitive.

        Returns:
            (bool): Returns True if the name matches the pattern.

        Note:
            The matching permits to use:
                - "*" any chars
                - "?" single char
                - "[abc]" character set

        Examples:
            >>> geo.name_matches("nametos*")
            >>> geo.name_matches("*pier*")
            >>> geo.name_matches("Abutment??", case_senstive=False)
        """
        if not case_sensitive:
            return fnmatchcase(self.name.casefold(), pattern.casefold())
        return fnmatchcase(self.name, pattern)
