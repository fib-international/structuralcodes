"""Abstract base class for geometries."""

from __future__ import annotations  # To have clean hints of ArrayLike in docs

import typing as t
from fnmatch import fnmatchcase

from ._material import Material


class Geometry:
    """Base class for a geometry object."""

    section_counter: t.ClassVar[int] = 0

    def __init__(
        self, name: t.Optional[str] = None, group_label: t.Optional[str] = None
    ) -> None:
        """Initializes a geometry object.

        The name and grouplabel serve for filtering in a compound object. By
        default it creates a new name each time.

        Arguments:
            name (Optional(str)): The name to be given to the object.
            group_label (Optional(str)): A label for grouping several objects.
        """
        if name is not None:
            self._name = name
        else:
            counter = Geometry.return_global_counter_and_increase()
            self._name = f'Geometry_{counter}'
        self._group_label = group_label

    @property
    def name(self):
        """Returns the name of the Geometry."""
        return self._name

    @property
    def group_label(self):
        """Returns the group_label fo the Geometry."""
        return self._group_label

    @classmethod
    def _increase_global_counter(cls):
        """Increases the global counter by one."""
        cls.section_counter += 1

    @classmethod
    def return_global_counter_and_increase(cls):
        """Returns the current counter and increases it by one."""
        counter = cls.section_counter
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
