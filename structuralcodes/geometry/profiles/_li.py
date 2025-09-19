"""L profiles."""

from shapely import (
    Polygon,
)

from ._base_profile import BaseProfile
from ._common_functions import (
    _create_L_section,
)


class LI(BaseProfile):
    """Simple class for representing a LI profile.

    European standard unequal side angle.

    """

    parameters = {
        'L120x80x8': {'h': 120.0, 'b': 80.0, 't': 8.0, 'r1': 11.0, 'r2': 5.5},
        'L120x80x10': {
            'h': 120.0,
            'b': 80.0,
            't': 10.0,
            'r1': 11.0,
            'r2': 5.5,
        },
        'L120x80x12': {
            'h': 120.0,
            'b': 80.0,
            't': 12.0,
            'r1': 11.0,
            'r2': 5.5,
        },
        'L150x75x9': {'h': 150.0, 'b': 75.0, 't': 9.0, 'r1': 12.0, 'r2': 6.0},
        'L150x75x10': {
            'h': 150.0,
            'b': 75.0,
            't': 10.0,
            'r1': 12.0,
            'r2': 6.0,
        },
        'L150x75x11': {
            'h': 150.0,
            'b': 75.0,
            't': 11.0,
            'r1': 12.0,
            'r2': 6.0,
        },
        'L150x75x12': {
            'h': 150.0,
            'b': 75.0,
            't': 12.0,
            'r1': 12.0,
            'r2': 6.0,
        },
        'L150x90x10': {
            'h': 150.0,
            'b': 90.0,
            't': 10.0,
            'r1': 12.0,
            'r2': 6.0,
        },
        'L150x90x11': {
            'h': 150.0,
            'b': 90.0,
            't': 11.0,
            'r1': 12.0,
            'r2': 6.0,
        },
        'L150x100x10': {
            'h': 150.0,
            'b': 100.0,
            't': 10.0,
            'r1': 12.0,
            'r2': 6.0,
        },
        'L150x100x12': {
            'h': 150.0,
            'b': 100.0,
            't': 12.0,
            'r1': 12.0,
            'r2': 6.0,
        },
        'L150x100x14': {
            'h': 150.0,
            'b': 100.0,
            't': 14.0,
            'r1': 12.0,
            'r2': 6.0,
        },
        'L200x100x10': {
            'h': 200.0,
            'b': 100.0,
            't': 10.0,
            'r1': 15.0,
            'r2': 7.5,
        },
        'L200x100x12': {
            'h': 200.0,
            'b': 100.0,
            't': 12.0,
            'r1': 15.0,
            'r2': 7.5,
        },
        'L200x100x14': {
            'h': 200.0,
            'b': 100.0,
            't': 14.0,
            'r1': 15.0,
            'r2': 7.5,
        },
    }

    @classmethod
    def get_polygon(cls, name: str) -> Polygon:
        """Returns a shapely polygon representing an L section."""
        parameters = cls.parameters.get(name)
        if parameters is None:
            raise ValueError(
                f"Profile '{name}' not found in L sections. "
                "Select a valid profile (available ones: "
                f"{cls.profiles()})"
            )
        parameters['t1'] = parameters['t']
        parameters['t2'] = parameters['t']
        return _create_L_section(
            **{
                key: parameters[key]
                for key in parameters
                if key in ['h', 'b', 't1', 't2', 'r1', 'r2']
            }
        )

    @classmethod
    def profiles(cls) -> list:
        """Returns a list containing all available profiles."""
        return list(cls.parameters.keys())

    def __init__(self, name: str) -> None:
        """Creates a new L object."""
        parameters = self.parameters.get(name)
        if parameters is None:
            raise ValueError(
                f"Profile '{name}' not found in LI sections. "
                "Select a valid profile (available ones: "
                f"{self.profiles()})"
            )
        super().__init__()
        self._h = parameters.get('h')
        self._b = parameters.get('b')
        self._t = parameters.get('t')
        self._r1 = parameters.get('r1')
        self._r2 = parameters.get('r2')
        self._polygon = _create_L_section(
            h=self._h,
            b=self._b,
            t1=self._t,
            t2=self._t,
            r1=self._r1,
            r2=self._r2,
        )

    @property
    def polygon(self) -> Polygon:
        """Returns shapely Polygon of section.

        Returns:
            Polygon: The represention of the L section.
        """
        return self._polygon

    @property
    def h(self) -> float:
        """Returns height of L section.

        Returns:
            float: Height h of L section.
        """
        return self._h

    @property
    def b(self) -> float:
        """Returns width of L section.

        Returns:
            float: Width b of L section.
        """
        return self._b

    @property
    def t(self) -> float:
        """Returns thickness of L section.

        Returns:
            float: Thickness t of L section.
        """
        return self._t

    @property
    def r1(self) -> float:
        """Returns fillet radius of L section.

        Returns:
            float: Fillet radius r1 of L section.
        """
        return self._r1

    @property
    def r2(self) -> float:
        """Returns fillet radius of L section.

        Returns:
            float: Fillet radius r2 of L section.
        """
        return self._r2
