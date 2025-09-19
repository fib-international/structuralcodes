"""L profiles."""

from shapely import (
    Polygon,
)

from ._base_profile import BaseProfile
from ._common_functions import (
    _create_L_section,
)


class L(BaseProfile):
    """Simple class for representing a L profile.

    European standard equal side angle L20x20 - L250x250.

    """

    parameters = {
        'L20x20x3': {'h': 20.0, 'b': 20.0, 't': 3.0, 'r1': 3.5, 'r2': 2.0},
        'L25x25x3': {'h': 25.0, 'b': 25.0, 't': 3.0, 'r1': 3.5, 'r2': 2.0},
        'L25x25x4': {'h': 25.0, 'b': 25.0, 't': 4.0, 'r1': 3.5, 'r2': 2.0},
        'L30x30x3': {'h': 30.0, 'b': 30.0, 't': 3.0, 'r1': 5.0, 'r2': 2.5},
        'L30x30x4': {'h': 30.0, 'b': 30.0, 't': 4.0, 'r1': 5.0, 'r2': 2.5},
        'L35x35x4': {'h': 35.0, 'b': 35.0, 't': 4.0, 'r1': 5.0, 'r2': 2.5},
        'L40x40x4': {'h': 40.0, 'b': 40.0, 't': 4.0, 'r1': 6.0, 'r2': 3.0},
        'L40x40x5': {'h': 40.0, 'b': 40.0, 't': 5.0, 'r1': 6.0, 'r2': 3.0},
        'L45x45x4.5': {'h': 45.0, 'b': 45.0, 't': 4.5, 'r1': 7.0, 'r2': 3.5},
        'L50x50x4': {'h': 50.0, 'b': 50.0, 't': 4.0, 'r1': 7.0, 'r2': 3.5},
        'L50x50x5': {'h': 50.0, 'b': 50.0, 't': 5.0, 'r1': 7.0, 'r2': 3.5},
        'L50x50x6': {'h': 50.0, 'b': 50.0, 't': 6.0, 'r1': 7.0, 'r2': 3.5},
        'L60x60x5': {'h': 60.0, 'b': 60.0, 't': 5.0, 'r1': 8.0, 'r2': 4.0},
        'L60x60x6': {'h': 60.0, 'b': 60.0, 't': 6.0, 'r1': 8.0, 'r2': 4.0},
        'L60x60x8': {'h': 60.0, 'b': 60.0, 't': 8.0, 'r1': 8.0, 'r2': 4.0},
        'L65x65x7': {'h': 65.0, 'b': 65.0, 't': 7.0, 'r1': 9.0, 'r2': 4.5},
        'L70x70x6': {'h': 70.0, 'b': 70.0, 't': 6.0, 'r1': 9.0, 'r2': 4.5},
        'L70x70x7': {'h': 70.0, 'b': 70.0, 't': 7.0, 'r1': 9.0, 'r2': 4.5},
        'L75x75x6': {'h': 75.0, 'b': 75.0, 't': 6.0, 'r1': 10.0, 'r2': 5.0},
        'L75x75x8': {'h': 75.0, 'b': 75.0, 't': 8.0, 'r1': 10.0, 'r2': 5.0},
        'L80x80x8': {'h': 80.0, 'b': 80.0, 't': 8.0, 'r1': 10.0, 'r2': 5.0},
        'L80x80x10': {'h': 80.0, 'b': 80.0, 't': 10.0, 'r1': 10.0, 'r2': 5.0},
        'L90x90x7': {'h': 90.0, 'b': 90.0, 't': 7.0, 'r1': 11.0, 'r2': 5.5},
        'L90x90x8': {'h': 90.0, 'b': 90.0, 't': 8.0, 'r1': 11.0, 'r2': 5.5},
        'L90x90x9': {'h': 90.0, 'b': 90.0, 't': 9.0, 'r1': 11.0, 'r2': 5.5},
        'L90x90x10': {'h': 90.0, 'b': 90.0, 't': 10.0, 'r1': 11.0, 'r2': 5.5},
        'L100x100x8': {
            'h': 100.0,
            'b': 100.0,
            't': 8.0,
            'r1': 12.0,
            'r2': 6.0,
        },
        'L100x100x10': {
            'h': 100.0,
            'b': 100.0,
            't': 10.0,
            'r1': 12.0,
            'r2': 6.0,
        },
        'L100x100x12': {
            'h': 100.0,
            'b': 100.0,
            't': 12.0,
            'r1': 12.0,
            'r2': 6.0,
        },
        'L110x110x10': {
            'h': 110.0,
            'b': 110.0,
            't': 10.0,
            'r1': 13.0,
            'r2': 6.5,
        },
        'L110x110x12': {
            'h': 110.0,
            'b': 110.0,
            't': 12.0,
            'r1': 13.0,
            'r2': 6.5,
        },
        'L120x120x10': {
            'h': 120.0,
            'b': 120.0,
            't': 10.0,
            'r1': 13.0,
            'r2': 6.5,
        },
        'L120x120x11': {
            'h': 120.0,
            'b': 120.0,
            't': 11.0,
            'r1': 13.0,
            'r2': 6.5,
        },
        'L120x120x12': {
            'h': 120.0,
            'b': 120.0,
            't': 12.0,
            'r1': 13.0,
            'r2': 6.5,
        },
        'L120x120x13': {
            'h': 120.0,
            'b': 120.0,
            't': 13.0,
            'r1': 13.0,
            'r2': 6.5,
        },
        'L120x120x15': {
            'h': 120.0,
            'b': 120.0,
            't': 15.0,
            'r1': 13.0,
            'r2': 6.5,
        },
        'L130x130x12': {
            'h': 130.0,
            'b': 130.0,
            't': 12.0,
            'r1': 14.0,
            'r2': 7.0,
        },
        'L150x150x10': {
            'h': 150.0,
            'b': 150.0,
            't': 10.0,
            'r1': 16.0,
            'r2': 8.0,
        },
        'L150x150x12': {
            'h': 150.0,
            'b': 150.0,
            't': 12.0,
            'r1': 16.0,
            'r2': 8.0,
        },
        'L150x150x14': {
            'h': 150.0,
            'b': 150.0,
            't': 14.0,
            'r1': 16.0,
            'r2': 8.0,
        },
        'L150x150x15': {
            'h': 150.0,
            'b': 150.0,
            't': 15.0,
            'r1': 16.0,
            'r2': 8.0,
        },
        'L150x150x18': {
            'h': 150.0,
            'b': 150.0,
            't': 18.0,
            'r1': 16.0,
            'r2': 8.0,
        },
        'L160x160x14': {
            'h': 160.0,
            'b': 160.0,
            't': 14.0,
            'r1': 17.0,
            'r2': 8.5,
        },
        'L160x160x15': {
            'h': 160.0,
            'b': 160.0,
            't': 15.0,
            'r1': 17.0,
            'r2': 8.5,
        },
        'L160x160x16': {
            'h': 160.0,
            'b': 160.0,
            't': 16.0,
            'r1': 17.0,
            'r2': 8.5,
        },
        'L160x160x17': {
            'h': 160.0,
            'b': 160.0,
            't': 17.0,
            'r1': 17.0,
            'r2': 8.5,
        },
        'L180x180x13': {
            'h': 180.0,
            'b': 180.0,
            't': 13.0,
            'r1': 18.0,
            'r2': 9.0,
        },
        'L180x180x14': {
            'h': 180.0,
            'b': 180.0,
            't': 14.0,
            'r1': 18.0,
            'r2': 9.0,
        },
        'L180x180x15': {
            'h': 180.0,
            'b': 180.0,
            't': 15.0,
            'r1': 18.0,
            'r2': 9.0,
        },
        'L180x180x16': {
            'h': 180.0,
            'b': 180.0,
            't': 16.0,
            'r1': 18.0,
            'r2': 9.0,
        },
        'L180x180x17': {
            'h': 180.0,
            'b': 180.0,
            't': 17.0,
            'r1': 18.0,
            'r2': 9.0,
        },
        'L180x180x18': {
            'h': 180.0,
            'b': 180.0,
            't': 18.0,
            'r1': 18.0,
            'r2': 9.0,
        },
        'L180x180x19': {
            'h': 180.0,
            'b': 180.0,
            't': 19.0,
            'r1': 18.0,
            'r2': 9.0,
        },
        'L180x180x20': {
            'h': 180.0,
            'b': 180.0,
            't': 20.0,
            'r1': 18.0,
            'r2': 9.0,
        },
        'L200x200x15': {
            'h': 200.0,
            'b': 200.0,
            't': 15.0,
            'r1': 18.0,
            'r2': 9.0,
        },
        'L200x200x16': {
            'h': 200.0,
            'b': 200.0,
            't': 16.0,
            'r1': 18.0,
            'r2': 9.0,
        },
        'L200x200x17': {
            'h': 200.0,
            'b': 200.0,
            't': 17.0,
            'r1': 18.0,
            'r2': 9.0,
        },
        'L200x200x18': {
            'h': 200.0,
            'b': 200.0,
            't': 18.0,
            'r1': 18.0,
            'r2': 9.0,
        },
        'L200x200x19': {
            'h': 200.0,
            'b': 200.0,
            't': 19.0,
            'r1': 18.0,
            'r2': 9.0,
        },
        'L200x200x20': {
            'h': 200.0,
            'b': 200.0,
            't': 20.0,
            'r1': 18.0,
            'r2': 9.0,
        },
        'L200x200x21': {
            'h': 200.0,
            'b': 200.0,
            't': 21.0,
            'r1': 18.0,
            'r2': 9.0,
        },
        'L200x200x22': {
            'h': 200.0,
            'b': 200.0,
            't': 22.0,
            'r1': 18.0,
            'r2': 9.0,
        },
        'L200x200x23': {
            'h': 200.0,
            'b': 200.0,
            't': 23.0,
            'r1': 18.0,
            'r2': 9.0,
        },
        'L200x200x24': {
            'h': 200.0,
            'b': 200.0,
            't': 24.0,
            'r1': 18.0,
            'r2': 9.0,
        },
        'L200x200x25': {
            'h': 200.0,
            'b': 200.0,
            't': 25.0,
            'r1': 18.0,
            'r2': 9.0,
        },
        'L200x200x26': {
            'h': 200.0,
            'b': 200.0,
            't': 26.0,
            'r1': 18.0,
            'r2': 9.0,
        },
        'L250x250x20': {
            'h': 250.0,
            'b': 250.0,
            't': 20.0,
            'r1': 18.0,
            'r2': 9.0,
        },
        'L250x250x21': {
            'h': 250.0,
            'b': 250.0,
            't': 21.0,
            'r1': 18.0,
            'r2': 9.0,
        },
        'L250x250x22': {
            'h': 250.0,
            'b': 250.0,
            't': 22.0,
            'r1': 18.0,
            'r2': 9.0,
        },
        'L250x250x23': {
            'h': 250.0,
            'b': 250.0,
            't': 23.0,
            'r1': 18.0,
            'r2': 9.0,
        },
        'L250x250x24': {
            'h': 250.0,
            'b': 250.0,
            't': 24.0,
            'r1': 18.0,
            'r2': 9.0,
        },
        'L250x250x25': {
            'h': 250.0,
            'b': 250.0,
            't': 25.0,
            'r1': 18.0,
            'r2': 9.0,
        },
        'L250x250x26': {
            'h': 250.0,
            'b': 250.0,
            't': 26.0,
            'r1': 18.0,
            'r2': 9.0,
        },
        'L250x250x27': {
            'h': 250.0,
            'b': 250.0,
            't': 27.0,
            'r1': 18.0,
            'r2': 9.0,
        },
        'L250x250x28': {
            'h': 250.0,
            'b': 250.0,
            't': 28.0,
            'r1': 18.0,
            'r2': 9.0,
        },
        'L250x250x35': {
            'h': 250.0,
            'b': 250.0,
            't': 35.0,
            'r1': 18.0,
            'r2': 9.0,
        },
        'L203x203x19': {
            'h': 203.0,
            'b': 203.0,
            't': 19.0,
            'r1': 8.0,
            'r2': 4.0,
        },
        'L203x203x22.2': {
            'h': 203.0,
            'b': 203.0,
            't': 22.2,
            'r1': 8.0,
            'r2': 4.0,
        },
        'L203x203x25.4': {
            'h': 203.0,
            'b': 203.0,
            't': 25.4,
            'r1': 8.0,
            'r2': 4.0,
        },
        'L203x203x28.6': {
            'h': 203.0,
            'b': 203.0,
            't': 28.6,
            'r1': 8.0,
            'r2': 4.0,
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
                f"Profile '{name}' not found in L sections. "
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
