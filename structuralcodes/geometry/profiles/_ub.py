"""UB profiles."""

from shapely import (
    Polygon,
)

from ._base_profile import BaseProfile
from ._common_functions import (
    _create_I_section,
)


class UB(BaseProfile):
    """Simple class for representing a UB profile.

    Universal Beams.
    """

    parameters = {
        'UB127x76x13': {'h': 127.0, 'b': 76.0, 'tw': 4.0, 'tf': 7.6, 'r': 8.0},
        'UB152x89x16': {'h': 152.0, 'b': 89.0, 'tw': 4.5, 'tf': 7.7, 'r': 8.0},
        'UB178x102x19': {
            'h': 178.0,
            'b': 101.0,
            'tw': 4.8,
            'tf': 7.9,
            'r': 8.0,
        },
        'UB203x102x23': {
            'h': 203.0,
            'b': 102.0,
            'tw': 5.4,
            'tf': 9.3,
            'r': 8.0,
        },
        'UB203x133x25': {
            'h': 203.0,
            'b': 133.0,
            'tw': 5.7,
            'tf': 7.8,
            'r': 8.0,
        },
        'UB203x133x30': {
            'h': 207.0,
            'b': 134.0,
            'tw': 6.4,
            'tf': 9.6,
            'r': 8.0,
        },
        'UB254x102x22': {
            'h': 254.0,
            'b': 102.0,
            'tw': 5.7,
            'tf': 6.8,
            'r': 8.0,
        },
        'UB254x102x25': {
            'h': 257.0,
            'b': 102.0,
            'tw': 6.0,
            'tf': 8.4,
            'r': 8.0,
        },
        'UB254x102x28': {
            'h': 260.0,
            'b': 102.0,
            'tw': 6.3,
            'tf': 10.0,
            'r': 8.0,
        },
        'UB254x146x31': {
            'h': 251.0,
            'b': 146.0,
            'tw': 6.0,
            'tf': 8.6,
            'r': 8.0,
        },
        'UB254x146x37': {
            'h': 256.0,
            'b': 146.0,
            'tw': 6.3,
            'tf': 10.9,
            'r': 8.0,
        },
        'UB254x146x43': {
            'h': 260.0,
            'b': 147.0,
            'tw': 7.2,
            'tf': 12.7,
            'r': 8.0,
        },
        'UB305x102x25': {
            'h': 305.0,
            'b': 102.0,
            'tw': 5.8,
            'tf': 7.0,
            'r': 8.0,
        },
        'UB305x102x28': {
            'h': 309.0,
            'b': 102.0,
            'tw': 6.0,
            'tf': 8.8,
            'r': 8.0,
        },
        'UB305x102x33': {
            'h': 313.0,
            'b': 102.0,
            'tw': 6.6,
            'tf': 10.8,
            'r': 8.0,
        },
        'UB305x127x37': {
            'h': 304.0,
            'b': 123.0,
            'tw': 7.1,
            'tf': 10.7,
            'r': 9.0,
        },
        'UB305x127x42': {
            'h': 307.0,
            'b': 124.0,
            'tw': 8.0,
            'tf': 12.1,
            'r': 9.0,
        },
        'UB305x127x48': {
            'h': 311.0,
            'b': 125.0,
            'tw': 9.0,
            'tf': 14.0,
            'r': 9.0,
        },
        'UB305x165x40': {
            'h': 303.0,
            'b': 165.0,
            'tw': 6.0,
            'tf': 10.2,
            'r': 9.0,
        },
        'UB305x165x46': {
            'h': 307.0,
            'b': 166.0,
            'tw': 6.7,
            'tf': 11.8,
            'r': 9.0,
        },
        'UB305x165x54': {
            'h': 310.0,
            'b': 167.0,
            'tw': 7.9,
            'tf': 13.7,
            'r': 9.0,
        },
        'UB356x171x45': {
            'h': 351.0,
            'b': 171.0,
            'tw': 7.0,
            'tf': 9.7,
            'r': 13.0,
        },
        'UB356x171x51': {
            'h': 355.0,
            'b': 172.0,
            'tw': 7.4,
            'tf': 11.5,
            'r': 13.0,
        },
        'UB356x171x57': {
            'h': 358.0,
            'b': 172.0,
            'tw': 8.1,
            'tf': 13.0,
            'r': 13.0,
        },
    }

    @classmethod
    def get_polygon(cls, name: str) -> Polygon:
        """Returns a shapely polygon representing a UB section."""
        parameters = cls.parameters.get(name)
        if parameters is None:
            raise ValueError(
                f"Profile '{name}' not found in UB sections. "
                "Select a valid profile (available ones: "
                f"{cls.profiles()})"
            )
        return _create_I_section(**parameters)

    @classmethod
    def profiles(cls) -> list:
        """Returns a list containing all available profiles."""
        return list(cls.parameters.keys())

    def __init__(self, name: str) -> None:
        """Creates a new UB object."""
        parameters = self.parameters.get(name)
        if parameters is None:
            raise ValueError(
                f"Profile '{name}' not found in UB sections. "
                "Select a valid profile (available ones: "
                f"{self.profiles()})"
            )
        super().__init__()
        self._h = parameters.get('h')
        self._b = parameters.get('b')
        self._tw = parameters.get('tw')
        self._tf = parameters.get('tf')
        self._r = parameters.get('r')
        self._polygon = _create_I_section(**parameters)

    @property
    def polygon(self) -> Polygon:
        """Returns shapely Polygon of section.

        Returns:
            Polygon: The represention of the UB section.
        """
        return self._polygon

    @property
    def h(self) -> float:
        """Returns height of UB section.

        Returns:
            float: Height h of UB section.
        """
        return self._h

    @property
    def b(self) -> float:
        """Returns width of UB section.

        Returns:
            float: Width b of UB section.
        """
        return self._b

    @property
    def tw(self) -> float:
        """Returns thickness of web of UB section.

        Returns:
            float: Web thickness tw of UB section.
        """
        return self._tw

    @property
    def tf(self) -> float:
        """Returns thickness of flange of UB section.

        Returns:
            float: Flange thickness tw of UB section.
        """
        return self._tf

    @property
    def r(self) -> float:
        """Returns fillet radius of UB section.

        Returns:
            float: Fillet radius r of UB section.
        """
        return self._r
