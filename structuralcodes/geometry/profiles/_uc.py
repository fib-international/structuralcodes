"""UC profiles."""

from shapely import (
    Polygon,
)

from ._base_profile import BaseProfile
from ._common_functions import (
    _create_I_section,
)


class UC(BaseProfile):
    """Simple class for representing a UC profile.

    Universal Columns.
    """

    parameters = {
        'UC152x152x23': {
            'h': 152.0,
            'b': 152.0,
            'tw': 5.8,
            'tf': 6.8,
            'r': 8.0,
        },
        'UC152x152x30': {
            'h': 158.0,
            'b': 153.0,
            'tw': 6.5,
            'tf': 9.4,
            'r': 8.0,
        },
        'UC152x152x37': {
            'h': 162.0,
            'b': 154.0,
            'tw': 8.0,
            'tf': 11.5,
            'r': 8.0,
        },
        'UC152x152x44': {
            'h': 166.0,
            'b': 156.0,
            'tw': 9.5,
            'tf': 13.6,
            'r': 8.0,
        },
        'UC152x152x51': {
            'h': 170.0,
            'b': 157.0,
            'tw': 11.0,
            'tf': 15.7,
            'r': 8.0,
        },
        'UC203x203x46': {
            'h': 203.0,
            'b': 204.0,
            'tw': 7.2,
            'tf': 11.0,
            'r': 13.0,
        },
        'UC203x203x52': {
            'h': 206.0,
            'b': 204.0,
            'tw': 7.9,
            'tf': 12.5,
            'r': 13.0,
        },
        'UC203x203x60': {
            'h': 210.0,
            'b': 206.0,
            'tw': 9.4,
            'tf': 14.2,
            'r': 13.0,
        },
        'UC203x203x71': {
            'h': 216.0,
            'b': 206.0,
            'tw': 10.0,
            'tf': 17.3,
            'r': 13.0,
        },
        'UC203x203x86': {
            'h': 222.0,
            'b': 209.0,
            'tw': 12.7,
            'tf': 20.5,
            'r': 13.0,
        },
        'UC203x203x100': {
            'h': 229.0,
            'b': 210.0,
            'tw': 14.5,
            'tf': 23.7,
            'r': 13.0,
        },
        'UC203x203x113': {
            'h': 235.0,
            'b': 212.0,
            'tw': 16.3,
            'tf': 26.9,
            'r': 13.0,
        },
        'UC203x203x127': {
            'h': 241.0,
            'b': 214.0,
            'tw': 18.1,
            'tf': 30.1,
            'r': 13.0,
        },
        'UC254x254x73': {
            'h': 254.0,
            'b': 255.0,
            'tw': 8.6,
            'tf': 14.2,
            'r': 20.0,
        },
        'UC254x254x81': {
            'h': 256.0,
            'b': 255.0,
            'tw': 9.4,
            'tf': 15.6,
            'r': 20.0,
        },
        'UC254x254x89': {
            'h': 260.0,
            'b': 256.0,
            'tw': 10.3,
            'tf': 17.3,
            'r': 20.0,
        },
        'UC254x254x101': {
            'h': 264.0,
            'b': 257.0,
            'tw': 11.9,
            'tf': 19.6,
            'r': 20.0,
        },
        'UC254x254x107': {
            'h': 267.0,
            'b': 259.0,
            'tw': 12.8,
            'tf': 20.5,
            'r': 20.0,
        },
        'UC254x254x115': {
            'h': 269.0,
            'b': 259.0,
            'tw': 13.5,
            'tf': 22.1,
            'r': 20.0,
        },
        'UC254x254x132': {
            'h': 276.0,
            'b': 261.0,
            'tw': 15.3,
            'tf': 25.3,
            'r': 20.0,
        },
        'UC254x254x149': {
            'h': 282.0,
            'b': 263.0,
            'tw': 17.3,
            'tf': 28.4,
            'r': 20.0,
        },
        'UC254x254x167': {
            'h': 289.0,
            'b': 265.0,
            'tw': 19.2,
            'tf': 31.7,
            'r': 20.0,
        },
        'UC305x305x97': {
            'h': 308.0,
            'b': 305.0,
            'tw': 9.9,
            'tf': 15.4,
            'r': 20.0,
        },
        'UC305x305x107': {
            'h': 311.0,
            'b': 306.0,
            'tw': 10.9,
            'tf': 17.0,
            'r': 20.0,
        },
    }

    @classmethod
    def get_polygon(cls, name: str) -> Polygon:
        """Returns a shapely polygon representing a UC section."""
        parameters = cls.parameters.get(name)
        if parameters is None:
            raise ValueError(
                f"Profile '{name}' not found in UC sections. "
                "Select a valid profile (available ones: "
                f"{cls.profiles()})"
            )
        return _create_I_section(**parameters)

    @classmethod
    def profiles(cls) -> list:
        """Returns a list containing all available profiles."""
        return list(cls.parameters.keys())

    def __init__(self, name: str) -> None:
        """Creates a new UC object."""
        parameters = self.parameters.get(name)
        if parameters is None:
            raise ValueError(
                f"Profile '{name}' not found in UC sections. "
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
            Polygon: The represention of the UC section.
        """
        return self._polygon

    @property
    def h(self) -> float:
        """Returns height of UC section.

        Returns:
            float: Height h of UC section.
        """
        return self._h

    @property
    def b(self) -> float:
        """Returns width of UC section.

        Returns:
            float: Width b of UC section.
        """
        return self._b

    @property
    def tw(self) -> float:
        """Returns thickness of web of UC section.

        Returns:
            float: Web thickness tw of UC section.
        """
        return self._tw

    @property
    def tf(self) -> float:
        """Returns thickness of flange of UC section.

        Returns:
            float: Flange thickness tw of UC section.
        """
        return self._tf

    @property
    def r(self) -> float:
        """Returns fillet radius of UC section.

        Returns:
            float: Fillet radius r of UC section.
        """
        return self._r
