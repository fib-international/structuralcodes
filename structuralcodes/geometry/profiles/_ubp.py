"""UBP profiles."""

from shapely import (
    Polygon,
)

from ._base_profile import BaseProfile
from ._common_functions import (
    _create_I_section,
)


class UBP(BaseProfile):
    """Simple class for representing a UBP profile.

    Universal Bearing Pile.
    """

    parameters = {
        'UBP203x203x45': {
            'h': 200.0,
            'b': 206.0,
            'tw': 9.5,
            'tf': 9.5,
            'r': 10.0,
        },
        'UBP203x203x54': {
            'h': 204.0,
            'b': 208.0,
            'tw': 11.3,
            'tf': 11.4,
            'r': 10.0,
        },
        'UBP254x254x63': {
            'h': 247.0,
            'b': 257.0,
            'tw': 10.6,
            'tf': 10.7,
            'r': 20.0,
        },
        'UBP254x254x71': {
            'h': 250.0,
            'b': 258.0,
            'tw': 12.0,
            'tf': 12.0,
            'r': 20.0,
        },
        'UBP254x254x85': {
            'h': 254.0,
            'b': 260.0,
            'tw': 14.4,
            'tf': 14.3,
            'r': 20.0,
        },
        'UBP305x305x79': {
            'h': 299.0,
            'b': 306.0,
            'tw': 11.0,
            'tf': 11.1,
            'r': 20.0,
        },
        'UBP305x305x88': {
            'h': 302.0,
            'b': 308.0,
            'tw': 12.4,
            'tf': 12.3,
            'r': 20.0,
        },
        'UBP305x305x95': {
            'h': 304.0,
            'b': 309.0,
            'tw': 13.3,
            'tf': 13.3,
            'r': 20.0,
        },
        'UBP305x305x110': {
            'h': 308.0,
            'b': 311.0,
            'tw': 15.3,
            'tf': 15.4,
            'r': 20.0,
        },
        'UBP305x305x126': {
            'h': 312.0,
            'b': 313.0,
            'tw': 17.5,
            'tf': 17.6,
            'r': 20.0,
        },
        'UBP305x305x149': {
            'h': 318.0,
            'b': 316.0,
            'tw': 20.6,
            'tf': 20.7,
            'r': 20.0,
        },
        'UBP305x305x186': {
            'h': 328.0,
            'b': 321.0,
            'tw': 25.5,
            'tf': 25.6,
            'r': 2.0,
        },
        'UBP305x305x223': {
            'h': 338.0,
            'b': 326.0,
            'tw': 30.3,
            'tf': 30.4,
            'r': 20.0,
        },
        'UBP356x368x109': {
            'h': 346.0,
            'b': 371.0,
            'tw': 12.8,
            'tf': 12.9,
            'r': 20.0,
        },
        'UBP356x368x133': {
            'h': 352.0,
            'b': 374.0,
            'tw': 15.6,
            'tf': 15.7,
            'r': 20.0,
        },
        'UBP356x368x152': {
            'h': 356.0,
            'b': 376.0,
            'tw': 17.8,
            'tf': 17.9,
            'r': 20.0,
        },
        'UBP356x368x174': {
            'h': 361.0,
            'b': 378.0,
            'tw': 20.3,
            'tf': 20.4,
            'r': 20.0,
        },
    }

    @classmethod
    def get_polygon(cls, name: str) -> Polygon:
        """Returns a shapely polygon representing a UBP section."""
        parameters = cls.parameters.get(name)
        if parameters is None:
            raise ValueError(
                f"Profile '{name}' not found in UBP sections. "
                "Select a valid profile (available ones: "
                f"{cls.profiles()})"
            )
        return _create_I_section(**parameters)

    @classmethod
    def profiles(cls) -> list:
        """Returns a list containing all available profiles."""
        return list(cls.parameters.keys())

    def __init__(self, name: str) -> None:
        """Creates a new UBP object."""
        parameters = self.parameters.get(name)
        if parameters is None:
            raise ValueError(
                f"Profile '{name}' not found in UBP sections. "
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
            Polygon: The represention of the UBP section.
        """
        return self._polygon

    @property
    def h(self) -> float:
        """Returns height of UBP section.

        Returns:
            float: Height h of UBP section.
        """
        return self._h

    @property
    def b(self) -> float:
        """Returns width of UBP section.

        Returns:
            float: Width b of UBP section.
        """
        return self._b

    @property
    def tw(self) -> float:
        """Returns thickness of web of UBP section.

        Returns:
            float: Web thickness tw of UBP section.
        """
        return self._tw

    @property
    def tf(self) -> float:
        """Returns thickness of flange of UBP section.

        Returns:
            float: Flange thickness tw of UBP section.
        """
        return self._tf

    @property
    def r(self) -> float:
        """Returns fillet radius of UBP section.

        Returns:
            float: Fillet radius r of UBP section.
        """
        return self._r
