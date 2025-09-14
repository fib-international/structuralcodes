"""HD profiles."""

from shapely import (
    Polygon,
)

from ._base_profile import BaseProfile
from ._common_functions import (
    _create_I_section,
)


class HD(BaseProfile):
    """Simple class for representing an HD profile.

    HD A, HD B, HD M 100-1000 in accordance with Standard Euronorm 53-62.
    """

    parameters = {
        'HD260x54.1': {
            'h': 244.0,
            'b': 260.0,
            'tw': 6.5,
            'tf': 9.5,
            'r': 24.0,
        },
        'HD260x68.2': {
            'h': 250.0,
            'b': 260.0,
            'tw': 7.5,
            'tf': 12.5,
            'r': 24.0,
        },
        'HD260x93.0': {
            'h': 260.0,
            'b': 260.0,
            'tw': 10.0,
            'tf': 17.5,
            'r': 24.0,
        },
        'HD260x114.0': {
            'h': 268.0,
            'b': 262.0,
            'tw': 12.5,
            'tf': 21.5,
            'r': 24.0,
        },
        'HD260x142.0': {
            'h': 278.0,
            'b': 265.0,
            'tw': 15.5,
            'tf': 26.5,
            'r': 24.0,
        },
        'HD260x172.0': {
            'h': 290.0,
            'b': 268.0,
            'tw': 18.0,
            'tf': 32.5,
            'r': 24.0,
        },
        'HD320x74.2': {
            'h': 301.0,
            'b': 300.0,
            'tw': 8.0,
            'tf': 11.0,
            'r': 27.0,
        },
        'HD320x97.6': {
            'h': 310.0,
            'b': 300.0,
            'tw': 9.0,
            'tf': 15.5,
            'r': 27.0,
        },
        'HD320x127.0': {
            'h': 320.0,
            'b': 300.0,
            'tw': 11.5,
            'tf': 20.5,
            'r': 27.0,
        },
        'HD320x158.0': {
            'h': 330.0,
            'b': 303.0,
            'tw': 14.5,
            'tf': 25.5,
            'r': 27.0,
        },
        'HD320x198.0': {
            'h': 343.0,
            'b': 306.0,
            'tw': 18.0,
            'tf': 32.0,
            'r': 27.0,
        },
        'HD320x245.0': {
            'h': 359.0,
            'b': 309.0,
            'tw': 21.0,
            'tf': 40.0,
            'r': 27.0,
        },
        'HD320x300.0': {
            'h': 375.0,
            'b': 313.0,
            'tw': 27.0,
            'tf': 48.0,
            'r': 27.0,
        },
        'HD360x134.0': {
            'h': 356.0,
            'b': 369.0,
            'tw': 11.2,
            'tf': 18.0,
            'r': 15.0,
        },
        'HD360x147.0': {
            'h': 360.0,
            'b': 370.0,
            'tw': 12.3,
            'tf': 19.8,
            'r': 15.0,
        },
        'HD360x162.0': {
            'h': 364.0,
            'b': 371.0,
            'tw': 13.3,
            'tf': 21.8,
            'r': 15.0,
        },
        'HD360x179.0': {
            'h': 368.0,
            'b': 373.0,
            'tw': 15.0,
            'tf': 23.9,
            'r': 15.0,
        },
        'HD360x196.0': {
            'h': 372.0,
            'b': 374.0,
            'tw': 16.4,
            'tf': 26.2,
            'r': 15.0,
        },
        'HD400x187.0': {
            'h': 368.0,
            'b': 391.0,
            'tw': 15.0,
            'tf': 24.0,
            'r': 15.0,
        },
        'HD400x216.0': {
            'h': 375.0,
            'b': 394.0,
            'tw': 17.3,
            'tf': 27.7,
            'r': 15.0,
        },
        'HD400x237.0': {
            'h': 380.0,
            'b': 395.0,
            'tw': 18.9,
            'tf': 30.2,
            'r': 15.0,
        },
        'HD400x262.0': {
            'h': 387.0,
            'b': 398.0,
            'tw': 21.1,
            'tf': 33.3,
            'r': 15.0,
        },
        'HD400x287.0': {
            'h': 393.0,
            'b': 399.0,
            'tw': 22.6,
            'tf': 36.6,
            'r': 15.0,
        },
        'HD400x314.0': {
            'h': 399.0,
            'b': 401.0,
            'tw': 24.9,
            'tf': 39.6,
            'r': 15.0,
        },
        'HD400x347.0': {
            'h': 407.0,
            'b': 404.0,
            'tw': 27.2,
            'tf': 43.7,
            'r': 15.0,
        },
        'HD400x382.0': {
            'h': 416.0,
            'b': 406.0,
            'tw': 29.8,
            'tf': 48.0,
            'r': 15.0,
        },
        'HD400x421.0': {
            'h': 425.0,
            'b': 409.0,
            'tw': 32.8,
            'tf': 52.6,
            'r': 15.0,
        },
        'HD400x463.0': {
            'h': 435.0,
            'b': 412.0,
            'tw': 35.8,
            'tf': 57.4,
            'r': 15.0,
        },
        'HD400x509.0': {
            'h': 446.0,
            'b': 416.0,
            'tw': 39.1,
            'tf': 62.7,
            'r': 15.0,
        },
        'HD400x551.0': {
            'h': 455.0,
            'b': 418.0,
            'tw': 42.0,
            'tf': 67.6,
            'r': 15.0,
        },
        'HD400x592.0': {
            'h': 465.0,
            'b': 421.0,
            'tw': 45.0,
            'tf': 72.3,
            'r': 15.0,
        },
        'HD400x634.0': {
            'h': 474.0,
            'b': 424.0,
            'tw': 47.6,
            'tf': 77.1,
            'r': 15.0,
        },
        'HD400x677.0': {
            'h': 483.0,
            'b': 428.0,
            'tw': 51.2,
            'tf': 81.5,
            'r': 15.0,
        },
        'HD400x744.0': {
            'h': 498.0,
            'b': 432.0,
            'tw': 55.6,
            'tf': 88.9,
            'r': 15.0,
        },
        'HD400x818.0': {
            'h': 514.0,
            'b': 437.0,
            'tw': 60.5,
            'tf': 97.0,
            'r': 15.0,
        },
        'HD400x900.0': {
            'h': 531.0,
            'b': 442.0,
            'tw': 65.9,
            'tf': 106.0,
            'r': 15.0,
        },
        'HD400x990.0': {
            'h': 550.0,
            'b': 448.0,
            'tw': 71.9,
            'tf': 115.0,
            'r': 15.0,
        },
        'HD400x1086.0': {
            'h': 569.0,
            'b': 454.0,
            'tw': 78.0,
            'tf': 125.0,
            'r': 15.0,
        },
    }

    @classmethod
    def get_polygon(cls, name: str) -> Polygon:
        """Returns a shapely polygon representing an HD section."""
        parameters = cls.parameters.get(name)
        if parameters is None:
            raise ValueError(
                f"Profile '{name}' not found in HD sections. "
                "Select a valid profile (available ones: "
                f"{cls.profiles()})"
            )
        return _create_I_section(**parameters)

    @classmethod
    def profiles(cls) -> list:
        """Returns a list containing all available profiles."""
        return list(cls.parameters.keys())

    def __init__(self, name: str) -> None:
        """Creates a new HD object."""
        parameters = self.parameters.get(name)
        if parameters is None:
            raise ValueError(
                f"Profile '{name}' not found in HD sections. "
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
            Polygon: THD represention of tHD HD section.
        """
        return self._polygon

    @property
    def h(self) -> float:
        """Returns Height of HD section.

        Returns:
            float: HHeight h of HD section.
        """
        return self._h

    @property
    def b(self) -> float:
        """Returns width of HD section.

        Returns:
            float: Width b of HD section.
        """
        return self._b

    @property
    def tw(self) -> float:
        """Returns thickness of web of HD section.

        Returns:
            float: Web thickness tw of HD section.
        """
        return self._tw

    @property
    def tf(self) -> float:
        """Returns thickness of flange of HD section.

        Returns:
            float: Flange thickness tw of HD section.
        """
        return self._tf

    @property
    def r(self) -> float:
        """Returns fillet radius of HD section.

        Returns:
            float: Fillet radius r of HD section.
        """
        return self._r
