"""IPE profiles."""

from shapely import (
    Polygon,
)

from ._base_profile import BaseProfile
from ._common_functions import (
    _create_I_section,
)


class IPE(BaseProfile):
    """Simple class for representing an IPE profile.

    IPE 80-600 in accordance with standard Euronorm 19-57.
    """

    parameters = {
        'IPE80': {'h': 80.0, 'b': 46.0, 'tw': 3.8, 'tf': 5.2, 'r': 5.0},
        'IPE100': {'h': 100.0, 'b': 55.0, 'tw': 4.1, 'tf': 5.7, 'r': 7.0},
        'IPE120': {'h': 120.0, 'b': 64.0, 'tw': 4.4, 'tf': 6.3, 'r': 7.0},
        'IPE140': {'h': 140.0, 'b': 73.0, 'tw': 4.7, 'tf': 6.9, 'r': 7.0},
        'IPE160': {'h': 160.0, 'b': 82.0, 'tw': 5.0, 'tf': 7.4, 'r': 9.0},
        'IPE180': {'h': 180.0, 'b': 91.0, 'tw': 5.3, 'tf': 8.0, 'r': 9.0},
        'IPE200': {'h': 200.0, 'b': 100.0, 'tw': 5.6, 'tf': 8.5, 'r': 12.0},
        'IPE220': {'h': 220.0, 'b': 110.0, 'tw': 5.9, 'tf': 9.2, 'r': 12.0},
        'IPE240': {'h': 240.0, 'b': 120.0, 'tw': 6.2, 'tf': 9.8, 'r': 15.0},
        'IPE270': {'h': 270.0, 'b': 135.0, 'tw': 6.6, 'tf': 10.2, 'r': 15.0},
        'IPE300': {'h': 300.0, 'b': 150.0, 'tw': 7.1, 'tf': 10.7, 'r': 15.0},
        'IPE330': {'h': 330.0, 'b': 160.0, 'tw': 7.5, 'tf': 11.5, 'r': 18.0},
        'IPE360': {'h': 360.0, 'b': 170.0, 'tw': 8.0, 'tf': 12.7, 'r': 18.0},
        'IPE400': {'h': 400.0, 'b': 180.0, 'tw': 8.6, 'tf': 13.5, 'r': 21.0},
        'IPE450': {'h': 450.0, 'b': 190.0, 'tw': 9.4, 'tf': 14.6, 'r': 21.0},
        'IPE500': {'h': 500.0, 'b': 200.0, 'tw': 10.2, 'tf': 16.0, 'r': 21.0},
        'IPE550': {'h': 550.0, 'b': 210.0, 'tw': 11.1, 'tf': 17.2, 'r': 24.0},
        'IPE600': {'h': 600.0, 'b': 220.0, 'tw': 12.0, 'tf': 19.0, 'r': 24.0},
    }

    @classmethod
    def get_polygon(cls, name: str) -> Polygon:
        """Returns a shapely polygon representing an IPE section."""
        if isinstance(name, (float, int)):
            name = f'IPE{int(name):0d}'
        parameters = cls.parameters.get(name)
        if parameters is None:
            raise ValueError(
                f"Profile '{name}' not found in IPE sections. "
                "Select a valid profile (available ones: "
                f"{cls.profiles()})"
            )
        return _create_I_section(**parameters)

    @classmethod
    def profiles(cls) -> list:
        """Returns a list containing all available profiles."""
        return list(cls.parameters.keys())

    def __init__(self, name: str) -> None:
        """Creates a new IPE object."""
        if isinstance(name, (float, int)):
            name = f'IPE{int(name):0d}'
        parameters = self.parameters.get(name)
        if parameters is None:
            raise ValueError(
                f"Profile '{name}' not found in IPE sections. "
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
            Polygon: The represention of the IPE section.
        """
        return self._polygon

    @property
    def h(self) -> float:
        """Returns height of IPE section.

        Returns:
            float: Height h of IPE section.
        """
        return self._h

    @property
    def b(self) -> float:
        """Returns width of IPE section.

        Returns:
            float: Width b of IPE section.
        """
        return self._b

    @property
    def tw(self) -> float:
        """Returns thickness of web of IPE section.

        Returns:
            float: Web thickness tw of IPE section.
        """
        return self._tw

    @property
    def tf(self) -> float:
        """Returns thickness of flange of IPE section.

        Returns:
            float: Flange thickness tw of IPE section.
        """
        return self._tf

    @property
    def r(self) -> float:
        """Returns fillet radius of IPE section.

        Returns:
            float: Fillet radius r of IPE section.
        """
        return self._r
