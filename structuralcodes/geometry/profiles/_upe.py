"""UPN profiles."""

from shapely import (
    Polygon,
)

from ._base_profile import BaseProfile
from ._common_functions import _create_parallel_U_section


class UPE(BaseProfile):
    """Simple class for representing an UPE profile.

    European Parallel flange channels UPE 80 - 400.

    Dimensions: EN 10 365:2017
    """

    parameters = {
        'UPE80': {'h': 80.0, 'b': 50.0, 'tw': 4.0, 'tf': 7.0, 'r': 10.0},
        'UPE100': {'h': 100.0, 'b': 55.0, 'tw': 4.5, 'tf': 7.5, 'r': 10.0},
        'UPE120': {'h': 120.0, 'b': 60.0, 'tw': 5.0, 'tf': 8.0, 'r': 12.0},
        'UPE140': {'h': 140.0, 'b': 65.0, 'tw': 5.0, 'tf': 9.0, 'r': 12.0},
        'UPE160': {'h': 160.0, 'b': 70.0, 'tw': 5.5, 'tf': 9.5, 'r': 12.0},
        'UPE180': {'h': 180.0, 'b': 75.0, 'tw': 5.5, 'tf': 10.5, 'r': 12.0},
        'UPE200': {'h': 200.0, 'b': 80.0, 'tw': 6.0, 'tf': 11.0, 'r': 13.0},
        'UPE220': {'h': 220.0, 'b': 85.0, 'tw': 6.5, 'tf': 12.0, 'r': 13.0},
        'UPE240': {'h': 240.0, 'b': 90.0, 'tw': 7.0, 'tf': 12.5, 'r': 15.0},
        'UPE270': {'h': 270.0, 'b': 95.0, 'tw': 7.5, 'tf': 13.5, 'r': 15.0},
        'UPE300': {'h': 300.0, 'b': 100.0, 'tw': 9.5, 'tf': 15.0, 'r': 15.0},
        'UPE330': {'h': 330.0, 'b': 105.0, 'tw': 11.0, 'tf': 16.0, 'r': 18.0},
        'UPE360': {'h': 360.0, 'b': 110.0, 'tw': 12.0, 'tf': 17.0, 'r': 18.0},
        'UPE400': {'h': 400.0, 'b': 115.0, 'tw': 13.5, 'tf': 18.0, 'r': 18.0},
    }

    @classmethod
    def get_polygon(cls, name: str) -> Polygon:
        """Returns a shapely polygon representing an UPN section."""
        if isinstance(name, (float, int)):
            name = f'UPN{int(name):0d}'
        parameters = cls.parameters.get(name)
        if parameters is None:
            raise ValueError(
                f"Profile '{name}' not found in UPN sections. "
                "Select a valid profile (available ones: "
                f"{cls.profiles()})"
            )
        return _create_parallel_U_section(**parameters)

    @classmethod
    def profiles(cls) -> list:
        """Returns a list containing all available profiles."""
        return list(cls.parameters.keys())

    def __init__(self, name: str) -> None:
        """Creates a new UPN object."""
        if isinstance(name, (float, int)):
            name = f'UPN{int(name):0d}'
        parameters = self.parameters.get(name)
        if parameters is None:
            raise ValueError(
                f"Profile '{name}' not found in UPN sections. "
                "Select a valid profile (available ones: "
                f"{self.profiles()})"
            )
        super().__init__()
        self._h = parameters.get('h')
        self._b = parameters.get('b')
        self._tw = parameters.get('tw')
        self._tf = parameters.get('tf')
        self._r = parameters.get('r')

        self._polygon = _create_parallel_U_section(
            h=self._h,
            b=self._b,
            tw=self._tw,
            tf=self._tf,
            r=self._r,
        )

    @property
    def polygon(self) -> Polygon:
        """Returns shapely Polygon of section.

        Returns:
            Polygon: The represention of the UPN section.
        """
        return self._polygon

    @property
    def h(self) -> float:
        """Returns height of UPN section.

        Returns:
            float: Height h of UPN section.
        """
        return self._h

    @property
    def b(self) -> float:
        """Returns width of UPN section.

        Returns:
            float: Width b of UPN section.
        """
        return self._b

    @property
    def tw(self) -> float:
        """Returns thickness of web of UPN section.

        Returns:
            float: Web thickness tw of UPN section.
        """
        return self._tw

    @property
    def tf(self) -> float:
        """Returns thickness of flange of UPN section.

        Returns:
            float: Flange thickness tw of UPN section.
        """
        return self._tf

    @property
    def r(self) -> float:
        """Returns fillet radius of UPN section.

        Returns:
            float: Fillet radius r1 of UPN section.
        """
        return self._r
