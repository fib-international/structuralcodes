"""UPN profiles."""

from shapely import (
    Polygon,
)

from ._base_profile import BaseProfile
from ._common_functions import (
    _create_taper_U_section,
)


class UPN(BaseProfile):
    """Simple class for representing an UPN profile.

    European standard channels UPN 50 - 400.

    Taper flange Channels.

    14% slope in flange.
    """

    parameters = {
        'UPN50': {
            'h': 50.0,
            'b': 38.0,
            'tw': 5.0,
            'tf': 7.0,
            'r1': 7.0,
            'r2': 4.0,
            'd': 21.0,
        },
        'UPN65': {
            'h': 65.0,
            'b': 42.0,
            'tw': 5.5,
            'tf': 7.5,
            'r1': 8.0,
            'r2': 4.0,
            'd': 34.0,
        },
        'UPN80': {
            'h': 80.0,
            'b': 45.0,
            'tw': 6.0,
            'tf': 8.0,
            'r1': 8.0,
            'r2': 4.0,
            'd': 47.0,
        },
        'UPN100': {
            'h': 100.0,
            'b': 50.0,
            'tw': 6.0,
            'tf': 8.5,
            'r1': 9.0,
            'r2': 5.0,
            'd': 64.0,
        },
        'UPN120': {
            'h': 120.0,
            'b': 55.0,
            'tw': 7.0,
            'tf': 9.0,
            'r1': 9.0,
            'r2': 5.0,
            'd': 82.0,
        },
        'UPN140': {
            'h': 140.0,
            'b': 60.0,
            'tw': 7.0,
            'tf': 10.0,
            'r1': 10.0,
            'r2': 5.0,
            'd': 98.0,
        },
        'UPN160': {
            'h': 160.0,
            'b': 65.0,
            'tw': 7.5,
            'tf': 10.5,
            'r1': 11.0,
            'r2': 6.0,
            'd': 115.0,
        },
        'UPN180': {
            'h': 180.0,
            'b': 70.0,
            'tw': 8.0,
            'tf': 11.0,
            'r1': 11.0,
            'r2': 6.0,
            'd': 133.0,
        },
        'UPN200': {
            'h': 200.0,
            'b': 75.0,
            'tw': 8.5,
            'tf': 11.5,
            'r1': 12.0,
            'r2': 6.0,
            'd': 151.0,
        },
        'UPN220': {
            'h': 220.0,
            'b': 80.0,
            'tw': 9.0,
            'tf': 12.5,
            'r1': 13.0,
            'r2': 7.0,
            'd': 167.0,
        },
        'UPN240': {
            'h': 240.0,
            'b': 85.0,
            'tw': 9.5,
            'tf': 13.0,
            'r1': 13.0,
            'r2': 7.0,
            'd': 184.0,
        },
        'UPN260': {
            'h': 260.0,
            'b': 90.0,
            'tw': 10.0,
            'tf': 14.0,
            'r1': 14.0,
            'r2': 7.0,
            'd': 200.0,
        },
        'UPN280': {
            'h': 280.0,
            'b': 95.0,
            'tw': 10.0,
            'tf': 15.0,
            'r1': 15.0,
            'r2': 8.0,
            'd': 216.0,
        },
        'UPN300': {
            'h': 300.0,
            'b': 100.0,
            'tw': 10.0,
            'tf': 16.0,
            'r1': 16.0,
            'r2': 8.0,
            'd': 232.0,
        },
        'UPN320': {
            'h': 320.0,
            'b': 100.0,
            'tw': 14.0,
            'tf': 17.5,
            'r1': 18.0,
            'r2': 9.0,
            'd': 246.0,
        },
        'UPN350': {
            'h': 350.0,
            'b': 100.0,
            'tw': 14.0,
            'tf': 16.0,
            'r1': 16.0,
            'r2': 8.0,
            'd': 282.0,
        },
        'UPN380': {
            'h': 380.0,
            'b': 102.0,
            'tw': 13.5,
            'tf': 16.0,
            'r1': 16.0,
            'r2': 8.0,
            'd': 313.0,
        },
        'UPN400': {
            'h': 400.0,
            'b': 110.0,
            'tw': 14.0,
            'tf': 18.0,
            'r1': 18.0,
            'r2': 9.0,
            'd': 324.0,
        },
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
        if parameters['h'] <= 300:
            parameters['slope'] = 0.08
            parameters['u'] = parameters['b'] / 2.0
        else:
            parameters['slope'] = 0.05
            parameters['u'] = (parameters['b'] - parameters['tw']) / 2.0
        return _create_taper_U_section(
            **{
                key: parameters[key]
                for key in parameters
                if key in ['h', 'b', 'tw', 'tf', 'r1', 'r2', 'slope', 'u']
            }
        )

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
        self._r1 = parameters.get('r1')
        self._r2 = parameters.get('r2')
        if self._h <= 300:
            self._flange_slope = 0.08
            self._u = self._b / 2.0
        else:
            self._flange_slope = 0.05
            self._u = (self._b - self._tw) / 2.0
        self._polygon = _create_taper_U_section(
            h=self._h,
            b=self._b,
            tw=self._tw,
            tf=self._tf,
            r1=self._r1,
            r2=self._r2,
            slope=self._flange_slope,
            u=self._u,
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
    def r1(self) -> float:
        """Returns fillet radius of UPN section.

        Returns:
            float: Fillet radius r1 of UPN section.
        """
        return self._r1

    @property
    def r2(self) -> float:
        """Returns fillet radius of UPN section.

        Returns:
            float: Fillet radius r2 of UPN section.
        """
        return self._r2
