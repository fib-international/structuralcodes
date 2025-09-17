"""IPN profiles."""

from shapely import (
    Polygon,
)

from ._base_profile import BaseProfile
from ._common_functions import (
    _create_taper_I_section,
)


class IPN(BaseProfile):
    """Simple class for representing an IPN profile.

    IPN in accordance with standard.

    14% slope in flange.
    """

    parameters = {
        'IPN80': {
            'h': 80.0,
            'b': 42.0,
            'tw': 3.9,
            'tf': 5.9,
            'r1': 3.9,
            'r2': 2.3,
            'd': 59.0,
        },
        'IPN100': {
            'h': 100.0,
            'b': 50.0,
            'tw': 4.5,
            'tf': 6.8,
            'r1': 4.5,
            'r2': 2.7,
            'd': 75.7,
        },
        'IPN120': {
            'h': 120.0,
            'b': 58.0,
            'tw': 5.1,
            'tf': 7.7,
            'r1': 5.1,
            'r2': 3.1,
            'd': 92.4,
        },
        'IPN140': {
            'h': 140.0,
            'b': 66.0,
            'tw': 5.7,
            'tf': 8.6,
            'r1': 5.7,
            'r2': 3.4,
            'd': 109.1,
        },
        'IPN160': {
            'h': 160.0,
            'b': 74.0,
            'tw': 6.3,
            'tf': 9.5,
            'r1': 6.3,
            'r2': 3.8,
            'd': 125.8,
        },
        'IPN180': {
            'h': 180.0,
            'b': 82.0,
            'tw': 6.9,
            'tf': 10.4,
            'r1': 6.9,
            'r2': 4.1,
            'd': 142.4,
        },
        'IPN200': {
            'h': 200.0,
            'b': 90.0,
            'tw': 7.5,
            'tf': 11.3,
            'r1': 7.5,
            'r2': 4.5,
            'd': 159.1,
        },
        'IPN220': {
            'h': 220.0,
            'b': 98.0,
            'tw': 8.1,
            'tf': 12.2,
            'r1': 8.1,
            'r2': 4.9,
            'd': 175.8,
        },
        'IPN240': {
            'h': 240.0,
            'b': 106.0,
            'tw': 8.7,
            'tf': 13.1,
            'r1': 8.7,
            'r2': 5.2,
            'd': 192.5,
        },
        'IPN260': {
            'h': 260.0,
            'b': 113.0,
            'tw': 9.4,
            'tf': 14.1,
            'r1': 9.4,
            'r2': 5.6,
            'd': 208.9,
        },
        'IPN280': {
            'h': 280.0,
            'b': 119.0,
            'tw': 10.1,
            'tf': 15.2,
            'r1': 10.1,
            'r2': 6.1,
            'd': 225.1,
        },
        'IPN300': {
            'h': 300.0,
            'b': 125.0,
            'tw': 10.8,
            'tf': 16.2,
            'r1': 10.8,
            'r2': 6.5,
            'd': 241.6,
        },
        'IPN320': {
            'h': 320.0,
            'b': 131.0,
            'tw': 11.5,
            'tf': 17.3,
            'r1': 11.5,
            'r2': 6.9,
            'd': 257.9,
        },
        'IPN340': {
            'h': 340.0,
            'b': 137.0,
            'tw': 12.2,
            'tf': 18.3,
            'r1': 12.2,
            'r2': 7.3,
            'd': 274.3,
        },
        'IPN360': {
            'h': 360.0,
            'b': 143.0,
            'tw': 13.0,
            'tf': 19.5,
            'r1': 13.0,
            'r2': 7.8,
            'd': 290.2,
        },
        'IPN380': {
            'h': 380.0,
            'b': 149.0,
            'tw': 13.7,
            'tf': 20.5,
            'r1': 13.7,
            'r2': 8.2,
            'd': 306.7,
        },
        'IPN400': {
            'h': 400.0,
            'b': 155.0,
            'tw': 14.4,
            'tf': 21.6,
            'r1': 14.4,
            'r2': 8.6,
            'd': 322.9,
        },
        'IPN450': {
            'h': 450.0,
            'b': 170.0,
            'tw': 16.2,
            'tf': 24.3,
            'r1': 16.2,
            'r2': 9.7,
            'd': 363.6,
        },
        'IPN500': {
            'h': 500.0,
            'b': 185.0,
            'tw': 18.0,
            'tf': 27.0,
            'r1': 18.0,
            'r2': 10.8,
            'd': 404.3,
        },
        'IPN550': {
            'h': 550.0,
            'b': 200.0,
            'tw': 19.0,
            'tf': 30.0,
            'r1': 19.0,
            'r2': 11.9,
            'd': 445.6,
        },
        'IPN600': {
            'h': 600.0,
            'b': 215.0,
            'tw': 21.6,
            'tf': 32.4,
            'r1': 21.6,
            'r2': 13.0,
            'd': 485.8,
        },
    }

    @classmethod
    def get_polygon(cls, name: str) -> Polygon:
        """Returns a shapely polygon representing an IPN section."""
        if isinstance(name, (float, int)):
            name = f'IPN{int(name):0d}'
        parameters = cls.parameters.get(name)
        if parameters is None:
            raise ValueError(
                f"Profile '{name}' not found in IPN sections. "
                "Select a valid profile (available ones: "
                f"{cls.profiles()})"
            )
        parameters['slope'] = 0.14
        return _create_taper_I_section(
            **{
                key: parameters[key]
                for key in parameters
                if key in ['h', 'b', 'tw', 'tf', 'r1', 'r2', 'slope']
            }
        )

    @classmethod
    def profiles(cls) -> list:
        """Returns a list containing all available profiles."""
        return list(cls.parameters.keys())

    def __init__(self, name: str) -> None:
        """Creates a new IPN object."""
        if isinstance(name, (float, int)):
            name = f'IPN{int(name):0d}'
        parameters = self.parameters.get(name)
        if parameters is None:
            raise ValueError(
                f"Profile '{name}' not found in IPN sections. "
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
        self._flange_slope = 0.14
        self._polygon = _create_taper_I_section(
            h=self._h,
            b=self._b,
            tw=self._tw,
            tf=self._tf,
            r1=self._r1,
            r2=self._r2,
            slope=self._flange_slope,
        )

    @property
    def polygon(self) -> Polygon:
        """Returns shapely Polygon of section.

        Returns:
            Polygon: The represention of the IPN section.
        """
        return self._polygon

    @property
    def h(self) -> float:
        """Returns height of IPN section.

        Returns:
            float: Height h of IPN section.
        """
        return self._h

    @property
    def b(self) -> float:
        """Returns width of IPN section.

        Returns:
            float: Width b of IPN section.
        """
        return self._b

    @property
    def tw(self) -> float:
        """Returns thickness of web of IPN section.

        Returns:
            float: Web thickness tw of IPN section.
        """
        return self._tw

    @property
    def tf(self) -> float:
        """Returns thickness of flange of IPN section.

        Returns:
            float: Flange thickness tw of IPN section.
        """
        return self._tf

    @property
    def r1(self) -> float:
        """Returns fillet radius of IPN section.

        Returns:
            float: Fillet radius r1 of IPN section.
        """
        return self._r1

    @property
    def r2(self) -> float:
        """Returns fillet radius of IPN section.

        Returns:
            float: Fillet radius r2 of IPN section.
        """
        return self._r2
