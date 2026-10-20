"""HE profiles."""

from shapely import (
    Polygon,
)

from ._base_profile import BaseProfile
from ._common_functions import (
    _create_I_section,
)


class HE(BaseProfile):
    """Simple class for representing an HE profile.

    HE A, HE B, HE M 100-1000 in accordance with Standard Euronorm 53-62.
    """

    parameters = {
        'HEA100': {'h': 96.0, 'b': 100.0, 'tw': 5.0, 'tf': 8.0, 'r': 12.0},
        'HEA120': {'h': 114.0, 'b': 120.0, 'tw': 5.0, 'tf': 8.0, 'r': 12.0},
        'HEA140': {'h': 133.0, 'b': 140.0, 'tw': 5.5, 'tf': 8.5, 'r': 12.0},
        'HEA160': {'h': 152.0, 'b': 160.0, 'tw': 6.0, 'tf': 9.0, 'r': 15.0},
        'HEA180': {'h': 171.0, 'b': 180.0, 'tw': 6.0, 'tf': 9.5, 'r': 15.0},
        'HEA200': {'h': 190.0, 'b': 200.0, 'tw': 6.5, 'tf': 10.0, 'r': 18.0},
        'HEA220': {'h': 210.0, 'b': 220.0, 'tw': 7.0, 'tf': 11.0, 'r': 18.0},
        'HEA240': {'h': 230.0, 'b': 240.0, 'tw': 7.5, 'tf': 12.0, 'r': 21.0},
        'HEA260': {'h': 250.0, 'b': 260.0, 'tw': 7.5, 'tf': 12.5, 'r': 24.0},
        'HEA280': {'h': 270.0, 'b': 280.0, 'tw': 8.0, 'tf': 13.0, 'r': 24.0},
        'HEA300': {'h': 290.0, 'b': 300.0, 'tw': 8.5, 'tf': 14.0, 'r': 27.0},
        'HEA320': {'h': 310.0, 'b': 300.0, 'tw': 9.0, 'tf': 15.5, 'r': 27.0},
        'HEA340': {'h': 330.0, 'b': 300.0, 'tw': 9.5, 'tf': 16.5, 'r': 27.0},
        'HEA360': {'h': 350.0, 'b': 300.0, 'tw': 10.0, 'tf': 17.5, 'r': 27.0},
        'HEA400': {'h': 390.0, 'b': 300.0, 'tw': 11.0, 'tf': 19.0, 'r': 27.0},
        'HEA450': {'h': 440.0, 'b': 300.0, 'tw': 11.5, 'tf': 21.0, 'r': 27.0},
        'HEA500': {'h': 490.0, 'b': 300.0, 'tw': 12.0, 'tf': 23.0, 'r': 27.0},
        'HEA550': {'h': 540.0, 'b': 300.0, 'tw': 12.5, 'tf': 24.0, 'r': 27.0},
        'HEA600': {'h': 590.0, 'b': 300.0, 'tw': 13.0, 'tf': 25.0, 'r': 27.0},
        'HEA650': {'h': 640.0, 'b': 300.0, 'tw': 13.5, 'tf': 26.0, 'r': 27.0},
        'HEA700': {'h': 690.0, 'b': 300.0, 'tw': 14.5, 'tf': 27.0, 'r': 27.0},
        'HEA800': {'h': 790.0, 'b': 300.0, 'tw': 15.0, 'tf': 28.0, 'r': 30.0},
        'HEA900': {'h': 890.0, 'b': 300.0, 'tw': 16.0, 'tf': 30.0, 'r': 30.0},
        'HEA1000': {'h': 990.0, 'b': 300.0, 'tw': 16.5, 'tf': 31.0, 'r': 30.0},
        'HEB100': {'h': 100.0, 'b': 100.0, 'tw': 6.0, 'tf': 10.0, 'r': 12.0},
        'HEB120': {'h': 120.0, 'b': 120.0, 'tw': 6.5, 'tf': 11.0, 'r': 12.0},
        'HEB140': {'h': 140.0, 'b': 140.0, 'tw': 7.0, 'tf': 12.0, 'r': 12.0},
        'HEB160': {'h': 160.0, 'b': 160.0, 'tw': 8.0, 'tf': 13.0, 'r': 15.0},
        'HEB180': {'h': 180.0, 'b': 180.0, 'tw': 8.5, 'tf': 14.0, 'r': 15.0},
        'HEB200': {'h': 200.0, 'b': 200.0, 'tw': 9.0, 'tf': 15.0, 'r': 18.0},
        'HEB220': {'h': 220.0, 'b': 220.0, 'tw': 9.5, 'tf': 16.0, 'r': 18.0},
        'HEB240': {'h': 240.0, 'b': 240.0, 'tw': 10.0, 'tf': 17.0, 'r': 21.0},
        'HEB260': {'h': 260.0, 'b': 260.0, 'tw': 10.0, 'tf': 17.5, 'r': 24.0},
        'HEB280': {'h': 280.0, 'b': 280.0, 'tw': 10.5, 'tf': 18.0, 'r': 24.0},
        'HEB300': {'h': 300.0, 'b': 300.0, 'tw': 11.0, 'tf': 19.0, 'r': 27.0},
        'HEB320': {'h': 320.0, 'b': 300.0, 'tw': 11.5, 'tf': 20.5, 'r': 27.0},
        'HEB340': {'h': 340.0, 'b': 300.0, 'tw': 12.0, 'tf': 21.5, 'r': 27.0},
        'HEB360': {'h': 360.0, 'b': 300.0, 'tw': 12.5, 'tf': 22.5, 'r': 27.0},
        'HEB400': {'h': 400.0, 'b': 300.0, 'tw': 13.5, 'tf': 24.0, 'r': 27.0},
        'HEB450': {'h': 450.0, 'b': 300.0, 'tw': 14.0, 'tf': 26.0, 'r': 27.0},
        'HEB500': {'h': 500.0, 'b': 300.0, 'tw': 14.5, 'tf': 28.0, 'r': 27.0},
        'HEB550': {'h': 550.0, 'b': 300.0, 'tw': 15.0, 'tf': 29.0, 'r': 27.0},
        'HEB600': {'h': 600.0, 'b': 300.0, 'tw': 15.5, 'tf': 30.0, 'r': 27.0},
        'HEB650': {'h': 650.0, 'b': 300.0, 'tw': 16.0, 'tf': 31.0, 'r': 27.0},
        'HEB700': {'h': 700.0, 'b': 300.0, 'tw': 17.0, 'tf': 32.0, 'r': 27.0},
        'HEB800': {'h': 800.0, 'b': 300.0, 'tw': 17.5, 'tf': 33.0, 'r': 30.0},
        'HEB900': {'h': 900.0, 'b': 300.0, 'tw': 18.5, 'tf': 35.0, 'r': 30.0},
        'HEB1000': {
            'h': 1000.0,
            'b': 300.0,
            'tw': 19.0,
            'tf': 36.0,
            'r': 30.0,
        },
        'HEM100': {'h': 120.0, 'b': 106.0, 'tw': 12.0, 'tf': 20.0, 'r': 12.0},
        'HEM120': {'h': 140.0, 'b': 126.0, 'tw': 12.5, 'tf': 21.0, 'r': 12.0},
        'HEM140': {'h': 160.0, 'b': 146.0, 'tw': 13.0, 'tf': 22.0, 'r': 12.0},
        'HEM160': {'h': 180.0, 'b': 166.0, 'tw': 14.0, 'tf': 23.0, 'r': 15.0},
        'HEM180': {'h': 200.0, 'b': 186.0, 'tw': 14.5, 'tf': 24.0, 'r': 15.0},
        'HEM200': {'h': 220.0, 'b': 206.0, 'tw': 15.0, 'tf': 25.0, 'r': 18.0},
        'HEM220': {'h': 240.0, 'b': 226.0, 'tw': 15.5, 'tf': 26.0, 'r': 18.0},
        'HEM240': {'h': 270.0, 'b': 248.0, 'tw': 18.0, 'tf': 32.0, 'r': 21.0},
        'HEM260': {'h': 290.0, 'b': 268.0, 'tw': 18.0, 'tf': 32.5, 'r': 24.0},
        'HEM280': {'h': 310.0, 'b': 288.0, 'tw': 18.5, 'tf': 33.0, 'r': 24.0},
        'HEM300': {'h': 340.0, 'b': 310.0, 'tw': 21.0, 'tf': 39.0, 'r': 27.0},
        'HEM320': {'h': 359.0, 'b': 309.0, 'tw': 21.0, 'tf': 40.0, 'r': 27.0},
        'HEM340': {'h': 377.0, 'b': 309.0, 'tw': 21.0, 'tf': 40.0, 'r': 27.0},
        'HEM360': {'h': 395.0, 'b': 308.0, 'tw': 21.0, 'tf': 40.0, 'r': 27.0},
        'HEM400': {'h': 432.0, 'b': 307.0, 'tw': 21.0, 'tf': 40.0, 'r': 27.0},
        'HEM450': {'h': 478.0, 'b': 307.0, 'tw': 21.0, 'tf': 40.0, 'r': 27.0},
        'HEM500': {'h': 524.0, 'b': 306.0, 'tw': 21.0, 'tf': 40.0, 'r': 27.0},
        'HEM550': {'h': 572.0, 'b': 306.0, 'tw': 21.0, 'tf': 40.0, 'r': 27.0},
        'HEM600': {'h': 620.0, 'b': 305.0, 'tw': 21.0, 'tf': 40.0, 'r': 27.0},
        'HEM650': {'h': 668.0, 'b': 305.0, 'tw': 21.0, 'tf': 40.0, 'r': 27.0},
        'HEM700': {'h': 716.0, 'b': 304.0, 'tw': 21.0, 'tf': 40.0, 'r': 27.0},
        'HEM800': {'h': 814.0, 'b': 303.0, 'tw': 21.0, 'tf': 40.0, 'r': 30.0},
        'HEM900': {'h': 910.0, 'b': 302.0, 'tw': 21.0, 'tf': 40.0, 'r': 30.0},
        'HEM1000': {
            'h': 1008.0,
            'b': 302.0,
            'tw': 21.0,
            'tf': 40.0,
            'r': 30.0,
        },
    }

    @classmethod
    def get_polygon(cls, name: str) -> Polygon:
        """Returns a shapely polygon representing an HE section."""
        parameters = cls.parameters.get(name)
        if parameters is None:
            raise ValueError(
                f"Profile '{name}' not found in HE sections. "
                "Select a valid profile (available ones: "
                f"{cls.profiles()})"
            )
        return _create_I_section(**parameters)

    @classmethod
    def profiles(cls) -> list:
        """Returns a list containing all available profiles."""
        return list(cls.parameters.keys())

    def __init__(self, name: str) -> None:
        """Creates a new HE object."""
        parameters = self.parameters.get(name)
        if parameters is None:
            raise ValueError(
                f"Profile '{name}' not found in HE sections. "
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
            Polygon: The represention of the HE section.
        """
        return self._polygon

    @property
    def h(self) -> float:
        """Returns height of HE section.

        Returns:
            float: Height h of HE section.
        """
        return self._h

    @property
    def b(self) -> float:
        """Returns width of HE section.

        Returns:
            float: Width b of HE section.
        """
        return self._b

    @property
    def tw(self) -> float:
        """Returns thickness of web of HE section.

        Returns:
            float: Web thickness tw of HE section.
        """
        return self._tw

    @property
    def tf(self) -> float:
        """Returns thickness of flange of HE section.

        Returns:
            float: Flange thickness tw of HE section.
        """
        return self._tf

    @property
    def r(self) -> float:
        """Returns fillet radius of HE section.

        Returns:
            float: Fillet radius r of HE section.
        """
        return self._r
