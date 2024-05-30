"""Utility functions for importing standard steel profiles."""

import numpy as np
from shapely import (
    LineString,
    Point,
    Polygon,
    get_geometry,
    polygonize,
    set_precision,
)
from shapely.affinity import scale, translate
from shapely.geometry.polygon import orient
from shapely.ops import linemerge, unary_union


class IPE:
    """Simple class for representing an IPE profile.

    IPE 80-600 in accordance with standard Euronorm 19-57
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
            Polygon: the represention of the IPE section.
        """
        return self._polygon

    @property
    def h(self) -> Polygon:
        """Returns height of IPE section.

        Returns:
            float: height h of IPE section.
        """
        return self._h

    @property
    def b(self) -> Polygon:
        """Returns width of IPE section.

        Returns:
            float: width b of IPE section.
        """
        return self._b

    @property
    def tw(self) -> Polygon:
        """Returns thickness of web of IPE section.

        Returns:
            float: web thickness tw of IPE section.
        """
        return self._tw

    @property
    def tf(self) -> Polygon:
        """Returns thickness of flange of IPE section.

        Returns:
            float: flange thickness tw of IPE section.
        """
        return self._tf

    @property
    def r(self) -> Polygon:
        """Returns fillet radius of IPE section.

        Returns:
            float: fillet radius r of IPE section.
        """
        return self._r


class HE:
    """Simple class for representing an HE profile.

    HE A, HE B, HE M 100-1000 in accordance with Standard Euronorm 53-62
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
            Polygon: the represention of the HE section.
        """
        return self._polygon

    @property
    def h(self) -> Polygon:
        """Returns height of HE section.

        Returns:
            float: height h of HE section.
        """
        return self._h

    @property
    def b(self) -> Polygon:
        """Returns width of HE section.

        Returns:
            float: width b of HE section.
        """
        return self._b

    @property
    def tw(self) -> Polygon:
        """Returns thickness of web of HE section.

        Returns:
            float: web thickness tw of HE section.
        """
        return self._tw

    @property
    def tf(self) -> Polygon:
        """Returns thickness of flange of HE section.

        Returns:
            float: flange thickness tw of HE section.
        """
        return self._tf

    @property
    def r(self) -> Polygon:
        """Returns fillet radius of HE section.

        Returns:
            float: fillet radius r of HE section.
        """
        return self._r


class UB:
    """Simple class for representing a UB profile.

    Universal Beams
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
            Polygon: the represention of the UB section.
        """
        return self._polygon

    @property
    def h(self) -> Polygon:
        """Returns height of UB section.

        Returns:
            float: height h of UB section.
        """
        return self._h

    @property
    def b(self) -> Polygon:
        """Returns width of UB section.

        Returns:
            float: width b of UB section.
        """
        return self._b

    @property
    def tw(self) -> Polygon:
        """Returns thickness of web of UB section.

        Returns:
            float: web thickness tw of UB section.
        """
        return self._tw

    @property
    def tf(self) -> Polygon:
        """Returns thickness of flange of UB section.

        Returns:
            float: flange thickness tw of UB section.
        """
        return self._tf

    @property
    def r(self) -> Polygon:
        """Returns fillet radius of UB section.

        Returns:
            float: fillet radius r of UB section.
        """
        return self._r


class UC:
    """Simple class for representing a UC profile.

    Universal Columns
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
            Polygon: the represention of the UC section.
        """
        return self._polygon

    @property
    def h(self) -> Polygon:
        """Returns height of UC section.

        Returns:
            float: height h of UC section.
        """
        return self._h

    @property
    def b(self) -> Polygon:
        """Returns width of UC section.

        Returns:
            float: width b of UC section.
        """
        return self._b

    @property
    def tw(self) -> Polygon:
        """Returns thickness of web of UC section.

        Returns:
            float: web thickness tw of UC section.
        """
        return self._tw

    @property
    def tf(self) -> Polygon:
        """Returns thickness of flange of UC section.

        Returns:
            float: flange thickness tw of UC section.
        """
        return self._tf

    @property
    def r(self) -> Polygon:
        """Returns fillet radius of UC section.

        Returns:
            float: fillet radius r of UC section.
        """
        return self._r


class UBP:
    """Simple class for representing a UBP profile.

    Universal Bearing Pile
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
            Polygon: the represention of the UBP section.
        """
        return self._polygon

    @property
    def h(self) -> Polygon:
        """Returns height of UBP section.

        Returns:
            float: height h of UBP section.
        """
        return self._h

    @property
    def b(self) -> Polygon:
        """Returns width of UBP section.

        Returns:
            float: width b of UBP section.
        """
        return self._b

    @property
    def tw(self) -> Polygon:
        """Returns thickness of web of UBP section.

        Returns:
            float: web thickness tw of UBP section.
        """
        return self._tw

    @property
    def tf(self) -> Polygon:
        """Returns thickness of flange of UBP section.

        Returns:
            float: flange thickness tw of UBP section.
        """
        return self._tf

    @property
    def r(self) -> Polygon:
        """Returns fillet radius of UBP section.

        Returns:
            float: fillet radius r of UBP section.
        """
        return self._r


class IPN:
    """Simple class for representing an IPN profile.

    IPN  in accordance with standard

    14% slope in flange
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
            Polygon: the represention of the IPN section.
        """
        return self._polygon

    @property
    def h(self) -> Polygon:
        """Returns height of IPN section.

        Returns:
            float: height h of IPN section.
        """
        return self._h

    @property
    def b(self) -> Polygon:
        """Returns width of IPN section.

        Returns:
            float: width b of IPN section.
        """
        return self._b

    @property
    def tw(self) -> Polygon:
        """Returns thickness of web of IPN section.

        Returns:
            float: web thickness tw of IPN section.
        """
        return self._tw

    @property
    def tf(self) -> Polygon:
        """Returns thickness of flange of IPN section.

        Returns:
            float: flange thickness tw of IPN section.
        """
        return self._tf

    @property
    def r1(self) -> Polygon:
        """Returns fillet radius of IPN section.

        Returns:
            float: fillet radius r1 of IPN section.
        """
        return self._r1

    @property
    def r2(self) -> Polygon:
        """Returns fillet radius of IPN section.

        Returns:
            float: fillet radius r2 of IPN section.
        """
        return self._r2


class UPN:
    """Simple class for representing an UPN profile.

    European standard channels UPN 50 - 400
    Taper flange Channels

    14% slope in flange
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
            name = f'IPN{int(name):0d}'
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
        """Creates a new IPN object."""
        if isinstance(name, (float, int)):
            name = f'IPN{int(name):0d}'
        parameters = self.parameters.get(name)
        if parameters is None:
            raise ValueError(
                f"Profile '{name}' not found in UPN sections. "
                "Select a valid profile (available ones: "
                f"{self.profiles()})"
            )
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
            Polygon: the represention of the UPN section.
        """
        return self._polygon

    @property
    def h(self) -> Polygon:
        """Returns height of UPN section.

        Returns:
            float: height h of UPN section.
        """
        return self._h

    @property
    def b(self) -> Polygon:
        """Returns width of UPN section.

        Returns:
            float: width b of UPN section.
        """
        return self._b

    @property
    def tw(self) -> Polygon:
        """Returns thickness of web of UPN section.

        Returns:
            float: web thickness tw of UPN section.
        """
        return self._tw

    @property
    def tf(self) -> Polygon:
        """Returns thickness of flange of UPN section.

        Returns:
            float: flange thickness tw of UPN section.
        """
        return self._tf

    @property
    def r1(self) -> Polygon:
        """Returns fillet radius of UPN section.

        Returns:
            float: fillet radius r1 of UPN section.
        """
        return self._r1

    @property
    def r2(self) -> Polygon:
        """Returns fillet radius of UPN section.

        Returns:
            float: fillet radius r2 of UPN section.
        """
        return self._r2


def _create_I_section(h: float, b: float, tw: float, tf: float, r: float):
    """Returns a polygon for a I section."""
    # top flange
    top_flange = Polygon(
        [
            (-b / 2, -h / 2),
            (b / 2, -h / 2),
            (b / 2, -h / 2 + tf),
            (-b / 2, -h / 2 + tf),
        ]
    )
    # bottom flange
    bottom_flange = translate(top_flange, xoff=0, yoff=h - tf)
    web = Polygon(
        [
            (-tw / 2, -h / 2 + tf),
            (tw / 2, -h / 2 + tf),
            (tw / 2, h / 2 - tf),
            (-tw / 2, h / 2 - tf),
        ]
    )
    # fillets
    p = Point([tw / 2 + r, -h / 2 + tf + r]).buffer(r)
    s = Polygon(
        [
            (tw / 2, -h / 2 + tf),
            (tw / 2 + r, -h / 2 + tf),
            (tw / 2 + r, -h / 2 + tf + r),
            (tw / 2, -h / 2 + tf + r),
        ]
    )
    fillet = s.difference(p)
    p = Point([-tw / 2 - r, -h / 2 + tf + r]).buffer(r)
    s = Polygon(
        [
            (-tw / 2 - r, -h / 2 + tf),
            (-tw / 2, -h / 2 + tf),
            (-tw / 2, -h / 2 + tf + r),
            (-tw / 2 - r, -h / 2 + tf + r),
        ]
    )
    fillet = s.difference(p).union(fillet)
    fillet = translate(
        scale(fillet, 1, -1), xoff=0, yoff=h - 2 * tf - r
    ).union(fillet)
    # Create the geometry
    geometries = [
        set_precision(geometry, grid_size=1e-13)
        for geometry in [fillet, top_flange, bottom_flange, web]
    ]
    return orient(unary_union(geometries), 1)


def _create_taper_I_section(
    h: float,
    b: float,
    tw: float,
    tf: float,
    r1: float,
    r2: float,
    slope: float,
) -> Polygon:
    """Returns a shapely polygon representing a Taper Flange I Section."""
    # Create first part of line
    ls = [
        set_precision(
            LineString([[0, -h / 2], [b / 2, -h / 2]]), grid_size=1e-13
        )
    ]
    # Create first fillet
    xy = np.array(
        [
            [b / 2, -h / 2],
            [b / 2, -h / 2 + tf - (b / 4 * slope)],
            [b / 4, -h / 2 + tf],
        ]
    )
    ls.append(
        set_precision(
            LineString(xy).offset_curve(r2).offset_curve(-r2), grid_size=1e-13
        )
    )
    # Create second fillet
    xy = np.array(
        [
            [b / 4, -h / 2 + tf],
            [tw / 2, -h / 2 + tf + (b / 4 - tw / 2) * slope],
            [tw / 2, 0],
        ]
    )
    ls.append(
        set_precision(
            LineString(xy).offset_curve(-r1).offset_curve(r1), grid_size=1e-13
        )
    )

    # Merge filleted
    merged_ls = linemerge(ls)

    # mirror the parts
    merged_ls = linemerge(
        [
            merged_ls,
            translate(scale(merged_ls, -1, 1), -b / 2, 0),
        ]
    )
    merged_ls = linemerge(
        [
            merged_ls,
            translate(scale(merged_ls, 1, -1), 0, h / 2),
        ]
    )

    # Create a polygon
    poly = polygonize([merged_ls])

    # Return the first and only polygon of this collection
    return orient(get_geometry(poly, 0), 1)


def _create_taper_U_section(
    h: float,
    b: float,
    tw: float,
    tf: float,
    r1: float,
    r2: float,
    slope: float,
    u: float,
) -> Polygon:
    """Returns a shapely polygon representing a Taper Flange U Section."""
    # Create first part of line
    ls = [
        set_precision(
            LineString([[0, h / 2], [0, 0], [b, 0]]), grid_size=1e-13
        )
    ]
    # Create first fillet
    xy = np.array([[b, 0], [b, tf - slope * u], [b - u, tf]])
    ls.append(
        set_precision(
            LineString(xy).offset_curve(r2).offset_curve(-r2), grid_size=1e-13
        )
    )
    # Create second fillet
    xy = np.array(
        [
            [b - u, tf],
            [tw, tf + slope * (b - u - tw)],
            [tw, h / 2],
        ]
    )
    ls.append(
        set_precision(
            LineString(xy).offset_curve(-r1).offset_curve(r1), grid_size=1e-13
        )
    )

    # Merge filleted
    merged_ls = linemerge(ls)

    # mirror the parts
    merged_ls = linemerge(
        [
            merged_ls,
            translate(scale(merged_ls, 1, -1), 0, h / 2),
        ]
    )

    # Create a polygon
    poly = polygonize([merged_ls])

    # Get the first and only polygon of this collection
    poly = get_geometry(poly, 0)
    # Translate it to the centroid when returning
    return orient(
        translate(poly, xoff=-poly.centroid.x, yoff=-poly.centroid.y), 1
    )