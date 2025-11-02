"""HP profiles."""

from shapely import (
    Polygon,
)

from ._base_profile import BaseProfile
from ._common_functions import (
    _create_I_section,
)


class HP(BaseProfile):
    """Simple class for representing an HP profile.

    HP profiles.
    """

    parameters = {
        'HP200x43': {'h': 200.0, 'b': 205.0, 'tw': 9.0, 'tf': 9.0, 'r': 10.0},
        'HP200x53': {
            'h': 204.0,
            'b': 207.0,
            'tw': 11.3,
            'tf': 11.3,
            'r': 10.0,
        },
        'HP220x57': {
            'h': 210.0,
            'b': 224.5,
            'tw': 11.0,
            'tf': 11.0,
            'r': 18.0,
        },
        'HP260x75': {
            'h': 249.0,
            'b': 265.0,
            'tw': 12.0,
            'tf': 12.0,
            'r': 24.0,
        },
        'HP260x87': {
            'h': 253.0,
            'b': 267.0,
            'tw': 14.0,
            'tf': 14.0,
            'r': 24.0,
        },
        'HP305x79': {
            'h': 299.3,
            'b': 306.4,
            'tw': 11.0,
            'tf': 11.1,
            'r': 15.0,
        },
        'HP305x88': {
            'h': 301.7,
            'b': 307.8,
            'tw': 12.4,
            'tf': 12.3,
            'r': 15.0,
        },
        'HP305x95': {
            'h': 303.7,
            'b': 308.7,
            'tw': 13.3,
            'tf': 13.3,
            'r': 15.0,
        },
        'HP305x110': {
            'h': 307.9,
            'b': 310.7,
            'tw': 15.3,
            'tf': 15.4,
            'r': 15.0,
        },
        'HP305x126': {
            'h': 312.3,
            'b': 312.9,
            'tw': 17.5,
            'tf': 17.6,
            'r': 15.0,
        },
        'HP305x149': {
            'h': 318.5,
            'b': 316.0,
            'tw': 20.6,
            'tf': 20.7,
            'r': 15.0,
        },
        'HP305x180': {
            'h': 326.7,
            'b': 319.7,
            'tw': 24.8,
            'tf': 24.8,
            'r': 15.0,
        },
        'HP305x186': {
            'h': 328.3,
            'b': 320.9,
            'tw': 25.5,
            'tf': 25.6,
            'r': 15.0,
        },
        'HP305x223': {
            'h': 337.9,
            'b': 325.7,
            'tw': 30.3,
            'tf': 30.4,
            'r': 15.0,
        },
        'HP320x88': {
            'h': 303.0,
            'b': 304.0,
            'tw': 12.0,
            'tf': 12.0,
            'r': 27.0,
        },
        'HP320x103': {
            'h': 307.0,
            'b': 306.0,
            'tw': 14.0,
            'tf': 14.0,
            'r': 27.0,
        },
        'HP320x117': {
            'h': 311.0,
            'b': 308.0,
            'tw': 16.0,
            'tf': 16.0,
            'r': 27.0,
        },
        'HP320x147': {
            'h': 319.0,
            'b': 312.0,
            'tw': 20.0,
            'tf': 20.0,
            'r': 27.0,
        },
        'HP320x184': {
            'h': 329.0,
            'b': 317.0,
            'tw': 25.0,
            'tf': 25.0,
            'r': 27.0,
        },
        'HP360x109': {
            'h': 346.4,
            'b': 371.0,
            'tw': 12.8,
            'tf': 12.9,
            'r': 15.0,
        },
        'HP360x133': {
            'h': 352.0,
            'b': 373.8,
            'tw': 15.6,
            'tf': 15.7,
            'r': 15.0,
        },
        'HP360x152': {
            'h': 356.4,
            'b': 376.0,
            'tw': 17.8,
            'tf': 17.9,
            'r': 15.0,
        },
        'HP360x174': {
            'h': 361.4,
            'b': 378.5,
            'tw': 20.3,
            'tf': 20.4,
            'r': 15.0,
        },
        'HP360x180': {
            'h': 362.9,
            'b': 378.8,
            'tw': 21.1,
            'tf': 21.1,
            'r': 15.0,
        },
        'HP400x122': {
            'h': 348.0,
            'b': 390.0,
            'tw': 14.0,
            'tf': 14.0,
            'r': 15.0,
        },
        'HP400x140': {
            'h': 352.0,
            'b': 392.0,
            'tw': 16.0,
            'tf': 16.0,
            'r': 15.0,
        },
        'HP400x158': {
            'h': 356.0,
            'b': 394.0,
            'tw': 18.0,
            'tf': 18.0,
            'r': 15.0,
        },
        'HP400x176': {
            'h': 360.0,
            'b': 396.0,
            'tw': 20.0,
            'tf': 20.0,
            'r': 15.0,
        },
        'HP400x194': {
            'h': 364.0,
            'b': 398.0,
            'tw': 22.0,
            'tf': 22.0,
            'r': 15.0,
        },
        'HP400x213': {
            'h': 368.0,
            'b': 400.0,
            'tw': 24.0,
            'tf': 24.0,
            'r': 15.0,
        },
        'HP400x231': {
            'h': 372.0,
            'b': 402.0,
            'tw': 26.0,
            'tf': 26.0,
            'r': 15.0,
        },
    }

    @classmethod
    def get_polygon(cls, name: str) -> Polygon:
        """Returns a shapely polygon representing an HP section."""
        parameters = cls.parameters.get(name)
        if parameters is None:
            raise ValueError(
                f"Profile '{name}' not found in HP sections. "
                "Select a valid profile (available ones: "
                f"{cls.profiles()})"
            )
        return _create_I_section(**parameters)

    @classmethod
    def profiles(cls) -> list:
        """Returns a list containing all available profiles."""
        return list(cls.parameters.keys())

    def __init__(self, name: str) -> None:
        """Creates a new HP object."""
        parameters = self.parameters.get(name)
        if parameters is None:
            raise ValueError(
                f"Profile '{name}' not found in HP sections. "
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
            Polygon: The represention of the HP section.
        """
        return self._polygon

    @property
    def h(self) -> float:
        """Returns height of HP section.

        Returns:
            float: Height h of HP section.
        """
        return self._h

    @property
    def b(self) -> float:
        """Returns width of HP section.

        Returns:
            float: Width b of HP section.
        """
        return self._b

    @property
    def tw(self) -> float:
        """Returns thickness of web of HP section.

        Returns:
            float: Web thickness tw of HP section.
        """
        return self._tw

    @property
    def tf(self) -> float:
        """Returns thickness of flange of HP section.

        Returns:
            float: Flange thickness tw of HP section.
        """
        return self._tf

    @property
    def r(self) -> float:
        """Returns fillet radius of HP section.

        Returns:
            float: Fillet radius r of HP section.
        """
        return self._r
