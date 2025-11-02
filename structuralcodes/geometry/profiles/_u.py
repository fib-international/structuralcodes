"""U profiles."""

from shapely import (
    Polygon,
)

from ._base_profile import BaseProfile
from ._common_functions import (
    _create_taper_U_section,
)


class U(BaseProfile):
    """Simple class for representing an U profile.

    EN 10279.

    Taper flange Channels.

    8% slope in flange.
    """

    parameters = {
        'U40x20': {
            'h': 40.0,
            'b': 20.0,
            'tw': 5.0,
            'tf': 5.5,
            'r1': 5.0,
            'r2': 2.5,
        },
        'U50x25': {
            'h': 50.0,
            'b': 25.0,
            'tw': 5.0,
            'tf': 6.0,
            'r1': 6.0,
            'r2': 3.0,
        },
        'U60x30': {
            'h': 60.0,
            'b': 30.0,
            'tw': 6.0,
            'tf': 6.0,
            'r1': 6.0,
            'r2': 3.0,
        },
        'U65x42': {
            'h': 65.0,
            'b': 42.0,
            'tw': 5.5,
            'tf': 7.5,
            'r1': 7.5,
            'r2': 4.0,
        },
    }

    @classmethod
    def get_polygon(cls, name: str) -> Polygon:
        """Returns a shapely polygon representing an U section."""
        parameters = cls.parameters.get(name)
        if parameters is None:
            raise ValueError(
                f"Profile '{name}' not found in U sections. "
                "Select a valid profile (available ones: "
                f"{cls.profiles()})"
            )
        parameters['slope'] = 0.08
        parameters['u'] = parameters['b'] / 2.0
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
        """Creates a new U object."""
        parameters = self.parameters.get(name)
        if parameters is None:
            raise ValueError(
                f"Profile '{name}' not found in U sections. "
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
        self._flange_slope = 0.08
        self._u = self._b / 2.0
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
            Polygon: The represention of the U section.
        """
        return self._polygon

    @property
    def h(self) -> float:
        """Returns height of U section.

        Returns:
            float: Height h of U section.
        """
        return self._h

    @property
    def b(self) -> float:
        """Returns width of U section.

        Returns:
            float: Width b of U section.
        """
        return self._b

    @property
    def tw(self) -> float:
        """Returns thickness of web of U section.

        Returns:
            float: Web thickness tw of U section.
        """
        return self._tw

    @property
    def tf(self) -> float:
        """Returns thickness of flange of U section.

        Returns:
            float: Flange thickness tw of U section.
        """
        return self._tf

    @property
    def r1(self) -> float:
        """Returns fillet radius of U section.

        Returns:
            float: Fillet radius r1 of U section.
        """
        return self._r1

    @property
    def r2(self) -> float:
        """Returns fillet radius of U section.

        Returns:
            float: Fillet radius r2 of U section.
        """
        return self._r2
