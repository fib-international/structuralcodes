"""A concrete implementation of a shell geometry."""

import typing as t

from ..core.base import Material
from ._geometry import Geometry


class ShellReinforcement(Geometry):
    """A class for representing reinforcement in a shell geometry."""

    _z: float
    _n_bars: float
    _cc_bars: float
    _diameter_bar: float
    _material: Material
    _phi: float

    def __init__(
        self,
        z: float,
        n_bars: float,
        cc_bars: float,
        diameter_bar: float,
        material: Material,
        phi: float,
        name: t.Optional[str] = None,
        group_label: t.Optional[str] = None,
    ) -> None:
        """Initialize a shell reinforcement."""
        super().__init__(name, group_label)

        self._z = z
        self._n_bars = n_bars
        self._cc_bars = cc_bars
        self._diameter_bar = diameter_bar
        self._material = material
        self._phi = phi

    @property
    def z(self) -> float:
        """Return the reinforcement position over the thickness."""
        return self._z

    @property
    def n_bars(self) -> float:
        """Return the number of bars per unit width."""
        return self._n_bars

    @property
    def cc_bars(self) -> float:
        """Return the spacing between bars."""
        return self._cc_bars

    @property
    def diameter_bar(self) -> float:
        """Return the bar diameter."""
        return self._diameter_bar

    @property
    def material(self) -> Material:
        """Return the material of the reinforcement."""
        return self._material

    @property
    def phi(self) -> float:
        """Return the orientation angle of the reinforcement."""
        return self._phi


class ShellGeometry(Geometry):
    """A class for a shell with a thickness and material."""

    _reinforcement: t.List[ShellReinforcement]

    def __init__(
        self,
        thickness: float,
        material: Material,
        name: t.Optional[str] = None,
        group_label: t.Optional[str] = None,
    ) -> None:
        """Initialize a shell geometry."""
        super().__init__(name=name, group_label=group_label)

        if thickness <= 0:
            raise ValueError('Shell thickness must be positive.')

        self._thickness = thickness
        self._material = material

        self._reinforcement = []

    @property
    def thickness(self) -> float:
        """Return the shell thickness."""
        return self._thickness

    @property
    def material(self) -> Material:
        """Return the material of the shell."""
        return self._material

    @property
    def reinforcement(self) -> t.List[ShellReinforcement]:
        """Return all reinforcement layers."""
        return self._reinforcement

    def add_reinforcement(
        self,
        reinforcement: t.Union[ShellReinforcement, t.List[ShellReinforcement]],
    ) -> None:
        """Add reinforcement to the shell geometry."""
        if isinstance(reinforcement, ShellReinforcement):
            self._reinforcement.append(reinforcement)
        elif isinstance(reinforcement, list):
            self._reinforcement.extend(reinforcement)
