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

    def __init__(self):
        """Initialize a shell reinforcement."""
        raise NotImplementedError


class ShellGeometry(Geometry):
    """A class for a shell with a thickness and material."""

    _reinforcement: t.List[ShellReinforcement]

    def __init__(
        self,
        thickness: float,
        material: Material,
        name: t.Optional[str] = None,
        group_label: t.Optional[str] = None,
    ):
        """Initialize a shell geometry."""
        super().__init__(name=name, group_label=group_label)
        self._thickness = thickness
        self._material = material

    def add_reinforcement(
        self,
        reinforcement: t.Union[ShellReinforcement, t.List[ShellReinforcement]],
    ):
        """Add reinforcement to the shell geometry."""
        if isinstance(reinforcement, ShellReinforcement):
            self._reinforcement.append(reinforcement)
        elif isinstance(reinforcement, list):
            self._reinforcement.extend(reinforcement)
