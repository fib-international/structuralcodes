"""The unit manager."""

import typing as t
from dataclasses import dataclass

from numpy.typing import ArrayLike

# Type annotations for units
_LENGTH_LITERAL = t.Literal['mm', 'm']
_FORCE_LITERAL = t.Literal['N', 'kN', 'MN']

# Unit conversion
_METER = 1.0
_NEWTON = 1.0

_UNITS = {
    'length': {
        'm': _METER,
        'mm': _METER * 1e-3,
    },
    'force': {
        'N': _NEWTON,
        'kN': _NEWTON * 1e3,
        'MN': _NEWTON * 1e6,
    },
}


@dataclass
class UnitSet:
    """A set of units that can be used in the UnitConverter."""

    length: _LENGTH_LITERAL = 'mm'
    force: _FORCE_LITERAL = 'N'

    @property
    def length_unit(self):
        return _UNITS['length'][self.length]

    @property
    def force_unit(self):
        return _UNITS['force'][self.force]

    @property
    def stress_unit(self):
        return self.force_unit / self.length_unit**2

    def __post_init__(self):
        """Validate the provided units."""
        for attr in ('length', 'force'):
            try:
                # Try to find the provided unit among the available
                _UNITS[attr][getattr(self, attr)]
            except KeyError:
                # The provided unit was not found
                # Try to see if there was perhaps a mix of upper and lower case
                # letters
                for key in _UNITS[attr]:
                    if key.lower() == getattr(self, attr).lower().strip():
                        setattr(self, attr, key)
                        break
                else:
                    raise ValueError(
                        f'{getattr(self, attr)} is not a valid {attr} unit. '
                        f'Use one of {", ".join(_UNITS[attr].keys())}.'
                    )


class UnitConverter:
    """A class responsible for converting between different sets of units."""

    _default_units: UnitSet
    _alternative_units: UnitSet

    def __init__(
        self,
        default_units: UnitSet,
        alternative_units: t.Optional[UnitSet] = None,
    ) -> None:
        self._default_units = default_units
        self._alternative_units = alternative_units

    def convert_stress_to_default(
        self, stress: t.Union[float, ArrayLike]
    ) -> t.Union[float, ArrayLike]:
        """Convert stresses from alternative units to default units."""
        return (
            stress
            * self.alternative_units.stress_unit
            / self.default_units.stress_unit
        )

    def convert_stress_from_default(
        self, stress: t.Union[float, ArrayLike]
    ) -> t.Union[float, ArrayLike]:
        """Convert stresses from default units to alternative units."""
        return (
            stress
            * self.default_units.stress_unit
            / self.alternative_units.stress_unit
        )

    def convert_length_to_default(
        self, length: t.Union[float, ArrayLike]
    ) -> t.Union[float, ArrayLike]:
        """Convert lengths from alternative units to default units."""
        return (
            length
            * self.alternative_units.length_unit
            / self.default_units.length_unit
        )

    def convert_length_from_default(
        self, length: t.Union[float, ArrayLike]
    ) -> t.Union[float, ArrayLike]:
        """Convert lengths from default units to alternative units."""
        return (
            length
            * self.default_units.length_unit
            / self.alternative_units.length_unit
        )

    def convert_force_to_default(
        self, force: t.Union[float, ArrayLike]
    ) -> t.Union[float, ArrayLike]:
        """Convert forces from alternative units to default units."""
        return (
            force
            * self.alternative_units.force_unit
            / self.default_units.force_unit
        )

    def convert_force_from_default(
        self, force: t.Union[float, ArrayLike]
    ) -> t.Union[float, ArrayLike]:
        """Convert forces from default units to alternative units."""
        return (
            force
            * self.default_units.force_unit
            / self.alternative_units.force_unit
        )

    @property
    def default_units(self) -> UnitSet:
        """The default set of units."""
        return self._default_units

    @property
    def alternative_units(self) -> UnitSet:
        """The alternative set of units."""
        return (
            self._alternative_units
            if self._alternative_units is not None
            else self._default_units
        )
