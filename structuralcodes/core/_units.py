"""Classes related to unit handling."""

import typing as t
from dataclasses import dataclass

from numpy.typing import ArrayLike

# Type annotations for units
_LENGTH_LITERAL = t.Literal['mm', 'm', 'inch', 'foot']
_FORCE_LITERAL = t.Literal['N', 'kN', 'MN']

# Unit conversion
_MILLIMETER = 1.0  # Only used as a reference and has no physical meaning
_NEWTON = 1e-3  # Only used as a reference and has no physical meaning

_UNITS = {
    'length': {
        'm': _MILLIMETER * 1e3,  # 1000 mm in one m
        'mm': _MILLIMETER,
        'inch': _MILLIMETER * 25.4,  # 25.4 mm in one inch
        'foot': _MILLIMETER * 304.8,  # 304.8 mm in one foot
    },
    'force': {
        'N': _NEWTON,
        'kN': _NEWTON * 1e3,  # 1000 N in one kN
        'MN': _NEWTON * 1e6,  # 1000000 N in one MN
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

    _from_units: UnitSet
    _to_units: UnitSet

    def __init__(
        self,
        from_units: UnitSet,
        to_units: t.Optional[UnitSet] = None,
    ) -> None:
        """Initialize a UnitConverter.

        Args:
            from_units (UnitSet): The set of units to convert forwards from or
                backwards to.
            to_units (Optional(UnitSet)): The set of units to convert forwards
                to or backwards from. If None, it is treated as equal to
                from_units, and no conversion happens.
        """
        self._from_units = from_units
        self._to_units = to_units

    def convert_stress_backwards(
        self, stress: t.Union[float, ArrayLike]
    ) -> t.Union[float, ArrayLike]:
        """Convert stress backwards."""
        if self.from_units == self.to_units:
            return stress
        return stress * self.to_units.stress_unit / self.from_units.stress_unit

    def convert_stress_forwards(
        self, stress: t.Union[float, ArrayLike]
    ) -> t.Union[float, ArrayLike]:
        """Convert stress forwards."""
        if self.from_units == self.to_units:
            return stress
        return stress * self.from_units.stress_unit / self.to_units.stress_unit

    def convert_length_backwards(
        self, length: t.Union[float, ArrayLike]
    ) -> t.Union[float, ArrayLike]:
        """Convert length backwards."""
        if self.from_units == self.to_units:
            return length
        return length * self.to_units.length_unit / self.from_units.length_unit

    def convert_length_forwards(
        self, length: t.Union[float, ArrayLike]
    ) -> t.Union[float, ArrayLike]:
        """Convert length forwards."""
        if self.from_units == self.to_units:
            return length
        return length * self.from_units.length_unit / self.to_units.length_unit

    def convert_force_backwards(
        self, force: t.Union[float, ArrayLike]
    ) -> t.Union[float, ArrayLike]:
        """Convert force backwards."""
        if self.from_units == self.to_units:
            return force
        return force * self.to_units.force_unit / self.from_units.force_unit

    def convert_force_forwards(
        self, force: t.Union[float, ArrayLike]
    ) -> t.Union[float, ArrayLike]:
        """Convert length forwards."""
        if self.from_units == self.to_units:
            return force
        return force * self.from_units.force_unit / self.to_units.force_unit

    @property
    def from_units(self) -> UnitSet:
        """The set of units to convert from, i.e. we convert forwards from this
        set of units.
        """
        return self._from_units

    @property
    def to_units(self) -> UnitSet:
        """The set of units to convert to, i.e. we convert forwards to this set
        of units.
        """
        return (
            self._to_units if self._to_units is not None else self._from_units
        )
