"""The fib Model Code 2010."""

import typing as t

from ._concrete_material_properties import Gf, fcm, fctkmax, fctkmin, fctm
from ._reinforcement_material_properties import fyd

__all__ = [
    'fcm',
    'fctm',
    'fctkmin',
    'fctkmax',
    'Gf',
    'fyd',
]

__title__: str = 'fib Model Code 2010'
__year__: str = '2013'
__materials__: t.Tuple[str] = ('concrete', 'reinforcement')
