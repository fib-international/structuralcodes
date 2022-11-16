"""The fib Model Code 2010"""
import typing as t
from ._concrete_material_properties import fcm, fctm, fctkmin, fctkmax, Gf, fcd

__all__ = [
    'fcm',
    'fctm',
    'fctkmin',
    'fctkmax',
    'Gf',
    'fcd',
]

__title__: str = 'fib Model Code 2010'
__year__: str = '2013'
__materials__: t.Tuple[str] = ('concrete',)
