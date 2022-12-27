"""EUROCODE 2 1992-1-1:2004"""
import typing as t

from ._crack_control import w_max

__all__ = ['w_max']

__title__: str = 'EUROCODE 2 1992-1-1'
__year__: str = '2004'
__materials__: t.Tuple[str] = ('concrete',)
