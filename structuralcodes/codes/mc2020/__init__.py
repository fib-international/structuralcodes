"""The fib Model Code 2020."""

import typing as t

__title__: str = 'fib Model Code 2020'
__year__: str = '2024'
__materials__: t.Tuple[str] = ('concrete', 'reinforcement')

from ._corrosion import kc, pcorr_fractile, pcorr_rep

__all__ = ['pcorr_rep', 'pcorr_fractile', 'kc']
