"""ACI 318-19: Building Code Requirements for Structural Concrete."""

import typing as t

from ._shear import lambda_s
from ._two_way_shear import alpha_s, vc_two_way

__all__ = [
    'alpha_s',
    'lambda_s',
    'vc_two_way',
]

__title__: str = 'ACI 318-19'
__year__: str = '2019'
__materials__: t.Tuple[str, ...] = ()
