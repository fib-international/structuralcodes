"""The fib Model Code 2010."""

import typing as t

from ._concrete_material_properties import (
    E_ci,
    Gf,
    eps_c1,
    eps_c2,
    eps_c3,
    eps_cu1,
    eps_cu2,
    eps_cu3,
    fcd,
    fcm,
    fctkmax,
    fctkmin,
    fctm,
    k_sargin,
    n_parabolic_rectangular,
)
from ._reinforcement_material_properties import (
    epsud,
    fyd,
    reinforcement_duct_props,
)

__all__ = [
    'fcm',
    'E_ci',
    'fctm',
    'fctkmin',
    'fctkmax',
    'fcd',
    'Gf',
    'eps_c1',
    'eps_cu1',
    'k_sargin',
    'eps_c2',
    'eps_cu2',
    'n_parabolic_rectangular',
    'eps_c3',
    'eps_cu3',
    'fyd',
    'epsud',
    'reinforcement_duct_props',
]

__title__: str = 'fib Model Code 2010'
__year__: str = '2013'
__materials__: t.Tuple[str] = ('concrete', 'reinforcement')
