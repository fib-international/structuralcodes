"""ACI 318-19."""

import typing as t

from ._concrete_material_properties import (
    Ec,
    alpha1,
    beta1,
    eps_cu,
    fct,
    fr,
    lambda_factor,
)
from ._reinforcement_material_properties import (
    Es,
    epsyd,
    fy_design,
    reinforcement_grade_props,
)

__all__ = [
    'Ec',
    'Es',
    'alpha1',
    'beta1',
    'eps_cu',
    'epsyd',
    'fct',
    'fr',
    'fy_design',
    'lambda_factor',
    'reinforcement_grade_props',
]

__title__: str = 'ACI 318-19'
__year__: str = '2019'
__materials__: t.Tuple[str] = ('concrete', 'reinforcement')
