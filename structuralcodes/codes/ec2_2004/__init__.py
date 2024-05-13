"""EUROCODE 2 1992-1-1:2004."""

import typing as t

from ._concrete_material_properties import (
    Ecm,
    eps_c1,
    eps_c2,
    eps_c3,
    eps_cu1,
    eps_cu2,
    eps_cu3,
    fcd,
    fcm,
    fctk_5,
    fctk_95,
    fctm,
    n_parabolic_rectangular,
)
from ._section_7_3_crack_control import (
    As_min,
    As_min_2,
    As_min_p,
    alpha_e,
    eps_sm_eps_cm,
    hc_eff,
    k,
    k1,
    k2,
    k3,
    k4,
    kc_flanges_area,
    kc_rect_area,
    kc_tension,
    kt,
    phi_eq,
    rho_p_eff,
    sr_max_close,
    sr_max_far,
    sr_max_theta,
    w_max,
    w_spacing,
    wk,
    xi1,
)

__all__ = [
    'As_min',
    'As_min_2',
    'As_min_p',
    'alpha_e',
    'eps_sm_eps_cm',
    'hc_eff',
    'k',
    'k1',
    'k2',
    'k3',
    'k4',
    'kc_flanges_area',
    'kc_rect_area',
    'kc_tension',
    'kt',
    'phi_eq',
    'rho_p_eff',
    'sr_max_close',
    'sr_max_far',
    'sr_max_theta',
    'w_max',
    'w_spacing',
    'wk',
    'xi1',
    'fcd',
    'fcm',
    'fctm',
    'fctk_5',
    'fctk_95',
    'Ecm',
    'eps_c1',
    'eps_cu1',
    'eps_c2',
    'eps_cu2',
    'n_parabolic_rectangular',
    'eps_c3',
    'eps_cu3',
]

__title__: str = 'EUROCODE 2 1992-1-1:2004'
__year__: str = '2004'
__materials__: t.Tuple[str] = ('concrete', 'reinforcement')
