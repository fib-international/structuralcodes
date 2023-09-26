"""EUROCODE 2 1992-1-1:2004."""
import typing as t

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
]

__title__: str = 'EUROCODE 2 1992-1-1'
__year__: str = '2004'
__materials__: t.Tuple[str] = ('concrete',)
