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
    k_sargin,
    n_parabolic_rectangular,
)
from ._reinforcement_material_properties import (
    epsud,
    fyd,
    reinforcement_duct_props,
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
from .annex_b_shrink_and_creep import (
    beta_c,
    beta_fcm,
    beta_H,
    beta_RH,
    eps_cd_0,
    eps_cs,
    phi,
    phi_RH,
    t0_adj,
)
from .shear import (
    Asw_max,
    VEdmax_unreinf,
    VRdc,
    VRdc_prin_stress,
    VRdmax,
    VRds,
)

__all__ = [
    'As_min',
    'As_min_2',
    'As_min_p',
    'Asw_max',
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
    'reinforcement_duct_props',
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
    'k_sargin',
    'eps_c2',
    'eps_cu2',
    'n_parabolic_rectangular',
    'eps_c3',
    'eps_cu3',
    'fyd',
    'epsud',
    'VEdmax_unreinf',
    'VRdc',
    'VRdc_prin_stress',
    'VRdmax',
    'VRds',
    'beta_c',
    'beta_fcm',
    'beta_H',
    'beta_RH',
    'eps_cd_0',
    'eps_cs',
    'phi',
    'phi_RH',
    't0_adj',
]

__title__: str = 'EUROCODE 2 1992-1-1:2004'
__year__: str = '2004'
__materials__: t.Tuple[str] = ('concrete', 'reinforcement')
