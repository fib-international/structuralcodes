"""EUROCODE 2 1992-1-1:2023."""

import typing as t

from ._annexB_time_dependent import alpha_c
from ._section5_materials import (
    A_phi_correction_exp,
    Ecm,
    Es,
    alpha_c_th,
    alpha_s_th,
    eps_c1,
    eps_c2,
    eps_cs_50y,
    eps_cu1,
    eps_cu2,
    eps_ud,
    eta_cc,
    fcd,
    fcm,
    fctd,
    fctk_5,
    fctk_95,
    fctm,
    fpd,
    fyd,
    hn,
    k_sargin,
    k_tc,
    k_tt,
    n_parabolic_rectangular,
    p_steel_stress_params,
    phi_50y_t0,
    phi_correction_factor,
    reinforcement_duct_props,
    sigma_p,
    sigma_s,
    weight_c,
    weight_s,
)
from ._section9_sls import (
    As_min_y,
    Ec_eff,
    delta_simpl,
    epssm_epscm,
    k_1_r,
    kfl,
    kh,
    srm_cal,
    wk_cal,
    wk_cal2,
)
from ._section_8_3_torsion import (
    VEd_i,
    tau_t_i,
    tau_t_rd,
    tau_t_rd_max,
    tau_t_rd_sl,
    tau_t_rd_sw,
)

__all__ = [
    'tau_t_rd_max',
    'tau_t_rd',
    'tau_t_rd_sl',
    'tau_t_rd_sw',
    'tau_t_i',
    'VEd_i',
    'A_phi_correction_exp',
    'alpha_c_th',
    'alpha_s_th',
    'Ecm',
    'eps_c1',
    'eps_cs_50y',
    'eps_cu1',
    'k_sargin',
    'eps_c2',
    'eps_cu2',
    'n_parabolic_rectangular',
    'eps_ud',
    'Es',
    'eta_cc',
    'fcd',
    'fcm',
    'fctd',
    'fctk_5',
    'fctk_95',
    'fctm',
    'fpd',
    'fyd',
    'hn',
    'k_tc',
    'k_tt',
    'p_steel_stress_params',
    'phi_50y_t0',
    'phi_correction_factor',
    'reinforcement_duct_props',
    'sigma_c',
    'sigma_p',
    'sigma_s',
    'weight_c',
    'weight_s',
    'alpha_c',
    'As_min_y',
    'delta_simpl',
    'Ec_eff',
    'epssm_epscm',
    'k_1_r',
    'kfl',
    'kh',
    'srm_cal',
    'wk_cal',
    'wk_cal2',
]

__title__: str = 'EUROCODE 2 1992-1-1:2023'
__year__: str = '2023'
__materials__: t.Tuple[str] = ('concrete', 'reinforcement')
