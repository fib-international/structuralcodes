"""The fib Model Code 2020."""

import typing as t

from ._concrete_material_properties import (
    GF,
    GF_T,
    RH_T,
    D_Ec,
    D_fc,
    D_fct,
    E_ci,
    E_ci_red_el,
    E_ci_T,
    E_ci_t,
    E_cm_hT_T,
    Ec1,
    El_ci,
    El_ci_T,
    GF_l,
    K_hT,
    RH_eq,
    S_Ec,
    S_fc,
    S_fct,
    alpha_bs,
    alpha_ds_1,
    alpha_ds_2,
    alpha_Ec,
    alpha_fc,
    alpha_fct,
    alpha_sT,
    alpha_T,
    beta_bc_fcm,
    beta_bc_t_t0,
    beta_bs_t,
    beta_c_sus,
    beta_cc,
    beta_dc_fcm,
    beta_dc_t0,
    beta_dc_t_t0,
    beta_ds_t_ts,
    beta_E,
    beta_h,
    beta_h_T,
    beta_lcc,
    beta_RH,
    beta_RH_c,
    beta_RH_s,
    beta_RH_T,
    beta_sT,
    beta_t0,
    eps_c1,
    eps_c_lim,
    eps_c_sigma,
    eps_c_T,
    eps_cbs_0,
    eps_cbs_t,
    eps_cc,
    eps_cds_0_fcm,
    eps_cds_t_ts,
    eps_cs_t_ts,
    eps_lc1,
    eps_lcs_t_ts,
    fck,
    fck_cube,
    fcm,
    fcm_hT_T,
    fcm_sus,
    fcm_T,
    fcm_t,
    fctk_max,
    fctk_min,
    fctk_sus,
    fctm,
    fctm_from_flexural,
    fctm_from_splitting,
    fctm_hT_T,
    fctm_T,
    flcm,
    flcm_T,
    flcm_t,
    flctk_max,
    flctk_min,
    flctm,
    gamma_t0,
    hT_Tref,
    multiaxial_stress_equation,
    nu_c,
    phi_bc_T,
    phi_bc_t_t0,
    phi_dc_T,
    phi_dc_t_t0,
    phi_l_t_t0,
    phi_sigma_t_t0,
    phi_t_t0,
    phi_t_t0_T,
    sC,
    sigma_c,
    sigma_crack_friction,
    sigma_ct_cracked,
    sigma_ct_uncracked,
    slC,
    t0_adj,
    t_T_maturity,
    tau_crack_friction,
)

__all__ = [
    'GF',
    'GF_T',
    'RH_T',
    'D_Ec',
    'D_fc',
    'D_fct',
    'E_ci',
    'E_ci_red_el',
    'E_ci_T',
    'E_ci_t',
    'E_cm_hT_T',
    'Ec1',
    'El_ci',
    'El_ci_T',
    'GF_l',
    'K_hT',
    'RH_eq',
    'S_Ec',
    'S_fc',
    'S_fct',
    'alpha_bs',
    'alpha_ds_1',
    'alpha_ds_2',
    'alpha_Ec',
    'alpha_fc',
    'alpha_fct',
    'alpha_sT',
    'alpha_T',
    'beta_bc_fcm',
    'beta_bc_t_t0',
    'beta_bs_t',
    'beta_c_sus',
    'beta_cc',
    'beta_dc_fcm',
    'beta_dc_t0',
    'beta_dc_t_t0',
    'beta_ds_t_ts',
    'beta_E',
    'beta_h',
    'beta_h_T',
    'beta_lcc',
    'beta_RH',
    'beta_RH_c',
    'beta_RH_s',
    'beta_RH_T',
    'beta_sT',
    'beta_t0',
    'eps_c1',
    'eps_c_lim',
    'eps_c_sigma',
    'eps_c_T',
    'eps_cbs_0',
    'eps_cbs_t',
    'eps_cc',
    'eps_cds_0_fcm',
    'eps_cds_t_ts',
    'eps_cs_t_ts',
    'eps_lc1',
    'eps_lcs_t_ts',
    'fck',
    'fck_cube',
    'fcm',
    'fcm_hT_T',
    'fcm_sus',
    'fcm_T',
    'fcm_t',
    'fctk_max',
    'fctk_min',
    'fctk_sus',
    'fctm',
    'fctm_from_flexural',
    'fctm_from_splitting',
    'fctm_hT_T',
    'fctm_T',
    'flcm',
    'flcm_T',
    'flcm_t',
    'flctk_max',
    'flctk_min',
    'flctm',
    'gamma_t0',
    'hT_Tref',
    'multiaxial_stress_equation',
    'nu_c',
    'phi_bc_T',
    'phi_bc_t_t0',
    'phi_dc_T',
    'phi_dc_t_t0',
    'phi_l_t_t0',
    'phi_sigma_t_t0',
    'phi_t_t0',
    'phi_t_t0_T',
    'sC',
    'sigma_c',
    'sigma_crack_friction',
    'sigma_ct_cracked',
    'sigma_ct_uncracked',
    'slC',
    't0_adj',
    't_T_maturity',
    'tau_crack_friction',
]

__title__: str = 'fib Model Code 2020'
__year__: str = '2024'
__materials__: t.Tuple[str] = ('concrete', 'reinforcement')
