"""fib MC2020 Chapter 30.5.2.4."""

import math
from typing import List, Literal, Optional

from scipy.optimize import newton


def wlim_rfc(exposure_class: Literal['X0', 'XC', 'XD', 'XS']) -> float:
    """Calculate the crack width limit for regular reinforced concrete (RC)
    based on the exposure class.

    fib Model Code 2020, Table (30.5-1)

    Args:
        exposure_class (str): Exposure class ('X0', 'XC', 'XD', 'XS').

    Returns:
        float: The crack width limit in mm.

    Raises:
        ValueError: If exposure_class is invalid.
    """
    if exposure_class in ('XC', 'XD', 'XS'):
        return 0.3

    return 0.4


def wlim_prestressed(
    exposure_class: Literal['X0', 'XC', 'XD', 'XS'],
    protection_level: Literal['PL1', 'PL2', 'PL3'],
) -> float:
    """Calculate the crack width limit (wlim) forprestressed concrete based on
    the exposure class.

    fib Model Code 2020, eq. (30.5-1)

    Args:
        exposure_class (str): Exposure class ('X0', 'XC', 'XD', 'XS').
        protection_level (str): Protection level required.

    Returns:
        float: The crack width limit en mm.

    Raises:
        ValueError: If exposure_class is invalid or required parameters are not
            provided.
    """
    if protection_level in ('PL2', 'PL3'):
        if exposure_class == 'X0':
            return 0.4
        return 0.3

    if exposure_class == 'X0':
        return 0.4
    if exposure_class == 'XC':
        return 0.2
    return 0.0


def wlim_frp() -> float:
    """Calculate the crack width limit (wlim) for FRP-reinforced concrete.

    fib Model Code 2020, (30.5.2.4.2)

    Returns:
        float: The crack width limit in mm.
    """
    return 0.7  # mm, relaxed limit for FRP-reinforced members


def wlim_fluid_tightness(self_healing: bool = False) -> float:
    """Calculate the crack width limit (wlim) for fluid-tight concrete.

    fib Model Code 2020, (30.5.2.4.2)

    Args:
        self_healing (bool): Whether self-healing of cracks is considered.

    Returns:
        float: The crack width limit in mm.
    """
    return 0.15 if not self_healing else 0.2  # mm, fluid-tight conditions


def Nr_crack(Ac: float, fctm: float, alpha_e: float, rho_s: float) -> float:
    """Calculate the cracking load Nr.

    fib Model Code 2020, eq. (30.5-2)

    Args:
        Ac (float): Area of the cross-section in mm2.
        fctm (float): Mean tensile strength of concrete in MPa.
        alpha_e (float): Coefficient that considers the ratio between the
            modulus of elasticity of steel and concrete.
        rho_s (float): Reinforcement ratio (As/Ac).

    Returns:
        float: The cracking load Nr in kN.

    Raises:
        ValueError: If any input is negative.
    """
    if Ac < 0 or fctm < 0 or alpha_e < 0 or rho_s < 0:
        raise ValueError('All input values must be non-negative.')

    Nr = Ac * fctm * (1 + alpha_e * rho_s)  # Nr in N, return in kN
    return Nr / 1000  # Convert N to kN


def eps_crack(
    Es: float,
    beta_TS: float,
    rho_s: float,
    fctm: float,
) -> float:
    """Calculate the mean strain during the crack formation stage.

    fib Model Code 2020, eq. (30.5-2)

    Args:
        Es (float): Modulus of elasticity of steel in MPa.
        beta_TS (float): Coefficient to assess the mean strain depending on the
            type of loding.
        rho_s (float): Reinforcement ratio (As/Ac).
        f_ctm (float): Mean tensile strength of concrete in MPa.

    Returns:
        float: Mean strain (dimensionless).

    Raises:
        ValueError: If any input is negative.
    """
    if Es < 0 or rho_s < 0 or fctm < 0:
        raise ValueError('All input values must be non-negative.')

    return (1 - beta_TS) * fctm / (Es * rho_s)


def k1_r(h: float, d: float, x: float) -> float:
    """Calculate the factor k1/r which accounts for the increase of crack width
    with cover due to curvature in bending.

    fib Model Code 2020, eq. (30.5-4)

    Args:
        h (float): Height of the section in mm.
        d (float): Effective depth of the section in mm.
        x (float): Depth of the neutral axis of the cracked section in mm.

    Returns:
        float: The factor k1/r.
    """
    if h < 0 or d < 0 or x < 0:
        raise ValueError('All input values must be non-negative.')

    return (h - x) / (d - x)


def wcal(
    k1_r: float,
    sr_max: float,
    eps_sm_eps_cm: float,
) -> float:
    """Calculate the design (or calculated) crack width wcal.

    fib Model Code 2020, eq. (30.5-3)

    Args:
        k1_r (float): Factor accounting for curvature in bending.
        sr_max (float): Maximum crack spacing in mm.
        eps_sm_eps_cm (float): relavive mean strain.

    Returns:
        float: The calculated crack width wcal in mm.

    Raises:
        ValueError: If any input is negative.
    """
    if k1_r < 0 or sr_max < 0 or eps_sm_eps_cm < 0:
        raise ValueError('All input values must be non-negative.')

    return k1_r * sr_max * eps_sm_eps_cm


def kfl(h: float, xg: float, hc_ef: float) -> float:
    """Calculate the factor kfl which accounts for stress distributions before
    cracking.

    fib Model Code 2020, eq. (30.5-6)

    Args:
        h (float): Height of the section in mm.
        xg (float): Depth of the neutral axis before cracking in mm.
        hc_ef (float): Effective height according to the relevant section in
            mm.

    Returns:
        float: The factor kfl.

    Raises:
        ValueError: If any input is negative or if h <= xg.
    """
    if h < 0 or xg < 0 or hc_ef < 0:
        raise ValueError('All input values must be non-negative.')

    return 1 / 2 * (1 + (h - xg - hc_ef) / (h - xg))


def sr_max(
    cracking_stage: Literal['stabilized', 'formation'],
    c: float,
    kfl: float,
    tau_bms: float,
    fctm: float,
    rho_s_eff: float,
    casting_condition: Literal['good', 'poor'],
    kc: float = 1.50,
    k_phi_rho: float = 0.25,
) -> float:
    """Calculate the maximum crack spacing sr,max using the advanced formula.

    fib Model Code 2020, eq. (30.5-5)

    Args:
        cracking_stage (str): Cracking stage, either 'stabilized' or
            'formation'.
        c (float): Maximum of the vertical and horizontal clear concrete
            covers of the superficial bar in mm.
        kfl (float): Factor to account for stress distributions before
            cracking.
        tau_bms (float): Mean bond stress in MPa.
        fctm (float): Concrete tensile stress in MPa.
        rho_s_eff (float): Effective reinforcement ratio.
        casting_condition (str): Casting condition, either 'good' or 'poor'.
        kc (float): Empirical parameter to quantify the influence of the
            concrete cover. Defaults to 1.50.
        k_phi_rho (float): Parameter quantifying the influence of bond tau_bms.
            Defaults to 0.25.

    Returns:
        float: The maximum crack spacing sr,max in mm.

    Raises:
        ValueError: If any input is negative.
    """
    if any(
        value < 0
        for value in [kc, c, k_phi_rho, kfl, tau_bms, fctm, rho_s_eff]
    ):
        raise ValueError('All input values must be non-negative.')

    beta_w = 1.7 if cracking_stage == 'stabilized' else 2.0
    kb = 0.9 if casting_condition == 'good' else 1.2

    return beta_w * (
        kc * c + k_phi_rho * kfl * kb * (fctm / (tau_bms * rho_s_eff))
    )


def Ac_ef_bar(rx: float, ry: float, phi: float, h: float, x: float) -> float:
    """Calculate the effective area of concrete in tension around a single bar.

    fib Model Code 2020, eq. (30.5-9)

    Args:
        rx (float): Cover to the center of the bar in the x direction in mm.
        ry (float): Cover to the center of the bar in the y direction in mm.
        phi (float): Diameter of the bar in mm.
        h (float): Height of the section in mm.
        x (float): Depth of the neutral axis in mm.

    Returns:
        float: The effective area of concrete in tension around a single bar
        (Ac,ef,bar) in mm2.

    Raises:
        ValueError: If any input is negative.
    """
    if any(value < 0 for value in [rx, ry, phi, h, x]):
        raise ValueError('All input values must be non-negative.')

    bc_ef = min(rx + 5 * phi, 10 * phi, 3.5 * rx)
    hc_ef = min(ry + 5 * phi, 10 * phi, 3.5 * ry, h - x)

    return bc_ef * hc_ef


def Ac_ef_group(
    rx: float,
    ry: float,
    phi: float,
    h: float,
    x: float,
    nl: int,
    sy: float,
    b: float,
) -> float:
    """Calculate the effective area of concrete in tension around a group of
    bars.

    fib Model Code 2020, eq. (30.5-10)

    Args:
        rx (float): Cover to the center of the bar in the x direction in mm.
        ry (float): Cover to the center of the bar in the y direction in mm.
        phi (float): Diameter of the bar in mm.
        h (float): Height of the section in mm.
        x (float): Depth of the neutral axis in mm.
        nl (int): Number of reinforcement layers.
        sy (float): Spacing in the y direction between bars in the tensile zone
            in mm.
        b (float): Width of the section in mm.

    Returns:
        float: The effective area of concrete in tension around a group of bars
        (Ac,ef,group) in mm2.

    Raises:
        ValueError: If any input is negative.
    """
    if any(value < 0 for value in [rx, ry, phi, h, x, sy, b]) or nl <= 0:
        raise ValueError(
            'All input values must be non-negative, and nl must be positive.'
        )

    hc_eff = min(ry + 5 * phi, 10 * phi, 3.5 * ry) + (nl - 1) * sy
    hc_ef = min(hc_eff, h - x)
    bc_ef = b

    return bc_ef * hc_ef


def rho_s_ef(As: float, Ac_ef: float) -> float:
    """Calculate the effective reinforcement ratio (rho_s,ef).

    fib Model Code 2020, eq. (30.5-8)

    Args:
        As (float): Total area of the bars in mm.
        Ac_ef (float): Effective area of concrete in tension around the bars in
            mm2.

    Returns:
        float: The effective reinforcement ratio (rho_s,ef).

    Raises:
        ValueError: If any input is negative.
    """
    if As < 0 or Ac_ef < 0:
        raise ValueError('All input values must be non-negative.')

    return As / Ac_ef


def sr_max_frc(
    cracking_stage: Literal['stabilized', 'formation'], srm: float
) -> float:
    """Computes the maximum crack spacing in FRC.

    fib Model Code 2020, eq. (30.5-##)

    Args:
        cracking_stage (str): Cracking stage, either 'stabilized' or
            'formation'.

    """
    if srm < 0:
        raise ValueError('All input values must be non-negative.')

    beta_w = 1.7 if cracking_stage == 'stabilized' else 2.0
    return beta_w * srm


def eps_sm_eps_cm(
    sigma_s: float, sigma_sr_ef: float, Es: float, beta_TS: float
) -> float:
    """Calculate the relative mean strain for an element subjected to direct
    loads or imposed strains where end restraint dominates.

    fib Model Code 2020, eq. (30.5-11)

    Args:
        sigma_s (float): Steel stress in the crack in MPa.
        sigma_sr_ef (float): Steel stress in a crack in the crack formation
            stage in MPa.
        Es (float): Modulus of elasticity of steel in MPa.
        beta_TS (float): Empirical coefficient from Table 30.5-2.

    Returns:
        float: Relative mean strain.

    Raises:
        ValueError: If any input is negative.
    """
    if sigma_s < 0:
        raise ValueError(f'sigma_s must not be negative. Got {sigma_s}')
    if sigma_sr_ef < 0:
        raise ValueError(
            f'sigma_sr_ef must not be negative. Got {sigma_sr_ef}'
        )
    if Es <= 0:
        raise ValueError(f'Es must be positive. Got {Es}')
    if beta_TS < 0 or beta_TS > 1:
        raise ValueError(f'beta_TS must be between 0 and 1. Got {beta_TS}')

    epsilon_sm_minus_cm = (sigma_s - beta_TS * sigma_sr_ef) / Es
    minimum_strain = sigma_s / Es * (1 - beta_TS)

    return max(epsilon_sm_minus_cm, minimum_strain)


def eps_sm_eps_cm_restrained(
    Rax: float, eps_free: float, beta_TS: float, fct_eff: float, Ec: float
) -> float:
    """Calculate the relative mean strain for an element subjected to
    restrained imposed strains and restrained at the edges.

    fib Model Code 2020, eq. (30.5-12)

    Args:
        Rax (float): Restraint factor.
        epsilon_free (float): Imposed strain which develops after the
            construction stage.
        beta_TS (float): Empirical coefficient from Table 30.5-2.
        fct_eff (float): Effective tensile strength of concrete in MPa.
        Ec (float): Modulus of elasticity of concrete in MPa.

    Returns:
        float: Relative mean strain.

    Raises:
        ValueError: If any input is negative.
    """
    if Rax < 0 or Rax > 1:
        raise ValueError(f'Rax must be between 0 and 1. Got {Rax}')
    if eps_free < 0:
        raise ValueError(f'epsilon_free must not be negative. Got {eps_free}')
    if beta_TS < 0 or beta_TS > 1:
        raise ValueError(f'beta_TS must be between 0 and 1. Got {beta_TS}')
    if fct_eff < 0:
        raise ValueError(f'fct_eff must not be negative. Got {fct_eff}')
    if Ec <= 0:
        raise ValueError(f'Ec must be positive. Got {Ec}')

    return Rax * eps_free - beta_TS * fct_eff / Ec


def sigma_sr_ef(fct_eff: float, rho_s_ef: float, alpha_e: float) -> float:
    """Calculate the steel stress in a crack during the crack formation stage.

    fib Model Code 2020, eq. (30.5-13)

    Args:
        fct_eff (float): Effective tensile strength of concrete in MPa.
        rho_s_ef (float): Effective reinforcement ratio.
        alpha_e (float): Modular ratio Es/Ec.

    Returns:
        float: Steel stress in the crack in MPa.

    Raises:
        ValueError: If any input is negative.
    """
    if fct_eff < 0:
        raise ValueError(f'fct_eff must not be negative. Got {fct_eff}')
    if rho_s_ef < 0:
        raise ValueError(f'rho_s_ef must not be negative. Got {rho_s_ef}')
    if alpha_e <= 0:
        raise ValueError(f'alpha_e must be positive. Got {alpha_e}')

    return fct_eff / rho_s_ef * (1 + alpha_e * rho_s_ef)


def Rax(eps_restr: float, eps_imp: float) -> float:
    """Calculate the restraint factor.

    fib Model Code 2020, eq. (30.5-14)

    Args:
        eps_restr (float): Strain which develops in the restrained element.
        eps_imp (float): Imposed strain (i.e., unrestrained shrinkage or
            temperature strain).

    Returns:
        float: Restraint factor.

    Raises:
        ValueError: If any input is negative or epsilon_imp is zero.
    """
    if eps_restr < 0:
        raise ValueError(
            f'epsilon_restr must not be negative. Got {eps_restr}'
        )
    if eps_imp <= 0:
        raise ValueError(f'epsilon_imp must be positive. Got {eps_imp}')

    return 1 - (eps_restr / eps_imp)


def tau_bms(
    fctm_t: float,
    load_type: Literal['short-term', 'long-term'],
    stage: Literal['crack_formation', 'stabilized_cracking'],
) -> float:
    """Calculate the bond stress tau_bms for deformed reinforcing bars based on
    the load type and stage.

    fib Model Code 2020, Table 30.5-2

    Args:
        fctm_t (float): Mean tensile strength of concrete at time t in MPa.
        load_type (str): Type of loading ('short-term' or 'long-term').
        stage (str): Cracking stage ('crack formation' or 'stabilized
            cracking').

    Returns:
        float: Bond stress tau_bms in MPa.

    Raises:
        ValueError: fctm_t is negative.
    """
    if fctm_t < 0:
        raise ValueError(f'fctm_t must not be negative. Got {fctm_t}')

    if load_type == 'short-term':
        return 1.8 * fctm_t
    if load_type == 'long-term' and stage == 'crack_formation':
        return 1.35 * fctm_t

    # if load_type == 'long-term' and stage == 'stabilized_cracking'
    return 1.8 * fctm_t


def beta_TS(
    load_type: Literal['short-term', 'long-term'],
    stage: Literal['crack_formation', 'stabilized_cracking'],
) -> float:
    """Determine the empirical coefficient beta_TS based on the load type and
    stage.

    fib Model Code 2020, Table 30.5-2

    Args:
        load_type (str): Type of loading ('short-term' or 'long-term').
        stage (str): Cracking stage ('crack_formation' or
            'stabilized_cracking').

    Returns:
        float: Empirical coefficient beta_TS.

    Raises:
        ValueError: If load_type or stage is invalid.
    """
    if load_type == 'short-term':
        return 0.6
    if load_type == 'long-term' and stage == 'crack_formation':
        return 0.6

    # if load_type == 'long-term' and stage == 'stabilized_cracking'
    return 0.4


def phi_eq(diameters: List[float], counts: List[int]) -> float:
    """Calculate the equivalent bar diameter phi_eq for a section with
    different bar diameters.

    fib Model Code 2020, eq. (30.5-15)

    Args:
        diameters (List[float]): List of bar diameters in mm.
        counts (List[int]): List of number of bars corresponding to each
            diameter.

    Returns:
        float: Equivalent bar diameter phi_eq in mm.

    Raises:
        ValueError: If the lists diameters and counts are not of the same
            length, or if any diameter or count is non-positive.
    """
    if len(diameters) != len(counts):
        raise ValueError(
            'diameters and counts must have the same '
            + f' length. Got {len(diameters)} and {len(counts)}'
        )
    if any(d <= 0 for d in diameters):
        raise ValueError(f'All diameters must be positive. Got {diameters}')
    if any(c <= 0 for c in counts):
        raise ValueError(f'All counts must be positive. Got {counts}')

    numerator = sum(n * d**2 for n, d in zip(counts, diameters))
    denominator = sum(n * d for n, d in zip(counts, diameters))

    if denominator == 0:
        raise ValueError(
            'Sum of counts multiplied by diameters must not be zero.'
        )

    return numerator / denominator


def sr_max_theta(sr_max_x: float, sr_max_y: float, theta: float) -> float:
    """Calculate the maximum crack spacing s_r,max,theta for orthogonally
    reinforced members.

    fib Model Code 2020, eq. (30.5-16)

    Args:
        sr_max_x (float): Characteristic crack spacing in the x direction in
            mm.
        sr_max_y (float): Characteristic crack spacing in the y direction in
            mm.
        theta (float): Angle between the reinforcement in the x direction and
            the direction of the principal tensile stress in degrees.

    Returns:
        float: Maximum crack spacing s_r,max,theta in mm.

    Raises:
        ValueError: If sr_max_x or sr_max_y are negative, or if theta is out of
            range.
    """
    if sr_max_x < 0:
        raise ValueError(f'sr_max_x must not be negative. Got {sr_max_x}')
    if sr_max_y < 0:
        raise ValueError(f'sr_max_y must not be negative. Got {sr_max_y}')
    if theta < 0 or theta > 90:
        raise ValueError(
            f'theta must be between 0 and 90 degrees. Got {theta}'
        )

    theta_rad = math.radians(theta)
    cos_theta = math.cos(theta_rad)
    sin_theta = math.sin(theta_rad)

    return 1 / (cos_theta / sr_max_x + sin_theta / sr_max_y)


def theta_reinf(
    sigma_x: float,
    sigma_y: float,
    tau_xy: float,
    rho_s_ef_x: float,
    rho_s_ef_y: float,
) -> float:
    """Solve for the angle theta between the reinforcement in the x direction
    and the direction of the principal tensile stress.

    fib Model Code 2020, eq. (30.5-17)

    Args:
        sigma_x (float): Normal stress in x direction in MPa.
        sigma_y (float): Normal stress in y direction in MPa.
        tau_xy (float): Shear stress in MPa.
        rho_s_ef_x (float): Effective reinforcement ratio in the x direction.
        rho_s_ef_y (float): Effective reinforcement ratio in the y direction.

    Returns:
        float: Angle theta in degrees.

    Raises:
        ValueError: If any input is negative.
    """
    if sigma_x < 0 or sigma_y < 0 or tau_xy < 0:
        raise ValueError(
            'Normal and shear stresses must not be negative. '
            + f'Got sigma_x={sigma_x}, sigma_y={sigma_y}, tau_xy={tau_xy}'
        )
    if rho_s_ef_x <= 0 or rho_s_ef_y <= 0:
        raise ValueError(
            'Reinforcement ratios must be positive. '
            + f'Got rho_s_ef_x={rho_s_ef_x}, rho_s_ef_y={rho_s_ef_y}'
        )

    def equation(tan_theta):
        term1 = tan_theta**4
        term2 = (sigma_x / tau_xy) * (tan_theta**3)
        term3 = -(sigma_y / tau_xy) * (rho_s_ef_x / rho_s_ef_y) * tan_theta
        term4 = -(rho_s_ef_x / rho_s_ef_y)
        return term1 + term2 + term3 + term4

    tan_theta_initial_guess = 1.0
    tan_theta_solution = newton(equation, tan_theta_initial_guess)
    theta_rad = math.atan(tan_theta_solution)

    return math.degrees(theta_rad)


def eps_sm_x_eps_cm_x(
    sigma_s_x: float, sigma_sr_ef_x: float, Es: float, beta_TS_x: float
) -> float:
    """Calculate the relative mean strain (eps_sm,x - eps_cm,x) in the x
    direction.

    fib Model Code 2020, eq. (30.5-20)

    Args:
        sigma_s_x (float): Reinforcing steel stress at the crack in x direction
            in MPa.
        sigma_sr_ef_x (float): Steel stress at the crack formation stage in x
            direction in MPa.
        Es (float): Modulus of elasticity of steel MPa.
        beta_TS_x (float): Empirical coefficient for x direction.

    Returns:
        float: Relative mean strain (eps_sm,x - eps_cm,x).

    Raises:
        ValueError: If any input is invalid.
    """
    if (
        sigma_s_x < 0
        or sigma_sr_ef_x < 0
        or Es <= 0
        or not (0 <= beta_TS_x <= 1)
    ):
        raise ValueError(
            'Invalid input values for relative '
            + 'mean strain calculation in x direction.'
        )

    return max(
        (sigma_s_x - beta_TS_x * sigma_sr_ef_x) / Es,
        sigma_s_x / Es * (1 - beta_TS_x),
    )


def eps_sm_y_eps_cm_y(
    sigma_s_y: float, sigma_sr_ef_y: float, Es: float, beta_TS_y: float
) -> float:
    """Calculate the relative mean strain (eps_sm,y - eps_cm,y) in the y
    direction.

    fib Model Code 2020, eq. (30.5-21)

    Args:
        sigma_s_y (float): Reinforcing steel stress at the crack in y direction
            in MPa.
        sigma_sr_ef_y (float): Steel stress at the crack formation stage in y
            direction in MPa.
        Es (float): Modulus of elasticity of steel in MPa.
        beta_TS_y (float): Empirical coefficient for y direction.

    Returns:
        float: Relative mean strain (eps_sm,y - eps_cm,y).

    Raises:
        ValueError: If any input is invalid.
    """
    if (
        sigma_s_y < 0
        or sigma_sr_ef_y < 0
        or Es <= 0
        or not (0 <= beta_TS_y <= 1)
    ):
        raise ValueError(
            'Invalid input values for relative mean '
            + 'strain calculation in y direction.'
        )

    return max(
        (sigma_s_y - beta_TS_y * sigma_sr_ef_y) / Es,
        sigma_s_y / Es * (1 - beta_TS_y),
    )


def eps_2(tau_xy: float, Ec: float, theta: float) -> float:
    """Calculate the principal compressive strain eps_2.

    fib Model Code 2020, eq. (30.5-19)

    Args:
        tau_xy (float): Shear stress in MPa.
        Ec (float): Modulus of elasticity of concrete in MPa.
        theta (float): Angle between the x direction reinforcement and the
            direction of principal tensile stress in degrees.

    Returns:
        float: Principal compressive strain ε2.

    Raises:
        ValueError: If any input is invalid.
    """
    if tau_xy < 0 or Ec <= 0 or not (0 <= theta <= 90):
        raise ValueError(
            'Invalid input values for principal '
            + 'compressive strain calculation.'
        )

    theta_rad = math.radians(theta)
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)

    return -abs(tau_xy) / (Ec * sin_theta * cos_theta)


def sigma_s_x(
    sigma_x: float, tau_xy: float, rho_s_ef_x: float, theta: float
) -> float:
    """Calculate the reinforcing steel stress sigma_s,x at the crack in the x
    direction.

    fib Model Code 2020, eq. (30.5-21)

    Args:
        sigma_x (float): Normal stress in x direction in MPa.
        tau_xy (float): Shear stress in MPa.
        rho_s_ef_x (float): Effective reinforcement ratio in the y direction.
        theta (float): Angle between the y direction reinforcement and the
            direction of principal tensile stress in degrees.

    Returns:
        float: Steel stress sigma_s,y in the x direction in MPa.

    Raises:
        ValueError: If any input is invalid.
    """
    if sigma_x < 0 or tau_xy < 0 or rho_s_ef_x <= 0 or not (0 <= theta <= 90):
        raise ValueError(
            'Invalid input values for steel stress calculation in y direction.'
        )

    theta_rad = math.radians(theta)
    return (sigma_x + tau_xy * math.tan(theta_rad)) / rho_s_ef_x


def sigma_s_y(
    sigma_y: float, tau_xy: float, rho_s_ef_y: float, theta: float
) -> float:
    """Calculate the reinforcing steel stress sigma_s,y at the crack in the y
    direction.

    fib Model Code 2020, eq. (30.5-21)

    Args:
        sigma_y (float): Normal stress in y direction in MPa.
        tau_xy (float): Shear stress in MPa.
        rho_s_ef_y (float): Effective reinforcement ratio in the y direction.
        theta (float): Angle between the y direction reinforcement and the
            direction of principal tensile stress in degrees.

    Returns:
        float: Steel stress sigma_s,y in the y direction in MPa.

    Raises:
        ValueError: If any input is invalid
    """
    if sigma_y < 0 or tau_xy < 0 or rho_s_ef_y <= 0 or not (0 <= theta <= 90):
        raise ValueError(
            'Invalid input values for steel stress calculation in y direction.'
        )

    theta_rad = math.radians(theta)
    return (sigma_y + tau_xy / math.tan(theta_rad)) / rho_s_ef_y


def beta_TS_x(
    beta_TS: float, sr_max_theta: float, theta: float, sr_max_x: float
) -> float:
    """Calculate the empirical coefficient beta_TS for the x direction.

    fib Model Code 2020, eqs. (30.5-20) and (30.5-21)

    Args:
        beta_TS (float): Base beta_TS coefficient.
        sr_max_theta (float): Maximum crack spacing in the theta direction in
            mm.
        theta (float): Angle between the direction of reinforcement and the
            direction of principal tensile stress in degrees.
        sr_max_x (float): Characteristic crack spacing in the given direction x
            in mm.

    Returns:
        float: Empirical coefficient beta_TS for the specified direction.

    Raises:
        ValueError: If any input is invalid.
    """
    if sr_max_theta <= 0 or sr_max_x <= 0 or not (0 <= theta <= 90):
        raise ValueError(
            'Invalid input values for empirical coefficient calculation.'
        )

    theta_rad = math.radians(theta)

    return beta_TS * (sr_max_theta / math.sin(theta_rad)) / sr_max_x


def beta_TS_y(
    beta_TS: float, sr_max_theta: float, theta: float, sr_max_y: float
) -> float:
    """Calculate the empirical coefficient beta_TS for the y direction.

    fib Model Code 2020, eqs. (30.5-20) and (30.5-21)

    Args:
        beta_TS (float): Base beta_TS coefficient.
        sr_max_theta (float): Maximum crack spacing in the theta direction in
            mm.
        theta (float): Angle between the direction of reinforcement and the
            direction of principal tensile stress in degrees.
        sr_max_y (float): Characteristic crack spacing in the given direction y
            in mm.

    Returns:
        float: Empirical coefficient beta_TS for the specified direction.

    Raises:
        ValueError: If any input is invalid.
    """
    if sr_max_theta < 0 or sr_max_y <= 0 or not (0 <= theta <= 90):
        raise ValueError(
            'Invalid input values for empirical coefficient calculation.'
        )

    theta_rad = math.radians(theta)

    return beta_TS * (sr_max_theta / math.cos(theta_rad)) / sr_max_y


def eps_sm_eps_cm_theta(
    eps_sm_x_eps_cm_x: float, eps_sm_y_eps_cm_y: float, eps_2: float
) -> float:
    """Computes the principal relative mean strain given the relative strain in
    the x-y directions.

    fib Model Code 2020, eqs. (30.5-18)

    Args:
        eps_sm_x_eps_cm_x (float): Relative strain in the x direction.
        eps_sm_y_eps_cm_y (float): Relative strain in the y direction.

    Returns:
        float: The principal mean strain.

    Raises:
        ValueError: If any input is invalid
    """
    if (
        eps_sm_x_eps_cm_x < 0
        or eps_sm_y_eps_cm_y < 0
        or eps_sm_y_eps_cm_y < eps_2
    ):
        raise ValueError(
            'Invalid input values for the relative strain calculation.'
        )
    return eps_sm_x_eps_cm_x - eps_sm_y_eps_cm_y - eps_2


def rho_s_p_ef(rho_s_ef: float, rho_p_ef: float, xi_1: float) -> float:
    """Calculate the equivalent reinforcement ratio rho_s+p,ef for combined
    reinforcement.

    fib Model Code 2020, eq. (30.5-22)

    Args:
        rho_s_ef (float): Effective reinforcement ratio for reinforcing steel.
        rho_p_ef (float): Effective reinforcement ratio for prestressing steel.
        xi_1 (float): Bond factor.

    Returns:
        float: Equivalent reinforcement ratio rho_s+p,ef.

    Raises:
        ValueError: If any input is invalid.
    """
    if rho_s_ef < 0 or rho_p_ef < 0 or xi_1 < 0:
        raise ValueError(
            'Invalid input values for equivalent reinforcement '
            + 'ratio calculation.'
        )

    return rho_s_ef + (xi_1**2) * rho_p_ef


def xi_1(tau_bmp_tau_bms: float, phi: float, phi_p_eq: float) -> float:
    """Calculate the bond factor xi_1 for prestressing steel.

    fib Model Code 2020, eq. (30.5-23)

    Args:
        tau_bmp_tau_bms (float): Ratio between tau_bmp and tau_bms.
        phi (float): Diameter of reinforcing steel in mm.
        phi_p_eq (float): Equivalent diameter of prestressing steel in mm.

    Returns:
        float: Bond factor xi_1.

    Raises:
        ValueError: If any input is invalid.
    """
    if tau_bmp_tau_bms < 0 or phi < 0 or phi_p_eq < 0:
        raise ValueError('Invalid input values for bond factor calculation.')

    return math.sqrt(tau_bmp_tau_bms * phi / phi_p_eq)


def phi_p_eq(Ap: float, up: float) -> float:
    """Calculate the equivalent diameter phi_p,eq of the prestressing steel.

    fib Model Code 2020, eq. (30.5-23)

    Args:
        Ap (float): Total area of prestressing steel in mm2.
        up (float): Total equivalent perimeter of prestressing steel in mm.

    Returns:
        float: Equivalent diameter phi_p,eq in mm.

    Raises:
        ValueError: If any input is invalid.
    """
    if Ap < 0 or up < 0:
        raise ValueError(
            'Invalid input values for equivalent diameter calculation.'
        )

    return 4 * Ap / up


def up_i(
    type_: Literal['bundle', '7-wire', '3-wire'],
    Ap_i: Optional[float] = None,
    phi_wire: Optional[float] = None,
) -> float:
    """Calculate the equivalent perimeter up,i of the prestressing steel.

    fib Model Code 2020, eq. for perimeter calculation

    Args:
        type_ (str): Type of prestressing steel ('bundle', '7-wire', '3-wire').
        Ap_i (float, optional): Area of the prestressing steel bundle or strand
            in mm2. Required for 'bundle'.
        phi_wire (float, optional): Diameter of a single wire in the strand in
            mm, required for '7-wire' and '3-wire' types.

    Returns:
        float: Equivalent perimeter up,i in mm.

    Raises:
        ValueError: If any input is invalid.
    """
    if type_ == 'bundle':
        if not Ap_i or Ap_i < 0:
            raise ValueError(f'Not valid value for Ap_i. Got {Ap_i}.')
        return 1.6 * math.pi * math.sqrt(Ap_i)

    if not phi_wire or phi_wire < 0:
        raise ValueError(f'Not valid value for phi_wire. Got {phi_wire}.')

    if type_ == '7-wire':
        return 1.75 * math.pi * phi_wire

    # If 3-wire strands
    return 1.20 * math.pi * phi_wire


def tau_bmp_tau_bms(
    surface_condition: Literal[
        'no bond', 'smooth wire', 'strand', 'indented wire', 'ribbed bar'
    ],
    tensioning_type: Literal['pretensioned', 'post-tensioned', 'no-bond'],
) -> float:
    """Retrieve the bond factor τbmp/τbms for different types of prestressing
    steel based on the surface condition and tensioning type.

    fib Model Code 2020, Table 30.5-3

    Args:
        surface_condition (str): Surface condition of the prestressing steel
            ('no bond', 'smooth wire', 'strand', 'indented wire',
            'ribbed bar').
        tensioning_type (str): Type of tensioning ('pretensioned' or
            'post-tensioned').

    Returns:
        float: Bond factor tau_bmp/tau_bms.

    Raises:
        ValueError: If the surface_condition or tensioning_type is invalid or
            if the combination is not supported.
    """
    # Define the bond factors from Table 30.5-3
    if tensioning_type == 'no-bond':
        return 0.0

    table = {
        'pretensioned': {
            'no bond': 0.0,
            'smooth wire': 0.4,
            'strand': 0.6,
            'indented wire': 0.8,
            'ribbed bar': None,  # Not applicable for pretensioned
        },
        'post-tensioned': {
            'no bond': 0.0,
            'smooth wire': 0.2,
            'strand': 0.4,
            'indented wire': 0.6,
            'ribbed bar': 1.0,
        },
    }

    # Validate input
    if tensioning_type not in table:
        raise ValueError(
            'Invalid tensioning type. Expected "pretensioned" or '
            + f'"post-tensioned", got {tensioning_type}.'
        )
    if surface_condition not in table[tensioning_type]:
        raise ValueError(
            f'Invalid surface condition. Got {surface_condition}.'
        )

    # Retrieve bond factor
    bond_factor = table[tensioning_type][surface_condition]
    if bond_factor is None:
        raise ValueError(
            f'The combination of surface condition "{surface_condition}" '
            + f'and tensioning type "{tensioning_type}" is not applicable.'
        )

    return bond_factor


def As_min_frc(
    fctm: float,
    fFtsm: float,
    Act: float,
    sigma_s: float,
    h: float,
    kc: float = 1.0,
) -> float:
    """Calculate the minimum reinforcement area As,min for crack control in
    FRC.

    fib Model Code 2020, eq. (30.5-##)

    Args:
        fctm (float): Average value of the tensile strength of the concrete
            matrix in MPa.
        fFtsm (float): Average value of the residual strength of FRC in MPa.
        Act (float): Tensile part of the concrete cross-section in mm2.
        sigma_s (float): Maximum tensile stress in the reinforcement in the
            cracked state in MPa.
        h (float): thickness of the flange or the web in mm.
        kc (float): Coefficient accounting for the stress distribution in the
            cross-section. Defaults to 1.0.

    Returns:
        float: Minimum reinforcement area As,min mm2. If the calculated value
        is negative, returns 0.

    Raises:
        ValueError: If any input is invalid or negative.
    """
    if fctm < 0 or fFtsm < 0 or Act < 0 or sigma_s < 0 or kc < 0 or h <= 0:
        raise ValueError(
            'Invalid input values for minimum reinforcement calculation.'
        )

    if h <= 300:
        k = 1.0
    elif h >= 800:
        k = 0.65
    else:
        k = 1.0 - ((h - 300) / (800 - 300)) * (1.0 - 0.65)

    As_min = Act * kc * k * (fctm - fFtsm) / sigma_s

    # If As_min is negative, return 0 as the reinforcement can be
    # fully provided by FRC
    return max(As_min, 0.0)


def max_phi_crack(
    rho_s_ef: float,
    ry: float,
    d: float,
    w_lim_cal: float,
    beta_w: float,
    kfl_simp: float,
    kb: float,
    k1_r_simpl: float,
    sigma_s: float,
    Es: float,
    c: float,
) -> float:
    """Calculate the maximum allowable bar diameter phi for simplified crack
    calculation.

    fib Model Code 2020, eq. (30.5-24)

    Args:
        rho_s_ef (float): Effective reinforcement ratio.
        ry (float): Radius of the reinforcement in mm.
        d (float): Effective depth of the cross-section in mm.
        w_lim_cal (float): Design crack width limit in mm.
        beta_w (float): Coefficient related to crack width.
        kfl_simpl (float): Coefficient related to the stress state of the
            section.
        kb (float): Coefficient accounting for bond conditions.
        k1_r_simpl (float): Simplified coefficient for stress redistribution.
        sigma_s (float): Tensile stress in the steel reinforcement in MPa.
        Es (float): Modulus of elasticity of steel in MPa.
        c (float): Concrete cover in mm.

    Returns:
        float: Maximum allowable bar diameter φ in mm.

    Raises:
        ValueError: If any input is invalid or negative.
    """
    if any(
        val < 0
        for val in [
            rho_s_ef,
            ry,
            d,
            w_lim_cal,
            beta_w,
            kfl_simp,
            kb,
            k1_r_simpl,
            sigma_s,
            Es,
            c,
        ]
    ):
        raise ValueError('All input values must be positive.')

    term1 = 2.1 * rho_s_ef / ((ry / d) * kfl_simp * kb)
    term2 = w_lim_cal / (beta_w * k1_r_simpl * 0.9 * sigma_s / Es) - 1.5 * c

    return term1 * term2


def max_st_crack(
    rho_s_ef: float,
    ry: float,
    d: float,
    w_lim_cal: float,
    beta_w: float,
    kfl_simpl: float,
    kb: float,
    k1_r_simpl: float,
    sigma_s: float,
    Es: float,
    c: float,
) -> float:
    """Calculate the maximum allowable bar spacing sl for simplified crack
    calculation.

    fib Model Code 2020, eq. (30.5-25)

    Args:
        rho_s_ef (float): Effective reinforcement ratio.
        ry (float): Radius of the reinforcement in mm.
        d (float): Effective depth of the cross-section in mm.
        w_lim_cal (float): Design crack width limit in mm.
        beta_w (float): Coefficient related to crack width.
        kfl_simpl (float): Coefficient related to the stress state of the
            section.
        kb (float): Coefficient accounting for bond conditions.
        k1_r_simpl (float): Simplified coefficient for stress redistribution.
        sigma_s (float): Tensile stress in the steel reinforcement in MPa.
        Es (float): Modulus of elasticity of steel in MPa.
        c (float): Concrete cover in mm.

    Returns:
        float: Maximum allowable bar spacing sl in mm.

    Raises:
        ValueError: If any input is invalid or negative.
    """
    if any(
        val < 0
        for val in [
            rho_s_ef,
            ry,
            d,
            w_lim_cal,
            beta_w,
            kfl_simpl,
            kb,
            k1_r_simpl,
            sigma_s,
            Es,
            c,
        ]
    ):
        raise ValueError('All input values must be positive.')

    term1 = 3.45 * rho_s_ef / ((ry**2 / d) * kfl_simpl**2 * kb**2)
    term2 = w_lim_cal / (beta_w * k1_r_simpl * 0.9 * sigma_s / Es) - 1.5 * c

    return term1 * (term2**2)


def kfl_simpl(
    ry: float, h: float, tension: Literal['one-face', 'two-faces']
) -> float:
    """Calculate the coefficient kfl_simpl for the stress state of the section.

    fib Model Code 2020, related to Eq. (30.5-24)

    Args:
        ry (float): Radius of the reinforcement in mm.
        h (float): Total depth of the cross-section in mm.
        tension (str): Wether the element has two sides or both sides in
            tension.

    Returns:
        float: Coefficient kfl_simpl.

    Raises:
        ValueError: If any input is invalid or negative.
    """
    if tension == 'two-faces':
        return 1.0

    if ry < 0 or h < 0:
        raise ValueError(
            'Radius of reinforcement and depth of '
            + 'the cross-section must be positive.'
        )

    return min(1 - 3.5 * ry / h, 1)


def k1_r_simpl(
    rho_s_p_ef: float, h: float, d: float, state: Literal['bending', 'tension']
) -> float:
    """Calculate the simplified coefficient k1/r,simpl for stress
    redistribution.

    fib Model Code 2020, related to Eq. (30.5-24)

    Args:
        rho_s_p_ef (float): Effective reinforcement ratio for combined
            reinforcement.
        h (float): Total depth of the cross-section in mm.
        d (float): Effective depth of the cross-section in mm.
        state (str): 'bending' or 'tension'.

    Returns:
        float: Coefficient k1/r,simpl.

    Raises:
        ValueError: If any input is invalid or negative.
    """
    if state == 'tension':
        return 1.0

    if rho_s_p_ef < 0 or h < 0 or d < 0:
        raise ValueError(
            'Effective reinforcement ratio must be non-negative, '
            + 'and depths must be positive.'
        )

    ratio = h / d

    k1_r_simpl = 25 * (ratio - 1) * rho_s_p_ef + 1.15 * ratio - 0.15

    return min(k1_r_simpl, 1.0)
