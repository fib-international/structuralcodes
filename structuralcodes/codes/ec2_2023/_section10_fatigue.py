"""Functions from Section 10 of EN 1992-1-1:2023."""

import math
from typing import Literal, Optional


def Ae(Ap: float, phi: float, phi_p: float, xi: float) -> float:
    """Calculate the equivalent area of reinforcement.

    This function calculates the equivalent area of reinforcement
        for the prestressing tendon considering the ratio of bond stress ξ.

    EN1992-1-1:2023 Eq. (10.2)

    Args:
        Ap (float): Area of prestressing tendon or tendons in mm2.
        phi (float): Largest diameter of reinforcing steel in mm.
        phi_p (float): Diameter or equivalent diameter of prestressing steel
            in mm.
        xi (float): Ratio of bond strength between bonded tendons and
            ribbed or indented reinforcing steel in concrete.

    Returns:
        float: Equivalent area of reinforcement in mm2.

    Raises:
        ValueError: If any input is negative.
    """
    if Ap < 0:
        raise ValueError(f'Ap must not be negative. Got {Ap}')
    if phi < 0:
        raise ValueError(f'phi must not be negative. Got {phi}')
    if phi_p < 0:
        raise ValueError(f'phi_p must not be negative. Got {phi_p}')
    if xi < 0:
        raise ValueError(f'xi must not be negative. Got {xi}')

    return Ap * math.sqrt(xi * phi / phi_p)


def phi_p_eq(
    Ap: float,
    type: Literal['bundle', 'single_7_wire', 'single_3_wire'],
    phi_wire: Optional[float] = None,
) -> float:
    """Calculate the equivalent diameter of prestressing steel.

    This function calculates the equivalent diameter of prestressing
        steel based on the type of steel and optionally the wire
            diameter if applicable.

    EN1992-1-1:2023 10.3(2)

    Args:
        Ap (float): Area of prestressing tendon or tendons in mm2.
        type (str): Type of prestressing steel
            ('bundle', 'single_7_wire', 'single_3_wire').
        phi_wire (float, optional): Wire diameter in mm, required for
            'single_7_wire' and 'single_3_wire'.

    Returns:
        float: Equivalent diameter of prestressing steel in mm.

    Raises:
        ValueError: If any input is negative or if required inputs
            are missing for specific steel types.
    """
    if Ap < 0:
        raise ValueError(f'Ap must not be negative. Got {Ap}')
    if phi_wire is not None and phi_wire < 0:
        raise ValueError(f'phi_wire must not be negative. Got {phi_wire}')

    if type == 'bundle':
        phi_p_eq = 1.60 * math.sqrt(Ap)
    elif type == 'single_7_wire':
        if phi_wire is None:
            raise ValueError('phi_wire is required for single_7_wire')
        phi_p_eq = 1.75 * phi_wire
    elif type == 'single_3_wire':
        if phi_wire is None:
            raise ValueError('phi_wire is required for single_3_wire')
        phi_p_eq = 1.20 * phi_wire

    return phi_p_eq


def cot_theta_fat(cot_theta_uls: float) -> float:
    """Calculate the cotangent inclination of the compressive
        struts for fatigue verification.

    This function calculates the inclination of the compressive struts θfat
        using the compression field inclination θ at ULS.

    EN1992-1-1:2023 Eq. (10.3)

    Args:
        cot_theta_uls (float): Cotangent of the inclination
            of the compressive field at ultimate limit state (ULS) in degrees.

    Returns:
        float: Inclination of the compressive struts for
            fatigue verification in degrees.

    Raises:
        ValueError: If theta_uls is negative.
    """
    if cot_theta_uls < 0:
        raise ValueError(
            f'cot_theta_uls must not be negative. Got {cot_theta_uls}'
        )

    cot_theta_fat = math.sqrt(cot_theta_uls)
    return max(1, cot_theta_fat)


def xi_bond(  # noqa: PLR0911, PLR0912
    fck: float,
    state: Literal['pre', 'post'],
    prestressing_steel_type: Literal['smooth', 'strand', 'indented', 'ribbed'],
) -> float:
    """Determine the bond strength ratio ξ based on the
        compressive strength of the concrete and the type
        of prestressing steel.

    EN1992-1-1:2023 Table 10.1

    Args:
        fck (float): Compressive strength of concrete in MPa.
        state (str): pre or post tensioned tendons.
        prestressing_steel_type (str): Type of prestressing steel
            ('smooth', 'strand', 'indented', 'ribbed').

    Returns:
        float: Bond strength ratio ξ.

    Raises:
        ValueError: If fck is negative or if
            prestressing_steel_type is invalid.
    """
    fck = abs(fck)

    # Pre-tensioned strands
    if state == 'pre' and prestressing_steel_type == 'smooth':
        raise ValueError(
            'Smoothed bars/wires in pre-tensioned do notallow bond strength. '
        )

    if state == 'pre':
        data = {
            'strand': 0.60,
            'indented': 0.70,
            'ribbed': 0.80,
        }
        return data[prestressing_steel_type]

    # Post-tensioned strands
    if prestressing_steel_type == 'smooth':
        if fck <= 50:
            return 0.30
        if fck >= 70:
            return 0.15
        return 0.30 - 0.15 * (fck - 50) / 20  # Linear interpolation

    if prestressing_steel_type == 'strand':
        if fck <= 50:
            return 0.60
        if fck >= 70:
            return 0.25
        return 0.60 - 0.35 * (fck - 50) / 20  # Linear interpolation

    if prestressing_steel_type == 'indented':
        if fck <= 50:
            return 0.70
        if fck >= 70:
            return 0.30
        return 0.70 - 0.40 * (fck - 50) / 20  # Linear interpolation

    # Ribbed
    if fck <= 50:
        return 0.80
    if fck >= 70:
        return 0.35
    return 0.80 - 0.45 * (fck - 50) / 20  # Linear interpolation


def fcd_fat(
    fck: float, gamma_c: float, beta_cc: float, eta_cc: float
) -> float:
    """Calculate the design fatigue strength of concrete.

    EN1992-1-1:2023 Eq. (10.5)

    Args:
        fck (float): Characteristic compressive strength of concrete in MPa.
        gamma_c (float): Partial safety factor for concrete. Must be positive.
        beta_cc (float): Coefficient for concrete strength at first
            load application t0.
        eta_cc (float): Coefficient for long-term effects on the
            compressive strength of concrete.

    Returns:
        float: Design fatigue strength of concrete fcd,fat in MPa.

    Raises:
        ValueError: If any input is invalid.
    """
    if fck < 0:
        raise ValueError(f'fck must not be negative. Got {fck}')
    if gamma_c <= 0:
        raise ValueError(f'gamma_c must be positive. Got {gamma_c}')
    if beta_cc <= 0:
        raise ValueError(f'beta_cc must be positive. Got {beta_cc}')
    if eta_cc <= 0:
        raise ValueError(f'eta_cc must be positive. Got {eta_cc}')

    # Calculate the design fatigue strength of concrete
    eta_cc_fat = min(0.85 * eta_cc, 0.8)
    return beta_cc * (fck / gamma_c) * eta_cc_fat


def fatigue_compression_check(
    sigma_cd_max: float, sigma_cd_min: float, fcd_fat: float
) -> bool:
    """Check the fatigue resistance of concrete under compression.

    EN1992-1-1:2023 Eq. (10.4)

    Args:
        sigma_cd_max (float): Maximum compressive stress at a fibre under the
            fatigue load combination in MPa.
        sigma_cd_min (float): Minimum compressive
            stress at the same fibre in MPa.
        fcd_fat (float): Design fatigue strength
            of concrete in MPa.

    Returns:
        bool: True if the concrete complies with the fatigue
            resistance criteria, False otherwise.

    Raises:
        ValueError: If fcd_fat is invalid.
    """
    if fcd_fat < 0:
        raise ValueError(f'fcd_fat must not be negative. Got {fcd_fat}')

    # Check the condition for fatigue resistance
    left_side = abs(sigma_cd_max) / fcd_fat
    right_side = 0.5 + 0.45 * abs(sigma_cd_min) / fcd_fat

    return left_side <= right_side <= 0.9


def reinf_shear_fatigue_check(
    sigma_cd_max: float,
    sigma_cd_min: float,
    fcd_fat: float,
    nu: float = 0.5,
) -> bool:
    """Checks the fatigue resistance of concrete under
        shear in the strut model.

    EN1992-1-1:2023 (10.6(1))

    Args:
        sigma_cd_max (float): Maximum compressive stress at a fibre under the
            fatigue load combination in MPa.
        sigma_cd_min (float): Minimum compressive
            stress at the same fibre in MPa.
        fcd_fat (float): Design fatigue strength
            of concrete in MPa.
        nu (float): reduction coefficient for resistance in the strut.
            Defaults to 0.5.

    Returns:
        bool: True if the condition is satisfied, False otherwise.

    Raises:
        ValueError: If fcd_fat value is negative.
    """
    fcd_fat_red = nu * fcd_fat
    return fatigue_compression_check(sigma_cd_max, sigma_cd_min, fcd_fat_red)


def no_reinf_shear_fatigue_check(
    tau_ed_max: float,
    tau_ed_min: float,
    tau_rd_c: float,
) -> bool:
    """Checks the fatigue resistance of concrete under shear.

    EN1992-1-1:2023 Eq. (10.6) and Eq. (10.7)

    Args:
        tau_ed_max (float): Design shear stress due to the maximum
            applied shear force in MPa.
        tau_ed_min (float)): Design shear stress due to the minimum applied
            shear in MPa.
        tau_rd_c (float): Design value for shear resistance stress without
            shear reinforcement in MPa.

    Returns:
        bool: True if the condition is satisfied, False otherwise.

    Raises:
        ValueError: If any input value is negative.
    """
    if tau_rd_c <= 0:
        raise ValueError(f'tau_rd_c must be positive. Got {tau_rd_c}')

    if tau_ed_min / tau_ed_max >= 0:
        return (
            (abs(tau_ed_max) / tau_rd_c)
            <= (0.5 + 0.45 * abs(tau_ed_min) / tau_rd_c)
            <= 0.90
        )
    return (abs(tau_ed_max) / tau_rd_c) <= (0.5 - abs(tau_ed_min) / tau_rd_c)


def cv_1_fat(surface: Literal['rough', 'very rough', 'keyed']) -> float:
    """Returns the cv_1_fat coefficient for computing
    verification at shear interfaces.

    EN1992-1-1:2023 Table (10.2)

    Args:
        surface (str): The surface type.
    """
    if surface == 'rough':
        return 0.075
    if surface == 'very rough':
        return 0.095

    # Keyed
    return 0.185


def mu_v_fat(
    surface: Literal['rough', 'very rough', 'keyed'],
) -> float:
    """Returns mu_v_fat coefficient for computing
    verification at shear interfaces.

    EN1992-1-1:2023 Table (10.2)

    Args:
        surface (str): The surface type.
    """
    if surface == 'rough':
        return 0.7
    if surface == 'very rough':
        return 0.9

    # Keyed
    return 0.9


def tau_Rdi_fat_no_reinf(
    cv_1_fat: float,
    mu_v_fat: float,
    sigma_n: float,
    fck: float,
    gamma_c: float,
    rho_i: float,
    fyd: float,
    alpha: float,
) -> float:
    """Calculates the fatigue shear strength at an interface without
        reinforcement.

    EN1992-1-1:2023 Eq. (10.8)

    Args:
        cv_1_fat (float): Coefficient for the interface roughness
            according to Table 10.2 (unitless).
        mu_v_fat (float): Coefficient for the interface roughness with
            reinforcement according to Table 10.2 (unitless).
        sigma_n (float): Normal stress at the interface in MPa.
        fck (float): Characteristic compressive strength of concrete in MPa.
        gamma_c (float): Partial safety factor for concrete (unitless).
        rho_i (float): Ratio of reinforcement (unitless).
        fyd (float): Design racteristic yield strength of reinforcement in MPa.
        alpha (float): Angle of inclination (degrees).

    Returns:
        float: Shear strength in MPa.

    Raises:
        ValueError: If any input value is negative or
            zero where it shouldn't be.
    """
    # Validate inputs
    if (
        cv_1_fat < 0
        or mu_v_fat < 0
        or sigma_n < 0
        or fck < 0
        or gamma_c <= 0
        or rho_i < 0
        or fyd <= 0
    ):
        raise ValueError(
            'Input values must be non-negative '
            + ' and gamma_c, f_yk must be positive.'
        )

    # Convert angle to radians
    alpha_rad = math.radians(alpha)

    # Calculate shear strength according to Eurocode formula
    tau_r_di = cv_1_fat * math.sqrt(fck) / gamma_c + mu_v_fat * abs(sigma_n)
    return min(
        tau_r_di, 0.3 * fck / gamma_c + rho_i * fyd * math.cos(alpha_rad)
    )


def delta_tau_Rdi_fat_reinf(
    mu_v_fat: float,
    sigma_n: float,
    rho: float,
    sigma_Rsk: float,
    gamma_s: float,
    alpha: float,
) -> float:
    """Calculates the shear strength at an interface with reinforcement.

    EN1992-1-1:2023 Eq. (10.9)

    Args:
        mu_v_fat (float): Coefficient for the interface roughness
            with reinforcement according to Table 10.2 (unitless).
        sigma_n (float): Normal stress at the interface in MPa.
        rho (float): Ratio of reinforcement (unitless).
        sigma_Rsk (float): value taken from Table E.1.
        gamma_s (float): Partial safety factor for reinforcement (unitless).
        alpha (float): Angle of inclination (degrees).

    Returns:
        float: Shear strength in MPa.

    Raises:
        ValueError: If any input value is negative or
            zero where it shouldn't be.
    """
    # Validate inputs
    if (
        mu_v_fat < 0
        or sigma_n < 0
        or rho < 0
        or sigma_Rsk <= 0
        or gamma_s <= 0
    ):
        raise ValueError(
            'Input values must be non-negative and '
            + ' f_s, gamma_s must be positive.'
        )

    # Convert angle to radians
    alpha_rad = math.radians(alpha)

    # Calculate shear strength according to Eurocode formula
    return mu_v_fat * abs(sigma_n) + rho * sigma_Rsk / (0.45 * gamma_s) * (
        mu_v_fat * math.sin(alpha_rad) + math.cos(alpha_rad)
    )
