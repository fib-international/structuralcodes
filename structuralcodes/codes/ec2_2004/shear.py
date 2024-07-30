"""Design rules according to EN 1992-1-1 regarding shear."""

import math

# General functions


# Part of Equation (6.2).
def _k(d: float) -> float:
    """Compute a correction factor.

    Defined in EN 1992-1-1 (2005), Eq. (6.2).

    Args:
        d (float): The effective depth of the cross-section in mm.

    Returns:
        float: Correction factor to account for the cross-sectional size on the
        shear resistance.
    """
    return min(1.0 + math.sqrt(200.0 / d), 2.0)


# Part of Equation (6.2).
def _rho_L(Asl: float, bw: float, d: float) -> float:
    """Compute the longitudinal reinforcement ratio.

    Defined in EN 1992-1-1 (2005), Eq. (6.2).

    Args:
        Asl (float): The cross-sectional area of the tensile reinforcement,
            anchored at least (lbd + d) beyond the considered cross-section, in
            mm2.
        bw (float): The smallest width of the cross-section in tension in mm.
        d (float): The effective depth of the cross-section in mm.

    Returns:
        float: The maximum allowable reinforcement ratio of the longitudinal
        reinforcement, unitless.
    """
    return min(Asl / (bw * d), 0.02)


# Part of Equation(6.2).
def _sigma_cp(NEd: float, Ac: float, fcd: float) -> float:
    """Calculate the average prestress stress in the cross-section.

    Defined in EN 1992-1-1 (2005), Eq. (6.2).

    Args:
        NEd (float): The normal force in the cross-section due to loading or
            prestress (NEd > 0 for compression) in N.
        Ac (float): The cross-sectional area of the concrete in mm2.
        fcd (float): The design compressive strength in MPa.

    Returns:
        float: The maximum allowable average prestress in the cross-section in
        MPa.
    """
    return min(NEd / Ac, 0.2 * fcd)


# Part of Equation (6.4)
def _alpha_l(L_x, L_pt2) -> float:
    """Compute the relative anchorage length for prestreched prestressing
    steel.

    Defined in EN 1992-1-1 (2005), Eq. (6.4).

    Args:
        L_x (float): Distance from the considered cross-section until the
            starting point of the transference length of the prestress steel.
        L_pt2 (float): Maximum value of the transference length of the
            prestress steel, according to Eq. (8.18).

    Returns:
        float: Fraction (relative anchorage length) for determining the amount
        of prestress that may be used when determining the shear resistance
        using Mohr's circle.
    """
    return min(L_x / L_pt2, 1.0)


# Equation (6.7N)
def _theta(theta: float, cot_min: float = 1.0, cot_max: float = 2.5) -> None:
    """Check if the provided angle theta is within the bounds provided by the
    code.

    EN 1992-1-1 (2005). Eq. (6.7N)

    Args:
        theta (float): The chosen angle of the compression strut in degrees.

    Keyword Args:
        cot_min (float): The minimum value for cot(theta). Default value is
            1.0. Different value might be provided in the National Annexes.
        cot_max (float): The maximum value for cot(theta). Default value is
            2.5. Different value might be provided in the National Annexes.

    Raises:
        ValueError if the chosen angle is not within the given bounds.
    """
    # Use round to allow less precise angles (i.e. 21.8 degrees instead
    # of 21.801..... degrees for cot_max = 2.5).
    theta_ = math.radians(theta)
    if (
        round(1.0 / math.tan(theta_), 2) < cot_min
        or round(1.0 / math.tan(theta_), 2) > cot_max
    ):
        raise ValueError(
            'Wrong value for theta is chosen. Theta has '
            f'to be chosen such that 1/tan(theta) lies between '
            f'{cot_min} and {cot_max}. This corresponds to an angle '
            f'between {round(math.degrees(math.atan(1/cot_min)),2)} '
            f'and {round(math.degrees(math.atan(1/cot_max)),2)} '
            f'degrees, respectively. Current angle is set at {theta}'
            ' degrees.'
        )


# Equation (6.11N)
def alpha_cw(Ned: float, Ac: float, fcd: float) -> float:
    """Calculate factor that affects the maximum shear resistance of the
    concrete based on the prestress.

    EN 1992-1-1 (2005). Eq. (6.11N)

    Args:
        NEd (float): The normal force in the cross-section due to loading or
            prestress (NEd > 0 for compression) in N.
        Ac (float): The cross-sectional area of the concrete in mm2.
        fcd (float): The design strength of the concrete in MPa.

    Returns:
        float: Factor that affects the maximum shear resistance of the concrete
        based on the level of prestress.

    Raises:
        ValueError: The applied prestress exceeds the concrete design strength.
    """
    # No function call for sigma_cp, value is allowed to be higher than
    # 0.2fcd.
    sigma_cp = Ned / Ac
    if sigma_cp <= 0.0:
        value = 1.0
    elif sigma_cp <= 0.25 * fcd:
        value = 1.0 + sigma_cp / fcd
    elif sigma_cp <= 0.5 * fcd:
        value = 1.25
    elif sigma_cp < fcd:
        value = 2.5 * (1 - sigma_cp / fcd)
    else:
        raise ValueError(
            f'sigma_cp/fcd={sigma_cp/fcd}. Prestress has to be smaller'
            ' than design compressive strength.'
        )
    return value


# Without shear reinforcement


# Equation (6.2 a + b)
def VRdc(
    fck: float,
    d: float,
    Asl: float,
    bw: float,
    NEd: float,
    Ac: float,
    fcd: float,
    k1: float = 0.15,
    gamma_c: float = 1.5,
) -> float:
    """Compute the design strength of the shear resistance.

    EN 1992-1-1 (2005), Eq. (6.2)

    Args:
        fck (float): The characteristic compressive strength in MPa.
        d (float): The effective depth of the cross-section in mm.
        Asl (float): The cross-sectional area of the tensile reinforcement,
            anchored atleast (lbd + d) beyond the considered cross-section, in
            mm2.
        bw (float): The smallest width of the cross-section in tension in mm.
        NEd (float): The normal force in the cross-section due to loading or
            prestress (NEd > 0 for compression) in N.
        Ac (float): The cross-sectional area of the concrete in mm2.
        fcd (float): The design compressive strength in MPa.

    Keyword Args:
        k1 (float): Factor used to include the effect of the normal stress
            into the shear resistance of the concrete. Default value = 0.15,
            value might differ between National Annexes.
        gamma_c (float): Partial factor for concrete. Default value = 1.5,
            value might differ between National Annexes.

    Returns:
        float: The concrete shear resistance in MPa.
    """
    CRdc = 0.18 / gamma_c
    return (
        max(
            CRdc * _k(d) * (100 * _rho_L(Asl, bw, d) * fck) ** (1.0 / 3.0)
            + k1 * _sigma_cp(NEd, Ac, fcd),  # VRdc
            vmin(fck, d) + k1 * _sigma_cp(NEd, Ac, fcd),  # VRdcmin
        )
        * bw
        * d
    )


# Equation (6.3N)
def vmin(fck: float, d: float) -> float:
    """Compute the minimum shear resistance of the concrete.

    EN 1992-1-1 (2005), Eq. (6.3)

    Args:
        fck (float): The characteristic compressive strength in MPa.
        d (float): The effective depth of the cross-section in mm.

    Returns:
        float: The minimal shear stress resistance of the concrete in MPa.
    """
    return 0.035 * _k(d) ** (3.0 / 2.0) * fck ** (1 / 2)


# Equation (6.4)
def VRdc_prin_stress(
    Iy: float,
    bw: float,
    S: float,
    fctd: float,
    NEd: float,
    Ac: float,
    L_x: float = None,
    L_pt2: float = None,
) -> float:
    """Calculate the shear resistance in uncracked, prestressed elements
    without shear reinforcement, value is determined via Mohr's circle.

    The maximal value of the principle tensile stress does no necessarily lay
    at the centre of gravity. If this is the ase the minimum value of the shear
    resistance and corresponding stress needs to be found at the relevant
    location.

    EN 1992-1-1 (2005), Eq. (6.4).

    Args:
        Iy (float): The second moment of area of the considered cross-section
            in mm4.
        bw (float): The width of the cross-section at the centre of gravity.
        S (float): The first moment of area of the considered cross-section of
            the part above the centre of gravity, and with respect to the
            centre of gravity in mm3.
        fctd (float): Design value of the tensile strength of the concrete.
        NEd (float): The normal force in the cross-section due to loading or
            prestress (NEd > 0 for compression) in N.
        Ac (float): The cross-sectional area of the concrete in mm2.

    Keyword Args:
        L_x (float): Distance from the considered cross-section until the
            starting point of the transference length of the prestress steel.
            This value should be provided when the prestressing steel is
            prestreched. Default value is None.
        L_pt2 (float): Maximum value of the transference length of the
            prestress steel, according to Eq. (8.18). This value should be
            provided when the prestressing steel is prestreched. Default value
            is None.

    Returns:
        float: The maximum allowable shear force in N for an uncracked,
        prestressed element without shear reinforcement, determined from
        maximum allowable principle stress.
    """
    # No function call for sigma_cp, value is allowed to be higher than
    # 0.2fcd.
    sigma_cp = NEd / Ac
    alpha_L = 1.0 if L_x is None or L_pt2 is None else _alpha_l(L_x, L_pt2)

    return Iy * bw / S * math.sqrt(fctd**2 + alpha_L * sigma_cp * fctd)


# Equation (6.5)
def VEdmax_unreinf(
    bw: float,
    d: float,
    fck: float,
    fcd: float,
) -> float:
    """Calculate the maximum allowable shear force for cross-sections without
    shear reinforcement.

    En 1992-1-1 (2005), Eq. (6.5).

    Args:
        bw (float): The smallest width of the cross-section in tension in mm.
        d (float): The effective depth of the cross-section in mm.
        fck (float): The characteristic compressive strength in MPa.
        fcd (float): The design compressive strength in MPa.

    Returns:
        float: The maximum allowable shear force in the cross-section in N.
        When a reduced shear force may be considered for the calculations, the
        unreduced shear force has to comply to this value.
    """
    return 0.5 * bw * d * v(fck) * fcd


# Equation (6.6N)
def v(fck: float) -> float:
    """Calculate a strength redcution factor for concrete cracked by shear
    forces.

    EN 1992-1-1 (2005), Eq. (6.6N)

    Args:
        fck (float): The characteristic compressive strength in MPa.

    Returns:
        float: A concrete reduction factor to account for concrete cracked by
        shear forces.
    """
    return 0.6 * (1 - fck / 250.0)


# Equation (6.10N)
def v1(fck: float) -> float:
    """Calculate a strength redcution factor for concrete cracked by shear
    forces.

    EN 1992-1-1 (2005), Eq. (6.10N)

    Args:
        fck (float): The characteristic compressive strength in MPa.

    Returns:
        float: A concrete reduction factor to account for concrete cracked by
        shear forces.
    """
    return 0.6 if fck <= 60 else max(0.9 - fck / 200.0, 0.5)


# With shear reinforcement


# Equation (6.8 & 6.13)
# For alpha == 90 degrees, Equation (6.13) reduces to Equation (6.8).
def VRds(
    Asw: float,
    s: float,
    z: float,
    theta: float,
    fyk: float,
    alpha: float = 90.0,
    gamma_s: float = 1.15,
) -> float:
    """Calculate the shear resistance of vertical shear reinforcement.

    EN 1992-1-1 (2005). Eq. (6.8)

    Args:
        Asw (float): the cross-sectional area of the shear reinforcement in
            mm2.
        s (float): The centre-to-centre distance of the shear reinforcement in
            mm.
        z (float): The inner lever arm of internal forces in mm.
        theta (float): The angle of the compression strut in degrees.
        fyk (float): The characteristic strength of the reinforcement steel in
            MPa.

    Keyword Args:
        alpha (float): The angle of the shear reinforcement with respect to the
            neutral axis in degrees. Default value = 90 degrees.
        gamma_s (float): Partial factor of the reinforcement steel. Default
            value = 1.15. Value might differ between National Annexes.

    Returns:
        float: The shear resistance of the shear reinforcement in N.

    Raises:
        ValueError: When theta < 21.8 degrees or theta > 45 degrees.
    """
    fywd = fyk / gamma_s
    _theta(theta)
    theta = math.radians(theta)
    alpha = math.radians(alpha)
    return (
        Asw
        / s
        * z
        * fywd
        * (1 / math.tan(theta) + 1.0 / math.tan(alpha))
        * math.sin(alpha)
    )


# Equation (6.9 & 6.14)
# For alpha == 90 degrees, Equation (6.14) reduces to Equation (6.9).
def VRdmax(
    bw: float,
    z: float,
    fck: float,
    theta: float,
    NEd: float,
    Ac: float,
    fcd: float,
    alpha: float = 90.0,
    limit_fyd: bool = False,
) -> float:
    """Calculate the maximum shear strength of the compression strut.

    EN 1992-1-1 (2005). Eq. (6.9)

    Args:
        bw (float): The smallest width of the cross-section in tension in mm.
        z (float): The inner lever arm of internal forces in mm.
        fck (float): The characteristic compressive strength in MPa.
        theta (float): The angle of the compression strut in degrees.
        NEd (float): The normal force in the cross-section due to loading or
            prestress (NEd > 0 for compression) in N.
        Ac (float): The cross-sectional area of the concrete in mm2.
        fcd (float): The design compressive strength in MPa.

    Keyword Args:
        alpha (float): The angle of the shear reinforcement with respect to the
            neutral axis in degrees. Default value = 90 degrees.
        limit_fyd (bool): Flag to indicate if the design yield stress is
            limited to 0.8 * fyk or not. This controls whether the stress
            reduction factor of concrete is given by Eq. (6.6) (False) or
            (6.10) (True).

    Returns:
        float: The shear strength of the shear reinforcement in N.

    Raises:
        ValueError: When theta < 21.8 degrees or theta > 45 degrees.
        ValueError: When sigma_cp > fcd.
    """
    _theta(theta)
    theta = math.radians(theta)
    alpha = math.radians(alpha)
    strength_reduction = v(fck) if not limit_fyd else v1(fck)
    return (
        alpha_cw(NEd, Ac, fcd)
        * bw
        * z
        * strength_reduction
        * fcd
        * (1.0 / math.tan(theta) + 1.0 / math.tan(alpha))
        / (1 + 1.0 / math.tan(theta) ** 2)
    )


# Equation (6.12 & 6.15)
# For alpha == 90 degrees, Equation (6.15) reduces to Equation (6.12).
def Asw_max(
    fcd: float,
    fck: float,
    bw: float,
    s: float,
    fywd: float,
    NEd: float,
    Ac: float,
    alpha: float = 90.0,
) -> float:
    """Calculate the maximum cross-sectional area of the shear reinforcement,
    based on the assumption 1/tan(theta) == 1.

    EN 1992-1-1 (2005). Eq. (6.13)

    Args:
        fcd (float): The design strength of the concrete in MPa.
        fck (float): The characteristic compressive strength in MPa.
        bw (float): The smallest width of the cross-section in tension in mm.
        s (float): The centre-to-centre distance of the shear reinforcement in
            mm.
        fwyd (float): The design strength of the shear reinforcement steel in
            MPa.
        NEd (float): The normal force in the cross-section due to loading or
            prestress (NEd > 0 for compression) in N.
        Ac (float): The cross-sectional area of the concrete in mm2.

    Keyword Args:
        alpha (float): The angle of the shear reinforcement with respect to the
            neutral axis in degrees. Default value = 90 degrees.

    Returns:
        float: The maximum allowable cross-sectional area of the shear
        reinforcement in mm2.

    Raises:
        ValueError: When sigma_cp > fcd.
    """
    alpha = math.radians(alpha)
    return (
        1.0
        / 2.0
        * alpha_cw(NEd, Ac, fcd)
        * v(fck)
        * fcd
        * bw
        * s
        / (fywd * math.sin(alpha))
    )
