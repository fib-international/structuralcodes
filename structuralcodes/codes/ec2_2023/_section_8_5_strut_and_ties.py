"""Functions from Section 8.5 of EN 1992-1-1:2023."""

import math


def sigma_cd_strut(F_cd: float, b_c: float, t: float) -> float:
    """Calculate the compressive stress in a concrete
        strut or compression field.

    EN1992-1-1:2023 Eq. (8.113)

    Args:
        F_cd (float): Compressive force of the strut in kN.
        b_c (float): Width of the strut at the considered location in mm.
        t (float): Thickness of the strut in mm.

    Returns:
        float: Compressive stress σ_cd in MPa.

    Raises:
        ValueError: If b_c or t are not within valid ranges.
    """
    if b_c <= 0:
        raise ValueError(f'b_c must be positive. Got {b_c}')
    if t <= 0:
        raise ValueError(f't must be positive. Got {t}')

    # Convert F_cd from kN to N and calculate σ_cd in MPa
    return abs(F_cd) * 1000 / (b_c * t)


def nu_strut(theta_cs: float) -> float:
    """Determine the strength reduction factor ν based on
        the smallest angle θ_cs.

    EN1992-1-1:2023 Eqs. (8.119)

    Args:
        theta_cs (float): Angle between the strut and the tie in degrees.
            Must be between 20 and 90.
        transverse_cracking (bool): Indicates if the region has
            transverse cracking. Defaults to True.

    Returns:
        float: Strength reduction factor ν.

    Raises:
        ValueError: If θ_cs is not within the valid range.
    """
    if not 20 <= theta_cs <= 90:
        raise ValueError(
            f'theta_cs must be between 20° and 90°. Got {theta_cs}'
        )

    theta_rad = math.radians(theta_cs)
    cot_theta = 1 / math.tan(theta_rad)
    return 1 / (1.11 + 0.22 * (cot_theta**2))


def nu_strut_no_crack() -> float:
    """Determine the strength reduction factor ν in
        areas without transverse cracking.

    EN1992-1-1:2023 Eqs. (8.120)

    Returns:
        float: Strength reduction factor ν.
    """
    return 1.0


def nu_refined(eps_1: float) -> float:
    """Calculate a more refined value for the strength reduction factor ν
    for compression fields in cracked zones based on principal tensile strain.

    EN1992-1-1:2023 Eq. (8.121)

    Args:
        epsilon_1 (float): Maximum principal tensile strain (dimensionless).

    Returns:
        float: Refined strength reduction factor ν,
            capped at a maximum value of 1.0.

    Raises:
        ValueError: If epsilon_1 is negative.

    """
    if eps_1 < 0:
        raise ValueError(f'epsilon_1 must be non-negative. Got {eps_1}')

    nu = 1 / (1.0 + 110 * eps_1)
    return min(nu, 1.0)


def FRd_tie(
    As: float,
    fyd: float,
    Ap: float = 0.0,
    fpd: float = 0.0,
    sigma_pd: float = 0.0,
) -> float:
    """Calculate the resistance of a tie.

    EN1992-1-1:2023 Eq. (8.122)

    Args:
        As (float): Cross-sectional area of non-prestressed
            reinforcement in mm2.
        fyd (float): Design yield strength of non-prestressed
            reinforcement in MPa.
        Ap (float, optional): Cross-sectional area of prestressed
            reinforcement in mm2. Default is 0.0.
        fpd (float, optional): Design strength of prestressed
            reinforcement in MPa. Default is 0.0.
        sigma_pd (float, optional): Stress in the prestressed reinforcement
            considered as an external action in MPa. Default is 0.0.

    Returns:
        float: Resistance of the tie in kN.

    Raises:
        ValueError: If any of the input values are negative.
    """
    # Input validation
    if As < 0:
        raise ValueError(f'As must not be negative. Got {As}')
    if fyd < 0:
        raise ValueError(f'fyd must not be negative. Got {fyd}')
    if Ap < 0:
        raise ValueError(f'Ap must not be negative. Got {Ap}')
    if fpd < 0:
        raise ValueError(f'fpd must not be negative. Got {fpd}')
    if sigma_pd < 0:
        raise ValueError(f'sigma_pd must not be negative. Got {sigma_pd}')

    # Calculate tie resistance
    FRd = As * fyd + Ap * (fpd - sigma_pd)
    return FRd / 1000


def Ftd_conc(
    Fd: float, a: float, b: float, H: float, near_edge: bool = False
) -> float:
    """Calculate the transverse reinforcement for concentrated
        forces spreading into a member using the strut-and-ties model.

    EN1992-1-1:2023 Eq. (8.123), Eq. (8.124), Eq. (8.125)

    Args:
        Fd (float): Design value of concentrated force in kN.
        a (float): Width of the concentrated force application area in mm.
        b (float): Width of the member in which force spreads in mm.
        H (float): Height of the member section in mm.
        near_edge (bool, optional): Indicates if the force is acting
            near an edge. Default is False.

    Returns:
        float: the transverse reinforcement force F_td in kN.

    Raises:
        ValueError: If any of the input values are negative.
    """
    # Input validation
    if Fd < 0:
        raise ValueError(f'Fd must not be negative. Got {Fd}')
    if a < 0:
        raise ValueError(f'a must not be negative. Got {a}')
    if b < 0:
        raise ValueError(f'b must not be negative. Got {b}')
    if H < 0:
        raise ValueError(f'H must not be negative. Got {H}')

    # Calculate tan(theta_cf) based on position and geometry
    tan_theta_cf = (1 - a / b) / 2 if b <= a + H / 2 else 0.5

    # For forces near an edge, tan_theta_cf is assumed to be at least 1/4
    if near_edge:
        tan_theta_cf = max(tan_theta_cf, 1 / 4)
        F_td = Fd * tan_theta_cf
    else:
        F_td = Fd / 2 * tan_theta_cf

    return F_td
