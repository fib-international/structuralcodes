"""Shear at the interface between concrete with different casting times."""

from math import cos, pi, sin


def tau_edi(beta: float, v_ed: float, z: float, b_i: float) -> float:
    """Shear at the interface between concrete cast at different times.

    fib Model Code 2010, eq. (7.3-49).

    Args:
        beta (float): The ratio of longitudinal force in the new concrete and
            the longitudinal force in either compression or tension zone.
        v_ed (float): The shear force at the interface in N.
        z (float): The inner lever arm of the composed section in mm.
        b_i (float): The width of the inerface in mm.

    Returns:
        float: The shear force that should be used at the intersection.
    """
    return (beta * v_ed) / (z * b_i)


def tau_rdi_without_reinforcement(
    c_a: float,
    f_ctd: float,
    mu: float,
    sigma_n: float,
    f_ck: float,
    f_cd: float,
) -> float:
    """Shear resistance without reinforcement at the intesection with different
    casting time.

    fib Model Code 2010, eq. (7.3-50).

    Args:
        c_a (float): The coefficient for adhesive bond (Table 7.3-1).
        f_ctd (float): The design axial tensile strength of concrete.
        mu (float): The friction coefficient.
        sigma_n (float): The lowest expected compressive stress from normal
            forces in MPa.
        f_ck (float): Characteristic strength in MPa.
        f_cd (float): The design value of cylinder compressive strength
            concrete in MPa.

    Returns:
        float: The shear resistance without reinforcement at the intesection
        with different casting time.
    """
    v = min(0.55 * (30 / f_ck) ** (1 / 3), 0.55)
    return min((c_a * f_ctd) + (mu * sigma_n), 0.5 * v * f_cd)


def tau_rdi_with_reinforcement(
    c_r: float,
    k1: float,
    k2: float,
    mu: float,
    ro: float,
    sigma_n: float,
    alpha: float,
    beta_c: float,
    f_ck: float,
    f_yd: float,
    f_cd: float,
) -> float:
    """Shear resistance with reinforcement or dowels at the intesection with
    different casting time.

    fib Model Code 2010, eq. (7.3-51).

    Args:
        c_r (float): Coefficient for aggregate interlock effects (Table 7.3-2).
        k1 (float): The interaction coefficient for tensile force activated in
            reinforcment (Table 7.3-2).
        k2 (float): The interaction coeffiction for flexural resistance (Table
            7.3-2).
        mu (float): The friction coefficient (Table 7.3-2).
        ro (float): The reinforcement ratio of reinforcing steel crossing the
            interface.
        sigma_n (float): The lowest expected compressive stress resulting from
            normal forces acting on the interface in MPa.
        alpha (float): The inclination of reinforcement crossing the interface
            (Table 7.3-14).
        beta_c (float): The coefficient for strength of compression strut
            (Table 7.3-2).
        f_ck (float): Characteristic strength in MPa.
        f_yd (float): design strength of reinforment steel in MPa.
        f_cd (float): The design value of cylinder compressive strength
            concrete.

    Returns:
        float: Shear resistance with reinforcement at intesection with
        different casting time.
    """
    v = min(0.55 * (30 / f_ck) ** (1 / 3), 0.55)
    return min(
        (c_r * f_ck ** (1 / 3))
        + (mu * sigma_n)
        + k1 * ro * f_yd * (mu * sin(alpha * pi / 180) + cos(alpha * pi / 180))
        + k2 * ro * (f_yd * f_cd) ** 0.5,
        beta_c * v * f_cd,
    )
