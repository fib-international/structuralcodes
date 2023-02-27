"Shear at the interface between concrete with different casting times"
from math import pi, sin, cos


def tau_edi(beta: float, v_ed: float, z: float, b_i: float):
    """Shear at the interface between cocrete cast at different times
    fib Model Code 2010, eq. (7.3-49)
    Args:
        beta (float): The ratio of longitudinal force in the new concrete and
        the longitudinal force in either compression or tension zone
        z (float): The inner lever arm of the composed section in mm
        b_i (float): The width of the inerface in mm
        v_ed (float): The shear force at the interface in N

    return:
        The shear force that should be used at the intersection"""
    return (beta * v_ed) / (z * b_i)


def tau_rdi_without_reinforceent(
    c_a: float,
    f_ctd: float,
    mu: float,
    sigma_n: float,
    f_ck: float,
    f_cd: float,
):
    """Shear resistance without reinforcement at the intesection with
    different casting time

    fib Model Code 2010, eq. (7.3-50)

    Args:
        c_a (float): The coefficient for adhesive bond (tabel 7.3-1)
        mu (float): The friction coefficient
        sigma_n (float): The loweat expected compressiv stress from
        normal forces in MPa
        f_ck (float): Characteristic strength in MPa
        f_cd (float): The design value of cylinder compressive
        strength concrete in MPa

    return:
        The shear resistance without reinforcement at the intesection with
    different casting time"""

    v = min(0.55 * (30 / f_ck) ** (1 / 3), 0.55)
    return min((c_a * f_ctd) + (mu * sigma_n), 0.5 * v * f_cd)


def tau_rdi_with_reinforcement(
    c_r: float,
    k1: float,
    k2: float,
    mu: float,
    ro: float,
    sigma_n: float,
    alfa: float,
    beta_c: float,
    f_ck: float,
    f_yd: float,
    f_cd: float,
):
    """Shear resistance with reinforcement or dowels at the intesection with
    different casting time

    fib Model Code 2010, eq. (7.3-51)

    Args:
        c_r (float): Coefficient for aggregate interlock effects (tabel 7.3-2)
        k1 (float): The interction coefficient for tensile
        force activated in reinforcment (tabel 7.3-2)
        k2 (float): The interction coeffiction for flexural
        resistance (tabel 7.3-2)
        mu (float): The friction coefficient (tabel 7.3-2)
        ro (float): The reinforcement ratio of reinforing steel
        crossing the interface
        sigma_n (float): The loweat expected compressiv stress resulting
        from normal forces acting on the interface in MPa
        alfa (float): The inclination of reinforcement crossing the
        interface (tabel 7.3-14)
        beta_c (float): The coefficient for strength of
        compresstion strut (tabel 7.3-2)
        f_ck (float): Characteristic strength in MPa
        f_yd (float): design strength of reinforment steel in MPa
        f_cd (float): The design value of cylinder compressive
        strength concrete

    return:
        Shear resistance with reinforcement at intesection with
        different casting time"""
    v = min(0.55 * (30 / f_ck) ** (1 / 3), 0.55)
    return min(
        (c_r * f_ck ** (1 / 3))
        + (mu * sigma_n)
        + k1
        * ro
        * f_yd
        * (
            ro * sin(alfa * pi / 180)
            + cos(alfa * pi / 180)
            + k2 * ro * (f_yd * f_cd) ** 0.5
        ),
        beta_c * v * f_cd,
    )
