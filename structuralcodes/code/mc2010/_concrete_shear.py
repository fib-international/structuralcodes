"""A collection of shear formulas for concrete"""
import typing as t
import warnings
from math import pi, tan, sin


def epsilon_x(
    E: float,
    As: float,
    Med: float,
    Ved: float,
    Ned: float,
    z: float,
    delta_e: float,
) -> float:
    """The maximum allowed shear resistance

    fib Model Code 2010, eq. (7.3-26) and (7.3-24)

    Args:
        E (float): The E-modulus to the material. in n/mm^2
        AS (Float): The cross-section area. in mm^2
        Med (Float): The moment working on the material in nmm^2
        Ved (float): The shear working on the material in N
        Ned: (float): The normal force working on the material in N
        z: (float): The length to the areasenter of cross-section in mm
        delta_E (float):
    Returns:
        float: The maximum allowed shear resisThe design shear resistance
        providance regarled by ss of approximatirrups"""
    return (1 / (2 * E * As)) * (
        (abs(Med) / z) + abs(Ved) + abs(Ned) * ((1 / 2) + (delta_e / z))
    )


def v_rd(
    approx_lvl_c: int,
    approx_lvl_s: int,
    reinforcment: bool,
    fck: float,
    z: float,
    bw: float,
    dg: float,
    E: float,
    As: float,
    Med: float,
    Ved: float,
    Ned: float,
    delta_e: float,
    alfa: float,
    gamma_c: float,
    asw: float,
    sw: float,
    fywd: float,
    theta: float,
) -> float:
    """Compute the shear resistance of a web or slab.

    fib Model Code 2010, Eq. (7.3-11)

    Args:
        vrdc (float): Design shear resistance of concrete.
        vrds (float): Design shear resistance of shear reinforcement.

    Returns:
        float: Design shear resistance
    """
    if not (approx_lvl_c == 1 or 2):
        warnings.warn("Not a valid approximation level")

    if not reinforcment:
        return v_rdc(
            approx_lvl_c,
            fck,
            z,
            bw,
            dg,
            E,
            As,
            Med,
            Ved,
            Ned,
            delta_e,
            alfa,
            gamma_c,
            )
    if reinforcment and approx_lvl_s == 3:
        return min(v_rdc(
            approx_lvl_c,
            fck,
            z,
            bw,
            dg,
            E,
            As,
            Med,
            Ved,
            Ned,
            delta_e,
            alfa,
            gamma_c) +
            v_rds(asw, sw, z, fywd, theta, alfa),
            v_rd_max(
            approx_lvl_s,
            fck,
            bw,
            theta,
            z,
            E,
            As,
            Med,
            Ved,
            Ned,
            delta_e,
            alfa,
            gamma_c,
            )
            )
    elif reinforcment and approx_lvl_s == 1 or 2:
        return min(v_rds(asw, sw, z, fywd, theta, alfa), v_rd_max(
            approx_lvl_s,
            fck,
            bw,
            theta,
            z,
            E,
            As,
            Med,
            Ved,
            Ned,
            delta_e,
            alfa,
            gamma_c,
            )
            )


def v_rdc(
    approx_lvl_c: int,
    approx_lvl_s: int,
    fck: float,
    z: float,
    bw: float,
    dg: float,
    E: float,
    As: float,
    Med: float,
    Ved: float,
    Ned: float,
    delta_e: float,
    alfa: float,
    gamma_c: float = 1.5,
) -> float:

    """The design shear resistance of a web or a slab without
    shear reinforcement.

    fib Model Code 2010, Eq. (7.3-17)

    Args:
        vck (float): The characteristic compressive strength in MPa.
        z (float): the effective shear depth.
        gamma_c: Material factor.
        bw:

    Returns:
        float: Design shear resistance without shear reinforcement
    """

    if approx_lvl_c == 1:
        return v_rdc_approx1(fck, z, bw, gamma_c)

    elif approx_lvl_c == 2:
        return v_rdc_approx2(
            fck, z, bw, dg, E, As, Med, Ved, Ned, delta_e, gamma_c
        )

    elif approx_lvl_s == 3:
        return v_rds_approx3(
            fck, z, bw, E, As, Med, Ved, Ned, delta_e, alfa, gamma_c
        )


def v_rdc_approx1(
    fck: float,
    z: float,
    bw: float,
    gamma_c: float = 1.5,
) -> float:
    """The design shear resistance of a web or a slab without
    shear reinforcement.

    fib Model Code 2010, Eq. (7.3-17)

    Args:
        vck (float): The characteristic compressive strength in MPa.
        z (float): the effective shear depth.
        gamma_c: Material factor.
        bw:

    Returns:
        float: Design shear resistance without shear reinforcement
    """
    fsqr = min(fck**0.5, 8)
    kv = 180 / (1000 + 1.25 * z)
    return (kv * fsqr * z * bw) / gamma_c


def v_rdc_approx2(
    fck: float,
    z: float,
    bw: float,
    dg: float,
    E: float,
    As: float,
    Med: float,
    Ved: float,
    Ned: float,
    delta_e: float,
    gamma_c: float = 1.5,
) -> float:
    """The design shear resistance of a web or a slab without
    shear reinforcement.

    fib Model Code 2010, Eq. (7.3-17)

    Args:
        vck (float): The characteristic compressive strength in MPa.
        z (float): the effective shear depth.
        gamma_c: Material factor.
        bw:

    Returns:
        float: Design shear resistance without shear reinforcement
    """
    fsqr = min(fck**0.5, 8)

    kdg = max(32 / (16 + dg), 0.75)
    kv = (0.4 / (1 + 1500 * epsilon_x(E, As, Med, Ved, Ned, z, delta_e))) * (
        1300 / (1000 + kdg * z)
    )
    return (kv * fsqr * z * bw) / gamma_c


def v_rds_approx3(  # tror dette egentlig er vrds3
    fck: float,
    z: float,
    bw: float,
    E: float,
    As: float,
    Med: float,
    Ved: float,
    Ned: float,
    delta_e: float,
    alfa: float,
    gamma_c: float = 1.5,
) -> float:

    """The design shear resistance of a web or a slab without
    shear reinforcement.

    fib Model Code 2010, Eq. (7.3-17)

    Args:
        vck (float): The characteristic compressive strength in MPa.
        z (float): the effective shear depth.
        gamma_c: Material factor.
        bw:

    Returns:
        float: Design shear resistance without shear reinforcement
    """
    fsqr = min(fck**0.5, 8)
    theta_min = 20 + 10000 * epsilon_x(E, As, Med, Ved, Ned, z, delta_e)
    kv = max(
        (0.4 / (1 + 1500 * epsilon_x(E, As, Med, Ved, Ned, z, delta_e))) *
        (
            1 - Ved / (
                v_rd_max(
                    fck,
                    bw,
                    theta_min,
                    z,
                    E,
                    As,
                    Med,
                    Ved,
                    Ned,
                    delta_e,
                    alfa,
                    gamma_c,
                )
            ),
            0,
        )
    )
    return (kv * fsqr * z * bw) / gamma_c


def v_rds(
    asw: float,
    sw: float,
    z: float,
    fywd: float,
    theta: float,
    alpha: t.Optional[float] = pi / 2,
) -> float:
    """fib Model Code 2010, Eq. (7.3-29)
    Args:
        asw (float):
        sw (float):
        gamma_c:
        bw:

    Returns:
        float:
    ion level."""
    if 45 < theta < 20:
        warnings.warn("Too high or too low compression field angel")
    return (
        (asw / sw)
        * z
        * fywd
        * ((1 / tan(theta*pi/180)) + (1 / tan(alpha*pi/180)))
        * sin(alpha*pi/180)
    )


def v_rd_max(
    approx_lvl_s: int,
    fck: float,
    bw: float,
    theta: float,
    z: float,
    E: float,
    As: float,
    Med: float,
    Ved: float,
    Ned: float,
    delta_e: float,
    alfa: float = 0,
    gamma_c: float = 1.5,
) -> float:
    """The maximum allowed shear resistance, when there is shear reinforcment

    fib Model Code 2010, eq. (7.3-26) and (7.3-24)

    Args:
        fck (float): The characteristic compressive strength in MPa.
        bw (float): The width.
        theta (float): The incline of the reinforment relative to the beam axis

    Returns:
        float: The maximum allowed shear resistance regardless of
        approximation level"""
    if approx_lvl_s == 1:
        return v_rd_max_approx1(fck, bw, theta, z, alfa, gamma_c)

    elif approx_lvl_s == 2:
        return v_rd_max_approx2(
            fck, bw, theta, z, E, As, Med, Ved, Ned, delta_e, alfa, gamma_c
        )

    elif approx_lvl_s == 3:
        return v_rd_max_approx3(
            fck, bw, theta, z, E, As, Med, Ved, Ned, delta_e, alfa, gamma_c
        )


def v_rd_max_approx1(
    fck: float,
    bw: float,
    theta: float,
    z: float,
    alfa: float = 0,
    gamma_c: float = 1.5,
) -> float:
    """The maximum allowed shear resistance, when there is shear reinforcment

    fib Model Code 2010, eq. (7.3-26) and (7.3-24)

    Args:
        fck (float): The characteristic compressive strength in MPa.
        bw (float): The width.
        theta (float): The incline of the reinforment relative to the beam axis

    Returns:
        float: The maximum allowed shear resistance regardless of
        approximation level"""
    nfc = min((30 / fck) ** (1 / 3), 1)

    return (
        0.55
        * nfc
        * (fck / gamma_c)
        * bw
        * z
        * (
            ((1 / tan(theta*pi/180)) + (1 / tan(alfa*pi/180)))
            / (1 + (1 / tan(theta*pi/180)) ** 2)
        )
    )


def v_rd_max_approx2(
    fck: float,
    bw: float,
    theta: float,
    z: float,
    E: float,
    As: float,
    Med: float,
    Ved: float,
    Ned: float,
    delta_e: float,
    alfa: float = 0,
    gamma_c: float = 1.5,
) -> float:
    """The maximum allowed shear resistance, when there is shear reinforcment

    fib Model Code 2010, eq. (7.3-26) and (7.3-24)

    Args:
        fck (float): The characteristic compressive strength in MPa.
        bw (float): The width.
        theta (float): The incline of the reinforment relative to the beam axis

    Returns:
        float: The maximum allowed shear resistance regardless of
        approximation level"""
    nfc = min((30 / fck) ** (1 / 3), 1)

    epsilon_1 = epsilon_x(E, As, Med, Ved, Ned, z, delta_e) + (
        epsilon_x(E, As, Med, Ved, Ned, z, delta_e) + 0.002
    ) * ((1 / tan(theta*pi/180)) ** 2)
    k_epsilon = 1 / (1.2 + 55 * epsilon_1)
    if k_epsilon > 0.65:
        k_epsilon = min(0.65)

        return (
            k_epsilon
            * nfc
            * (fck / gamma_c)
            * bw
            * z
            * (
                ((1 / tan(theta*pi/180)) + (1 / tan(alfa*pi/180)))
                / (1 + (1 / tan(theta*pi/180)) ** 2)
            )
        )


def v_rd_max_approx3(
    fck: float,
    bw: float,
    theta: float,
    z: float,
    E: float,
    As: float,
    Med: float,
    Ved: float,
    Ned: float,
    delta_e: float,
    alfa: float = 0,
    gamma_c: float = 1.5,
) -> float:
    """The maximum allowed shear resistance, when there is shear reinforcment

    fib Model Code 2010, eq. (7.3-26) and (7.3-24)

    Args:
        fck (float): The characteristic compressive strength in MPa.
        bw (float): The width.
        theta (float): The incline of the reinforment relative to the beam axis

    Returns:
        float: The maximum allowed shear resistance regardless of
        approximation level"""
    nfc = min((30 / fck) ** (1 / 3), 1)

    epsilon_1 = epsilon_x(E, As, Med, Ved, Ned, z, delta_e) + (
        epsilon_x(E, As, Med, Ved, Ned, z, delta_e) + 0.002
    ) * ((1 / tan(theta*pi/180)) ** 2)
    k_epsilon = min(1 / (1.2 + 55 * epsilon_1), 0.65)

    theta_min = 20 + 10000 * epsilon_x
    return (
        k_epsilon
        * nfc
        * (fck / gamma_c)
        * bw
        * z
        * (
            ((1 / tan(theta_min*pi/180)) + (1 / tan(alfa*pi/180)))
            / (1 + (1 / tan(theta_min*pi/180)) ** 2)
        )
    )
