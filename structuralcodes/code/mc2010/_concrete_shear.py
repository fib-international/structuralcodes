"""A collection of shear formulas for concrete"""
import typing as t
import math


def vrd(
    fck: float,
    z: float,
    bw: float,
    gamma_c: float,
    asw: float,
    sw: float,
    fywd: float,
    theta: float,
    dg: float,
    Approx_lvl: float,
    alfa: float,
    ved: float,
    E: float,
    As: float,
    Med: float,
    Ved: float,
    Ned: float,
    deltaE: float
) -> float:
    """Compute the shear resistance of a web or slab.

    fib Model Code 2010, Eq. (7.3-11)

    Args:
        vrdc (float): Design shear resistance of concrete.
        vrds (float): Design shear resistance of shear reinforcement.

    Returns:
        float: Design shear resistance
    """

    return abs(
        vrdc(fck, z, bw, dg, Approx_lvl,
            epsilonx(E, As, Med, Ved, Ned, z, deltaE), alfa, ved, gamma_c)
    ) + abs(vrds(asw, sw, z, bw, fywd, theta))


def vrdc(
    fck: float,
    z: float,
    bw: float,
    dg: float,
    Approx_lvl: int,
    epsilonx: float,
    alfa: float,
    ved: float,
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

    if Approx_lvl == 1:
        kv = 180 / (1000 + 1.25 * z)
        return (kv * fsqr * z * bw) / gamma_c

    elif Approx_lvl == 2:
        kdg = max(32 / (16 + dg), 0.75)
        kv = (0.4 / (1 + 1500 * epsilonx)) * (1300 / (1000 + kdg * z))
        return (kv * fsqr * z * bw) / gamma_c
    elif Approx_lvl == 3:
        theta_min = 20 + 10000 * epsilonx
        kv = max((0.4 / (1 + 1500 * epsilonx)) * (1 - ved
                / (
                    vrdmax(
                        fck,
                        bw,
                        Approx_lvl,
                        theta_min,
                        z,
                        epsilonx,
                        alfa,
                        gamma_c,
                    )
                ),
                0,
            )
        )


def vrds(
    asw: float,
    sw: float,
    z: float,
    fywd: float,
    theta: float,
    alpha: t.optional[float] = 90.0,
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
    return (
        (asw / sw)
        * z
        * fywd
        * ((1 / math.tan(theta)) + (1 / math.tan(alpha)))
        * math.sin(alpha)
    )


def vrdmax(
    fck: float,
    bw: float,
    Approx_lvl: int,
    theta: float,
    z: float,
    epsilon_x: float,
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
    nfc = (30 / fck) ** (1 / 3)
    if nfc > 1:
        nfc = 1

    if Approx_lvl == 1:
        return (
            0.55
            * nfc
            * (fck / gamma_c)
            * bw
            * z
            * (
                ((1 / math.tan(theta)) + (1 / math.tan(alfa)))
                / (1 + (1 / math.tan(theta)) ** 2)
            )
        )

    elif Approx_lvl == 2:
        epsilon_1 = epsilon_x + (epsilon_x + 0.002) * (
            (1 / math.tan(theta)) ** 2
        )
        k_epsilon = 1 / (1.2 + 55 * epsilon_1)
        if k_epsilon > 0.65:
            k_epsilon = 0.65

        return (
            k_epsilon
            * nfc
            * (fck / gamma_c)
            * bw
            * z
            * (
                ((1 / math.tan(theta)) + (1 / math.tan(alfa)))
                / (1 + (1 / math.tan(theta)) ** 2)
            )
        )

    elif Approx_lvl == 3:
        epsilon_1 = epsilon_x + (epsilon_x + 0.002) * (
            (1 / math.tan(theta)) ** 2
        )
        k_epsilon = 1 / (1.2 + 55 * epsilon_1)
        if k_epsilon > 0.65:
            k_epsilon = 0.65

        theta_min = 20 + 10000 * epsilonx
        return (
            k_epsilon
            * nfc
            * (fck / gamma_c)
            * bw
            * z
            * (
                ((1 / math.tan(theta_min)) + (1 / math.tan(alfa)))
                / (1 + (1 / math.tan(theta_min)) ** 2)
            )
        )


def epsilonx(
    E: float, As: float, Med: float, Ved: float, Ned: float, z: float, deltaE
) -> float:
    """The maximum allowed shear resistance

    fib Model Code 2010, eq. (7.3-26) and (7.3-24)

    Args:
        fck (float): The characteristic compressive strength in MPa.
        bw (float): The width.
        theta (float): The incline of the reinforment relative to the beam axis

    Returns:
        float: The maximum allowed shear resisThe design shear resistance providance regarled by ss of
        approximatirrups"""
    return (1 / 2 * E * As) * (Med / z) + Ved + Ned * ((1 / 2) + (deltaE / z))
