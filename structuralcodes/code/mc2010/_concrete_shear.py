"""A collection of shear formulas for concrete"""
import typing as t
import math


def epsilon_x(
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
    return (1 / (2 * E * As)) * ((Med / z) + Ved + Ned * ((1 / 2) + (deltaE / z)))


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
    approx_lvl: float,
    alfa: float,
    E: float,
    As: float,
    Med: float,
    Ved: float,
    Ned: float,
    deltaE: float,
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
        vrdc(
            approx_lvl,
            fck,
            z,
            bw,
            dg,
            E,
            As,
            Med,
            Ved,
            Ned,
            deltaE,
            alfa,
            gamma_c,
        )
    ) + abs(vrds(asw, sw, z, bw, fywd, theta))


def vrdc(
    approx_lvl: int,
    fck: float,
    z: float,
    bw: float,
    dg: float,
    E: float,
    As: float,
    Med: float,
    Ved: float,
    Ned: float,
    deltaE: float,
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

    if approx_lvl == 1:
        return vrdc_approx1(fck, z, bw, gamma_c)

    elif approx_lvl == 2:
        return vrdc_approx2(fck, z, bw, dg, E, As, Med, Ved, Ned, deltaE, gamma_c)

    elif approx_lvl == 3:
        return vrdc_approx3(fck, z, bw, E, As, Med, Ved, Ned, deltaE, alfa, gamma_c)


def vrdc_approx1(
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


def vrdc_approx2(
    fck: float,
    z: float,
    bw: float,
    dg: float,
    E: float,
    As: float,
    Med: float,
    Ved: float,
    Ned: float,
    deltaE: float,
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
    kv = (0.4 / (1 + 1500 * epsilon_x(E, As, Med, Ved, Ned, z, deltaE))) * (
        1300 / (1000 + kdg * z)
    )
    return (kv * fsqr * z * bw) / gamma_c


def vrdc_approx3(
    fck: float,
    z: float,
    bw: float,
    E: float,
    As: float,
    Med: float,
    Ved: float,
    Ned: float,
    deltaE: float,
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
    theta_min = 20 + 10000 * epsilon_x(E, As, Med, Ved, Ned, z, deltaE)
    kv = max(
        (0.4 / (1 + 1500 * epsilon_x(E, As, Med, Ved, Ned, z, deltaE))) *
        (
            1 - Ved / (
                vrdmax(
                    fck,
                    bw,
                    theta_min,
                    z,
                    E,
                    As,
                    Med,
                    Ved,
                    Ned,
                    deltaE,
                    alfa,
                    gamma_c,
                )
            ),
            0,
        )
    )
    return (kv * fsqr * z * bw) / gamma_c


def vrds(
    asw: float,
    sw: float,
    z: float,
    fywd: float,
    theta: float,
    alpha: t.Optional[float] = math.pi / 2,
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
    approx_lvl: int,
    fck: float,
    bw: float,
    theta: float,
    z: float,
    E: float,
    As: float,
    Med: float,
    Ved: float,
    Ned: float,
    deltaE: float,
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
    if approx_lvl == 1:
        return vrdmax_approx1(fck, bw, theta, z, alfa, gamma_c)

    elif approx_lvl == 2:
        return vrdmax_approx2(fck, bw, theta, z, E, As, Med, Ved, Ned, deltaE, alfa, gamma_c)

    elif approx_lvl == 3:
        return vrdmax_approx3(fck, bw, theta, z, E, As, Med, Ved, Ned, deltaE, alfa, gamma_c)


def vrdmax_approx1(
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
            ((1 / math.tan(theta)) + (1 / math.tan(alfa)))
            / (1 + (1 / math.tan(theta)) ** 2)
        )
    )


def vrdmax_approx2(
    fck: float,
    bw: float,
    theta: float,
    z: float,
    E: float,
    As: float,
    Med: float,
    Ved: float,
    Ned: float,
    deltaE: float,
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

    epsilon_1 = epsilon_x(E, As, Med, Ved, Ned, z, deltaE) + (
        epsilon_x(E, As, Med, Ved, Ned, z, deltaE) + 0.002
    ) * ((1 / math.tan(theta)) ** 2)
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
                ((1 / math.tan(theta)) + (1 / math.tan(alfa)))
                / (1 + (1 / math.tan(theta)) ** 2)
            )
        )


def vrdmax_approx3(
    fck: float,
    bw: float,
    theta: float,
    z: float,
    E: float,
    As: float,
    Med: float,
    Ved: float,
    Ned: float,
    deltaE: float,
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

    epsilon_1 = epsilon_x(E, As, Med, Ved, Ned, z, deltaE) + (
        epsilon_x(E, As, Med, Ved, Ned, z, deltaE) + 0.002
    ) * ((1 / math.tan(theta)) ** 2)
    k_epsilon = 1 / (1.2 + 55 * epsilon_1)
    if k_epsilon > 0.65:
        k_epsilon = 0.65

    theta_min = 20 + 10000 * epsilon_x
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
