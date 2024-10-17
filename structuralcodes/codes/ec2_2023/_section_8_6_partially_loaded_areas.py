"""Functions from Section 8.6 of EN 1992-1-1:2023."""

import math


def sigma_Rd_u(
    fcd: float,
    a0: float,
    b0: float,
    a1: float,
    b: float,
    ea: float = 0,
    eb: float = 0,
    nu_part: float = 3.0,
) -> float:
    """Calculate the design resistance for partially loaded areas without
    horizontal force components.

    EN1992-1-1:2023 Eq. (8.126).

    Args:
        fcd (float): Design value of concrete compressive strength in MPa.
        a0 (float): Length of the loaded area in the direction perpendicular to
            the closest edge in mm.
        b0 (float): Width of the loaded area in mm.
        a1 (float): Length of the load introduction block parallel to a0 in mm.
        b (float): Total width of the supporting area in mm.
        ea (float): Eccentricity of the applied load parallel to a0 in mm.
            Default is 0.
        eb (float): Eccentricity of the applied load parallel to b0 (mm).
            Default is 0.
        nu_part (float): Confinement factor. Default is 3.0.

    Returns:
        float: Design resistance in MPa.

    Raises:
        ValueError: If any of the input dimensions or resistances are negative.
    """
    if fcd < 0:
        raise ValueError(f'fcd must not be negative. Got {fcd}')
    if a0 < 0:
        raise ValueError(f'a0 must not be negative. Got {a0}')
    if b0 < 0:
        raise ValueError(f'b0 must not be negative. Got {b0}')
    if a1 < 0:
        raise ValueError(f'a1 must not be negative. Got {a1}')
    if b < 0:
        raise ValueError(f'b must not be negative. Got {b}')
    if ea < 0:
        raise ValueError(f'ea must not be negative. Got {ea}')
    if eb < 0:
        raise ValueError(f'eb must not be negative. Got {eb}')

    # Eccentrically reduced loaded area
    a0_red = a0 - 2 * ea
    b0_red = b0 - 2 * eb
    Ac0_red = a0_red * b0_red  # in mm²

    # Contributing concrete area
    Ac1 = a1 * min(b0 + (a1 - a0), b)  # in mm²

    # Design resistance
    sigma_Rdu = fcd * math.sqrt(Ac1 / Ac0_red)
    return min(sigma_Rdu, nu_part * fcd)
