## AASHTO LRFD Functions for Deflections

import math


def Mu(
    fc_prime: float,
    fy: float,
    rho: float,
    b: float,
    h: float,
    d: float,
    L: float,
    phi: float,
) -> float:
    """Determines the design ultimate moment.

    AASHTO LRFD 2024 10th Edition, Eq (5.6.3.2a-1)

    Args:
        fc_prime (float): compressive strength of concrete in (ksi)
        fy (float): yield strength of the steel reinforcement in (ksi)
        rho (float): reinforcement ratio
        b (float): width of the cross section in (in)
        h (float): height of the cross section in (in)
        d (float): effective depth of cross section in (in)
        L (float): span in (in)
        phi (flaot): design factor

    Returns:
        Design moment in (k-in)

    Raises:
        ValueError: If fc_prime is less than 0
        ValueError: If fy is less than 0
        ValueError: If rho is less than 0
        ValueError: If b is less than 0
        ValueError: If h is less than 0
        ValueError: If d is less than 0
        ValueError: If L is less than 0
        ValueError: If phi is less than 0
    """
    if fc_prime < 0:
        raise ValueError(f'fc_prime={fc_prime} fc_prime cannot be less than 0')
    if fy < 0:
        raise ValueError(f'fy={fy} fy cannot be less than 0')
    if rho < 0:
        raise ValueError(f'rho={rho} rho cannot be less than 0')
    if b < 0:
        raise ValueError(f'b={b} b cannot be less than 0')
    if h < 0:
        raise ValueError(f'h={h} h cannot be less than 0')
    if d < 0:
        raise ValueError(f'd={d} d cannot be less than 0')
    if L < 0:
        raise ValueError(f'L={L} L cannot be less than 0')
    if phi < 0:
        raise ValueError(f'phi={phi} phi cannot be less than 0')

    As = rho * b * d
    x = As * fy / (fc_prime * b * 0.8)

    return phi * (As * fy * (d - 0.4 * x))


def Mcr(
    fc_prime: float,
    b: float,
    h: float,
    d: float,
    L: float,
) -> float:
    """Determines the cracking moment.

    AASHTO LRFD 2024 10th Edition, Eq (5.6.3.2a-1)

    Args:
        fc_prime (float): compressive strength of concrete in (ksi)
        b (float): width of the cross section in (in)
        h (float): height of the cross section in (in)
        d (float): effective depth of cross section in (in)
        L (float): span in (in)

    Returns:
        Cracking moment in (k-in)

    Raises:
        ValueError: If fc_prime is less than 0
        ValueError: If b is less than 0
        ValueError: If h is less than 0
        ValueError: If d is less than 0
        ValueError: If L is less than 0
    """
    if fc_prime < 0:
        raise ValueError(f'fc_prime={fc_prime} fc_prime cannot be less than 0')
    if b < 0:
        raise ValueError(f'b={b} b cannot be less than 0')
    if h < 0:
        raise ValueError(f'h={h} h cannot be less than 0')
    if d < 0:
        raise ValueError(f'd={d} d cannot be less than 0')
    if L < 0:
        raise ValueError(f'L={L} L cannot be less than 0')

    fr = 0.24 * math.sqrt(fc_prime)
    Ig = b * (h**3) / 12
    y_bar = h / 2
    yt = d - y_bar

    return fr * Ig / yt


def Ie(
    fc_prime: float,
    fy: float,
    rho: float,
    b: float,
    h: float,
    d: float,
    L: float,
    phi: float,
    Es: float,
    Ec: float,
) -> float:
    """Determines the effective moment of inertia.

    AASHTO LRFD 2024 10th Edition

    Args:
        fc_prime (float): compressive strength of concrete in (ksi)
        fy (float): yield strength of the steel reinforcement in (ksi)
        rho (float): reinforcement ratio
        b (float): width of the cross section in (in)
        h (float): height of the cross section in (in)
        d (float): effective depth of cross section in (in)
        L (float): span in (in)
        phi (float): design factor
        Es (float): modulus of elasticity of steel (ksi)
        Ec (float): modulus of elasticity of concrete (ksi)

    Returns:
        Effective moment of inertia (in^4)

    Raises:
        ValueError: If fc_prime is less than 0
        ValueError: If fy is less than 0
        ValueError: If rho is less than 0
        ValueError: If b is less than 0
        ValueError: If h is less than 0
        ValueError: If d is less than 0
        ValueError: If L is less than 0
        ValueError: If phi is less than 0
        ValueError: If Es is less than 0
        ValueError: If Ec is less than 0
    """
    if fc_prime < 0:
        raise ValueError(f'fc_prime={fc_prime} fc_prime cannot be less than 0')
    if fy < 0:
        raise ValueError(f'fy={fy} fy cannot be less than 0')
    if rho < 0:
        raise ValueError(f'rho={rho} rho cannot be less than 0')
    if b < 0:
        raise ValueError(f'b={b} b cannot be less than 0')
    if h < 0:
        raise ValueError(f'h={h} h cannot be less than 0')
    if d < 0:
        raise ValueError(f'd={d} d cannot be less than 0')
    if L < 0:
        raise ValueError(f'L={L} L cannot be less than 0')
    if phi < 0:
        raise ValueError(f'phi={phi} phi cannot be less than 0')
    if Es < 0:
        raise ValueError(f'Es={Es} Es cannot be less than 0')
    if Ec < 0:
        raise ValueError(f'Ec={Ec} Es cannot be less than 0')

    As = rho * b * d
    n = Es / Ec
    Ac = n * As
    Ig = b * (h**3) / 12
    Ma = Mu(fc_prime, fy, rho, b, h, d, L, phi)
    M_cr = Mcr(fc_prime, b, h, d, L)
    two_thirds_Mcr = (2 / 3) * M_cr
    x = As * fy / (fc_prime * b * 0.8)
    Icr = ((b * x**3) / 3) + Ac * (d - x) ** 2

    if Ma < two_thirds_Mcr:
        return Ig

    if Ma > two_thirds_Mcr:
        return Icr / (1 - ((two_thirds_Mcr / Ma) ** 2) * (1 - (Icr / Ig)))


def delta(
    qqp: float,
    L: float,
    E: float,
    Ieff: float,
) -> float:
    """Determines the deflection of the beam.

    Args:
        qqp (float): quasi permanent load in (k/in)
        L (float): span (in)
        E (float): modulus of elasticity of concrete (ksi)
        Ieff (float): effective moment of inertia (in^4)

    Returns:
        The deflection in (in)

    Raises:
        ValueError: If L is less than 0
        ValueError: If E is less than 0
        ValueError: If Ieff is less than 0
    """
    if L < 0:
        raise ValueError(f'L={L} cannot be less than 0')
    if E < 0:
        raise ValueError(f'E={E} cannot be less than 0')
    if Ieff < 0:
        raise ValueError(f'Ieff={Ieff} cannot be less than 0')

    return (5 * qqp * L**4) / (384 * E * Ieff)
