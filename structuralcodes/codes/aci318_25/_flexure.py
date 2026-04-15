"""Flexural strength functions according to ACI 318-25, Ch. 22.2-22.3."""

import math

# ---------------------------------------------------------------------------
# Equilibrium helpers
# ---------------------------------------------------------------------------


def stress_block_depth_sr(
    As: float,
    fy: float,
    fc: float,
    b: float,
) -> float:
    """Depth of equivalent rectangular stress block for a singly-reinforced
    section.

    ACI 318-25, Sec. 22.2.2.4.1.  Derived from horizontal force equilibrium:
    ``a = As * fy / (0.85 * fc * b)``.

    Args:
        As (float): Area of tension reinforcement in mm².
        fy (float): Specified yield strength of reinforcement in MPa.
        fc (float): Specified compressive strength of concrete in MPa.
        b (float): Width of compression face in mm.

    Returns:
        float: Stress-block depth *a* in mm.
    """
    return As * fy / (0.85 * fc * b)


def stress_block_depth_dr(
    As: float,
    As_prime: float,
    fy: float,
    fy_prime: float,
    fc: float,
    b: float,
) -> float:
    """Depth of equivalent rectangular stress block for a doubly-reinforced
    section.

    ACI 318-25, Sec. 22.2.2.4.1.  Force equilibrium with compression steel:
    ``a = (As*fy - As'*fy') / (0.85 * fc * b)``.

    Args:
        As (float): Area of tension reinforcement in mm².
        As_prime (float): Area of compression reinforcement in mm².
        fy (float): Yield strength of tension reinforcement in MPa.
        fy_prime (float): Yield (or stress) of compression reinforcement
            in MPa.
        fc (float): Specified compressive strength of concrete in MPa.
        b (float): Width of compression face in mm.

    Returns:
        float: Stress-block depth *a* in mm.
    """
    return (As * fy - As_prime * fy_prime) / (0.85 * fc * b)


def neutral_axis_depth(a: float, beta1: float) -> float:
    """Depth to neutral axis from the compression face.

    ACI 318-25, Sec. 22.2.2.4.1:  ``c = a / beta1``.

    Args:
        a (float): Depth of equivalent rectangular stress block in mm.
        beta1 (float): Stress-block factor (dimensionless).

    Returns:
        float: Neutral-axis depth *c* in mm.
    """
    return a / beta1


def eps_t_from_c(
    c: float,
    dt: float,
    eps_cu: float = 0.003,
) -> float:
    """Net tensile strain in the extreme tension reinforcement.

    ACI 318-25, Sec. 22.2.2.1.  Linear strain compatibility:
    ``eps_t = eps_cu * (dt - c) / c``.

    Args:
        c (float): Neutral-axis depth in mm.
        dt (float): Distance from extreme compression fibre to extreme tension
            reinforcement in mm.
        eps_cu (float): Maximum usable concrete compressive strain. Defaults to
            0.003 per Sec. 22.2.2.1.

    Returns:
        float: Net tensile strain eps_t (dimensionless, positive in tension).
    """
    return eps_cu * (dt - c) / c


def eps_s_prime(
    c: float,
    d_prime: float,
    eps_cu: float = 0.003,
) -> float:
    """Strain in compression reinforcement.

    ACI 318-25, Sec. 22.2.2.1.  Linear strain compatibility:
    ``eps_s' = eps_cu * (c - d') / c``.

    Args:
        c (float): Neutral-axis depth in mm.
        d_prime (float): Distance from extreme compression fibre to centroid of
            compression reinforcement in mm.
        eps_cu (float): Maximum usable concrete compressive strain. Defaults to
            0.003 per Sec. 22.2.2.1.

    Returns:
        float: Compression steel strain (dimensionless, positive in
            compression).
    """
    return eps_cu * (c - d_prime) / c


# ---------------------------------------------------------------------------
# Nominal moment strength
# ---------------------------------------------------------------------------


def Mn_singly_reinforced(
    As: float,
    fy: float,
    fc: float,
    b: float,
    d: float,
) -> float:
    """Nominal flexural strength of a singly-reinforced rectangular section.

    ACI 318-25, Sec. 22.3.2.1:
    ``Mn = As * fy * (d - a/2)``
    where *a* is obtained from :func:`stress_block_depth_sr`.

    Args:
        As (float): Area of tension reinforcement in mm².
        fy (float): Specified yield strength of reinforcement in MPa.
        fc (float): Specified compressive strength of concrete in MPa.
        b (float): Width of compression face in mm.
        d (float): Distance from extreme compression fibre to centroid of
            tension reinforcement in mm.

    Returns:
        float: Nominal moment strength *Mn* in N·mm.
    """
    a = stress_block_depth_sr(As, fy, fc, b)
    return As * fy * (d - a / 2.0)


def Mn_doubly_reinforced(
    As: float,
    As_prime: float,
    fy: float,
    fy_prime: float,
    fc: float,
    b: float,
    d: float,
    d_prime: float,
) -> float:
    """Nominal flexural strength of a doubly-reinforced rectangular section.

    ACI 318-25, Sec. 22.3.2.1.  The caller is responsible for verifying that
    compression steel has yielded
    (``fy_prime <= Es * eps_s_prime(c, d_prime)``)
    before passing *fy_prime* as the compression-steel stress.

    ``Mn = (As*fy - As'*fy') * (d - a/2) + As'*fy' * (d - d')``

    Args:
        As (float): Area of tension reinforcement in mm².
        As_prime (float): Area of compression reinforcement in mm².
        fy (float): Yield strength of tension reinforcement in MPa.
        fy_prime (float): Yield (or actual) stress of compression reinforcement
            in MPa.  Caller verifies compression steel yields via
            :func:`eps_s_prime`.
        fc (float): Specified compressive strength of concrete in MPa.
        b (float): Width of compression face in mm.
        d (float): Distance from extreme compression fibre to centroid of
            tension reinforcement in mm.
        d_prime (float): Distance from extreme compression fibre to centroid of
            compression reinforcement in mm.

    Returns:
        float: Nominal moment strength *Mn* in N·mm.
    """
    a = stress_block_depth_dr(As, As_prime, fy, fy_prime, fc, b)
    return (As * fy - As_prime * fy_prime) * (
        d - a / 2.0
    ) + As_prime * fy_prime * (d - d_prime)


# ---------------------------------------------------------------------------
# Reinforcement limits
# ---------------------------------------------------------------------------


def As_min_slab(
    fy: float,
    b: float,
    h: float,
) -> float:
    """Minimum flexural reinforcement area for a slab.

    ACI 318-25, Sec. 7.6.1.1, referring to Sec. 24.4.3.2.
    The minimum steel ratio depends on *fy* in psi:

    - *fy* <= 50 000 psi : rho_min = 0.0020
    - *fy* <= 60 000 psi : rho_min = 0.0018
    - *fy*  > 60 000 psi : rho_min = max(0.0014, 0.0018 * 60 000 / fy_psi)

    Args:
        fy (float): Specified yield strength of reinforcement in MPa.
        b (float): Width of slab strip in mm.
        h (float): Overall slab thickness in mm.

    Returns:
        float: Minimum steel area *As_min* in mm².
    """
    fy_psi = fy * 145.038  # 1 MPa = 145.038 psi
    if fy_psi <= 50000.0:
        ratio = 0.0020
    elif fy_psi <= 60000.0:
        ratio = 0.0018
    else:
        ratio = max(0.0014, 0.0018 * 60000.0 / fy_psi)
    return ratio * b * h


def As_min_beam(
    fc: float,
    fy: float,
    bw: float,
    d: float,
) -> float:
    """Minimum flexural reinforcement area for a beam.

    ACI 318-25, Sec. 9.6.1.2:
    ``As_min = max(0.25*sqrt(fc)/fy, 1.4/fy) * bw * d``.

    Args:
        fc (float): Specified compressive strength of concrete in MPa.
        fy (float): Specified yield strength of reinforcement in MPa.
        bw (float): Web width of beam in mm.
        d (float): Distance from extreme compression fibre to centroid of
            tension reinforcement in mm.

    Returns:
        float: Minimum steel area *As_min* in mm².
    """
    ratio = max(0.25 * math.sqrt(fc) / fy, 1.4 / fy)
    return ratio * bw * d


def As_max_check(
    eps_t: float,
    fy: float,
    Es: float = 200000.0,
) -> bool:
    """Check whether the section is tension-controlled (maximum steel check).

    ACI 318-25, Sec. 21.2.2.  Returns ``True`` when the net tensile strain
    satisfies the tension-controlled limit:
    ``eps_t >= fy/Es + 0.003``.

    Args:
        eps_t (float): Net tensile strain at extreme tension steel
            (dimensionless).
        fy (float): Specified yield strength of reinforcement in MPa.
        Es (float): Modulus of elasticity of reinforcement in MPa. Defaults to
            200 000 MPa per Sec. 20.2.2.2.

    Returns:
        bool: ``True`` if tension-controlled; ``False`` otherwise.
    """
    return eps_t >= fy / Es + 0.003


# ---------------------------------------------------------------------------
# Design helper
# ---------------------------------------------------------------------------


def As_required(
    Mu: float,
    phi: float,
    fy: float,
    fc: float,
    b: float,
    d: float,
) -> float:
    """Required tension reinforcement area for a given factored moment.

    Solves the quadratic that results from setting
    ``Mu = phi * As * fy * (d - As*fy / (1.7*fc*b))``.

    Expanding and rearranging:
    ``(fy²/(1.7*fc*b)) * As² - fy*d * As + Mu/phi = 0``

    The physically meaningful (smaller) root is returned.

    Args:
        Mu (float): Factored design moment in N·mm.
        phi (float): Strength reduction factor (dimensionless).
        fy (float): Specified yield strength of reinforcement in MPa.
        fc (float): Specified compressive strength of concrete in MPa.
        b (float): Width of compression face in mm.
        d (float): Effective depth in mm.

    Returns:
        float: Required steel area *As* in mm².

    Raises:
        ValueError: If the discriminant is negative (section is undersized for
            the given moment).
    """
    a_coeff = fy**2 / (1.7 * fc * b)
    b_coeff = -fy * d
    c_coeff = Mu / phi
    discriminant = b_coeff**2 - 4.0 * a_coeff * c_coeff
    if discriminant < 0.0:
        raise ValueError(
            f'Discriminant is negative ({discriminant:.6g}): '
            'section is undersized for the given factored moment.'
        )
    return (-b_coeff - math.sqrt(discriminant)) / (2.0 * a_coeff)
