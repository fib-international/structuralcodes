"""Strength reduction factors according to ACI 318-25, Chapter 21."""

import typing as t


def phi_shear() -> float:
    """Strength reduction factor for shear.

    ACI 318-25, Table 21.2.1(b).

    Returns:
        float: Strength reduction factor phi = 0.75.
    """
    return 0.75


def phi_torsion() -> float:
    """Strength reduction factor for torsion.

    ACI 318-25, Table 21.2.1(c).

    Returns:
        float: Strength reduction factor phi = 0.75.
    """
    return 0.75


def phi_bearing() -> float:
    """Strength reduction factor for bearing on concrete.

    ACI 318-25, Table 21.2.1(d).

    Returns:
        float: Strength reduction factor phi = 0.65.
    """
    return 0.65


def phi_flexure(
    eps_t: float,
    fy: float,
    Es: float = 200000.0,
    transverse: t.Literal['spiral', 'other'] = 'other',
) -> float:
    """Strength reduction factor for flexure and axial loads.

    ACI 318-25, Table 21.2.2. The factor depends on the net tensile strain
    eps_t at the extreme tension steel layer, the yield strain eps_ty = fy/Es,
    and the type of transverse reinforcement.

    Zones:
    - Compression-controlled (eps_t <= eps_ty): phi = 0.65 (other) or 0.75
      (spiral).
    - Tension-controlled (eps_t >= eps_ty + 0.003): phi = 0.90.
    - Transition: linearly interpolated between the two limits.
      - other:   phi = 0.65 + 0.25 * (eps_t - eps_ty) / 0.003
      - spiral:  phi = 0.75 + 0.15 * (eps_t - eps_ty) / 0.003

    Args:
        eps_t (float): Net tensile strain at extreme tension steel layer
            (dimensionless). Positive in tension.
        fy (float): Specified yield strength of reinforcement in MPa.
            Must be > 0.
        Es (float): Modulus of elasticity of reinforcement in MPa.
            Defaults to 200000.0 per Section 20.2.2.2.
        transverse (Literal['spiral', 'other']): Type of transverse
            reinforcement. 'spiral' gives a higher phi in the
            compression-controlled and transition zones.

    Returns:
        float: Strength reduction factor phi.

    Raises:
        ValueError: If fy <= 0 or Es <= 0 or transverse is not valid.
    """
    if fy <= 0:
        raise ValueError(f'fy must be positive, got {fy}')
    if Es <= 0:
        raise ValueError(f'Es must be positive, got {Es}')
    if transverse not in ('spiral', 'other'):
        raise ValueError(
            f"transverse must be 'spiral' or 'other', got {transverse!r}"
        )

    eps_ty = fy / Es
    eps_tension = eps_ty + 0.003

    if eps_t <= eps_ty:
        # Compression-controlled zone
        return 0.75 if transverse == 'spiral' else 0.65

    if eps_t >= eps_tension:
        # Tension-controlled zone
        return 0.90

    # Transition zone: linear interpolation
    ratio = (eps_t - eps_ty) / 0.003
    if transverse == 'spiral':
        return 0.75 + 0.15 * ratio
    return 0.65 + 0.25 * ratio


def section_classification(
    eps_t: float,
    fy: float,
    Es: float = 200000.0,
) -> t.Literal['tension-controlled', 'transition', 'compression-controlled']:
    """Classify a cross-section based on the net tensile strain.

    ACI 318-25, Table 21.2.2. Uses eps_ty = fy/Es as the yield strain.

    Args:
        eps_t (float): Net tensile strain at extreme tension steel layer
            (dimensionless). Positive in tension.
        fy (float): Specified yield strength of reinforcement in MPa.
            Must be > 0.
        Es (float): Modulus of elasticity of reinforcement in MPa.
            Defaults to 200000.0 per Section 20.2.2.2.

    Returns:
        Literal['tension-controlled', 'transition', 'compression-controlled']:
        Section classification string.

    Raises:
        ValueError: If fy <= 0 or Es <= 0.
    """
    if fy <= 0:
        raise ValueError(f'fy must be positive, got {fy}')
    if Es <= 0:
        raise ValueError(f'Es must be positive, got {Es}')

    eps_ty = fy / Es
    eps_tension = eps_ty + 0.003

    if eps_t <= eps_ty:
        return 'compression-controlled'
    if eps_t >= eps_tension:
        return 'tension-controlled'
    return 'transition'
