"""Reinforcement material properties according to ACI 318-25, Chapter 20."""

REINFORCEMENT_GRADES = {
    '40': {'fy': 280.0, 'fu': 420.0},
    '60': {'fy': 420.0, 'fu': 550.0},
    '80': {'fy': 550.0, 'fu': 690.0},
    '100': {'fy': 690.0, 'fu': 860.0},
}


def Es() -> float:
    """Modulus of elasticity of non-prestressed reinforcement.

    ACI 318-25, Section 20.2.2.2.

    Returns:
        float: Modulus of elasticity Es in MPa, equal to 200000.0.
    """
    return 200000.0


def fy_design(fy: float, phi: float = 1.0) -> float:
    """Design yield strength of reinforcement.

    ACI 318-25 applies strength reduction factors (phi) at the member
    level, not the material level. This function scales fy by phi when
    a value other than the default is supplied, allowing a consistent
    interface with codes that do apply material-level factors.

    Args:
        fy (float): Specified yield strength of reinforcement in MPa.
            Must be > 0.
        phi (float): Strength reduction factor. Must be in (0, 1].
            Defaults to 1.0 (no reduction).

    Returns:
        float: Design yield strength phi * fy in MPa.

    Raises:
        ValueError: If fy <= 0 or phi is not in (0, 1].
    """
    if fy <= 0:
        raise ValueError(f'fy must be positive, got {fy}')
    if not (0 < phi <= 1.0):
        raise ValueError(f'phi must be in (0, 1], got {phi}')
    return phi * fy


def epsyd(fy: float, _Es: float = 200000.0) -> float:
    """Design yield strain of reinforcement.

    ACI 318-25, Section 20.2.2.2 (derived from Es and fy).

    Args:
        fy (float): Specified yield strength of reinforcement in MPa.
            Must be > 0.
        _Es (float): Modulus of elasticity of reinforcement in MPa.
            Defaults to 200000.0 per Section 20.2.2.2.

    Returns:
        float: Yield strain fy / Es (dimensionless).

    Raises:
        ValueError: If fy <= 0.
    """
    if fy <= 0:
        raise ValueError(f'fy must be positive, got {fy}')
    return fy / _Es


def reinforcement_grade_props(grade: str) -> dict:
    """Yield and tensile strength for standard reinforcement grades.

    ACI 318-25, Table 20.2.1.3.

    Args:
        grade (str): Reinforcement grade designation. Must be one of
            '40', '60', '80', or '100'.

    Returns:
        dict: Dictionary with keys 'fy' (yield strength, MPa) and
        'fu' (tensile strength, MPa).

    Raises:
        ValueError: If grade is not a recognised designation.
    """
    props = REINFORCEMENT_GRADES.get(grade)
    if props is None:
        valid = ', '.join(f"'{k}'" for k in REINFORCEMENT_GRADES)
        raise ValueError(
            f"Unknown reinforcement grade '{grade}'. "
            f'Valid grades are: {valid}.'
        )
    return props
