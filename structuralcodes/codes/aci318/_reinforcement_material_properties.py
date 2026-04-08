"""Reinforcement material properties according to ACI 318-19."""

from __future__ import annotations

import typing as t

REINFORCEMENT_GRADES = {
    '40': {'fy': 280.0, 'fu': 420.0},
    '60': {'fy': 420.0, 'fu': 550.0},
    '80': {'fy': 550.0, 'fu': 690.0},
    '100': {'fy': 690.0, 'fu': 860.0},
}


def Es() -> float:
    """The modulus of elasticity of reinforcement.

    ACI 318-19, Section 20.2.2.2.

    Returns:
        float: The modulus of elasticity in MPa.
    """
    return 200000.0


def fy_design(fy: float, phi: float = 1.0) -> float:
    """The design yield strength of reinforcement.

    ACI 318-19 applies strength reduction factors (phi) at the
    member capacity level, not the material level. The default
    phi=1.0 returns the unreduced yield strength, which is the
    standard ACI convention for material properties.

    Args:
        fy (float): The specified yield strength in MPa.

    Keyword Args:
        phi (float): Optional strength reduction factor.
            Default is 1.0 (no reduction).

    Returns:
        float: The design yield strength in MPa.

    Raises:
        ValueError: If fy is not positive.
        ValueError: If phi is not in (0, 1].
    """
    if fy <= 0:
        raise ValueError(f'fy={fy} must be positive')
    if phi <= 0 or phi > 1.0:
        raise ValueError(f'phi={phi} must be in the range (0, 1]')
    return phi * fy


def epsyd(fy: float, _Es: float = 200000.0) -> float:
    """The yield strain of reinforcement.

    Args:
        fy (float): The specified yield strength in MPa.

    Keyword Args:
        _Es (float): The modulus of elasticity in MPa.
            Default is 200000 MPa.

    Returns:
        float: The yield strain (dimensionless).

    Raises:
        ValueError: If fy is not positive.
    """
    if fy <= 0:
        raise ValueError(f'fy={fy} must be positive')
    return fy / _Es


def reinforcement_grade_props(
    grade: t.Literal['40', '60', '80', '100'],
) -> t.Dict[str, float]:
    """Return the minimum specified properties for a reinforcement grade.

    ACI 318-19, Table 20.2.2.4a (SI equivalents).

    Args:
        grade (str): The ASTM reinforcement grade designation.
            One of '40', '60', '80', or '100'.

    Returns:
        Dict[str, float]: A dict with keys 'fy' (yield strength in
        MPa) and 'fu' (ultimate strength in MPa).

    Raises:
        ValueError: If the grade is not recognized.
    """
    props = REINFORCEMENT_GRADES.get(str(grade))
    if props is None:
        raise ValueError(
            f'Unknown reinforcement grade: {grade}. '
            f'Valid grades: {list(REINFORCEMENT_GRADES.keys())}'
        )
    return dict(props)
