"""Reinforcement material properties according to ACI 318-19."""

from __future__ import annotations

import typing as t

# ACI 318-19, Table 20.2.1.3(a), SI equivalents for ASTM A615
# reinforcement grades recognized by the code.
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


def fy_design(fy: float) -> float:
    """The design yield strength of reinforcement.

    ACI 318-19 applies strength reduction factors (phi) at the
    member capacity level, not the material level. This function
    therefore returns the unreduced yield strength.

    Args:
        fy (float): The specified yield strength in MPa.

    Returns:
        float: The design yield strength in MPa.

    Raises:
        ValueError: If fy is not positive.
    """
    if fy <= 0:
        raise ValueError(f'fy={fy} must be positive')
    return fy


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

    ACI 318-19, Table 20.2.1.3(a) (SI equivalents).

    Args:
        grade (str): The ASTM reinforcement grade designation.
            One of '40', '60', '80', or '100'.

    Returns:
        Dict[str, float]: A dict with keys 'fy' (yield strength in
        MPa) and 'fu' (ultimate strength in MPa).

    Note:
        ACI 318-19 Table 20.2.1.3(a) does not provide a single
        grade-level ultimate strain. Applicable elongation requirements
        depend on the reinforcement specification and bar size, so they
        are not returned here.

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
