"""Concrete material properties according to ACI 318-25, Chapter 19."""

import math

LAMBDA_FACTORS = {
    'normalweight': 1.0,
    'sand-lightweight': 0.85,
    'all-lightweight': 0.75,
}


def Ec(fc: float, wc: float = 2320.0) -> float:
    """Modulus of elasticity of concrete.

    ACI 318-25, Table 19.2.2.1.

    Args:
        fc (float): Specified compressive strength of concrete in MPa.
            Must be > 0.
        wc (float): Unit weight of concrete in kg/m3. Must be in the
            range [1440, 2560]. Defaults to 2320.0 (normalweight).

    Returns:
        float: Modulus of elasticity Ec in MPa.

    Raises:
        ValueError: If fc <= 0 or wc is outside [1440, 2560].
    """
    if fc <= 0:
        raise ValueError(f'fc must be positive, got {fc}')
    if not (1440 <= wc <= 2560):
        raise ValueError(
            f'wc must be in the range [1440, 2560] kg/m3, got {wc}'
        )
    return (wc**1.5) * 0.043 * math.sqrt(fc)


def fr(fc: float, lambda_s: float = 1.0) -> float:
    """Modulus of rupture of concrete.

    ACI 318-25, Eq. 19.2.3.1.

    Args:
        fc (float): Specified compressive strength of concrete in MPa.
            Must be > 0.
        lambda_s (float): Lightweight modification factor. Must be in
            (0, 1]. Defaults to 1.0 (normalweight).

    Returns:
        float: Modulus of rupture fr in MPa.

    Raises:
        ValueError: If fc <= 0 or lambda_s is not in (0, 1].
    """
    if fc <= 0:
        raise ValueError(f'fc must be positive, got {fc}')
    if not (0 < lambda_s <= 1.0):
        raise ValueError(
            f'lambda_s must be in (0, 1], got {lambda_s}'
        )
    return 0.62 * lambda_s * math.sqrt(fc)


def beta1(fc: float) -> float:
    """Whitney stress block depth factor.

    ACI 318-25, Table 22.2.2.4.3.

    Args:
        fc (float): Specified compressive strength of concrete in MPa.
            Must be > 0.

    Returns:
        float: Stress block depth factor beta1 (dimensionless).

    Raises:
        ValueError: If fc <= 0.
    """
    if fc <= 0:
        raise ValueError(f'fc must be positive, got {fc}')
    if fc <= 28:
        return 0.85
    if fc >= 55:
        return 0.65
    return max(0.65, 0.85 - 0.05 * (fc - 28) / 7)


def eps_cu() -> float:
    """Maximum usable compressive strain in concrete.

    ACI 318-25, Section 22.2.2.1.

    Returns:
        float: Ultimate concrete strain (dimensionless), equal to 0.003.
    """
    return 0.003


def alpha1() -> float:
    """Stress block intensity factor.

    ACI 318-25, Section 22.2.2.4.1.

    Returns:
        float: Stress block intensity factor alpha1 (dimensionless),
        equal to 0.85.
    """
    return 0.85


def fct(fc: float, lambda_s: float = 1.0) -> float:
    """Splitting tensile strength of concrete.

    ACI 318-25, Section 19.2.4.3.

    Args:
        fc (float): Specified compressive strength of concrete in MPa.
            Must be > 0.
        lambda_s (float): Lightweight modification factor. Defaults to
            1.0 (normalweight).

    Returns:
        float: Splitting tensile strength fct in MPa.

    Raises:
        ValueError: If fc <= 0.
    """
    if fc <= 0:
        raise ValueError(f'fc must be positive, got {fc}')
    return 0.56 * lambda_s * math.sqrt(fc)


def lambda_factor(concrete_type: str) -> float:
    """Lightweight modification factor lambda for concrete.

    ACI 318-25, Table 19.2.4.2.

    Args:
        concrete_type (str): Type of concrete. Must be one of
            'normalweight', 'sand-lightweight', or 'all-lightweight'.

    Returns:
        float: Lightweight modification factor lambda (dimensionless).

    Raises:
        ValueError: If concrete_type is not a recognised value.
    """
    value = LAMBDA_FACTORS.get(concrete_type)
    if value is None:
        valid = ', '.join(f"'{k}'" for k in LAMBDA_FACTORS)
        raise ValueError(
            f"Unknown concrete type '{concrete_type}'. "
            f'Valid types are: {valid}.'
        )
    return value
