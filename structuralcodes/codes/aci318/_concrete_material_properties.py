"""Concrete material properties according to ACI 318-19."""

from __future__ import annotations

import math
import typing as t

LAMBDA_FACTORS = {
    'normalweight': 1.0,
    'sand-lightweight': 0.85,
    'all-lightweight': 0.75,
}


def Ec(fc: float, wc: float = 2320.0) -> float:
    """The modulus of elasticity of concrete.

    ACI 318-19, Table 19.2.2.1.

    Args:
        fc (float): The specified compressive strength of concrete in
            MPa.

    Keyword Args:
        wc (float): The unit weight of concrete in kg/m3.
            Default is 2320 kg/m3 (normalweight concrete).

    Returns:
        float: The modulus of elasticity in MPa.

    Raises:
        ValueError: If fc is not positive.
        ValueError: If wc is outside the range 1440-2560 kg/m3.
    """
    if fc <= 0:
        raise ValueError(f'fc={fc} must be positive')
    if wc < 1440 or wc > 2560:
        raise ValueError(f'wc={wc} must be between 1440 and 2560 kg/m3')
    return wc**1.5 * 0.043 * math.sqrt(fc)


def fr(fc: float, lambda_s: float = 1.0) -> float:
    """The modulus of rupture of concrete.

    ACI 318-19, Eq. 19.2.3.1.

    Args:
        fc (float): The specified compressive strength of concrete in
            MPa.

    Keyword Args:
        lambda_s (float): The modification factor for lightweight
            concrete. Default is 1.0 (normalweight).

    Returns:
        float: The modulus of rupture in MPa.

    Raises:
        ValueError: If fc is not positive.
        ValueError: If lambda_s is not in (0, 1].
    """
    if fc <= 0:
        raise ValueError(f'fc={fc} must be positive')
    if lambda_s <= 0 or lambda_s > 1.0:
        raise ValueError(f'lambda_s={lambda_s} must be in the range (0, 1]')
    return 0.62 * lambda_s * math.sqrt(fc)


def beta1(fc: float) -> float:
    """The Whitney stress block depth factor.

    ACI 318-19, Table 22.2.2.4.3.

    Args:
        fc (float): The specified compressive strength of concrete in
            MPa.

    Returns:
        float: The stress block depth factor (dimensionless).

    Raises:
        ValueError: If fc is not positive.
    """
    if fc <= 0:
        raise ValueError(f'fc={fc} must be positive')
    if fc <= 28:
        return 0.85
    if fc >= 55:
        return 0.65
    return 0.85 - 0.05 * (fc - 28) / 7


def eps_cu() -> float:
    """The maximum usable strain at the extreme concrete compression fiber.

    ACI 318-19, Section 22.2.2.1.

    Returns:
        float: The ultimate concrete strain (dimensionless).
    """
    return 0.003


def lambda_factor(
    concrete_type: t.Literal[
        'normalweight', 'sand-lightweight', 'all-lightweight'
    ],
) -> float:
    """The modification factor for lightweight concrete.

    ACI 318-19, Table 19.2.4.2.

    Args:
        concrete_type (str): The concrete type. One of
            'normalweight', 'sand-lightweight', or 'all-lightweight'.

    Returns:
        float: The lightweight modification factor (dimensionless).

    Raises:
        ValueError: If concrete_type is not recognized.
    """
    result = LAMBDA_FACTORS.get(concrete_type.lower())
    if result is None:
        raise ValueError(
            f'Unknown concrete type: {concrete_type}. '
            f'Valid types: {list(LAMBDA_FACTORS.keys())}'
        )
    return result


def fct(fc: float, lambda_s: float = 1.0) -> float:
    """The approximate splitting tensile strength of concrete.

    ACI 318-19, Section 19.2.4.3.

    Args:
        fc (float): The specified compressive strength of concrete in
            MPa.

    Keyword Args:
        lambda_s (float): The modification factor for lightweight
            concrete. Default is 1.0 (normalweight).

    Returns:
        float: The splitting tensile strength in MPa.

    Raises:
        ValueError: If fc is not positive.
    """
    if fc <= 0:
        raise ValueError(f'fc={fc} must be positive')
    return 0.56 * lambda_s * math.sqrt(fc)


def alpha1() -> float:
    """The ratio of equivalent rectangular stress block intensity.

    ACI 318-19, Section 22.2.2.4.1.

    Returns:
        float: The stress block intensity factor (dimensionless).
    """
    return 0.85
