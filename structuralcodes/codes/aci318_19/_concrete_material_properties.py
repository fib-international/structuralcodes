"""Concrete material properties according to ACI 318-19."""

from __future__ import annotations

import math
import typing as t

from ._unit_conversions import LB_PER_CUBIC_FOOT_TO_KG_PER_CUBIC_METER

LIGHTWEIGHT_LAMBDA_LIMIT = 100.0 * LB_PER_CUBIC_FOOT_TO_KG_PER_CUBIC_METER
NORMALWEIGHT_LAMBDA_LIMIT = 135.0 * LB_PER_CUBIC_FOOT_TO_KG_PER_CUBIC_METER
MINIMUM_STRESS_BLOCK_FC = 17.0


def Ec(fc: float, wc: t.Optional[float] = None) -> float:
    """The modulus of elasticity of concrete.

    ACI 318-19, Table 19.2.2.1.

    Args:
        fc (float): The specified compressive strength of concrete in
            MPa.

    Keyword Args:
        wc (float, optional): The equilibrium density of concrete in
            kg/m3. If omitted, the normalweight expression
            ``4700 * sqrt(fc)`` is used.

    Returns:
        float: The modulus of elasticity in MPa.

    Raises:
        ValueError: If fc is not positive.
        ValueError: If wc is given outside the range 1440-2560 kg/m3.

    Note:
        wc is the ACI equilibrium density for Eq. 19.2.2.1, not the
        base material density stored on a concrete object. Omitting wc
        selects the ACI normalweight expression.
    """
    if fc <= 0:
        raise ValueError(f'fc={fc} must be positive')
    if wc is None:
        return 4700.0 * math.sqrt(fc)
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


def lambda_factor(wc: float) -> float:
    """The modification factor for lightweight concrete.

    ACI 318-19, Table 19.2.4.1(a).

    Args:
        wc (float): The equilibrium density of concrete in kg/m3.

    Returns:
        float: The lightweight modification factor (dimensionless).

    Raises:
        ValueError: If wc is not positive.
    """
    if wc <= 0:
        raise ValueError(f'wc={wc} must be positive')
    if wc <= LIGHTWEIGHT_LAMBDA_LIMIT:
        return 0.75
    if wc >= NORMALWEIGHT_LAMBDA_LIMIT:
        return 1.0
    return min(1.0, 0.0075 * wc / LB_PER_CUBIC_FOOT_TO_KG_PER_CUBIC_METER)


def eps_cu() -> float:
    """The maximum usable strain at the extreme concrete compression fiber.

    ACI 318-19, Section 22.2.2.1.

    Returns:
        float: The ultimate concrete strain (dimensionless).
    """
    return 0.003


def eps_c0() -> float:
    """The default strain at peak concrete compression.

    Returns:
        float: The peak concrete strain (dimensionless).

    Note:
        ACI 318-19 does not prescribe a complete concrete stress-strain
        curve. The value 0.002 is the library default peak-strain
        parameter for the ACI parabolic and bilinear constitutive laws.
    """
    return 0.002


def alpha1() -> float:
    """The ratio of equivalent rectangular stress block intensity.

    ACI 318-19, Section 22.2.2.4.1.

    Returns:
        float: The stress block intensity factor (dimensionless).
    """
    return 0.85


def beta1(fc: float) -> float:
    """The Whitney stress block depth factor.

    ACI 318-19, Table 22.2.2.4.3.

    Args:
        fc (float): The specified compressive strength of concrete in
            MPa.

    Returns:
        float: The stress block depth factor (dimensionless).

    Raises:
        ValueError: If fc is below 17 MPa.
    """
    if fc < MINIMUM_STRESS_BLOCK_FC:
        raise ValueError(
            f'fc={fc} must be at least {MINIMUM_STRESS_BLOCK_FC} MPa'
        )
    if fc <= 28:
        return 0.85
    if fc <= 55:
        return 0.85 - 0.20 / 27 * (fc - 28)
    return 0.65
