"""Functions from Section 8.3 of EN 1992-1-1:2023."""

from typing import List


def tau_t_i(TEd: float, Ak: float, teff_i: float) -> float:
    """Calculate the torsional shear stress in a wall element.

    EN1992-1-1:2023 Eq. (8.79)

    Args:
        TEd (float): Torsional moment applied to the section in kNm
        Ak (float): Area enclosed by the center-lines of the connecting walls,
            including inner hollow areas in mm2.
        teff_i (float): Effective wall thickness. It may be taken as A/u,
            but should not be taken as less than twice the distance
            between the outer concrete surface and the center of the
            longitudinal reinforcement in mm.

    Returns:
        float: Torsional shear stress in i-wall in MPa.

    Raises:
        ValueError: If any of the input values Ak or teff_i are not positive.
    """
    if Ak <= 0:
        raise ValueError(f'Ak must be positive. Got {Ak}')
    if teff_i <= 0:
        raise ValueError(f'teff_i must be positive. Got {teff_i}')

    return TEd * 1e6 / (2 * Ak * teff_i)


def VEd_i(tau_t_i: float, teff_i: float, zi: float) -> float:
    """Calculate the shear force in a wall element due to torsion.

    EN1992-1-1:2023 Eq. (8.80)

    Args:
        tau_t_i (float): Torsional shear stress in i-wall in MPa.
        teff_i (float): Effective wall thickness in mm.
        zi (float): Lever arm of i-wall element in mm.

    Returns:
        float: Shear force in wall element i due to torsion in kN.

    Raises:
        ValueError: If any of the input values are not positive.
    """
    if teff_i <= 0:
        raise ValueError(f'teff_i must be positive. Got {teff_i}')
    if zi <= 0:
        raise ValueError(f'zi must be positive. Got {zi}')

    return tau_t_i * teff_i * zi / 1000


def tau_t_rd_sw(
    Asw: float,
    fywd: float,
    teff: float,
    s: float,
    cot_theta: float,
    cot_theta_min: float,
) -> float:
    """Calculate the torsional capacity governed by
        yielding of the shear reinforcement.

    EN1992-1-1:2023 Eq. (8.82), (8.85)

    Args:
        Asw (float): Cross-sectional area of the shear reinforcement in mm2.
        fywd (float): Design yield stress of the shear reinforcement in MPa.
        teff (float): Effective wall thickness in mm.
        s (float): Spacing between the shear reinforcement in mm.
        theta (float): Cotangent of the angle of compression field with respect
            to the longitudinal axis.
        cot_theta_min (float): limit value for the cotangent.

    Returns:
        float: Torsional capacity in MPa.

    Raises:
        ValueError: If any of the input values are non-positive.
    """
    if Asw < 0:
        raise ValueError(f'Asw should not be negative. Got {Asw}')
    if fywd < 0:
        raise ValueError(f'fywd should not be negative. Got {fywd}')
    if teff < 0:
        raise ValueError(f'teff should not be negative. Got {teff}')
    cot_theta = max(1 / cot_theta_min, min(cot_theta, cot_theta_min))
    return cot_theta * Asw / (teff * s) * fywd


def tau_t_rd_sl(
    Asl: List[float],
    fyd: List[float],
    teff: float,
    uk: float,
    cot_theta: float,
    cot_theta_min: float,
) -> float:
    """Calculate the torsional capacity governed by yielding of
        the longitudinal reinforcement.

    EN1992-1-1:2023 Eq. (8.83), (8.85)

    Args:
        Asl (float): List of cross-sectional areas of the longitudinal
            reinforcement in mm2.
        fyd (float): List of design yield stresses of the longitudinal
            reinforcement in MPa.
        teff (float): Effective wall thickness in mm.
        uk (float): Perimeter of the area in mm.
        cot_theta (float): Cotangent of the angle of compression
            field with respect to the longitudinal axis.
        cot_theta_min (float): limit value for the cotangent.

    Returns:
        float: Torsional capacity in MPa.

    Raises:
        ValueError: If any of the input values are non-positive.
    """
    if len(Asl) != len(fyd):
        raise ValueError('Length of Asl and fyd should be the same.')
    for a in Asl:
        if a < 0:
            raise ValueError(f'Asl should not be negative. Got {a}')
    for f in fyd:
        if f < 0:
            raise ValueError(f'fyd should not be negative. Got {f}')
    if uk < 0:
        raise ValueError(f'uk should not be negative. Got {uk}')
    if teff < 0:
        raise ValueError(f'teff should not be negative. Got {teff}')

    sum_r = 0
    n = len(Asl)
    for i in range(n):
        sum_r += Asl[i] * fyd[i]

    cot_theta = max(1 / cot_theta_min, min(cot_theta, cot_theta_min))
    return sum_r / (teff * uk * cot_theta)


def tau_t_rd_max(
    nu: float, fcd: float, cot_theta: float, cot_theta_min: float
) -> float:
    """Calculate the torsional capacity governed by crushing
        of the compression field in concrete.

    EN1992-1-1:2023 Eq. (8.84), (8.85)

    Args:
        nu (float): Coefficient as determined by the formulae in Annex G.
        fcd (float): Design value of concrete compressive strength in MPa.
        cot_theta (float): Cotangent Angle of compression field with respect
            to the longitudinal axis.
        cot_theta_min (float): limit value for the cotangent.

    Returns:
        float: Torsional capacity in MPa.

    Raises:
        ValueError: If any of the input values are non-positive.
    """
    if nu < 0:
        raise ValueError(f'nu should not be negative. Got {nu}')
    if fcd < 0:
        raise ValueError(f'fcd should not be negative. Got {fcd}')

    cot_theta = max(1 / cot_theta_min, min(cot_theta, cot_theta_min))
    tan_theta = 1 / cot_theta
    return nu * fcd / (cot_theta + tan_theta)


def tau_t_rd(
    tau_t_rd_sw: float, tau_t_rd_sl: float, tau_t_rd_max: float
) -> float:
    """Calculate the design torsional capacity for a
        single cell or thin-walled section.

    EN1992-1-1:2023 Eq. (8.81)

    Args:
        tau_t_rd_sw (float): torsional capacity governed by yielding of the
            shear reinforcement in MPa.
        tau_t_rd_sl (float): torsional capacity governed by yielding of the
            longitudinal reinforcement in MPa.
        tau_t_rd_max (float): torsional capacity governed by crushing of the
            compression field in concrete in MPa.

    Returns:
        float: Design torsional capacity in MPa.

    Raises:
        ValueError: If any of the input values are non-positive.
    """
    if tau_t_rd_sw < 0:
        raise ValueError(
            f'tau_t_rd_sw should not be negative. Got {tau_t_rd_sw}'
        )
    if tau_t_rd_sl < 0:
        raise ValueError(
            f'tau_t_rd_sl should not be negative. Got {tau_t_rd_sl}'
        )
    if tau_t_rd_max < 0:
        raise ValueError(
            f'tau_t_rd_sl should not be negative. Got {tau_t_rd_max}'
        )

    return min(tau_t_rd_sw, tau_t_rd_sl, tau_t_rd_max)
