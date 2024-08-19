"""Functions from Section 14 of EN 1992-1-1:2023."""

import math
from typing import Literal


def fcd_pl(kc_pl: float, fcd: float) -> float:
    """Calculate the design compressive strength of plain concrete.

    EN1992-1-1:2023 Eq. (14.1)

    Args:
        kc_pl (float): Coefficient for plain concrete (dimensionless).
        fcd (float): Design compressive strength of concrete in MPa.

    Returns:
        float: Design compressive strength of plain concrete in MPa.

    Raises:
        ValueError: If kc_pl or fcd are negative.
    """
    if kc_pl < 0:
        raise ValueError(f'kc_pl must not be negative. Got {kc_pl}')
    if fcd < 0:
        raise ValueError(f'fcd must not be negative. Got {fcd}')

    return kc_pl * fcd


def fctd_pl(kt_pl: float, fctd: float) -> float:
    """Calculate the design tensile strength of plain concrete.

    EN1992-1-1:2023 Eq. (14.2)

    Args:
        kt_pl (float): Coefficient for plain concrete (dimensionless).
        fctd (float): Design tensile strength of concrete in MPa.

    Returns:
        float: Design tensile strength of plain concrete in MPa.

    Raises:
        ValueError: If kt_pl or fctd are negative.
    """
    if kt_pl < 0:
        raise ValueError(f'kt_pl must not be negative. Got {kt_pl}')
    if fctd < 0:
        raise ValueError(f'fctd must not be negative. Got {fctd}')

    return kt_pl * fctd


def Nrd_pl(
    fcd_pl: float,
    b: float,
    h: float,
    e: float,
) -> float:
    """Calculate the axial resistance of a rectangular
        cross-section with a uniaxial eccentricity.

    EN1992-1-1:2023 Eq. (14.3)

    Args:
        fcd_pl (float): Design compressive
            strength of plain concrete in MPa.
        b (float): Width of the cross-section in mm.
        h (float): Height of the cross-section in mm.
        e (float): Eccentricity in mm.

    Returns:
        float: Axial resistance of the cross-section in kN.

    Raises:
        ValueError: If any input value is negative or if e > h/2.
    """
    if fcd_pl < 0:
        raise ValueError(f'fcd_pl must not be negative. Got {fcd_pl}')
    if b < 0:
        raise ValueError(f'b must not be negative. Got {b}')
    if h < 0:
        raise ValueError(f'h must not be negative. Got {h}')
    if e < 0 or e > h / 2:
        raise ValueError(f'e must be between 0 and h/2. Got {e}')

    return fcd_pl * b * h * (1 - 2 * e / h) / 1000


def sigma_cp_pl(NEd: float, Acc: float) -> float:
    """Calculate the design compressive stress.

    EN1992-1-1:2023 Eq. (14.4)

    Args:
        NEd (float): Normal compressive force in kN.
        Acc (float): Compressive area in mm2.

    Returns:
        float: Design compressive stress in MPa.

    Raises:
        ValueError: If Acc is negative.
    """
    if Acc < 0:
        raise ValueError(f'Acc must not be negative. Got {Acc}')

    return abs(NEd) * 1000 / Acc


def tau_cp_pl(VEd: float, Acc: float) -> float:
    """Calculate the design shear stress for rectangular sections.

    EN1992-1-1:2023 Eq. (14.5)

    Args:
        VEd (float): Shear force kN.
        Acc (float): Compressive area mm2.

    Returns:
        float: Design shear stress MPa.

    Raises:
        ValueError: If Acc is negative.
    """
    if Acc < 0:
        raise ValueError(f'Acc must not be negative. Got {Acc}')

    return 1.5 * VEd * 1000 / Acc


def sigma_c_lim_pl(fcd_pl: float, fctd_pl: float) -> float:
    """Calculate the limiting compressive stress.

    EN1992-1-1:2023 Eq. (14.8)

    Args:
        fcd_pl (float): Design compressive strength
            of plain concrete in MPa.
        fctd_pl (float): Design tensile strength of
            plain concrete in MPa.

    Returns:
        float: Limiting compressive stress in MPa.
    """
    return fcd_pl - 2 * math.sqrt(fctd_pl * (fctd_pl + fcd_pl))


def tau_Rd_pl(
    fctd_pl: float,
    sigma_cp_pl: float,
    sigma_c_lim_pl: float,
) -> float:
    """Calculate the plain concrete design strength in shear.

    EN1992-1-1:2023 Eq. (14.6) and Eq. (14.7)

    Args:
        fctd_pl (float): Design tensile strength of plain concrete in MPa.
        sigma_cp_pl (float): Design compressive stress in MPa.
        sigma_c_lim_pl (float): Limiting compressive stress in MPa.

    Returns:
        float: Plain concrete design strength in shear in MPa.

    Raise:
        ValueError: if any of the inputs is negative.
    """
    if fctd_pl < 0:
        raise ValueError(f'fctd_pl must be positive. Got {fctd_pl}')
    if sigma_cp_pl < 0:
        raise ValueError(f'sigma_cp_pl must be positive. Got {sigma_cp_pl}')
    if sigma_c_lim_pl < 0:
        raise ValueError(
            f'sigma_c_lim_pl must be positive. Got {sigma_c_lim_pl}'
        )

    if sigma_cp_pl <= sigma_c_lim_pl:
        return math.sqrt(fctd_pl**2 + sigma_cp_pl * fctd_pl)

    return math.sqrt(
        fctd_pl**2
        + (sigma_cp_pl * fctd_pl)
        - ((sigma_cp_pl - sigma_c_lim_pl) / 2) ** 2
    )


def beta_Eul(b: float, lw: float, num_sides: Literal[3, 4]) -> float:
    """Calculate the Euler-coefficient based on the support conditions.

    EN1992-1-1:2023 Table (14.1)

    Args:
        b (float): Length of the wall in mm.
        lw (float): Clear height of the wall in mm.
        num_sides (int): Number of sides with lateral bearing (3 or 4).

    Returns:
        float: Euler-coefficient (dimensionless).

    Raises:
        ValueError: If any input value is negative.
    """
    if b < 0:
        raise ValueError(f'b must not be negative. Got {b}')
    if lw < 0:
        raise ValueError(f'lw must not be negative. Got {lw}')
    if num_sides not in [3, 4]:
        raise ValueError(f'num_sides must be 3 or 4. Got {num_sides}')

    if num_sides == 3:
        return 1 / (1 + (lw / (3 * b)) ** 2)
    if b >= lw:
        return 1 / (1 + (lw / b) ** 2)
    return b / (2 * lw)


def l0(beta_Eul: float, lw: float, factor: float = 1.0) -> float:
    """Calculate the effective length of a column or wall.

    EN1992-1-1:2023 Eq. (14.9)

    Args:
        beta_Eul (float): Euler-coefficient (dimensionless).
        lw (Union[int, float]): Clear height of the wall in mm.
        factor (float, optional): Factor for reducing the
            effective length (default is 1.0).

    Returns:
        float: Effective length of the column or wall in mm.

    Raises:
        ValueError: If any input value is negative or if
            factor is out of range.
    """
    if beta_Eul < 0:
        raise ValueError(f'beta_Eul must not be negative. Got {beta_Eul}')
    if lw < 0:
        raise ValueError(f'lw must not be negative. Got {lw}')
    if not (0 < factor <= 1):
        raise ValueError(f'factor must be in the range [0, 1]. Got {factor}')

    return beta_Eul * lw * factor


def phi(
    e_tot: float,
    l0: float,
    h: float,
    fcd_pl: float,
    phi_eff: float,
) -> float:
    """Calculate the factor phi considering eccentricity,
        second order effects, and other factors.

    EN1992-1-1:2023 Eq. (14.11)

    Args:
        e_tot (float): Total eccentricity in mm.
        l0 (float): Effective length of the column or wall in mm.
        h (float): Height of the cross-section in mm.
        fcd_pl (float): Design compressive
            strength of plain concrete in MPa.
        phi_eff (float): Effective creep coefficient (dimensionless).

    Returns:
        float: Factor Φ (dimensionless).

    Raises:
        ValueError: If any input value is negative.
    """
    if e_tot < 0:
        raise ValueError(f'e_tot must not be negative. Got {e_tot}')
    if l0 < 0:
        raise ValueError(f'l0 must not be negative. Got {l0}')
    if h < 0:
        raise ValueError(f'h must not be negative. Got {h}')
    if fcd_pl < 0:
        raise ValueError(f'fcd_pl must not be negative. Got {fcd_pl}')
    if phi_eff < 0:
        raise ValueError(f'phi_eff must not be negative. Got {phi_eff}')

    term1 = 2.1 + 0.02 * (l0 / h)
    term2 = 0.9 + 6 * (e_tot / h)
    term3 = (0.8 + phi_eff) / 1000
    term4 = (fcd_pl / 20) ** 0.6

    num = 1 - (term1 * e_tot / h)
    den = (1 + (l0 / h) ** 2) * term2 * term3 * term4

    return num / den


def e_tot(e0: float, ei: float) -> float:
    """Calculate the total eccentricity.

    EN1992-1-1:2023 Eq. (14.12)


    Args:
        e0 (float): First order eccentricity in mm.
        ei (float): Additional eccentricity covering the effects
            of geometrical imperfections in mm.

    Returns:
        float: Total eccentricity in mm.
    """
    if e0 < 0:
        raise ValueError(f'e0 must not be negative. Got {e0}')
    if ei < 0:
        raise ValueError(f'ei must not be negative. Got {ei}')

    return e0 + ei


def NRd_pl_simp(
    b: float,
    h: float,
    fcd_pl: float,
    phi: float,
) -> float:
    """Calculate the design resistance in terms of
        axial force for a braced wall or column in plain concrete.

    EN1992-1-1:2023 Eq. (14.10)

    Args:
        b (float): Width of the cross-section in mm.
        h (float): Height of the cross-section in mm.
        fcd_pl (float): Design compressive strength of plain concrete in MPa.
        phi (float): Factor considering eccentricity and
            second order effects (dimensionless).

    Returns:
        float: Design resistance in kN.

    Raises:
        ValueError: If any input value is negative.
    """
    if b < 0:
        raise ValueError(f'b must not be negative. Got {b}')
    if h < 0:
        raise ValueError(f'h must not be negative. Got {h}')
    if fcd_pl < 0:
        raise ValueError(f'fcd_pl must not be negative. Got {fcd_pl}')
    if phi < 0:
        raise ValueError(f'phi must not be negative. Got {phi}')

    return b * h * fcd_pl / 1000 * phi


def min_tw() -> bool:
    """Minimum requirement for overall thickness of a cast in-situ
        plain concrete wall.

    EN1992-1-1:2023 14.6.1 (1)

    Returns:
        bool: the minimum thickness in mm.
    """
    return 120.0


def min_footing_depth_ratio(
    sigma_gd: float,
    fctd_pl: float,
) -> bool:
    """Minimum footing depth for an axially loaded strip and pad footing.

    EN1992-1-1:2023 Eq. (14.13)

    Args:
        sigma_gd (float): The design value of ground pressure in MPa.
        fctd_pl (float): The design value of the tensile
            strength of plain concrete in MPa.

    Returns:
        float: the minimum hf/af ratio for a pad footing where hf is
            the footing depth and af is the distance from
            the footing edge to the column or wall face

    Raises:
        ValueError: If any of the input values are negative.
    """
    if sigma_gd < 0:
        raise ValueError(
            f'sigma_gd (ground pressure) must not be negative. Got {sigma_gd}'
        )
    if fctd_pl < 0:
        raise ValueError(
            f'fctd_pl (tensile strength) must not be negative. Got {fctd_pl}'
        )

    return math.sqrt(3 * sigma_gd / fctd_pl) / 0.85


def min_footing_depth_ratio_simp() -> bool:
    """Minimum footing depth for an axially loaded strip and pad footing
        with the simplified method.

    EN1992-1-1:2023 Eq. (14.14)

    Returns:
        float: the minimum hf/af ratio for a pad footing where hf is
            the footing depth and af is the distance from
            the footing edge to the column or wall face
            with the simplified method.
    """
    return 2.0
