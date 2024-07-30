from typing import Literal

from scipy.interpolate import interp1d


def MEd_min(h: float, NEd: float) -> float:
    """Minimum eccentricity for effects of imperfections.

    EN1992-1-1:2023 Eq.(8.1)

    Computes the minimum moment for a determined h-height section
    for taking into consideration geometric imperfections unless
    second order effects are used.

    Args:
        h (float): height of the element in mm.
        Ned (float): axial force in kN.

    Returns:
        float: minimum design moment in kN·m
    """
    ed_min = max(h / 30, 20) / 1000
    return NEd * ed_min


def NRd0(Ac: float, fcd: float, As: float, fyd: float) -> float:
    """Design value of axial resistance in compression.

    EN1992-1-1:2023 Eq. (8.3)

    Computes the design value of axial resistance under compression
    without accompanying moments.

    Args:
        Ac (float): concrete area in mm2
        fcd (float): compressive design resistance of concrete in MPa.
            If confined concrete, then replace by fcd,c (8.15)
        As (float): reinforcement area in mm2
        fyd (float): yield tensile resistance of steel in MPa

    Returns:
        float: axial resistance in compression in kN

    Raise:
        ValueError: if any of Ac, fcd, As or Ayd is less than 0.
    """
    if Ac < 0:
        raise ValueError('Concrete area cannot be negative')
    if fcd < 0:
        raise ValueError('Concrete resistance cannot be negative')
    if As < 0:
        raise ValueError('Steel area cannot be negative')
    if fyd < 0:
        raise ValueError('Steel resistance cannot be negative')

    return (Ac * fcd + As * fyd) / 1000


def biaxial_resistant_ratio(
    MEdz_MRdz: float,
    MEdy_MRdy: float,
    Ned_NRd: float,
    section_type: Literal['rectangular', 'circular', 'elliptical'],
) -> float:
    """Computes the resistant ratio in biaxial bending.

    EN1992-1-1:2023 Eq. (8.2)

    In the absence of an accurate cross-section design for biaxial beding
    this criterion may be used.

    Args:
        MEdz_MRdz (float): ratio between the design bending moment and the
            resistance in the Z-axis.
        MEdy_MRdy (float): ratio between the design bending moment and the
            resistance in the Y-axis.
        Ned_NRd (float): ratio between the design axial force and the
            axial compressive resistance.
        section_type (str): the section geometry type.

    Returns:
        float: the resistance ratio (non-dimensional).
    """
    if section_type in ('elliptical', 'circular'):
        an = 2.0
    elif Ned_NRd <= 0.1:
        an = 1.0
    elif Ned_NRd >= 1.0:
        an = 2.0
    else:
        an = interp1d([0.1, 0.7, 2.0], [1.0, 1.5, 2.0])(Ned_NRd)

    return abs(MEdz_MRdz) ** an + abs(MEdy_MRdy) ** an


def sigma_cd(fcd: float, eps_c: float) -> float:
    """Computes the stress distribution in the compression zones.

    EN1992-1-1:2023 Eq. (8.4)

    Computes the scress distribution in the cmpressiopn zones (compressive
    shown as positive).

    Args:
        fcd (float): compressive design resistance of concrete (MPa)
        eps_c (float): strain value of concrete (non dimensional)

    Returns:
        float: concrete stress in MPa

    Raises:
        ValueError: if strain greater than eps_c_u=0.0035
    """
    if eps_c <= 0:
        return 0.0
    if eps_c <= 0.002:
        return fcd * (1 - (1 - eps_c / 0.002) ** 2)
    if eps_c <= 0.0035:
        return fcd

    raise ValueError('Strain cannot be greater than eps_c_u=0.0035')
