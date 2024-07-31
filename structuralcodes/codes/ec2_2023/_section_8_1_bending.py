from typing import List, Literal

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


def delta_fcd_confined(sigma_c2d: float, f_cd: float, ddg: float) -> float:
    """Calculate the compressive strength increase (Δf_cd) due to a
       transverse compressive stress.

    EN1992-1-1:2023 Eq. (8.9 and 8.10)

    Parameters:
        sigma_c2d (float): Transverse compressive stress in MPa
        f_cd (float): Compressive design strength in MPa
        ddg (float): Maximum aggregate size in mm

    Returns:
        float: Compressive strength increase in MPa
    """
    if sigma_c2d < 0:
        raise ValueError(
            f'signma_c2d must be non-negative. Got {sigma_c2d} instead.'
        )
    if f_cd < 0:
        raise ValueError(
            f'f_cd must be non-negative. Got {sigma_c2d} instead.'
        )
    if ddg < 0:
        raise ValueError(f'ddg must be non-negative. Got {sigma_c2d} instead.')

    if sigma_c2d <= 0.6 * f_cd:
        delta_f_cd = 4 * sigma_c2d
    else:
        delta_f_cd = 3.5 * sigma_c2d ** (3 / 4) * f_cd ** (1 / 4)

    if ddg < 32:
        delta_f_cd *= ddg / 32

    return delta_f_cd


def confinement_sigma_c2d_circular_square(
    A_s_conf: float, f_yd: float, b_cs: float, s: float
) -> float:
    """Calculate the confinement stress (σc2d) for circular and square members
    with single confinement reinforcement.

    EN1992-1-1:2023 Eq (8.11)

    Args:
        A_s_conf (float): Cross-sectional area of one leg of confinement
            reinforcement in mm2
        f_yd (float): Yield strength of reinforcement in MPa
        b_cs (float): Width of the confinement core in mm
        s (float): Spacing of confinement reinforcement in mm

    Returns:
        float: Confinement stress in MPa
    """
    if A_s_conf < 0:
        raise ValueError(
            f'A_s_conf must be non-negative. Got {A_s_conf} instead.'
        )
    if f_yd < 0:
        raise ValueError(f'f_yd must be non-negative. Got {f_yd} instead.')
    if b_cs < 0:
        raise ValueError(f'b_cs must be non-negative. Got {b_cs} instead.')
    if s < 0:
        raise ValueError(f's must be non-negative. Got {s} instead.')

    return 2 * A_s_conf * f_yd / (b_cs * s)


def confinement_sigma_c2d_rectangular(
    A_s_conf: float, f_yd: float, b_csx: float, b_csy: float, s: float
) -> float:
    """Calculate the confinement stress (σc2d)
        for rectangular members with single confinement reinforcement.

    EN1992-1-1:2023 Eq. (8.12)

    Args:
        A_s_conf (float): Cross-sectional area of one leg of confinement
        reinforcement
        f_yd (float): Yield strength of reinforcement
        b_csx (float): Width of the confinement core in x direction
        b_csy (float): Width of the confinement core in y direction
        s (float): Spacing of confinement reinforcement

    Returns:
        float: Confinement stress
    """
    if A_s_conf < 0:
        raise ValueError(
            f'A_s_conf must be non-negative. Got {A_s_conf} instead.'
        )
    if f_yd < 0:
        raise ValueError(f'f_yd must be non-negative. Got {f_yd} instead.')
    if b_csx < 0:
        raise ValueError(f'b_csx must be non-negative. Got {b_csx} instead.')
    if b_csy < 0:
        raise ValueError(f'b_csy must be non-negative. Got {b_csy} instead.')
    if s < 0:
        raise ValueError(f's must be non-negative. Got {s} instead.')

    return 2 * A_s_conf * f_yd / (max(b_csx, b_csy) * s)


def confinement_sigma_c2d_multiple(
    A_s_confx: List[float],
    A_s_confy: List[float],
    f_yd: float,
    b_csx: float,
    b_csy: float,
    s: float,
) -> float:
    """Calculate the confinement stress (σc2d) for members
        with multiple confinement reinforcement.

    EN1992-1-1:2023 Eq. (8.13)

    Parameters:
        A_s_confx (List[float]): Cross-sectional areas of confinement
            reinforcement in x direction in mm2
        A_s_confy (List[float]): Cross-sectional areas of confinement
            reinforcement in y direction in mm2
        f_yd (float): Yield strength of reinforcement in MPa
        b_csx (float): Width of the confinement core in x direction in mm
        b_csy (float): Width of the confinement core in y direction in mm
        s (float): Spacing of confinement reinforcement in mm

    Returns:
        float: Confinement stress in mm
    """
    for value in A_s_confx:
        if value < 0:
            raise ValueError(
                f'A_s_confx must be non-negative. Got {value} instead.'
            )
    for value in A_s_confy:
        if value < 0:
            raise ValueError(
                f'A_s_confy must be non-negative. Got {value} instead.'
            )
    if f_yd < 0:
        raise ValueError(f'f_yd must be non-negative. Got {f_yd} instead.')
    if b_csx < 0:
        raise ValueError(f'b_csx must be non-negative. Got {b_csx} instead.')
    if b_csy < 0:
        raise ValueError(f'b_csy must be non-negative. Got {b_csy} instead.')
    if s < 0:
        raise ValueError(f's must be non-negative. Got {s} instead.')

    return min((sum(A_s_confx) / b_csy), (sum(A_s_confy) / b_csx)) * f_yd / s


def confinement_sigma_c2d_compression_zones(
    A_s_confx: List[float],
    A_s_confy: List[float],
    f_yd: float,
    b_csy: float,
    x_cs: float,
    s: float,
) -> float:
    """Calculate the confinement stress (σc2d) for compression zones.

    EN1992-1-1:2023 Eq. (8.14)

    Parameters:
        A_s_confx (List[float]): Cross-sectional area of confinement
            reinforcement in x direction in mm2
        A_s_confy (List[float]): Cross-sectional area of confinement
            reinforcement in y direction in mm2
        f_yd (float): Yield strength of reinforcement in MPa
        b_csy (float): Width of the confinement core in y direction in mm
        x_cs (float): Width of the confinement core in x direction in mm
        s (float): Spacing of confinement reinforcement in mm

    Returns:
        float: Confinement stress in mm2
    """
    for value in A_s_confx:
        if value < 0:
            raise ValueError(
                f'A_s_confx must be non-negative. Got {value} instead.'
            )
    for value in A_s_confy:
        if value < 0:
            raise ValueError(
                f'A_s_confy must be non-negative. Got {value} instead.'
            )
    if f_yd < 0:
        raise ValueError(f'f_yd must be non-negative. Got {f_yd} instead.')
    if x_cs < 0:
        raise ValueError(f'x_cs must be non-negative. Got {x_cs} instead.')
    if b_csy < 0:
        raise ValueError(f'b_csy must be non-negative. Got {b_csy} instead.')
    if s < 0:
        raise ValueError(f's must be non-negative. Got {s} instead.')

    return min((sum(A_s_confx) / b_csy), (sum(A_s_confy) / x_cs)) * f_yd / s
