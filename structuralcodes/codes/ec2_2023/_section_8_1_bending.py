from typing import Literal

import numpy as np
from numpy.typing import ArrayLike


def MEd_min(h: float, NEd: float) -> float:
    """Minimum eccentricity for effects of imperfections.

    EN1992-1-1:2023 Eq.(8.1).

    Computes the minimum moment for a determined h-height section for taking
    into consideration geometric imperfections unless second order effects are
    used.

    Args:
        h (float): Height of the element in mm.
        Ned (float): Axial force in kN.

    Returns:
        float: minimum design moment in kNm.
    """
    ed_min = max(h / 30, 20) / 1000
    return NEd * ed_min


def e_min(h: float) -> float:
    """Compute the minimum eccentricity for geometric imperfections.

    EN1992-1-1:2023 Eq.(8.1).

    Args:
        h (float): Height of the element in mm.

    Returns:
        float: Minimum eccentricity in m.
    """
    return max(h / 30, 20) / 1000


def NRd0(Ac: float, fcd: float, As: float, fyd: float) -> float:
    """Design value of axial resistance in compression.

    EN1992-1-1:2023 Eq. (8.3).

    Computes the design value of axial resistance under compression without
    accompanying moments.

    Args:
        Ac (float): Concrete area in mm2.
        fcd (float): Compressive design resistance of concrete in MPa. If
            confined concrete, then replace by fcd,c (8.15).
        As (float): Reinforcement area in mm2.
        fyd (float): Yield tensile resistance of steel in MPa.

    Returns:
        float: Axial resistance in compression in kN.

    Raises:
        ValueError: If any of Ac, fcd, As or Ayd is less than 0.
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

    EN1992-1-1:2023 Eq. (8.2).

    In the absence of an accurate cross-section design for biaxial beding this
    criterion may be used.

    Args:
        MEdz_MRdz (float): Ratio between the design bending moment and the
            resistance in the Z-axis.
        MEdy_MRdy (float): Ratio between the design bending moment and the
            resistance in the Y-axis.
        Ned_NRd (float): Ratio between the design axial force and the axial
            compressive resistance.
        section_type (str): The section geometry type.

    Returns:
        float: The resistance ratio (non-dimensional).
    """
    if section_type in ('elliptical', 'circular'):
        an = 2.0
    else:
        an = np.interp(
            Ned_NRd,
            xp=[0.1, 0.7, 2.0],
            fp=[1.0, 1.5, 2.0],
        )

    return abs(MEdz_MRdz) ** an + abs(MEdy_MRdy) ** an


def sigma_cd(fcd: float, eps_c: float) -> float:
    """Computes the stress distribution in the compression zones.

    EN1992-1-1:2023 Eq. (8.4).

    Computes the scress distribution in the compression zones (compressive
    shown as positive).

    Args:
        fcd (float): Compressive design resistance of concrete (MPa).
        eps_c (float): Strain value of concrete (non dimensional).

    Returns:
        float: Concrete stress in MPa.

    Raises:
        ValueError: If strain greater than eps_c_u=0.0035.
    """
    if eps_c <= 0:
        return 0.0
    if eps_c <= 0.002:
        return fcd * (1 - (1 - eps_c / 0.002) ** 2)
    if eps_c <= 0.0035:
        return fcd

    raise ValueError('Strain cannot be greater than eps_c_u=0.0035')


def delta_fcd_confined(sigma_c2d: float, f_cd: float, ddg: float) -> float:
    """Calculate the compressive strength increase (delta_f_cd) due to a
    transverse compressive stress.

    EN1992-1-1:2023 Eq. (8.9 and 8.10).

    Args:
        sigma_c2d (float): the absolute alue of the minum principal
            transverse compressive stress in MPa.
        f_cd (float): Compressive design strength in MPa.
        ddg (float): Maximum aggregate size in mm.

    Returns:
        float: Compressive strength increase in MPa.
    """
    if sigma_c2d < 0:
        raise ValueError(
            f'sigma_c2d must be non-negative. Got {sigma_c2d} instead.'
        )
    if f_cd < 0:
        raise ValueError(f'f_cd must be non-negative. Got {f_cd} instead.')
    if ddg < 0:
        raise ValueError(f'ddg must be non-negative. Got {ddg} instead.')

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
    """Calculate the confinement stress (sigma_c2d) for circular and square
    members with single confinement reinforcement.

    EN1992-1-1:2023 Eq (8.11).

    Args:
        A_s_conf (float): Cross-sectional area of one leg of confinement
            reinforcement in mm2.
        f_yd (float): Yield design strength of reinforcement in MPa.
        b_cs (float): Width of the confinement core in mm.
        s (float): Spacing of confinement reinforcement in mm.

    Returns:
        float: Confinement stress in MPa.
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
    """Calculate the confinement stress (sigma_c2d) for rectangular members
    with single confinement reinforcement.

    EN1992-1-1:2023 Eq. (8.12).

    Args:
        A_s_conf (float): Cross-sectional area of one leg of confinement
            reinforcement.
        f_yd (float): Yield strength of reinforcement.
        b_csx (float): Width of the confinement core in x direction.
        b_csy (float): Width of the confinement core in y direction.
        s (float): Spacing of confinement reinforcement.

    Returns:
        float: Confinement stress.
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
    A_s_confx: ArrayLike,
    A_s_confy: ArrayLike,
    f_yd: float,
    b_csx: float,
    b_csy: float,
    s: float,
) -> float:
    """Calculate the confinement stress (sigma_c2d) for members with multiple
    confinement reinforcement.

    EN1992-1-1:2023 Eq. (8.13).

    Args:
        A_s_confx (ArrayLike): Cross-sectional areas of confinement
            reinforcement in x direction in mm2.
        A_s_confy (ArrayLike): Cross-sectional areas of confinement
            reinforcement in y direction in mm2.
        f_yd (float): Yield strength of reinforcement in MPa.
        b_csx (float): Width of the confinement core in x direction in mm.
        b_csy (float): Width of the confinement core in y direction in mm.
        s (float): Spacing of confinement reinforcement in mm.

    Returns:
        float: Confinement stress in mm.
    """
    A_s_confx = np.atleast_1d(A_s_confx)
    A_s_confy = np.atleast_1d(A_s_confy)

    if np.any(A_s_confx < 0):
        raise ValueError(f'A_s_confx must be non-negative. Got {A_s_confx}')
    if np.any(A_s_confy < 0):
        raise ValueError(f'A_s_confy must be non-negative. Got {A_s_confy}')

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
    A_s_confx: ArrayLike,
    A_s_confy: ArrayLike,
    f_yd: float,
    b_csy: float,
    x_cs: float,
    s: float,
) -> float:
    """Calculate the confinement stress (sigma_c2d) for compression zones.

    EN1992-1-1:2023 Eq. (8.14).

    Args:
        A_s_confx (ArrayLike): Cross-sectional area of confinement
            reinforcement in x direction in mm2.
        A_s_confy (ArrayLike): Cross-sectional area of confinement
            reinforcement in y direction in mm2.
        f_yd (float): Yield strength of reinforcement in MPa.
        b_csy (float): Width of the confinement core in mm.
        x_cs (float): Width of the confinement core in mm.
        s (float): Spacing of confinement reinforcement in mm.

    Returns:
        float: Confinement stress in mm2.
    """
    A_s_confx = np.atleast_1d(A_s_confx)
    A_s_confy = np.atleast_1d(A_s_confy)

    if np.any(A_s_confx < 0):
        raise ValueError(f'A_s_confx must be non-negative. Got {A_s_confx}')
    if np.any(A_s_confy < 0):
        raise ValueError(f'A_s_confy must be non-negative. Got {A_s_confy}')

    if f_yd < 0:
        raise ValueError(f'f_yd must be non-negative. Got {f_yd} instead.')
    if x_cs < 0:
        raise ValueError(f'x_cs must be non-negative. Got {x_cs} instead.')
    if b_csy < 0:
        raise ValueError(f'b_csy must be non-negative. Got {b_csy} instead.')
    if s < 0:
        raise ValueError(f's must be non-negative. Got {s} instead.')

    return min((sum(A_s_confx) / b_csy), (sum(A_s_confy) / x_cs)) * f_yd / s


def fcd_c(
    fcd: float, kconf_b: float, kconf_s: float, delta_fcd: float
) -> float:
    """Calculate the average concrete strength increase in the confined areas.

    EN1992-1-1:2023 Eq. (8.15).

    Args:
        fcd (float): The concrete design strength in MPa.
        kconf_b (float): Effectiveness factor for the shape of the compression
            zone and confinement reinforcement (non-dimensional).
        kconf_s (float): Effectiveness factor for the spacing of the
            confinement reinforcement (non-dimensional).
        delta_fcd (float): The increase in concrete design strength due to
            confinement in MPa.

    Returns:
        float: The average concrete strength increase in MPa.
    """
    if fcd < 0:
        raise ValueError(f'fcd must be non-negative. Got {fcd} instead.')
    if kconf_b < 0:
        raise ValueError(
            f'kconf_b must be non-negative. Got {kconf_b} instead.'
        )
    if kconf_s < 0:
        raise ValueError(
            f'kconf_s must be non-negative. Got {kconf_s} instead.'
        )
    if delta_fcd < 0:
        raise ValueError(
            f'delta_fcd must be non-negative. Got {delta_fcd} instead.'
        )

    return fcd + kconf_b * kconf_s * delta_fcd


def kconf_b_square_single(bcs: float, b: float) -> float:
    """Calculate the kconf_b effectiveness factor for square members in
    compression with single confinement reinforcement.

    EN1992-1-1:2023 Table (8.1).

    Args:
        bcs (float): Width of the confined section in mm.
        b (float): Width of the unconfined section in mm.

    Returns:
        float: The effectiveness factor kconf_b (non-dimensional).
    """
    if bcs < 0:
        raise ValueError(f'bcs must be non-negative. Got {bcs} instead.')
    if b < 0:
        raise ValueError(f'b must be non-negative. Got {b} instead.')

    return (1 / 3) * (bcs / b) ** 2


def kconf_s_square_single(s: float, bcs: float) -> float:
    """Calculate the kconf_s effectiveness factor for square members in
    compression with single confinement reinforcement.

    EN1992-1-1:2023 Table (8.1).

    Args:
        s (float): Spacing of the confinement reinforcement in mm.
        bcs (float): Width of the confined section in mm.

    Returns:
        float: The effectiveness factor kconf_s (non-dimensional).
    """
    if s < 0:
        raise ValueError(f's must be non-negative. Got {s} instead.')
    if bcs < 0:
        raise ValueError(f'bcs must be non-negative. Got {bcs} instead.')

    return (max(1 - s / (2 * bcs), 0)) ** 2


def kconf_b_circular(bcs: float, b: float) -> float:
    """Calculate the kconf_b effectiveness factor for circular members in
    compression with circular confinement reinforcement.

    EN1992-1-1:2023 Table (8.1).

    Args:
        bcs (float): Diameter of the confined section in mm.
        b (float): Diameter of the unconfined section in mm.

    Returns:
        float: The effectiveness factor kconf_b (non-dimensional).
    """
    if bcs < 0:
        raise ValueError(f'bcs must be non-negative. Got {bcs} instead.')
    if b < 0:
        raise ValueError(f'b must be non-negative. Got {b} instead.')

    return (bcs / b) ** 2


def kconf_b_multiple(
    bcsx: float, bcsy: float, b_i: ArrayLike, bx: float, by: float
) -> float:
    """Calculate the kconf_b effectiveness factor for square and rectangular
    members in compression with multiple confinement reinforcement.

    EN1992-1-1:2023 Table (8.1).

    Args:
        bcsx (float): Width of the confined section in x-direction in mm.
        bcsy (float): Width of the confined section in y-direction in mm.
        b_i (ArrayLike): Distances between bends of straight segments in mm.
        bx (float): Width of the unconfined section in x-direction in mm.
        by (float): Width of the unconfined section in y-direction in mm.

    Returns:
        float: The effectiveness factor kconf_b.
    """
    b_i = np.atleast_1d(b_i)
    if bcsx < 0:
        raise ValueError(f'bcsx must be non-negative. Got {bcsx} instead.')
    if bcsy < 0:
        raise ValueError(f'bcsy must be non-negative. Got {bcsy} instead.')
    if np.any(b_i < 0):
        raise ValueError(f'b_i must be non-negative. Got {b_i} instead.')
    if bx < 0:
        raise ValueError(f'bx must be non-negative. Got {bx} instead.')
    if by < 0:
        raise ValueError(f'by must be non-negative. Got {by} instead.')

    sq_sum = np.sum(b_i**2)

    return (bcsx * bcsy - (1 / 6) * sq_sum) / (bx * by)


def kconf_s_multiple(s: float, bcsx: float, bcsy: float) -> float:
    """Calculate the kconf_s effectiveness factor for square and rectangular
    members in compression with multiple confinement reinforcement.

    EN1992-1-1:2023 Table (8.1).

    Args:
        s (float): Spacing of the confinement reinforcement in mm.
        bcsx (float): Width of the confined section in x-direction in mm.
        bcsy (float): Width of the confined section in y-direction in mm.

    Returns:
        float: The effectiveness factor kconf_s.
    """
    if s < 0:
        raise ValueError(f's must be non-negative. Got {s} instead.')
    if bcsx < 0:
        raise ValueError(f'bcsx must be non-negative. Got {bcsx} instead.')
    if bcsy < 0:
        raise ValueError(f'bcsy must be non-negative. Got {bcsy} instead.')

    return max((1 - s / (2 * bcsx)), 0) * max((1 - s / (2 * bcsy)), 0)


def kconf_b_bending(Ac_conf: float, Acc: float, b_i: ArrayLike) -> float:
    """Calculate the kconf_b effectiveness factor for compression zones due to
    bending and axial force.

    EN1992-1-1:2023 Table (8.1).

    Args:
        Ac_conf (float): Confined area within the centrelines of the
            confinement reinforcement and the neutral axis in mm2.
        Acc (float): Compressive area mm2.
        b_i (ArrayLike): Distances between bends of straight segments in mm.

    Returns:
        float: The effectiveness factor kconf_b.
    """
    if Ac_conf < 0:
        raise ValueError(
            f'Ac_conf must be non-negative. Got {Ac_conf} instead.'
        )
    if Acc < 0:
        raise ValueError(f'Acc must be non-negative. Got {Acc} instead.')
    b_i = np.atleast_1d(b_i)
    if np.any(b_i < 0):
        raise ValueError(f'b_i must be non-negative. Got {b_i} instead.')

    sq_sum = np.sum(b_i**2)

    return (Ac_conf - (1 / 6) * sq_sum) / Acc


def kconf_s_bending(s: float, xcs: float, bcsx: float, bcsy: float) -> float:
    """Calculate the kconf_s effectiveness factor for compression zones
        under bending and axial force.

    Based on EN1992-1-1:2023 Table (8.1).

    Args:
        s (float): Spacing of the confinement reinforcement.
        xcs (float): Shortest perpendicular distance from the confined
            reinforcement to the neutral axis.
        bcsx (float): Width of the total confined section.
        bcsy (float): Width of the confined section
            delimited by the neutral axis and confinement reinforcement.

    Returns:
        float: The effectiveness factor kconf_s (dimensionless).
    """
    if s < 0:
        raise ValueError(f's must be non-negative. Got {s} instead.')
    if xcs < 0:
        raise ValueError(f'xcs must be non-negative. Got {xcs} instead.')
    if bcsx < 0:
        raise ValueError(f'bcsx must be non-negative. Got {bcsx} instead.')
    if bcsy < 0:
        raise ValueError(f'bcsy must be non-negative. Got {bcsy} instead.')

    xcs = min(xcs, bcsx / 2)
    return max((1 - s / (4 * xcs)), 0) * max((1 - s / (2 * bcsy)), 0)


def epsc2_c(epsc2: float, delta_fcd: float, fcd: float) -> float:
    """Calculate the confined concrete strain limit at maximum stress.

    EN1992-1-1:2023 Eq. (8.16).

    Args:
        epsc2 (float): The strain limit at maximum stress for unconfined
            concrete.
        delta_fcd (float): The increase in concrete design strength due to
            confinement in MPa.
        fcd (float): The concrete design strength in MPa.

    Returns:
        float: The confined concrete strain limit at maximum stress.
    """
    if epsc2 < 0:
        raise ValueError(f'epsc2 must be non-negative. Got {epsc2} instead.')
    if delta_fcd < 0:
        raise ValueError(
            f'delta_fcd must be non-negative. Got {delta_fcd} instead.'
        )
    if fcd < 0:
        raise ValueError(f'fcd must be non-negative. Got {fcd} instead.')

    return epsc2 * (1 + 5 * delta_fcd / fcd)


def epscu_c(epscu: float, sigma_c2d: float, fcd: float) -> float:
    """Calculate the confined concrete ultimate strain limit.

    EN1992-1-1:2023 Eq. (8.17).

    Args:
        epscu (float): The ultimate strain limit for unconfined concrete.
        sigma_c2d (float): The stress at maximum strain for confined concrete
            in MPa.
        fcd (float): The concrete design strength in MPa.

    Returns:
        float: The confined concrete ultimate strain limit.
    """
    if epscu < 0:
        raise ValueError(f'epscu must be non-negative. Got {epscu} instead.')
    if sigma_c2d < 0:
        raise ValueError(
            f'sigma_c2d must be non-negative. Got {sigma_c2d} instead.'
        )
    if fcd < 0:
        raise ValueError(f'fcd must be non-negative. Got {fcd} instead.')

    return epscu + 0.2 * sigma_c2d / fcd
