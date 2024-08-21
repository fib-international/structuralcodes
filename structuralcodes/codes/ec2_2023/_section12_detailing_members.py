"""Functions from Section 12 of EN 1992-1-1:2023."""

import math
from typing import Literal


def As_min(Ac: float, fctm: float, fyk: float) -> float:
    """Calculate the minimum reinforcement area
        for members subjected to pure tension.

    EN1992-1-1:2023 Eq. (12.2)

    Args:
        Ac (float): The concrete cross-sectional area in mm2.
        fctm (float): The mean value of axial tensile
            strength of the concrete in MPa.
        fyk (float): The characteristic yield
            strength of the reinforcement in MPa.

    Returns:
        float: The minimum reinforcement area As,min in mm2.

    Raises:
        ValueError: If any input value is negative.
    """
    if Ac < 0:
        raise ValueError(
            'Ac (concrete cross-sectional area) must not be negative. '
            + f'Got {Ac}'
        )
    if fctm < 0:
        raise ValueError(
            'fctm (tensile strength of concrete) must not be negative.'
            + f' Got {fctm}'
        )
    if fyk < 0:
        raise ValueError(
            'fyk (yield strength of reinforcement) must not be negative. '
            + f'Got {fyk}'
        )

    return Ac * fctm / fyk


def As_w_min(
    fck: float,
    fyk: float,
    alpha: float,
    bw: float,
    s: float,
    ductility_class: Literal['A', 'B', 'C'],
) -> float:
    """Calculate the minimum shear reinforcement.

    EN1992-1-1:2023 Eq. (12.4)

    Args:
        fck (float): The characteristic compressive
            strength of concrete in MPa.
        fyk (float): The characteristic yield strength
            of reinforcement in MPa.
        alpha (float): The angle between shear
            reinforcement and the longitudinal axis in degrees.
        bw (float): The width of the
            web of the member in mm.
        s (float): The spacing of the shear reinforcement
            along the longitudinal axis in mm.
        ductility_class (str): The ductility class
            of the reinforcement ('A', 'B', or 'C').

    Returns:
        float: The minimum shear reinforcement in mm2.

    Raises:
        ValueError: If any input value is negative.
    """
    if fck < 0:
        raise ValueError(
            'fck (compressive strength of concrete) '
            + f'must not be negative. Got {fck}'
        )
    if fyk < 0:
        raise ValueError(
            'f_yk (yield strength of reinforcement) '
            + f'must not be negative. Got {fyk}'
        )
    if alpha < 0 or alpha > 90:
        raise ValueError(
            f'alpha (angle) must be between 0 and 90 degrees. Got {alpha}'
        )
    if bw < 0:
        raise ValueError(f'b_w (width of web) must not be negative. Got {bw}')
    if s < 0:
        raise ValueError(
            f's (spacing of shear reinforcement) must not be negative. Got {s}'
        )

    alpha_rad = math.radians(alpha)
    rho_w_min = 0.08 * math.sqrt(fck) / fyk
    As_w_min = rho_w_min * s * bw * math.sin(alpha_rad)

    if ductility_class == 'B':
        As_w_min *= 0.9
    elif ductility_class == 'C':
        As_w_min *= 0.8

    return As_w_min


def MRd_min(
    Med: float,
    ductility_class: Literal['A', 'B', 'C'],
) -> float:
    """Calculate the minimum moment resistance.

    EN1992-1-1:2023 Eq. (12.3)

    Args:
        Med (float): The design moment in kNm.
        ductility_class (str): The ductility class of the reinforcement
            (A, B, or C).

    Returns:
        float: The minimum moment resistance M_rd,min in kN·m.
    """
    kdc_values = {'A': 1.3, 'B': 1.1, 'C': 1.0}
    kdc = kdc_values[ductility_class]

    return kdc * Med


def As_min_bot(As_min_req_span: float) -> float:
    """Calculate the minimum bottom reinforcement at inner and end supports
    taking account unforeseen postive moments.

    EN1992-1-1:2023 Table (12.1)

    Args:
        As_min_req_span (float): the minimum reinforcement area requirement in
            the span of the beam in mm2.

    Returns:
        float: The minimum reinforcement area in mm2
    """
    if As_min_req_span < 0:
        raise ValueError(
            f'As_min_req_span must not be negative. Got {As_min_req_span}'
        )

    return 0.25 * As_min_req_span


def sl_max(d: float, alpha: float) -> float:
    """Calculate maximum longitudinal spacing of shear assemblies/stirrups.

    EN1992-1-1:2023 Table (12.1)

    Args:
        d (float): Effective depth in mm.
        alpha (float): Angle of shear reinforcement in degrees.

    Returns:
        float: Maximum spacing sl,max in mm.

    Raises:
        ValueError: If d is negative or alpha is out of range.
    """
    if d < 0:
        raise ValueError(f'd must not be negative. Got {d}')
    if not (0 <= alpha <= 90):
        raise ValueError(
            f'alpha must be between 0 and 90 degrees. Got {alpha}'
        )

    rad_alpha = math.radians(alpha)
    cot_alpha = 1.0 / math.tan(rad_alpha)

    return 0.75 * d * (1 + cot_alpha)


def sbu_max(d: float, alpha: float) -> float:
    """Calculate maximum longitudinal spacing of bent-up bars.

    EN1992-1-1:2023 Table (12.1)

    Args:
        d (float): Effective depth in mm.
        alpha (float): Angle of shear reinforcement in degrees.

    Returns:
        float: Maximum spacing sbu,max in mm.

    Raises:
        ValueError: If d is negative or alpha is out of range.
    """
    if d < 0:
        raise ValueError(f'd must not be negative. Got {d}')
    if not (0 <= alpha <= 90):
        raise ValueError(
            f'alpha must be between 0 and 90 degrees. Got {alpha}'
        )

    rad_alpha = math.radians(alpha)
    cot_alpha = 1.0 / math.tan(rad_alpha)

    return 0.6 * d * (1 + cot_alpha)


def str_max(d: float) -> float:
    """Calculate maximum transverse spacing of shear legs.

    EN1992-1-1:2023 Table (12.1)

    Args:
        d (float): Effective depth in mm.

    Returns:
        float: Maximum transverse spacing str,max ion mm.

    Raises:
        ValueError: If d is negative.
    """
    if d < 0:
        raise ValueError(f'd must not be negative. Got {d}')

    return min(0.75 * d, 600)


def rho_w_stir_min(rho_w_req: float) -> float:
    """Calculate the minimum ratio of shear reinforcement in the
        form of stirrups with respect to the required reinforcement ratio.

    EN1992-1-1:2023 Table (12.1)

    Args:
        rho_w_req (float): Required reinforcement ratio.

    Returns:
        float: Minimum reinforcement ratio ρw,stir.

    Raises:
        ValueError: If rho_w_req is negative.
    """
    if rho_w_req < 0:
        raise ValueError(f'rho_w_req must not be negative. Got {rho_w_req}')

    return 0.5 * rho_w_req


def rho_w_tors_min(rho_w_req: float) -> float:
    """Calculate the minimum ratio of torsion reinforcement in the
        form of closed stirrups with respect to the
        required reinforcement ratio.

    EN1992-1-1:2023 Table (12.1)

    Args:
        rho_w_req (float): Required reinforcement ratio.

    Returns:
        float: Minimum torsion reinforcement ratio ρw,stir.

    Raises:
        ValueError: If rho_w_req is negative.
    """
    if rho_w_req < 0:
        raise ValueError(f'rho_w_req must not be negative. Got {rho_w_req}')

    return 0.2 * rho_w_req


def s_stir_max(u: float, b: float, h: float) -> float:
    """Calculate the maximum spacing for torsion assemblies/stirrups.

    EN1992-1-1:2023 Table (12.1)

    Args:
        u (float): Perimeter of the section in mm.
        b (float): Width of the section in mm.
        h (float): Height of the section in mm.

    Returns:
        float: Maximum spacing for torsion assemblies/stirrups sstir,max in mm.

    Raises:
        ValueError: If any of the inputs are negative.
    """
    if u < 0:
        raise ValueError(f'u must not be negative. Got {u}')
    if b < 0:
        raise ValueError(f'b must not be negative. Got {b}')
    if h < 0:
        raise ValueError(f'h must not be negative. Got {h}')

    return min(u / 8, b, h)


def sl_surf_max() -> float:
    """Calculate the minimum spacing of longitudinal surface reinforcement.

    EN1992-1-1:2023 Table (12.1)

    Returns:
        float: Minimum spacing in mm.
    """
    return 300.0


def al_with_shear_reinforcement(z: float, theta: float, alpha: float) -> float:
    """Calculate the shift distance al for members
        with shear reinforcement using the given parameters.

    EN1992-1-1:2023 Eq. (12.5)

    Args:
        z (float): Internal lever arm in mm.
        theta (float): Angle of compression strut in degrees.
        alpha (float): Angle of shear reinforcement in degrees.

    Returns:
        float: Shift distance al in mm.

    Raises:
        ValueError: If z, theta, or alpha are out of their valid ranges.
    """
    if z < 0:
        raise ValueError(f'z must not be negative. Got {z}')
    if not (0 <= theta <= alpha):
        raise ValueError(
            f'theta must be between 0 and 90 degrees. Got {theta}'
        )

    rad_theta = math.radians(theta)
    rad_alpha = math.radians(alpha)
    cot_theta = 1 / math.tan(rad_theta)
    cot_alpha = 1 / math.tan(rad_alpha)

    return z * (cot_theta - cot_alpha) / 2


def al_without_shear_reinforcement(d: float) -> float:
    """Calculate the shift distance al for members without shear reinforcement.

    EN1992-1-1:2023 Eq. (12.6)

    Args:
        d (float): Effective depth in mm.

    Returns:
        float: Shift distance al in mm.

    Raises:
        ValueError: If d is negative.

    """
    if d < 0:
        raise ValueError(f'd must not be negative. Got {d}')
    return d


def min_anchorage_length_bent_up_bar(
    lbd: float, zone: Literal['tension', 'compression']
) -> float:
    """Calculate the anchorage length of a bent-up bar
      contributing to shear resistance.

    EN1992-1-1:2023 12.3.2(4)

    Args:
        lbd (float): Basic anchorage length in mm.
        zone (str): The zone in which the bar is anchored,
            either "tension" or "compression".

    Returns:
        float: Required anchorage length ion mm.

    Raises:
        ValueError: If `lbd` is negative.
    """
    if lbd < 0:
        raise ValueError(f'lbd must not be negative. Got {lbd}')

    if zone == 'tension':
        return 1.3 * lbd
    return 0.7 * lbd


def bottom_reinforcement_extension(phi: float) -> float:
    """Calculate the minimum extension of bottom reinforcement
        at intermediate supports.

    EN1992-1-1:2023 12.3.2(5)

    Args:
        phi (float): Diameter of the reinforcement bar in mm.

    Returns:
        float: Minimum extension length in mm.

    Raises:
        ValueError: If phi is negative.
    """
    if phi < 0:
        raise ValueError(f'phi must not be negative. Got {phi}')
    return 10 * phi


def As_slab_min_secondary_reinforcement(As_req_span: float) -> float:
    """Calculate the minimum secondary reinforcement in slabs.

    EN1992-1-1:2023 Table (12.2)

    Args:
        As_req_span (float): Required reinforcement for positive
            bending moments at the span in mm2.

    Returns:
        float: Minimum secondary reinforcement area in mm2.

    Raises:
        ValueError: If As_req_span is negative.

    """
    if As_req_span < 0:
        raise ValueError(
            f'As_req_span must not be negative. Got {As_req_span}'
        )

    return 0.2 * As_req_span


def As_slab_min_bottom_reinforcement_inner_supports(
    As_req_span: float,
) -> float:
    """Calculate the minimum longitudinal bottom
        reinforcement at inner supports of slabs.

    EN1992-1-1:2023 Table (12.2)

    Args:
        As_req_span (float): Required reinforcement for positive
            bending moments at the span in mm2.

    Returns:
        float: Minimum bottom reinforcement area in mm2.

    Raises:
        ValueError: If As_req_span is negative.
    """
    if As_req_span < 0:
        raise ValueError(
            f'As_req_span must not be negative. Got {As_req_span}'
        )

    return 0.25 * As_req_span


def As_slab_min_bottom_reinforcement_end_supports(As_req_span: float) -> float:
    """Calculate the minimum longitudinal
         bottom reinforcement at end supports of slabs.

    EN1992-1-1:2023 Table (12.2)

    Args:
        As_req_span (float): Required reinforcement
            for positive bending moments at the span in mm2.

    Returns:
        float: Minimum bottom reinforcement area in mm2.

    Raises:
        ValueError: If As_req_span is negative.
    """
    if As_req_span < 0:
        raise ValueError(
            f'As_req_span must not be negative. Got {As_req_span}'
        )

    return 0.25 * As_req_span


def As_slab_min_top_reinforcement_end_supports(
    As_req_span: float, As_min: float
) -> float:
    """Calculate the minimum top reinforcement at
        end supports in buildings where unintentional restraint may occur.

    EN1992-1-1:2023 Table (12.2)

    Args:
        As_req_span (float): Required reinforcement for
            positive bending moments at the span in mm2.
        As_min (float): Minimum area of reinforcement
            according to 12.2(2) in mm2.

    Returns:
        float: Minimum top reinforcement area in mm2.

    Raises:
        ValueError: If As_req_span or As_min are negative
    )
    """
    if As_req_span < 0:
        raise ValueError(
            f'As_req_span must not be negative. Got {As_req_span}'
        )
    if As_min < 0:
        raise ValueError(f'As_min must not be negative. Got {As_min}')

    return max(0.25 * As_req_span, As_min)


def s_slab_max(h: float) -> float:
    """Calculate the maximum spacing of bars for concrete in tension in slabs.

    EN1992-1-1:2023 Table (12.2)

    Args:
        h (float): Thickness of the slab in mm.

    Returns:
        float: Maximum spacing of bars in mm.

    Raises:
        ValueError: If h is negative.
    """
    if h < 0:
        raise ValueError(f'h must not be negative. Got {h}')

    return min(3 * h, 400)


def sl_max_slab(d: float, alpha: float) -> float:
    """Calculate the maximum longitudinal spacing of
        shear assemblies/stirrups in slabs.

    EN1992-1-1:2023 Table (12.2)

    Args:
        d (float): Effective depth of the slab in mm.
        alpha (float): Angle of shear reinforcement in degrees.

    Returns:
        float: Maximum longitudinal spacing sl,max in mm.

    Raises:
        ValueError: If d is negative or alpha is out of range.

    """
    if d < 0:
        raise ValueError(f'd must not be negative. Got {d}')
    if not (0 <= alpha <= 90):
        raise ValueError(
            f'alpha must be between 0 and 90 degrees. Got {alpha}'
        )

    rad_alpha = math.radians(alpha)
    cot_alpha = 1.0 / math.tan(rad_alpha)

    return 0.75 * d * (1 + cot_alpha)


def sbu_max_slab(d: float) -> float:
    """Calculate the maximum longitudinal spacing of bent-up bars in slabs.

    EN1992-1-1:2023 Table (12.2)

    Args:
        d (float): Effective depth of the slab in mm.

    Returns:
        float: Maximum longitudinal spacing of bent-up bars sbu,max in mm.

    Raises:
        ValueError: If d is negative.
    """
    if d < 0:
        raise ValueError(f'd must not be negative. Got {d}')

    return d


def str_max_slab(d: float) -> float:
    """Calculate the maximum transverse spacing of shear legs in slabs.

    EN1992-1-1:2023 Table (12.2)

    Args:
        d (float): Effective depth of the slab in mm.

    Returns:
        float: Maximum transverse spacing of shear legs str,max in mm.

    Raises:
        ValueError: If d is negative.
    """
    if d < 0:
        raise ValueError(f'd must not be negative. Got {d}')

    return 1.5 * d


def min_slab_reinforcement_along_free_edge(h: float, lbd: float) -> bool:
    """Get the slab reinforcement minimum length along free edge.

    EN1992-1-1:2023 Fig. (12.4)

    Args:
        h (float): Height of the free edge in mm.
        lbd (float): Base anchor length in mm.

    Returns:
        bool: The mimimun reinforcement length in mm.

    Raises:
        ValueError: If as_top or as_bottom is negative.
    """
    if h < 0:
        raise ValueError(f'd must not be negative. Got {h}')
    if lbd < 0:
        raise ValueError(f'd must not be negative. Got {lbd}')

    return max(2 * h, lbd)


def check_slab_depth_for_shear_reinforcement(
    h: float,
    reinforcement_type: Literal[
        'stirrups', 'links', 'headed_bars', 'bent_up_bars'
    ],
) -> bool:
    """Check if the slab depth meets the minimum
        requirement for the specified type of shear reinforcement.

    EN1992-1-1:2023 Section 12.4.2(2)

    Args:
        h (float): Depth of the slab in mm.
        reinforcement_type (str): Type of reinforcement
            ('stirrups', 'links', 'headed_bars', 'bent_up_bars').

    Returns:
        bool: True if the slab depth meets the requirement, False otherwise.

    Raises:
        ValueError: If h is negative.
    """
    if h < 0:
        raise ValueError(f'd must not be negative. Got {h}')

    if reinforcement_type in ('stirrups', 'links', 'headed_bars'):
        return h >= 200
    # If bent_up_bars
    return h >= 160


def check_slab_shear_stress_conditions(
    tau_ed: float, fcd: float, alignment: float
) -> bool:
    """Check if the shear reinforcement
        condition is met for a slab based on the stress and alignment.
    Interpolates the allowable shear stress based on the alignment angle.

    EN1992-1-1:2023 12.4.2(5)

    Args:
        tau_ed (float): Design shear stress in MPa.
        fcd (float): Design compressive strength of concrete in MPa.
        alignment (float): Alignment of shear reinforcement
            in degrees (0 for vertical, 45 for 45°).

    Returns:
        bool: True if the shear stress condition is satisfied, False otherwise.

    Raises:
        ValueError: If tau_ed or fcd is negative or
            if alignment is out of range.

    """
    if tau_ed < 0:
        raise ValueError(f'tau_ed must not be negative. Got {tau_ed}')
    if fcd < 0:
        raise ValueError(f'tau_ed must not be negative. Got {fcd}')

    if not (0 <= alignment <= 45):
        raise ValueError(
            f'alignment must be between 0 and 45 degrees. Got {alignment}'
        )

    # Linear interpolation between the vertical and 45-degree values
    tau_ed_limit = (0.08 + (alignment / 45) * 0.08) * fcd

    return tau_ed <= tau_ed_limit


def slab_column_max_leg_shear_reinforcement_diameter(
    d: float,
    reinforcement_type: Literal[
        'single_leg',
        'open_stirrups',
        'closed_stirrups',
        'bent_up_bars',
        'headed_bars',
    ],
) -> float:
    """Calculate the maximum effective diameter of shear reinforcement
        based on the type and effective depth.

    EN1992-1-1:2023 Eq. (12.7), (12.8), (12.9)

    Args:
        d (float): Effective depth of the slab in mm.
        reinforcement_type (str): Type of reinforcement
            ('single_leg', 'open_stirrups', 'closed_stirrups',
            'bent_up_bars', 'headed_bars').

    Returns:
        float: Maximum effective diameter of shear reinforcement in mm.

    Raises:
        ValueError: If the reinforcement_type is not
            recognized or if d is negative.
    """
    if d < 0:
        raise ValueError(f'd must not be negative. Got {d}')

    if reinforcement_type in ('single_leg', 'open_stirrups'):
        return 10 * math.sqrt(d / 200)
    if reinforcement_type in (
        'closed_stirrups',
        'bars_with_similar_anchorage',
    ):
        return 11 * math.sqrt(d / 200)
    # if ('bent_up_bars', 'headed_bars')
    return 16 * math.sqrt(d / 200)


def slab_column_tangential_spacing_limit(
    dv: float,
    distance_to_column_edge: float,
    slab_type: Literal['column_edge', 'flat_slab', 'column_base'],
) -> float:
    """Calculate the tangential spacing limit for shear
        reinforcement based on the distance to the column edge.

    EN1992-1-1:2023 Section 12.5.1(2)

    Args:
        dv (float): Effective depth of the slab in mm.
        distance_to_column_edge (float): Distance from
            the shear reinforcement to the column edge in mm.
        slab_type (str): Type of slab
            ('column_edge', 'flat_slab', 'column_base').

    Returns:
        float: Maximum tangential spacing in mm.

    Raises:
        ValueError: If dv or distance_to_column_edge is negative,
            or if slab_type is not recognized.
    """
    if dv < 0:
        raise ValueError(f'dv must not be negative. Got {dv}')
    if distance_to_column_edge < 0:
        raise ValueError(
            'distance_to_column_edge must not be negative. '
            + f'Got {distance_to_column_edge}'
        )
    if slab_type == 'column_edge':
        return 1.5 * dv
    if slab_type == 'flat_slab':
        return 0.75 * dv

    # If column_base
    return 0.5 * dv


def Vrd_int_flat_slab(
    As_int: float, fyd: float, ductility_class: Literal['B', 'C'], VEd: float
) -> float:
    """Calculate the integrity reinforcement
        resistance for progressive collapse.

    EN1992-1-1:2023 Eq. (12.10)

    Args:
        As_int (float): Sum of the cross-sections of all
            integrity reinforcement bars crossing a column edge in mm2.
        fyd (float): Yield strength of the integrity reinforcement in MPa.
        kint (float): Ductility class (B or C).
        VEd (float): Design shear in kN

    Returns:
        float: Integrity reinforcement resistance VRd,int in kN.

    Raises:
        ValueError: If any input value is negative.
    """
    if As_int < 0:
        raise ValueError(f'As_int must not be negative. Got {As_int}')
    if fyd < 0:
        raise ValueError(f'fyd must not be negative. Got {fyd}')
    kint = 0.37 if ductility_class == 'B' else 0.49

    return max(As_int * fyd * kint / 1000, abs(VEd))


def Vrd_w_int_flat_slab(
    rho_w: float, fywd: float, b0_5: float, dv: float, VEd: float
) -> float:
    """Calculate the integrity reinforcement resistance
        for slabs with shear reinforcement.

    EN1992-1-1:2023 Eq. (12.11)

    Args:
        rho_w (float): Ratio of shear reinforcement.
        fywd (float): Yield strength of the shear reinforcement in MPa.
        b0_5 (float): Length of the control perimeter in mm.
        dv (float): Effective depth of the slab in mm.
        VEd (float): Design shear value in kN.

    Returns:
        float: Integrity reinforcement resistance VRd,w,int in kN.

    Raises:
        ValueError: If any input value is negative.
    """
    if rho_w < 0:
        raise ValueError(f'rho_w must not be negative. Got {rho_w}')
    if fywd < 0:
        raise ValueError(f'fywd must not be negative. Got {fywd}')
    if b0_5 < 0:
        raise ValueError(f'b0_5 must not be negative. Got {b0_5}')
    if dv < 0:
        raise ValueError(f'dv must not be negative. Got {dv}')

    return min(rho_w * fywd * b0_5 * dv / 1000, abs(VEd))  # Convert to kN


def Vrd_hog_flat_slab(
    nhog: int, fck: float, gamma_c: float, phi: float, s: float, c: float
) -> float:
    """Calculate the contribution of hogging reinforcement
        to the resistance against progressive collapse.

    EN1992-1-1:2023 Eq. (12.12)

    Args:
        nhog (int): Number of hogging bars crossing the
            control perimeter and fully anchored.
        fck (float): Characteristic compressive
            strength of concrete in MPa.
        gamma_c (float): Partial factor for concrete
            for accidental design situation.
        phi (float): Diameter of the
            hogging reinforcement in mm.
        s (float): Spacing of the hogging reinforcement in mm.
        c (float): Cover of the hogging
            reinforcement in mm.

    Returns:
        float: Hogging reinforcement resistance VRd,hog (kN).

    Raises:
        ValueError: If any input value is negative.
    """
    if nhog < 0:
        raise ValueError(f'nhog must not be negative. Got {nhog}')
    if fck < 0:
        raise ValueError(f'fck must not be negative. Got {fck}')
    if gamma_c < 0:
        raise ValueError(f'gamma_c must not be negative. Got {gamma_c}')
    if phi < 0:
        raise ValueError(f'phi must not be negative. Got {phi}')
    if s < 0:
        raise ValueError(f's must not be negative. Got {s}')
    if c < 0:
        raise ValueError(f'c must not be negative. Got {c}')

    bef_hog = min(s - phi, 6 * phi, 4 * c)
    return (
        nhog * (math.sqrt(fck) / gamma_c) * bef_hog * phi / 1000
    )  # Convert to kN


def As_column_min(NEd: float, fyd: float, Ac: float) -> float:
    """Calculate the minimum amount of longitudinal
        reinforcement for robustness and to avoid
        compressive yielding due to creep and shrinkage.

    EN1992-1-1:2023 Table 12.3

    Args:
        NEd (float): Design axial force in kN.
        fyd (float): Design yield strength of the reinforcement in MPa.
        Ac (float): Cross-sectional area of the column in mm2.

    Returns:
        float: Minimum area of longitudinal reinforcement in mm2.

    Raises:
        ValueError: If any input value is negative.
    """
    NEd = abs(NEd)
    if fyd < 0:
        raise ValueError(f'fyd must not be negative. Got {fyd}')
    if Ac < 0:
        raise ValueError(f'Ac must not be negative. Got {Ac}')

    return max(0.1 * NEd * 1000 / fyd, 0.002 * Ac)


def s_max_poly_col(h: float, b: float) -> float:
    """Calculate the maximum longitudinal spacing of
        transverse reinforcement for polygonal cross-sections.

    EN1992-1-1:2023 Table 12.3

    Args:
        h (float): Height of the column in mm.
        b (float): Width of the column in mm.

    Returns:
        float: Maximum longitudinal spacing in mm.

    Raises:
        ValueError: If any input value is negative.
    """
    if h < 0:
        raise ValueError(f'h must not be negative. Got {h}')
    if b < 0:
        raise ValueError(f'b must not be negative. Got {b}')

    return min(h, b, 400)


def s_max_circular_col(n_bars: int, diameter: float) -> float:
    """Calculate the maximum longitudinal spacing of
        transverse reinforcement for circular cross-sections.

    EN1992-1-1:2023 Table 12.3

    Args:
        n_bars (int): Number of longitudinal bars.
        diameter (float): Diameter of the circular column in mm.

    Returns:
        float: Maximum longitudinal spacing in mm.

    Raises:
        ValueError: If any input value is negative.
    """
    if n_bars < 0:
        raise ValueError(f'n_bars must not be negative. Got {n_bars}')
    if diameter < 0:
        raise ValueError(f'diameter must not be negative. Got {diameter}')

    return min(diameter * math.pi / n_bars, 400)


def s_max_col_int(
    phi_l_max: float, h: float, b: float, long_bars_res: bool
) -> float:
    """Calculate the maximum spacing of transverse
        reinforcement for columns in the intermediate region.

    EN1992-1-1:2023 Table 12.3

    Args:
        phi_l_max (float): Maximum diameter of longitudinal bars in mm.
        h (float): height of the column in mm.
        b (float): width of the column in mm.
        long_bars_res (bool): True if the longitudinal bars account
            for column resistance. Otherwise False.

    Returns:
        float: Maximum spacing of transverse reinforcement in mm.

    Raises:
        ValueError: If any input value is negative.
    """
    if phi_l_max < 0:
        raise ValueError(f'phi_l_max must not be negative. Got {phi_l_max}')
    if h < 0:
        raise ValueError(f'h must not be negative. Got {h}')
    if b < 0:
        raise ValueError(f'b must not be negative. Got {b}')
    if long_bars_res:
        return min(h, b, 400)

    return min(20 * phi_l_max, h, b, 300)


def s_max_col_end(s_max_col: float) -> float:
    """Calculate the maximum spacing of transverse
        reinforcement for columns in the end region.

    EN1992-1-1:2023 Table 12.3

    Args:
        s_max_col (float): Maximum spacing of transverse reinforcement in mm.

    Returns:
        float: Maximum spacing of transverse reinforcement in mm.

    Raises:
        ValueError: If s_max_col is negative.
    """
    if s_max_col < 0:
        raise ValueError(f's_max_col must not be negative. Got {s_max_col}')
    return 0.6 * s_max_col


def As_wall_min_v(
    Ac: float,
    fctm: float,
    fyk: float,
    design_case: Literal['in_plane_stress', 'compression_bending'],
) -> float:
    """Calculate the minimum amount of vertical
        reinforcement for walls and deep beams.

    EN1992-1-1:2023 Table 12.4

    Args:
        Ac (float): Cross-sectional area of the wall or deep beam .
        fctm (float): Mean tensile strength of concrete in MPa.
        fyk (float): Characteristic yield strength of reinforcement in MPa.
        design_case (str): Design case
            ('in_plane_stress' or 'compression_bending').

    Returns:
        float: Minimum area of vertical reinforcement in mm2.

    Raises:
        ValueError: If any input value is negative or
            if design_case is not recognized.
    """
    if Ac < 0:
        raise ValueError(f'Ac must not be negative. Got {Ac}')
    if fctm < 0:
        raise ValueError(f'fctm must not be negative. Got {fctm}')
    if fyk < 0:
        raise ValueError(f'fyk must not be negative. Got {fyk}')

    if design_case == 'in_plane_stress':
        return 0.25 * Ac * fctm / fyk

    # If compression_bending
    return 0.001 * Ac


def As_wall_min_h(
    Ac: float,
    As_v: float,
    fctm: float,
    fyk: float,
    design_case: Literal['in_plane_stress', 'compression_bending'],
) -> float:
    """Calculate the minimum amount of horizontal reinforcement
        for walls and deep beams.

    EN1992-1-1:2023 Table 12.4

    Args:
        Ac (float): Cross-sectional area of the wall or deep beam in mm2.
        As_v (float): Area of vertical reinforcement in mm2.
        fctm (float): Mean tensile strength of concrete in MPa.
        fyk (float): Characteristic yield strength of reinforcement in MPa.
        design_case (str): Design case
            ('in_plane_stress' or 'compression_bending').

    Returns:
        float: Minimum area of horizontal reinforcement in mm2.

    Raises:
        ValueError: If any input value is negative or if
            design_case is not recognized.
    """
    if Ac < 0:
        raise ValueError(f'Ac must not be negative. Got {Ac}')
    if As_v < 0:
        raise ValueError(f'As_v must not be negative. Got {As_v}')
    if fctm < 0:
        raise ValueError(f'fctm must not be negative. Got {fctm}')
    if fyk < 0:
        raise ValueError(f'fyk must not be negative. Got {fyk}')

    if design_case == 'in_plane_stress':
        return 0.25 * Ac * fctm / fyk
    # If compression_bending
    return 0.25 * As_v


def s_max_wall_v(h: float) -> float:
    """Calculate the maximum spacing of vertical reinforcement in walls.

    EN1992-1-1:2023 Table 12.4

    Args:
        h (float): Thickness of the wall in mm.

    Returns:
        float: Maximum spacing of vertical reinforcement in mm.

    Raises:
        ValueError: If h is negative.
    """
    if h < 0:
        raise ValueError(f'h must not be negative. Got {h}')

    return min(3 * h, 400)


def s_max_wall_h() -> float:
    """Calculate the maximum spacing of horizontal reinforcement
        in walls.

    EN1992-1-1:2023 Table 12.4

    Returns:
        float: Maximum spacing of horizontal reinforcement in mm.
    """
    return 400


def min_di_support_and_joint(
    ch_i: float,
    delta_a: float,
    loop_type: Literal['horizontal_loops', 'vertical_bent'],
    r_i: float = 0,
) -> float:
    """Calculate the nominal length di of a simple support or bearing
        based on the type of loop and reinforcement.

    EN1992-1-1:2023 12.10(5)

    Args:
        ch_i (float): Nominal cover (horizontal or vertical) in mm.
        delta_a (float): Allowance for construction deviations in mm.
        r_i (float, optional): Mandrel radius of bend of longitudinal
            reinforcement in mm. Required if loop_type is "vertical_bent".
        loop_type (str): Type of loop ('horizontal_loops', 'vertical_bent').

    Returns:
        float: Nominal length di in mm.

    Raises:
        ValueError: If any input value is negative.
    """
    if ch_i < 0:
        raise ValueError(f'ch_i must not be negative. Got {ch_i}')
    if delta_a < 0:
        raise ValueError(f'delta_a must not be negative. Got {delta_a}')
    if r_i < 0:
        raise ValueError(f'r_i must not be negative. Got {r_i}')

    if loop_type == 'horizontal_loops':
        return ch_i + delta_a

    # If vertical_bent
    return ch_i + delta_a + r_i
