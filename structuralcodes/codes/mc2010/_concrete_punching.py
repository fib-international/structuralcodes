"""Covers punching in Model code 2010, 7.3.5.1 to 7.3.5.4."""

import typing as t
import warnings
from math import cos, pi, sin


def b_0(v_ed: float, v_prep_d_max: float) -> float:
    """Gives the general output for b_0, shear-resisting control perimeter.

    fib Model Code 2010, Eq. (7.3-57).

    Args:
        V_ed (float): The acting shear force from the columns.
        v_prep_d_max (float): The maximum shear force per unit length
            perpendicular to the basic control parameter (Figure 7.3-24).

    Returns:
        float: The shear-resisting control perimeter, b_0.
    """
    return v_ed / v_prep_d_max


def b_s(
    l_x: float,
    l_y: float,
) -> float:
    """The width of the support strip for calculating m_ed.

    fib Model Code 2010, Eq. (7.3-76).

    Args:
        l_x (float): The width in x direction that the collumn carries.
        l_y (float): The width in y direction that the collumn carries.

    Returns:
        float: The width of the support strip for calculating m_ed.
    """
    r_sx = 0.22 * l_x  # see  MC2010 7.3.5.3
    r_sy = 0.22 * l_y  # see  MC2010 7.3.5.3
    l_min = min(l_x, l_y)
    return min(1.5 * (r_sx * r_sy) ** 0.5, l_min)


def m_ed(
    v_ed: float,
    e_u: float,
    b_s: float,
    inner: bool,
    edge_par: bool,
    edge_per: bool,
    corner: bool,
) -> float:
    """The average bending moment acting in the support strip.

    fib Model Code 2010, Eq. (7.3-71), (7.3-72), (7.3-73) and (7.3-74).

    Args:
        v_ed (float): The acting shear force from the columns.
        e_u (float): Refers to the eccentricity of the resultant of shear
            forces with respect to the centroid.
        b_s (float): The width of the support strip for calculating m_ed.
        inner (bool): Is true only if the column is a inner column.
        edge_par (bool): Is true only if the column is a edge column with
            tension reinforcement parallel to the edge.
        edge_per (bool): Is true only if the column is a edge column with
            tension reinforcement perpendicular to the edge.
        corner (bool): Is true only if the column is a corner column.

    Returns:
        float: The bending moment acting in the support strip regardless of the
        position of the column.
    """
    if inner:
        return v_ed * ((1 / 8) + abs(e_u) / (2 * b_s))
    if edge_par:
        return max(v_ed * ((1 / 8) + e_u / (2 * b_s)), v_ed / 4)
    if edge_per:
        return v_ed * ((1 / 8) + e_u / (b_s))
    if corner:
        return max(v_ed * ((1 / 8) + e_u / (b_s)), v_ed / 2)
    raise ValueError('Placement is not defined, only one needs to be True')


def b_sr(
    l_x: float,
    l_y: float,
) -> float:
    """The width of the support strip in the radial direction.

    fib Model Code 2010, 7.3.5.3.

    Args:
        l_x (float): The width in x direction that the collumn carries.
        l_y (float): The width in y direction that the collumn carries.

    Returns:
        float: The width of the support strip in the radial direction.
    """
    return min(l_x, l_y)


def r_s(
    l_x: float,
    l_y: float,
    x_direction: bool,
    is_level_three_approximation: bool = False,
    column_edge_or_corner: bool = False,
    b_sr: t.Optional[float] = None,
) -> float:
    """The position where the radial bending moment is zero with respect to the
    support axis.

    fib Model Code 2010, 7.3.5.3 and Eq. (7.3-78) for Level III of
    Approximation.

    Args:
        l_x (float): The width in x direction that the collumn carries.
        l_y (float): The width in y direction that the collumn carries.
        x_direction (bool): True if the radial bending moment is zero in the x
            direction, False if it is in the y direction.
        column_edge_or_corner (bool): True if the column is an edge or corner
            column, False if it is an inner column. To be used for Level III of
            Approximation.
        b_sr (float): The width of the support strip in the radial direction.

    Returns:
        float: The position where the radial bending moment is zero with
        respect to the support axis.
    """
    r_s = 0.22 * l_x if x_direction is True else 0.22 * l_y

    if column_edge_or_corner and is_level_three_approximation:
        if b_sr is None:
            raise ValueError(
                'b_sr is not defined for Level 3 of Approximation'
            )
        return max(r_s, 0.67 * b_sr)
    return r_s


def psi_punching_level_one(
    l_x: float,
    l_y: float,
    f_yd: float,
    d_eff: float,
    e_s: float,
) -> float:
    """The psi value for the punching level one.

    fib Model Code 2010, Eq. (7.3-70).

    Args:
        r_s (float): The position where the radial bending moment is zero with
            respect to the support axis.
        f_yd (float): Design strength of reinforment steel in MPa.
        d_eff (float): The mean value of the effective depth in mm.
        e_s (float): The E_modulus for steel in MPa.

    Returns:
        float: The psi value for the punching level one.
    """
    r_s = max(0.22 * l_x, 0.22 * l_y)  # MC2010 7.3.5.3 - Unique for (7.3-70)
    return 1.5 * r_s * f_yd / (d_eff * e_s)


def psi_punching_level_two(
    r_s: float,
    f_yd: float,
    d_eff: float,
    e_s: float,
    m_ed: float,
    m_rd: float,
    m_Pd: float = 0,
) -> float:
    """The psi value for the punching level two.

    fib Model Code 2010, Eq. (7.3-75) and (7.3-77).

    Args:
        r_s (float): The position where the radial bending moment is zero with
            respect to the support axis.
        f_yd (float): Design strength of reinforment steel in MPa.
        d_eff (float): The mean value of the effective depth in mm.
        e_s (float): The E_modulus for steel in MPa.
        m_ed (float): The average bending moment acting in the support strip.
        m_rd (float): The design average strength per unit length in MPa.
        m_Pd (float): Optional to cover Eq. (7.3-77) for prestressed slabs.
            The average decompression moment over the width of the support
            strip (b_s) due to prestressing.

    Returns:
        float: The psi value for the punching level two.
    """
    return (1.5 * r_s * f_yd / (d_eff * e_s)) * (
        (m_ed - m_Pd) / (m_rd - m_Pd)
    ) ** 1.5


def psi_punching_level_three(
    psi_punching_level_two: float,
    is_uncracked_model: bool = False,
    is_moment_from_uncracked_model: bool = False,
) -> float:
    """The psi value for the punching level three.

    fib Model Code 2010, Level III of Approximation Eq. (7.3-75) and Eq.
    (7.3-77) with coefficient 1.2 instead of 1.5 under specific conditions.

    Args:
        psi_punching_level_two (float): The psi value for the punching level 2.
        is_uncracked_model (bool): True if r_s is calculated using a linear
            elastic (uncracked) model.
        is_moment_from_uncracked_model (bool): True if m_sd is calculated from
            a linear elastic (uncracked) model as the average value of the
            moment for design of the flexural reinforcement over the width of
            the support strip (b_s).

    Returns:
        float: The psi value for the punching level three.
    """
    if is_uncracked_model and is_moment_from_uncracked_model:
        return 1.2 / 1.5 * psi_punching_level_two
    return psi_punching_level_two


def psi_punching(
    psi_punching_level_one: float,
    psi_punching_level_two: float,
    psi_punching_level_three: float,
    approx_lvl_p: float,
) -> float:
    """The rotation of the slab around the supported area.

    fib Model Code 2010, Clause 7.3.5.4.

    Args:
        psi_punching_level_one (float): The psi value for the punching level 1.
        psi_punching_level_two (float): The psi value for the punching level 2.
        psi_punching_level_three (float): The psi value for the punching level
            3.
        approx_lvl_p (float): The approx level for punching.

    Returns:
        float: psi for the chosen approx level in punching.
    """
    if approx_lvl_p == 1:
        return psi_punching_level_one
    if approx_lvl_p == 2:
        return psi_punching_level_two
    if approx_lvl_p == 3:
        return psi_punching_level_three
    raise ValueError('Approximation level is not defined')


def k_dg(
    d_g: float,
) -> float:
    """Calculate k_dg factor for punching resistance.

    fib Model Code 2010, Eq. (7.3-62).

    Args:
        d_g (float): Maximum size of aggregate.

    Returns:
        float: k_dg factor.
    """
    return max(32 / (16 + d_g), 0.75)


def k_psi(
    k_dg: float,
    d_eff: float,
    psi_punching: float,
) -> float:
    """Calculate k_psi factor for punching resistance.

    fib Model Code 2010, Eq. (7.3-63).

    Args:
        k_dg (float): k_dg factor.
        d_eff (float): The mean value of the effective depth in mm.
        psi_punching (float): psi value from psi_punching.

    Returns:
        float: k_psi factor.
    """
    return min(1 / (1.5 + 0.9 * k_dg * d_eff * psi_punching), 0.6)


def v_rdc_punching(
    k_psi_val: float,
    b_0: float,
    d_v: float,
    f_ck: float,
    gamma_c: float = 1.5,
) -> float:
    """Punching resistance from the concrete.

    fib Model Code 2010, Eq. (7.3-61).

    Args:
        k_psi_val (float): k_psi value from k_psi.
        b_0 (float): The shear-resisting control perimeter from b_0.
        d_v (float): The effective depth considering support in mm.
        f_ck (float): Characteristic strength in MPa.
        gamma_c: Safety factor for concrete.

    Returns:
        float: v_rdc for punching with the right approx level.
    """
    return k_psi_val * b_0 * d_v * (f_ck**0.5) / gamma_c


def f_ywd(
    f_ywk: float,
    gamma_s: float,
) -> float:
    """Calculate f_ywd for punching resistance.

    fib Model Code 2010, Eq. (7.3-64).

    Args:
        f_ywk (float): Characteristic yield strength of the shear reinforcement
            in MPa.
        gamma_s (float): Safety factor for reinforcement.

    Returns:
        float: f_ywd for punching resistance.
    """
    return f_ywk / gamma_s


def sigma_swd(
    e_s: float,
    psi_punching: float,
    alpha: float,
    f_bd: float,
    d_eff: float,
    f_ywd: float,
    phi_w: float,
) -> float:
    """Calculate sigma_swd for punching resistance.

    fib Model Code 2010, Eq. (7.3-65).

    Args:
        e_s (float): The E_modulus for steel in MPa.
        psi_punching (float): psi value from psi_punching.
        alpha (float): Inclination of the stirrups in degrees.
        f_bd (float): The design bond strength in MPa.
        d_eff (float): The mean value of the effective depth in mm.
        f_ywd (float): Design yield strength of the shear reinforcement in MPa.
        phi_w (float): The diameter of the shear reinforcement.

    Returns:
        float: sigma_swd.
    """
    return min(
        (e_s * psi_punching / 6)
        * (sin(alpha * pi / 180) + cos(alpha * pi / 180))
        * (sin(alpha * pi / 180) + f_bd * d_eff / (f_ywd * phi_w)),
        f_ywd,
    )


def v_rds_punching(
    f_ywd: float,
    e_u: float,
    b_u: float,
    alpha: float,
    sigma_swd: float,
    a_sw: float,
    v_ed: float,
) -> float:
    """The punching resistance from shear reinforcement.

    fib Model Code 2010, Eq. (7.3-64).

    Args:
        f_ywd (float): Design yield strength of the shear reinforcement in MPa.
        e_u (float): The ecentrisity of the result of shear forces with respect
            to the centroid (Figure 7.3-27b).
        b_u (float): The diamter of a circle with same surface as the region
            inside the basic control perimeter (Figure 7.3-27b).
        alpha (float): Inclination of the stirrups in degrees.
        sigma_swd (float): sigma_swd from sigma_swd.
        a_sw (float): The area of the shear reinforcement in mm^2.
        v_ed (float): The acting shear force from the columns.

    Returns:
        float: Punching resistance that comes from reinforcement.
    """
    k_e = 1 / (1 + e_u / b_u)
    if (a_sw * k_e * f_ywd) < 0.5 * v_ed:
        warnings.warn(
            'Consider increasing punching shear reinforcement for sufficient '
            'deformation capacity'
        )
    return a_sw * k_e * sigma_swd * sin(alpha * pi / 180)


def v_rd_max_punching(
    d_v: float,
    f_ck: float,
    d_head: bool,
    stirrups_compression: bool,
    b0_val: float,
    k_psi_val: float,
    gamma_c: float = 1.5,
) -> float:
    """Finds the maximum value you can have for v_rd_punching.

    fib Model Code 2010, Eq. (7.3-69).

    Args:
        d_v (float): The effective depth considering support in mm.
        f_ck (float): Characteristic strength in MPa.
        d_head (bool): True if diameter of heads is three times larger than.
        stirrups_compression: (bool): Stirrups with sufficient length at
            compression face, and bent on tension face.
        b0_val (float): The shear-resisting control perimeter from b_0.
        k_psi_val (float): k_psi value from k_psi.
        gamma_c (float): Safety factor for concrete.

    Return:
        float: The maximum allowed punching resistance.
    """
    if d_head:
        k_sys = 2.8
    elif stirrups_compression:
        k_sys = 2.4
    else:
        k_sys = 2

    base_resistance = b0_val * d_v * (f_ck**0.5 / gamma_c)

    return min(k_sys * k_psi_val * base_resistance, base_resistance)


def v_rd_punching(v_rd_c: float, v_rd_s: float, v_rd_max: float) -> float:
    """The total resistance for punching.

    fib Model Code 2010, Eq. (7.3-60).

    Args:
        v_rd_c: Concrete contribution to punching resistance.
        v_rd_s: Shear reinforcement contribution to punching resistance.
        v_rd_max: Maximum punching resistance.

    Returns:
        float: Total punching resistance as min(v_rd_c + v_rd_s, v_rd_max).
    """
    return min(v_rd_c + v_rd_s, v_rd_max)
