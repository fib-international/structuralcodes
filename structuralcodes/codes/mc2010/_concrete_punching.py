"""Covers punching in Model code 2010, 7.3.5.1 to 7.3.5.4."""

import warnings
from math import cos, pi, sin


def b_0(v_ed: float, v_prep_d_max: float) -> float:
    """Gives the general output for b_0, shear-resisting control perimeter.

    fib Model Code 2010, eq. (7.3-57).

    Args:
        V_ed (float): The acting shear force from the columns.
        v_prep_d_max (float): The maximum shear force per unit length
            perpendicular to the basic control parameter (Figure 7.3-24).

    Returns:
        float: The shear-resisting control perimeter, b_0.
    """
    return v_ed / v_prep_d_max


def m_ed(
    v_ed: float,
    e_u: float,
    l_x: float,
    l_y: float,
    inner: bool,
    edge_par: bool,
    edge_per: bool,
    corner: bool,
) -> float:
    """The average bending moment acting in the support strip.

    fib Model Code 2010, eq. (7.3-76), (7.3-71), (7.3-72), (7.3-73)
    and (7.3-74).

    Args:
        v_ed (float): The acting shear force from the columns.
        e_u (float): Refers to the eccentricity of the resultant of shear
            forces with respect to the centroid.
        l_x (float): The width in x direction that the collumn carries.
        l_y (float): The width in y direction that the collumn carries.
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
    r_sx = 0.22 * l_x
    r_sy = 0.22 * l_y
    l_min = min(l_x, l_y)
    b_s = min(1.5 * (r_sx * r_sy) ** 0.5, l_min)
    if inner:
        return v_ed * ((1 / 8) + abs(e_u) / (2 * b_s))
    if edge_par:
        return max(v_ed * ((1 / 8) + e_u / (2 * b_s)), v_ed / 4)
    if edge_per:
        return v_ed * ((1 / 8) + e_u / (b_s))
    if corner:
        return max(v_ed * ((1 / 8) + e_u / (b_s)), v_ed / 2)
    raise ValueError('Placement is not defined, only one needs to be True')


def psi_punching(
    l_x: float,
    l_y: float,
    f_yd: float,
    d: float,
    e_s: float,
    approx_lvl_p: float,
    v_ed: float,
    e_u: float,
    inner: bool,
    edge_par: bool,
    edge_per: bool,
    corner: bool,
    m_rd: float,
    x_direction: bool,
) -> float:
    """The rotation of the slab around the supported area.

    fib Model Code 2010, eq. (7.3-70), (7.3-75) and (7.3-77).

    Args:
        l_x (float): The distance between two columns in x direction.
        l_y (float): The distance between two columns in y direction.
        f_yd (float): Design strength of reinforment steel in MPa.
        d (float): The mean value of the effective depth in mm.
        e_s (float): The E_modulus for steel in MPa.
        approx_lvl_p (float): The approx level for punching.
        v_ed (float): The acting shear force from the columns.
        e_u (float): Refers to the eccentricity of the resultant of shear
            forces with respect to the centroid.
        inner (bool): Is true only if the column is a inner column.
        edge_par (bool): Is true only if the column is a edge column with
            tension reinforcement parallel to the edge.
        edge_per (bool): Is true only if the column is a edge column with
            tension reinforcement perpendicular to the edge.
        corner (bool): Is true only if the column is a corner column.
        m_rd (float): The design average strength per unit length in MPa.
        m_pd: (float): The average decompresstion moment due to prestressing
            in MPa.

    Returns:
        float: psi for the chosen approx level in punching.
    """
    r_s = max(0.22 * l_x, 0.22 * l_y)
    if approx_lvl_p == 1:
        psi = 1.5 * r_s * f_yd / (d * e_s)

    elif approx_lvl_p == 2:
        r_s = 0.22 * l_x if x_direction is True else 0.22 * l_y
        psi = (1.5 * r_s * f_yd / (d * e_s)) * (
            (m_ed(v_ed, e_u, l_x, l_y, inner, edge_par, edge_per, corner))
            / (m_rd)
        ) ** 1.5

    return psi


def v_rdc_punching(
    l_x: float,
    l_y: float,
    f_yd: float,
    d: float,
    e_s: float,
    approx_lvl_p: float,
    dg: float,
    f_ck: float,
    d_v: float,
    v_ed: float,
    e_u: float,
    inner: bool,
    edge_par: bool,
    edge_per: bool,
    corner: bool,
    m_rd: float,
    m_pd: float,
    v_prep_d_max: float,
    gamma_c: float = 1.5,
) -> float:
    """Punching resistance from the concrete.

    fib Model Code 2010, eq. (7.3-61), (7.3-62) and (7.3-63).

    Args:
        l_x (float): The distance between two columns in x direction.
        l_y (float): The distance between two columns in y direction.
        f_yd (float): Design strength of reinforment steel in MPa.
        d (float): The mean value of the effective depth in mm.
        e_s (float): The E_modulus for steel in MPa.
        approx_lvl_p (float): The approx level for punching.
        dg (float): Maximum size of aggregate.
        f_ck (float): Characteristic strength in MPa.
        d_v (float): The effective depth considering support in mm.
        v_ed (float): The acting shear force from the columns.
        e_u (float): Refers to the eccentricity of the resultant of shear
            forces with respect to the centroid.
        inner (bool): Is true only if the column is a inner column.
        edge_par (bool): Is true only if the column is a edge column with
            tension reinforcement parallel to the edge.
        edge_per (bool): Is true only if the column is a edge column with
            tension reinforcement perpendicular to the edge.
        corner (bool): Is true only if the column is a corner column.
        m_rd (float): The design average strength per unit length in MPa.
        m_pd: (float): The average decompresstion moment due to prestressing
            in MPa.
        v_prep_d_max (float): The maximum shear force per unit length
            perpendicular to the basic control parameter (Figure 7.3-24).
        gamma_c: Safety factor for concrete.

    Returns:
        float: v_rdc for punching with the right approx level.
    """
    k_dg = max(32 / (16 + dg), 0.75)
    k_psi = min(
        1
        / (
            1.5
            + 0.9
            * k_dg
            * d
            * psi_punching(
                l_x,
                l_y,
                f_yd,
                d,
                e_s,
                approx_lvl_p,
                v_ed,
                e_u,
                inner,
                edge_par,
                edge_per,
                corner,
                m_rd,
                m_pd,
            )
        ),
        0.6,
    )
    return k_psi * b_0(v_ed, v_prep_d_max) * d_v * (f_ck**0.5) / gamma_c


def v_rds_punching(
    e_u: float,
    b_u: float,
    l_x: float,
    l_y: float,
    f_yd: float,
    d: float,
    e_s: float,
    approx_lvl_p: float,
    v_ed: float,
    inner: bool,
    edge_par: bool,
    edge_per: bool,
    corner: bool,
    m_rd: float,
    m_pd: float,
    alpha: float,
    f_bd: float,
    f_ywk: float,
    phi_w: float,
    a_sw: float,
    gamma_s: float,
):
    """The punching resistance from shear reinforcement.

    fib Model Code 2010, eq. (7.3-64) and (7.3-65).

    Args:
        e_u (float): The ecentrisity of the result of shear forces
            with respect to the centroid (Figure 7.3-27b).
        b_u (float): The diamter of a circle with same surface as the region
            inside the basic control perimeter (Figure 7.3-27b).
        l_x (float): The distance between two columns in x direction.
        l_y (float): The distance between two columns in y direction.
        f_yd (float): Design strength of reinforment steel in MPa.
        d (float): The mean value of the effective depth in mm.
        e_s (float): The E_modulus for steel in MPa.
        approx_lvl_p (float): The approx level for punching.
        v_ed (float): The acting shear force from the columns.
        inner (bool): Is true only if the column is a inner column.
        edge_par (bool): Is true only if the column is a edge column with
            tension reinforcement parallel to the edge.
        edge_per (bool): Is true only if the column is a edge column with
            tension reinforcement perpendicular to the edge.
        corner (bool): Is true only if the column is a corner column.
        m_rd (float): The design average strength per unit length in MPa.
        m_pd: (float): The average decompresstion moment due to prestressing in
            MPa.
        alpha (float): Inclination of the stirrups in degrees.
        f_bd (float): The design bond strength in MPa.
        f_ywk (float): Characteristic yield strength of the shear reinforcement
            in MPa.
        phi_w (float): The diameter of the shear reinforcement.
        a_sw (float): The area of the shear reinforcement in mm^2.
        gamma_s (float): Safety factor for reinforcement.

    Returns:
        float: Punching resistance that comes from reinforcement.
    """
    f_ywd = f_ywk / gamma_s
    k_e = 1 / (1 + e_u / b_u)
    sigma_swd = min(
        (
            e_s
            * psi_punching(
                l_x,
                l_y,
                f_yd,
                d,
                e_s,
                approx_lvl_p,
                v_ed,
                e_u,
                inner,
                edge_par,
                edge_per,
                corner,
                m_rd,
                m_pd,
            )
            / 6
        )
        * (sin(alpha * pi / 180) + cos(alpha * pi / 180))
        * (sin(alpha * pi / 180) + f_bd * d / (f_ywd * phi_w)),
        f_ywd,
    )

    if (a_sw * k_e * f_ywd) < 0.5 * v_ed:
        warnings.warn(
            """In order to ensure sufficent deformation capacity,
            consider increasing the amount of punching shear reinforcement"""
        )
    return a_sw * k_e * sigma_swd * sin(alpha * pi / 180)


def v_rd_max_punching(
    l_x: float,
    l_y: float,
    f_yd: float,
    d: float,
    e_s: float,
    approx_lvl_p: float,
    v_ed: float,
    e_u: float,
    inner: bool,
    edge_par: bool,
    edge_per: bool,
    dg: float,
    corner: bool,
    m_rd: float,
    m_pd: float,
    v_prep_d_max: float,
    d_v: float,
    f_ck: float,
    d_head: bool,
    stirrups_compression: bool,
    gamma_c: float = 1.5,
) -> float:
    """Finds the maximum value you can have for v_rd_punching.

    fib Model Code 2010, eq. (7.3-68) and (7.3-69).

    Args:
        l_x (float): The distance between two columns in x direction.
        l_y (float): The distance between two columns in y direction.
        f_yd (float): Design strength of reinforment steel in MPa.
        d (float): The mean value of the effective depth in mm.
        e_s (float): The E_modulus for steel in MPa.
        approx_lvl_p (float): The approx level for punching.
        v_ed (float): The acting shear force from the columns.
        e_u (float): Refers to the eccentricity of the resultant of shear
            forces with respect to the centroid.
        inner (bool): Is true only if the column is a inner column.
        edge_par (bool): Is true only if the column is a edge column with
            tension reinforcement parallel to the edge.
        edge_per (bool): Is true only if the column is a edge column with
            tension reinforcement perpendicular to the edge.
        dg (float): Maximum size of aggregate.
        corner (bool): Is true only if the column is a corner column.
        m_rd (float): The design average strength per unit length in MPa.
        m_pd: (float): The average decompresstion moment due to prestressing in
            MPa.
        v_prep_d_max (float): The maximum shear force per unit length
            perpendiculerer to the basic control parameter (Figure 7.3-24).
        d_v (float): The effective depth considering support in mm.
        f_ck (float): Characteristic strength in MPa.
        d_head (bool): True if diameter of heads is three times larger than.
        stirrups_compression: (bool): Stirrups with sufficient length at
            compression face, and bent on tension face.
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

    k_dg = max(32 / (16 + dg), 0.75)
    k_psi = min(
        1
        / (
            1.5
            + 0.9
            * k_dg
            * d
            * psi_punching(
                l_x,
                l_y,
                f_yd,
                d,
                e_s,
                approx_lvl_p,
                v_ed,
                e_u,
                inner,
                edge_par,
                edge_per,
                corner,
                m_rd,
                m_pd,
            )
        ),
        0.6,
    )
    return min(
        (k_sys * k_psi * b_0(v_ed, v_prep_d_max) * d_v * f_ck**0.5) / gamma_c,
        (b_0(v_ed, v_prep_d_max) * d_v * f_ck**0.5) / gamma_c,
    )


def v_rd_punching(
    e_u: float,
    b_u: float,
    l_x: float,
    l_y: float,
    f_yd: float,
    d: float,
    e_s: float,
    approx_lvl_p: float,
    v_ed: float,
    inner: bool,
    edge_par: bool,
    edge_per: bool,
    corner: bool,
    m_rd: float,
    m_pd: float,
    alpha: float,
    f_bd: float,
    f_ywk: float,
    phi_w: float,
    a_sw: float,
    dg: float,
    f_ck: float,
    d_v: float,
    v_prep_d_max: float,
    d_head: bool,
    stirrups_compression: bool,
    gamma_c: float = 1.5,
    gamma_s: float = 1.15,
) -> float:
    """The total resistance for punching, both Vrd,c and Vrd,s.

    fib Model Code 2010, eq. (7.3-60).

    Args:
        e_u (float): The ecentrisity of the result of shear forces with respect
            to the centroid (Figure 7.3-27b).
        b_u (float): The diamter of a circle with same surface as the region
            inside the basic control perimeter (Figure 7.3-27b).
        l_x (float): The distance between two columns in x direction.
        l_y (float): The distance between two columns in y direction.
        f_yd (float): Design strength of reinforment steel in MPa.
        d (float): The mean value of the effective depth in mm.
        e_s (float): The E_modulus for steel in MPa.
        approx_lvl_p (float): The approx level for punching.
        v_ed (float): The acting shear force from the columns.
        inner (bool): Is true only if the column is a inner column.
        edge_par (bool): Is true only if the column is a edge column with
            tension reinforcement parallel to the edge.
        edge_per (bool): Is true only if the column is a edge column with
            tension reinforcement perpendicular to the edge.
        corner (bool): Is true only if the column is a corner column.
        m_rd (float): The design average strength per unit length in MPa.
        m_pd: (float): The average decompresstion moment due to prestressing in
            MPa.
        alpha (float): Inclination of the stirrups in degrees.
        f_bd (float): The design bond strength in MPa.
        f_ywk (float): Characteristic yield strength of the shear reinforcement
            in MPa.
        phi_w (float): The diameter of the shear reinforcement.
        a_sw (float): The area of the shear reinforcement in mm^2.
        dg (float): Maximum size of aggregate.
        f_ck (float): Characteristic strength in MPa.
        d_v (float): The effective depth considering support in mm.
        v_prep_d_max (float): The maximum shear force per unit length
            perpendicular to the basic control parameter (Figure 7.3-24).

    Return:
        float: The maximum allowed punching resistance, regardless of values
        from v_rdc and v_rds.
    """
    return min(
        v_rdc_punching(
            l_x,
            l_y,
            f_yd,
            d,
            e_s,
            approx_lvl_p,
            dg,
            f_ck,
            d_v,
            v_ed,
            e_u,
            inner,
            edge_par,
            edge_per,
            corner,
            m_rd,
            v_prep_d_max,
            gamma_c,
        )
        + v_rds_punching(
            e_u,
            b_u,
            l_x,
            l_y,
            f_yd,
            d,
            e_s,
            approx_lvl_p,
            v_ed,
            inner,
            edge_par,
            edge_per,
            corner,
            m_rd,
            m_pd,
            alpha,
            f_bd,
            f_ywk,
            phi_w,
            a_sw,
            gamma_s,
        ),
        v_rd_max_punching(
            l_x,
            l_y,
            f_yd,
            d,
            e_s,
            approx_lvl_p,
            v_ed,
            e_u,
            inner,
            edge_par,
            edge_per,
            dg,
            corner,
            m_rd,
            m_pd,
            v_prep_d_max,
            d_v,
            f_ck,
            d_head,
            stirrups_compression,
            gamma_c,
        ),
    )
