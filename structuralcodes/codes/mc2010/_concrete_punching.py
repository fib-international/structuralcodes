"Covers punching in Model code 2010, 7.3.5.1 to 7.3.5.4"
import warnings
from math import pi, sin, cos


def b_0(v_ed: float, v_prep_d_max: float) -> float:
    """Gives the general output for b_0, shear-resisting control perimeter.

    fib Model Code 2010, eq. (7.3-57)
    Args:
        V_ed (float): The acting shear force from the columns
        v_prep_d_max (float): The maximum shear force per unit length
        perpendiculerer to the basic control parameter (Figure 7.3-24)

    Return:
        The shear-resisting control perimeter, b_0"""

    return v_ed / v_prep_d_max


def m_ed(
    v_ed: float,
    e_u: float,
    l_x: float,
    l_y: float,
    l_min: float,
    inner: bool,
    edge_par1: bool,
    edge_per2: bool,
    corner: bool,
) -> float:
    """The average bending moment acting in the support strip.

    fib Model Code 2010, eq. (7.3-76), (7.3-71), (7.3-72), (7.3-73)
    and (7.3-74)

    Args:
        v_ed (float): The acting shear force from the columns
        e_u (float): Refers to the eccentricity of the resultant of shear
        forces with respect to the centroid
        r_sx (float): Denotes the position where the radial bending moment is
        zero with respect to support axis in x direction
        r_sy (float): Denotes the position where the radial bending moment is
        zero with respect to support axis in x direction
        l_min (float): The shorter side of the the L_x and L_y
        inner (bool): Is true only if the column is a inner column
        edge_par (bool): Is true only if the column is a edge column with
        tention reinforcement paralell to the edge
        edge_per (bool): Is true only if the column is a edge column with
        tention reinforcement perpendicular to the edge
        corner (bool): Is true only if the column is a corner column

    return:
        The bending moment acting in the support strip regardless of the
        position of the column
    """
    r_sx = 0.22 * l_x
    r_sy = 0.22 * l_y
    b_s = min(1.5 * (r_sx * r_sy) ** 0.5, l_min)
    if inner:
        return v_ed * ((1 / 8) + e_u / (2 * b_s))
    if edge_par1:
        return max(v_ed * ((1 / 8) + e_u / (2 * b_s)), v_ed / 4)
    if edge_per2:
        return v_ed * ((1 / 8) + e_u / (b_s))
    if corner:
        return max(v_ed * ((1 / 8) + e_u / (b_s)), v_ed / 2)
    raise ValueError("the placement is not defined, one needs to be True")


def psi_punching(
    l_x: float,
    l_y: float,
    f_yd: float,
    d: float,
    e_s: float,
    approx_lvl_p: float,
    v_ed: float,
    e_u: float,
    l_min: float,
    inner: bool,
    edge_par: bool,
    edge_per: bool,
    corner: bool,
    m_rd: float,
    m_pd: float,
) -> float:
    """The rotation of the slab around the supported area

    fib Model Code 2010, eq. (7.3-70), (7.3-75) and (7.3-77)
    args:
        r_s (float): Denotes the position where the radial bending moment is
        zero with respect to support axis
        l_x (float): The distance between two columns in x direction
        l_y (float): The distance between two columns in y direction
        f_yd (float): Design strength of reinforment steel in MPa
        d (float): The mean value of the effective depth in mm
        e_s (float): The E_s-modulus for steel in Mpa
        approx_lvl_p (float): The approx level for punching
        v_ed (float): The acting shear force from the columns
        e_u (float): Refers to the eccentricity of the resultant of shear
        forces with respect to the centroid
        r_sx (float): Denotes the position where the radial bending moment is
        zero with respect to support axis in x direction
        r_sy (float): Denotes the position where the radial bending moment is
        zero with respect to support axis in x direction
        l_min (float): The shorter side of the the L_x and L_y
        inner (bool): Is true only if the column is a inner column
        edge_par (bool): Is true only if the column is a edge column with
        tention reinforcement paralell to the edge
        edge_per (bool): Is true only if the column is a edge column with
        tention reinforcement perpendicular to the edge
        corner (bool): Is true only if the column is a corner column
        m_rd (float): The design average strength per unit length in MPa
        m_pd: (float): The average decompresstion moment due to prestressing
        in MPa

    return:
        psi for the choosen approx level in punching"""

    r_s = max(0.22 * l_x, 0.22 * l_y)

    if approx_lvl_p == 1:
        if not 0.5 < l_x / l_y < 2:
            warnings.warn("Reconsider maximum r_s value")
        psi = 1.5 * r_s * f_yd / (d * e_s)

    elif approx_lvl_p in (2, 3):
        if 0.5 < l_x / l_y < 2:
            warnings.warn("Reconsider maximum r_s value")
        psi = (1.5 * r_s * f_yd / (d * e_s)) * (
            (
                m_ed(
                    v_ed,
                    e_u,
                    l_x,
                    l_y,
                    l_min,
                    inner,
                    edge_par,
                    edge_per,
                    corner,
                )
                - m_pd
            )
            / (m_rd - m_pd)
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
    l_min: float,
    inner: bool,
    edge_par: bool,
    edge_per: bool,
    corner: bool,
    m_rd: float,
    m_pd: float,
    v_prep_d_max: float,
    gamma_c: float = 1.5,
) -> float:
    """Punching resistance from the concrete
    fib Model Code 2010, eq. (7.3-61), (7.3-62) and (7.3-63)
    args:
        l_x (float): The distance between two columns in x direction
        l_y (float): The distance between two columns in y direction
        f_yd (float): Design strength of reinforment steel in MPa
        d (float): The mean value of the effective depth in mm
        e_s (float): The E_s-modulus for steel in Mpa
        approx_lvl_p (float): The approx level for punching
        dg (float): Maximum size of aggregate
        f_ck (float): Characteristic strength in MPa
        d_v (float): The effective depth considering support in mm
        v_ed (float): The acting shear force from the columns
        e_u (float): Refers to the eccentricity of the resultant of shear
        forces with respect to the centroid
        r_sx (float): Denotes the position where the radial bending moment is
        zero with respect to support axis in x direction
        r_sy (float): Denotes the position where the radial bending moment is
        zero with respect to support axis in x direction
        l_min (float): The shorter side of the the L_x and L_y
        inner (bool): Is true only if the column is a inner column
        edge_par (bool): Is true only if the column is a edge column with
        tention reinforcement paralell to the edge
        edge_per (bool): Is true only if the column is a edge column with
        tention reinforcement perpendicular to the edge
        corner (bool): Is true only if the column is a corner column
        m_rd (float): The design average strength per unit length in MPa
        m_pd: (float): The average decompresstion moment due to prestressing
        in MPa
    return:
        v_rdc for punching with the right approx level"""

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
                l_min,
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
    result = k_psi * b_0(v_ed, v_prep_d_max) * d_v * (f_ck**0.5) / gamma_c
    return result


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
    l_min: float,
    inner: bool,
    edge_par: bool,
    edge_per: bool,
    corner: bool,
    m_rd: float,
    m_pd: float,
    alfa: float,
    f_bd: float,
    f_ywd: float,
    kam_w: float,
    a_sw: float,
):
    """The punching resistance from shear reinforcement
     fib Model Code 2010, eq. (7.3-64) and (7.3-65)
    Args:
        e_u (float): The ecentrisity of the result of shear forces
        with respect to the centroid (Figure 7.3-27b)
        b_u (float): The diamter of a circle with same surface as the
        region inside the basic control perimeter (Figure 7.3-27b)
        r_s (float): Denotes the position where the radial bending moment is
        zero with respect to support axis
        l_x (float): The distance between two columns in x direction
        l_y (float): The distance between two columns in y direction
        f_yd (float): Design strength of reinforment steel in MPa
        d (float): The mean value of the effective depth in mm
        e_s (float): The E_s-modulus for steel in Mpa
        approx_lvl_p (float): The approx level for punching
        v_ed (float): The acting shear force from the columns
        r_sx (float): Denotes the position where the radial bending moment is
        zero with respect to support axis in x direction
        r_sy (float): Denotes the position where the radial bending moment is
        zero with respect to support axis in x direction
        l_min (float): The shorter side of the the L_x and L_y
        inner (bool): Is true only if the column is a inner column
        edge_par (bool): Is true only if the column is a edge column with
        tention reinforcement paralell to the edge
        edge_per (bool): Is true only if the column is a edge column with
        tention reinforcement perpendicular to the edge
        corner (bool): Is true only if the column is a corner column
        m_rd (float): The design average strength per unit length in MPa
        m_pd: (float): The average decompresstion moment due to prestressing
        in MPa
        alfa (float): Inclination of the stirrups in degrees
        f_bd (float): The design bond strength in MPa
        f_ywd (float): Design yield strength of the shear reinforcement in Mpa
        kam_w (float): The diameter of the shear reinforcement
        a_sw (float): The area of the shear reinforcement in mm^2


    return: Punching resistance that comes from reinforcement
    """
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
                l_min,
                inner,
                edge_par,
                edge_per,
                corner,
                m_rd,
                m_pd,
            )
            / 6
        )
        * (sin(alfa * pi / 180) + cos(alfa * pi / 180))
        * (sin(alfa * pi / 180) + f_bd * d / (f_ywd * kam_w)),
        f_ywd,
    )

    if (a_sw * k_e * sigma_swd * sin(alfa * pi / 180)) < 0.5 * v_ed:
        warnings.warn(
            """In order to ensure sufficent deformation capacity,
                      the shear resistance in punching most increase"""
        )
    result = a_sw * k_e * sigma_swd * sin(alfa * pi / 180)
    return result


def v_rd_max_punching(
    l_x: float,
    l_y: float,
    f_yd: float,
    d: float,
    e_s: float,
    approx_lvl_p: float,
    v_ed: float,
    e_u: float,
    l_min: float,
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
):
    """Finds the maximum value you can have for v_rd
    fib Model Code 2010, eq. (7.3-68) and (7.3-69)
    Args:
        l_x (float): The distance between two columns in x direction
        l_y (float): The distance between two columns in y direction
        f_yd (float): Design strength of reinforment steel in MPa
        d (float): The mean value of the effective depth in mm
        e_s (float): The E_s-modulus for steel in Mpa
        approx_lvl_p (float): The approx level for punching
        dg (float): Maximum size of aggregate
        d_v (float): The effective depth considering support in mm
        v_ed (float): The acting shear force from the columns
        e_u (float): Refers to the eccentricity of the resultant of shear
        forces with respect to the centroid
        l_min (float): The shorter side of the the L_x and L_y
        inner (bool): Is true only if the column is a inner column
        edge_par (bool): Is true only if the column is a edge column with
        tention reinforcement paralell to the edge
        edge_per (bool): Is true only if the column is a edge column with
        tention reinforcement perpendicular to the edge
        corner (bool): Is true only if the column is a corner column
        m_rd (float): The design average strength per unit length in MPa
        m_pd: (float): The average decompresstion moment due to prestressing
        in MPa
        stirrups_compression: (bool): Stirrups with sufficient length at
        compression face, and bent on tension face

    Return
        The maximum allowed
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
                l_min,
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
    result = min(
        (k_sys * k_psi * b_0(v_ed, v_prep_d_max) * d_v * f_ck**0.5)
        / gamma_c,
        (b_0(v_ed, v_prep_d_max) * d_v * f_ck**0.5) / gamma_c,
    )
    return result


def v_rd_punching(
    e_u,
    b_u,
    l_x: float,
    l_y: float,
    f_yd: float,
    d: float,
    e_s: float,
    approx_lvl_p: float,
    v_ed: float,
    l_min: float,
    inner: bool,
    edge_par: bool,
    edge_per: bool,
    corner: bool,
    m_rd: float,
    m_pd: float,
    alfa: float,
    f_bd: float,
    f_ywd: float,
    kam_w: float,
    a_sw: float,
    dg: float,
    f_ck: float,
    d_v: float,
    v_prep_d_max: float,
    d_head: bool,
    stirrups_compression: bool,
    gamma_c: float = 1.5,
):
    """The total resistance for punching, both Vrd,c and Vrd,s
     fib Model Code 2010, eq. (7.3-60)
        Args:
        e_u (float): The ecentrisity of the result of shear forces
        with respect to the centroid (Figure 7.3-27b)
        b_u (float): The diamter of a circle with same surface as the
        region inside the basic control perimeter (Figure 7.3-27b)
        r_s (float): Denotes the position where the radial bending moment is
        zero with respect to support axis
        l_x (float): The distance between two columns in x direction
        l_y (float): The distance between two columns in y direction
        f_yd (float): Design strength of reinforment steel in MPa
        d (float): The mean value of the effective depth in mm
        e_s (float): The E_s-modulus for steel in Mpa
        approx_lvl_p (float): The approx level for punching
        v_ed (float): The acting shear force from the columns
        r_sx (float): Denotes the position where the radial bending moment is
        zero with respect to support axis in x direction
        r_sy (float): Denotes the position where the radial bending moment is
        zero with respect to support axis in x direction
        l_min (float): The shorter side of the the L_x and L_y
        inner (bool): Is true only if the column is a inner column
        edge_par (bool): Is true only if the column is a edge column with
        tention reinforcement paralell to the edge
        edge_per (bool): Is true only if the column is a edge column with
        tention reinforcement perpendicular to the edge
        corner (bool): Is true only if the column is a corner column
        m_rd (float): The design average strength per unit length in MPa
        m_pd: (float): The average decompresstion moment due to prestressing
        in MPa
        alfa (float): Inclination of the stirrups in degrees
        f_bd (float): The design bond strength in MPa
        f_ywd (float): Design yield strength of the shear reinforcement in Mpa
        kam_w (float): The diameter of the shear reinforcement
        a_sw (float): The area of the shear reinforcement in mm^2

    return: The maximum allowed punching resistance, regardless of
    values from v_rdc and v_rds"""

    result = min(
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
            l_min,
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
            l_min,
            inner,
            edge_par,
            edge_per,
            corner,
            m_rd,
            m_pd,
            alfa,
            f_bd,
            f_ywd,
            kam_w,
            a_sw,
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
            l_min,
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
    return result
