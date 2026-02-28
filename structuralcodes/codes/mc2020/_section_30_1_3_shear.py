"""Functions from Section 30.1.3 of fib Model Code 2020."""

import math
import typing as t


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                           Section 30.1.3.1                          #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


def unity_check_V(V_Rd: float, V_Ed: float) -> float:
    """Perform a unity check on shear.

    fib Model Code 2020: Eqs. (30.1-7) and (30.1-53).

    Args:
        V_Rd (float): Design shear resistance.
        V_Ed (float): Design shear force.

    Returns:
        float: The unity check (unitless).
    """
    return abs(V_Ed) / V_Rd


def k_dir(a_v: float, d: float) -> float:
    """Calculate the reduction factor k_dir for direct strut or arch action.
    N.B. Equation 30.1-8 is changed w.r.t. how it is written in the online
    documentation (k_dir = max(.) instead of k_dir = min(.).

    fib Model Code 2020: Eq. (30.1-8).

    Args:
        a_v (float): Effective shear span.
        d (float): Effective depth.

    Returns:
        float: The coefficient k_dir (unitless).
    """
    if a_v > 2 * d:
        k_dir = 1.0
    elif a_v < d:
        k_dir = 0.5
    else:
        k_dir = max(a_v / (2 * d), 0.5)

    return k_dir


def z_v(z_s: float, z_p: float, A_s: float, A_p: float, h: float) -> float:
    """Calculate the effective shear depth z_v for members containing both
    ordinary steel reinforcement and prestressed tendons. z_v may also be
    assumed as 0.9d.

    fib Model Code 2020: Eq. (30.1-9).

    Args:
        z_s (float): Distance between centerline compressive chord and
        reinforcement axis.
        z_p (float): Distance between centerline compressive chord and tendon
        axis.
        A_s (float): Area of reinforcement.
        A_p (float): Area of prestressing.
        h (float): Height of the member.

    Returns:
        float: The effective shear depth z_v.
    """
    z_v = (z_s**2 * A_s + z_p**2 * A_p) / (z_s * A_s + z_p * A_p)

    return max(z_v, 0.72 * h)


def eps_x(
    M_Ed: float,
    V_Ed: float,
    N_Ed: float,
    E_s: float,
    A_s: float,
    z_v: float,
    cot_theta: float,
    delta_e: float,
    E_c: float = 0.0,
    A_c_ten: float = 0.0,
) -> float:
    """Calculate the longitudinal strain at the mid-depth of the effective
    shear depth (epsilon_x).

    fib Model Code 2020: Eq. (30.1-10).

    N.B. Strictly speaking Eq. (30.1-10) requires non-negative values as
    outcome; however, in the description it seems that negative values are
    allowed.

    Args:
        M_Ed (float): Absolute value bending moment.
        V_Ed (float): Absolute value shear force.
        N_Ed (float): Axial force. Positive value for tension and negative
        value for compression.
        E_s (float): Modulus of elasticity of steel reinforcement.
        A_s (float): Area of reinforcement.
        z_v (float): Effective shear depth.
        cot_theta (float): Cotangent inclination of compression field in
        radians.
        delta_e (float): Distance between centre of gravity and mid-depth of
        the effective shear depth.
        E_c (float): Modulus of elasticity of concrete.
        A_c_ten (float): Area of tension chord due to bending.

    Returns:
        float: epsilon_x (strain).
    """

    force_rebar = (
        M_Ed / z_v + 0.5 * V_Ed * cot_theta + N_Ed * (0.5 - delta_e / z_v)
    )

    if force_rebar < 0.0:
        eps_x = force_rebar / (2.0 * (E_c * A_c_ten + E_s * A_s))
    else:
        eps_x = force_rebar / (2.0 * E_s * A_s)

    return max(eps_x, 0.0)


def eps_x_bond(
    M_Ed: float,
    V_Ed: float,
    N_Ed: float,
    E_s: float,
    E_p: float,
    A_s: float,
    A_p: float,
    z_v: float,
    z_s: float,
    z_p: float,
    cot_theta: float,
    e_p: float,
) -> float:
    """Calculate the longitudinal strain at the mid-depth of the effective
    shear depth (epsilon_x) for prestressed members with bonded tendons.

    fib Model Code 2020: Eq. (30.1-11).

    Args:
        M_Ed (float): Absolute value bending moment.
        V_Ed (float): Absolute value shear force.
        N_Ed (float): Axial force. Positive value for tension and negative
        value for compression.
        E_s (float): Modulus of elasticity of steel reinforcement.
        E_p (float): Modulus of elasticity of prestressing reinforcement.
        A_s (float): Area of reinforcement.
        A_p (float): Area of prestressing.
        z_v (float): Effective shear depth.
        z_s (float): Distance between centreline compressive chord and
        reinforcement axis.
        z_p (float): Distance between centreline compressive chord and tendon
        axis.
        cot_theta (float): Cotangent inclination of compression field in
        radians.
        e_p (float): Eccentricity prestressing force.

    Returns:
        float: epsilon_x (strain).

    Raises:
        ValueError: If any of the input parameters are negative.
    """
    epsilon_x = (
        M_Ed / z_v + 0.5 * V_Ed * cot_theta + N_Ed * (z_p - e_p) / z_v
    ) / (2 * ((z_s / z_v) * E_s * A_s + (z_p / z_v) * E_p * A_p))

    return max(epsilon_x, 0)


def M_Ed(M_Ed0: float, M_Pd: float) -> float:
    """Calculate the design value of the applied bending moment (M_Ed).

    fib Model Code 2020: Eq. (30.1-12a).

    Args:
        M_Ed0 (float): Design value of applied moment resulting from
        sectional analysis.
        M_Pd (float): Design bending moment due to prestressing.

    Returns:
        float: The design value of the applied bending moment M_Ed.
    """
    return M_Ed0 + M_Pd


def N_Ed(N_Ed0: float, F_p: float, delta_p: float) -> float:
    """Calculate the design value of the applied axial force (N_Ed).

    fib Model Code 2020: Eq. (30.1-12b).

    Args:
        N_Ed0 (float): Design value of applied axial force resulting from
        sectional analysis.
        F_p (float): Prestressing force.
        delta_p (float): Angle between tendon and the neutral axis of the
        member in radians.

    Returns:
        float: The design value of the applied axial force N_Ed.
    """
    return N_Ed0 - F_p * math.cos(delta_p)


def V_Ed(V_Ed0: float, F_p: float, delta_p: float) -> float:
    """Calculate the design value of the shear force (V_Ed).

    fib Model Code 2020: Eq. (30.1-12c).

    Args:
        V_Ed0 (float): Design value of shear force resulting from sectional
        analysis.
        F_p (float): Prestressing force.
        delta_p (float): Angle between tendon and the neutral axis of the
        member in radians.

    Returns:
        float: The design value of the shear force V_Ed.
    """
    return V_Ed0 - F_p * math.sin(delta_p)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                           Section 30.1.3.2                          #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
def V_Rdc_1(
    k_v: float, f_ck: float, gamma_C: float, b_w: float, z_v: float
) -> float:
    """Calculate the design shear resistance of a web or slab without
        shear reinforcement.

    fib Model Code 2020: Eqs. (30.1-13), (30.1-23), (30.1-55).

    Args:
        k_v (float): Multiplication factor.
        f_ck (float): Characteristic compressive strength of concrete
            in MPa. In case z_v >= 800 mm sqrt(f_ck) is max 8 MPa.
        gamma_C (float): Partial factor for concrete (unitless).
        b_w (float): Effective width in mm.
        z_v (float): Effective shear depth in mm.

    Returns:
        float: The design value of the shear resistance without shear
            reinforcement in N.
    """
    if z_v >= 800 and math.sqrt(f_ck) > 8:
        f_ck = 8**2
    return k_v * math.sqrt(f_ck) / gamma_C * b_w * z_v


def k_v_LOI_1(z_v: float, f_yd: float, E_s: float) -> float:
    """Multiplcation factor for calculating concrete shear resistance.

    fib Model Code 2020: Eq. (30.1-15).

    Args:
        z_v (float): Effective shear depth in mm.
        f_yd (float): Design value of reinforcement yield strength in
            MPa.
        E_s (float): Young's modulus of reinforcement steel in MPa.

    Returns:
        float: Multiplcation factor for calculating concrete shear
            resistance (unitless).
    """
    return 0.4 / (1 + 750 * f_yd / E_s) * 1300 / (1000 + 1.25 * z_v)


def k_v_LOI_2(z_v: float, k_dg: float, eps_x: float) -> float:
    """Multiplcation factor for calculating concrete shear resistance.

    fib Model Code 2020: Eq. (30.1-16).

    Args:
        z_v (float): Effective shear depth in mm.
        k_dg (float): multiplication factor to account for roughness of
            critical shear crack.
        eps_x (float): Longitudinal strain at mid-depth of the
            effective shear depth (unitless).

    Returns:
        float: Multiplcation factor for calculating concrete shear
            resistance (unitless).
    """
    return 0.4 / (1 + 1500 * eps_x) * 1300 / (1000 * k_dg * z_v)


def k_dg(
    d_g: float,
    z_v: float,
    f_ck: float,
    lightweight_aggregate: bool = False,
) -> float:
    """Calculate multiplication factor to account for roughness of
        critical shear crack.

    fib Model Code 2020: Eq. (30.1-17).

    Args:
        d_g (float): Aggregate size in mm.
        z_v (float): Effective shear depth in mm.
        f_ck (float): Characteristic compressive strength of concrete
            in MPa.
        lightweight_aggregate (bool): Concrete based on lightweight
            aggregates is used. Default value = False.

    Returns:
        float: Multiplication factor to account for roughness of
            critical shear crack.
    """
    if z_v >= 800 and lightweight_aggregate:
        return 2
    elif z_v >= 800 and f_ck >= 70:
        return 2
    return max(32 / (16 + d_g), 0.75)


def V_Rdc_2(
    a: float,
    d: float,
    f_ck: float,
    rho_l: float,
    k_dg: float,
    z: float,
    gamma_C: float,
    b_w: float,
) -> float:
    """Calculate the design shear resistance of a web or slab without
        shear reinforcement. This equation can be derived from Equations
        (30.1-13) to (30.1-17).

    fib Model Code 2020: Eq. (30.1-18).

    Args:
        a (float): Shear span in mm. Can be assumed as M_Ed / V_Ed + 1.
        d (float): Effective depth in mm.
        f_ck (float): Characteristic compressive strength of concrete
            in MPa.
        rho_l (float): Longitudinal reinforcement ratio (unitless).
        k_dg (float): Multiplication factor to account for roughness of
            critical shear crack (unitless).
        z (float): Internal lever arm in mm.
        gamma_C (float): Partial factor for concrete (unitless).
        b_w (float): Effective width in mm.

    Returns:
        float: Design shear resistance of a web or slab without shear
            reinforcement in N.
    """
    return (
        (
            math.sqrt(
                1 + 7.8 * a / d * math.sqrt(f_ck) / (rho_l * (1000 * k_dg * z))
            )
            - 1
        )
        * 120
        * rho_l
        / (a / d)
        * 1
        / gamma_C
        * b_w
        * z
    )


def V_Rdc_FprEC2(
    gamma_V: float,
    rho_l: float,
    E_s: float,
    f_ck: float,
    d_dg: float,
    a_cs: float,
    d: float,
    b_w: float,
    z: float,
) -> float:
    """Calculate the design shear resistance of a web or slab without
        shear reinforcement. This is based on the Critical Shear Crack
        Theory.

    fib Model Code 2020: Eq. (30.1-19).

    Args:
        gamma_V (float): Partial safety factor for shear (unitless).
        rho_l (float): Longitudinal reinforcement ratio (unitless).
        E_s (float): Young's modulus of reinforcement steel in MPa.
        f_ck (float): Characteristic compressive strength of concrete
            in MPa.
        d_dg (float): Parameter related to aggregate size related to d_dg.
        a_cs (float): Shear span in mm.
        d (float): Effective depth in mm.
        b_w (float): Effective width in mm.
        z (float): Internal lever arm in mm.

    Returns:
        Design shear resistance of a web or slab without shear
            reinforcement according to Critical Shear Crack Theory in
            N.
    """
    return (
        0.083
        / gamma_V
        * (rho_l * E_s * f_ck * d_dg / math.sqrt(a_cs * d)) ** (1 / 3)
        * b_w
        * z
    )


def d_dg(k_dg: float) -> float:
    """Parameter related to aggregate size d_g. Used when
        calculating V_Rdc according to the Critical Shear Crack Theory.

    Args:
        k_dg (float): Multiplication factor to account for roughness of
            critical shear crack (unitless).

    Returns:
        float: Parameter related to aggregate size d_g in
            mm.
    """
    return 32 / k_dg


def tan_theta_v(v_y: float, v_x: float) -> float:
    """Calculate the principle direction of the shear force in a planar
        member (slab and shell).

    fib Model Code 2020: Eq. (30.1-20).

    Args:
        v_y (float): Shear stress resultant per unit width in shell
            elements, y-direction.
        v_x (float): Shear stress resultant per unit width in shell
            elements, x-direction.

    Returns:
        float: The principle direction of the shear force in a planar
        member (slab and shell) (unitless).
    """
    return v_y / v_x


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                           Section 30.1.3.3                          #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
def rho_w(f_ck: float, f_yk: float) -> float:
    """Calculate demand for minimum shear reinforcement according to Section
    30.1.3.3.

    fib ModelCode 2020: Eq. (30.1-21).

    Args:
        f_ck (float): Characteristic cylinder compressive strength of
        concrete in MPa.
        f_yk (float): Characteristic yield stress of reinforcement in MPa.

    Returns:
        float: Demand for minimum shear reinforcement according to Section
        30.1.3.3.
    """
    return 0.08 * f_ck**0.5 / f_yk


def V_Rd(V_Rd_c: float, V_Rd_s: float, V_Rd_max: float) -> float:
    """Calculate the design shear resistance.

    fib ModelCode 2020: Eq. (30.1-22).

    Args:
        V_Rd_c (float): Design shear resistance attributed to concrete.
        V_Rd_s (float): Design shear resistance provided by shear reinforcement.
        V_Rd_max (float): Maximum shear resistance related to crushing of
        concrete carrying the compression field.

    Returns:
        float: The design shear resistance.
    """
    return min(V_Rd_c + V_Rd_s, V_Rd_max)


def V_Rd_s_orthog(
    A_sw: float,
    s_w: float,
    z_v: float,
    f_ywd: float,
    cot_theta: float,
) -> float:
    """Calculate the design shear resistance provided by shear reinforcement,
    with the shear reinforcement orthogonal to
    the beam axis.

    fib Model Code 2020: Eq. (30.1-24).

    Args:
        A_sw (float): Area of shear reinforcement.
        s_w (float): Stirrup spacing.
        z_v (float): Effective shear depth.
        f_ywd (float): Design yield strength of the shear reinforcement.
        cot_theta (float): Cotangent of theta, which describes the
        inclination of the compressive stress field
        depending on the LoA.

    Returns:
        float: The design shear resistance.

    """
    return A_sw / s_w * z_v * f_ywd * cot_theta


def V_Rd_max_orthog(
    k_eps: float,
    f_cd: float,
    b_w: float,
    z_v: float,
    cot_theta: float,
    tan_theta: float,
) -> float:
    """Calculate the maximum shear resistance related to crushing of concrete
    carrying the compression field, with the
    shear reinforcement orthogonal to the beam axis.

    fib Model Code 2020: Eq. (30.1-25).

    Args:
        k_eps (float): Strength reduction factor of concrete carrying the
        compression field which depends on LoA.
        f_cd (float): Design cylinder compressive strength of concrete.
        b_w (float): Web width.
        z_v (float): Effective shear depth.
        cot_theta (float): Cotangent of theta which describes the inclination
        of the compressive stress field
        depending on the LoA.
        tan_theta (float): Tangent of theta which describes the inclination
        of the compressive stress field
        depending on the LoA.

    Returns:
        float: The maximum shear resistance related to crushing of concrete
        carrying the compression field.

    Raises:
        ValueError: if alpha > 90 degrees.
    """
    return k_eps * f_cd * b_w * z_v / (cot_theta + tan_theta)


def V_Rd_s_nonorthog(
    A_sw: float,
    s_w: float,
    z_v: float,
    f_ywd: float,
    cot_theta: float,
    alpha: float,
) -> float:
    """Calculate the design shear resistance provided by shear reinforcement,
    with the shear reinforcement at a
    non-orthogonal angle with respect to the beam main axis.

    fib Model Code 2020: Eq. (30.1-26).

    Args:
        A_sw (float): Area of shear reinforcement.
        s_w (float): Stirrup spacing.
        z_v (float): Effective shear depth.
        f_ywd (float): Design yield strength of the shear reinforcement.
        cot_theta (float): Cotangent of theta, which describes the
        inclination of the compressive stress field
        depending on the LoA.
        alpha (float): Angle between the shear reinforcement w.r.t.
        longitudinal axis, with alpha > 90 degrees when
        shear reinforcement is inclined in the same direction as the
        compression field.

    Returns:
        float: The design shear resistance.

    Raises:
        ValueError: if alpha > 90 degrees.
    """
    if alpha < 90:
        return (
            A_sw
            / s_w
            * z_v
            * f_ywd
            * (cot_theta + 1 / math.tan(alpha))
            * math.sin(alpha)
        )
    else:
        raise ValueError(
            'Alpha > 90 degrees should be avoided due to potentially '
            'unacceptable cracking at SLS and stronger '
            'concrete softening at ULS.'
        )


def V_Rd_max_nonorthog(
    k_eps: float,
    f_cd: float,
    b_w: float,
    z_v: float,
    cot_theta: float,
    alpha: float,
) -> float:
    """Calculate the maximum shear resistance related to crushing of concrete
    carrying the compression field, with the
    shear reinforcement at a non-orthogonal angle with respect to the beam
    main axis.

    fib Model Code 2020: Eq. (30.1-27).

    Args:
        k_eps (float): Strength reduction factor of concrete carrying the
        compression field which depends on LoA.
        f_cd (float): Design cylinder compressive strength of concrete.
        b_w (float): Web width.
        z_v (float): Effective shear depth.
        cot_theta (float): Cotangent of theta which describes the inclination
        of the compressive stress field
        depending on the LoA.
        alpha (float): Angle between the shear reinforcement w.r.t.
        longitudinal axis, with alpha > 90 degrees when
        shear reinforcement is inclined in the same direction as the
        compression field.

    Returns:
        float: The maximum shear resistance related to crushing of concrete
        carrying the compression field.

    Raises:
        ValueError: if alpha > 90 degrees.
    """
    if alpha < 90:
        return (
            k_eps
            * f_cd
            * b_w
            * z_v
            * (cot_theta + 1 / math.tan(alpha))
            / (1 + cot_theta**2)
        )
    else:
        raise ValueError(
            'Alpha > 90 degrees should be avoided due to potentially '
            'unacceptable cracking at SLS and stronger '
            'concrete softening at ULS.'
        )


def b_w_nom(
    b_w: float,
    k_D: float,
    sum_phi_D: float,
) -> float:
    """Calculate nominal value of the web width.

     fib Model Code 2020: Eq. (30.1-28).

     Args:
        b_w (float): Web width.
        k_D (float): Parameter depending on material of duct and
        grouted/ungrouted. Suggested design values are: 0.5
        for grouted steel duct, 0.8 for grouted plastic duct, 1.2 for
        ungrouted duct.
        sum_phi_D (float): sum of the diameters of the ducts at the most
        unfavorable prestressing tendon configuration.

    Returns:
        float: The nominal value of the web width.
    """
    return b_w - k_D * sum_phi_D


def delta_F_td_orthog(
    V_Ed: float,
    cot_theta: float,
) -> float:
    """Additional force due to shear by the longitudinal reinforcement to be
    resisted by the flexural tensile chord,
    with the shear reinforcement orthogonal to the beam axis.

    fib Model Code 2020: Eq. (30.1-29).

    Args:
        V_Ed (float): Design shear force.
        cot_theta (float): Cotangent of theta which describes the inclination
        of the compressive stress field
        depending on the LoA.

    Returns:
        float: The additional force due to shear by the longitudinal
        reinforcement to be resisted by the flexural
        tensile chord.

    Raises:
        ValueError: if alpha > 90 degrees.
    """
    return V_Ed / 2 * cot_theta


def delta_F_td_nonorthog(V_Ed: float, cot_theta: float, alpha: float) -> float:
    """Additional force due to shear by the longitudinal reinforcement to be
    resisted by the flexural tensile chord,
    with the shear reinforcement at a non-orthogonal angle with respect to
    the beam main axis.

    fib Model Code 2020: Eq. (30.1-30).

    Args:
        V_Ed (float): Design shear force.
        cot_theta (float): Cotangent of theta which describes the inclination
        of the compressive stress field
        depending on the LoA.
        alpha (float): Angle between the shear reinforcement w.r.t.
        longitudinal axis, with alpha > 90 degrees when
        shear reinforcement is inclined in the same direction as the
        compression field.

    Returns:
        float: The additional force due to shear by the longitudinal
        reinforcement to be resisted by the flexural
        tensile chord.

    Raises:
        ValueError: if alpha > 90 degrees.
    """
    if alpha < 90:
        return V_Ed / 2 * (cot_theta - 1 / math.tan(alpha))
    else:
        raise ValueError(
            'Alpha > 90 degrees should be avoided due to potentially '
            'unacceptable cracking at SLS and stronger '
            'concrete softening at ULS.'
        )


def V_Ed_LoA3(V_Ed: float, V_Rd_c: float) -> float:
    """Calculate the design shear force to be used in Eq. (30.1-29) and Eq. (
    30.1-30) in the LoA III approach.

    fib Model Code 2020: Eq. (30.1-31).

    Args:
        V_Ed (float): Design shear force.
        V_Rd_c (float): Design shear resistance attributed to concrete.

    Returns:
        float: LoA III approach design shear force to be used in Eq. (
        30.1-29) and Eq. (30.1-30).
    """
    return V_Ed + V_Rd_c


def cot_theta_min_limits_LoA1_2a(cot_theta: float, cot_theta_min: float):
    """Checks the limits of cot_theta for LoA I and IIa.

    fib Model Code 2020: Eq. (30.1-32).

    Args:
        cot_theta (float): Compressive stress field inclination
        cot_theta_min (float): Minimum compressive stress field inclination
        according to fib Model Code 2020 Table
        30.1-1.

    Returns:
        bool: True if cot_theta fulfills condition of Eq. (30.1-32) and False
        if not.
    """
    if 1 <= cot_theta <= cot_theta_min:
        return True
    else:
        return False


def cot_theta_min_LoA1(load_condition) -> float:
    """Determines minimum inclination of the compressive stress field in the
    LoA I approach.

    fib Model Code 2020: Eq. (30.1-34), Eq. (30.1-35) and Eq. (30.1-36).

    Args:
        load_condition (str): Load condition as defined in Eq. (30.1-34),
        Eq. (30.1-35) and Eq. (30.1-36).

    Returns:
        float: Minimum inclination of the compressive stress field in the LoA
        I approach.
    """
    if load_condition == 'axial_compression':
        return 2.2
    if load_condition == 'rc_member':
        return 1.7
    if load_condition == 'axial_tension':
        return 1.2


def k_eps_LoA1():
    return 0.55


def cot_theta_simul_yield_crush(
    k_e: float,
    f_cd: float,
    rho_w: float,
    f_ywd: float,
    cot_theta_min,
) -> float:
    """Calculate the cot_theta where simultaneous yielding of the shear
    reinforcement and failure of the inclined
    compression field is expected.

    fib Model Code 2020: Eq. (30.1-38).

    Args:
        k_e (float): Strength reduction factor for V_Rd_max according to
        Table 30.1-1.
        f_cd (float): Design cylinder compressive strength of concrete.
        rho_w (float): Demand for minimum shear reinforcement according to
        Section 30.1.3.3.
        f_ywd (float): Design yield strength of the shear reinforcement.
        cot_theta_min (float): Minimum inclination according to Table 30.1-1.

    Returns:
        float: cot_theta where simultaneous yielding of the shear
        reinforcement and failure of the inclined
               compression field is expected.
    """
    cot_theta = ((k_e * f_cd) / (rho_w * f_ywd)) ** 0.5

    return max(1.0, min(cot_theta, cot_theta_min))


def V_Rd_LoA1_LoA2a(
    V_Rd_s: float,
    V_Rd_max: float,
) -> float:
    """Calculate the design shear resistance with the LoA I and LoA IIa
    approach.

    fib Model Code 2020: Eq. (30.1-33) and Eq. (30.1-).

    Args:
        V_Rd_s (float): Design shear resistance provided by shear reinforcement.
        V_Rd_max (float): Design maximum shear resistance related to crushing
        of concrete.

    Returns:
        float: LoA1 and LoA2a approach:Design shear resistance.
    """
    return min(V_Rd_s, V_Rd_max)


def cot_theta_min_LoA2a_A(
    eps_x: float,
) -> float:
    """Calculate the LoA 2a approach minimum inclination of the compressive
    stress field.

    fib Model Code 2020: Eq. (30.1-40).

    Args:
        eps_x (float): Longitudinal strain at mid-depth of the effective
        shear depth core layer as shown in Fig. 30.1-8.

    Returns:
        float: LoA IIa approach minimum inclination of the compressive stress
        field.
    """
    return 1 / math.tan(math.radians(20) + 4000 * eps_x)


def cot_theta_min_LoA2a_B(
    eps_x: float,
) -> float:
    """Calculate the LoA IIa approach reduced minimum inclination of the
    compressive stress field for ductility class
    (?), provided that cracking at SLS is verified explicitly.

    fib Model Code 2020: Eq. (30.1-41).

    Args:
        eps_x (float): Longitudinal strain at mid-depth of the effective
        shear depth core layer as shown in Fig. 30.1-8.

    Returns:
        float: LoA IIa approach minimum inclination of the compressive stress
        field.
    """
    return 1 / math.tan(math.radians(15) + 3000 * eps_x)


def cot_theta_min_LoA2a_CD(
    eps_x: float,
) -> float:
    """Calculate LoA IIa approach reduced minimum inclination of the
    compressive stress field for ductility classes C
    and D, provided that cracking at SLS is verified explicitly.

    fib Model Code 2020: Eq. (30.1-42).

    Args:
        eps_x (float): Longitudinal strain at mid-depth of the effective
        shear depth core layer as shown in Fig. 30.1-8.

    Returns:
        float: LoA IIa approach reduced minimum inclination of the
        compressive stress field for ductility classes C and
        D, provided that cracking at SLS is verified explicitly.
    """
    return 1 / math.tan(math.radians(13) + 2500 * eps_x)


def k_eps_LoAIIa(eps_1: float) -> float:
    """Strength reduction factor for V_Rd_max in LoA IIa approach.

    fib Model Code 2020: Eq. (30.1-43).

    Args:
        eps_1 (str): Strain calculated with Eq. (30.1-44).

    Returns:
        float: Strength reduction factor for V_Rd_max in LoA IIa approach.
    """
    return min(1.0, 1.2 + 60 * eps_1)


def eps_1(eps_x: float, cot_theta: float) -> float:
    """Strain value to calculate k_eps in the LoA IIa approach.

    fib Model Code 2020: Eq. (30.1-44).

    Args:
        eps_x (float): Longitudinal strain at mid-depth of the effective
        shear depth core layer as shown in Fig. 30.1-8.
        cot_theta (float): Compressive stress field inclination

    Returns:
        float: Strain value to calculate k_eps in the LoA IIa approach.
    """
    return eps_x + (eps_x + 0.001) * cot_theta**2


def V_Rd_LoA2b(
    V_Rd_s: float,
    V_Rd_c: float,
    V_Rd_max: float,
) -> float:
    """Calculate the design shear resistance with the LoA IIb approach.

    fib Model Code 2020: Eq. (30.1-45).

    Args:
        V_Rd_s (float): Design shear resistance provided by shear reinforcement.
        V_Rd_s (float): Design shear resistance provided by concrete.
        V_Rd_max (float): Design maximum shear resistance related to crushing
        of concrete.

    Returns:
        float: Design shear resistance in LoA IIb approach.
    """
    return min(V_Rd_s + V_Rd_c, V_Rd_max)


def cot_theta_LoA2b_nocrush(eps_x: float) -> float:
    """Calculate cot_theta in case V_Rd is not governed by V_Rd_max, but by
    V_Rd_c + V_Rd_s, in the LoA IIb approach.

    fib Model Code 2020: Eq. (30.1-46).

    Args:
        eps_x (float): Longitudinal strain at mid-depth of the effective
        shear depth core layer as shown in Fig. 30.1-8.

    Returns:
        float: Cot_theta in case V_Rd is not governed by V_Rd_max, but by
        V_Rd_c + V_Rd_s, in the LoA IIb approach.
    """
    return 1 / math.tan(math.radians(29) + 7000 * eps_x)


def k_v_LoA2b(eps_x: float) -> float:
    """Calculate factor k_v in the design shear resistance attributed to
    concrete in LoA IIb approach.

    fib Model Code 2020: Eq. (30.1-47).

    Args:
        eps_x (float): Longitudinal strain at mid-depth of the effective
        shear depth core layer as shown in Fig. 30.1-8.

    Returns:
        float: k_v in the design sheat resistance V_Rd_c in LoA IIb approach.
    """
    return 0.4 / (1 + 1500 * eps_x)


def V_Rd_direct(
    k_eps: float,
    f_cd: float,
    cot_theta: float,
    cot_beta: float,
    b_w: float,
    z_v: float,
    A_sw: float,
    s_w: float,
    f_ywd: float,
    V_Rd_max: float,
) -> float:
    """Calculate the design shear resistance in case loads are carried
    directly to a support through strut or arch
    action.

    fib Model Code 2020: Eq. (30.1-48).

    Args:
        eps_x (float): Longitudinal strain at mid-depth of the effective
        shear depth core layer as shown in Fig. 30.1-8.
        f_cd (float): Design cylinder compressive strength of concrete.
        cot_theta (float): Cotangent of theta which describes the inclination
        of the compressive stress field
        depending on the LoA.
        cot_beta (float): Cotangent of beta as shown in Fig. 30.1-9
        b_w (float): Web width.
        z_v (float): Effective shear depth.
        A_sw (float): Area of shear reinforcement.
        s_w (float): Stirrup spacing.
        f_ywd (float): Design yield strength of the shear reinforcement.
        V_Rd_max (float): Design maximum shear resistance related to crushing
        of concrete.

    Returns:
        float: The design shear resistance in case loads are carried directly
        to a support through strut or arch
               action.
    """
    return min(
        V_Rd_max,
        (
            k_eps
            * f_cd
            * (cot_theta - cot_beta)
            / (1 + cot_theta**2)
            * b_w
            * z_v
            + A_sw / s_w * z_v * f_ywd * cot_beta
        ),
    )


def delta_M_Ed(
    V_Ed: float,
    A_sw: float,
    s_w: float,
    z_v: float,
    f_ywd: float,
    cot_theta: float,
    a: float,
    x: float,
) -> float:
    """Calculate the increment moment due to shear in case loads are carried
    directly to a support through strut or arch
        action.

    fib Model Code 2020: Eq. (30.1-49).

    Args:
        V_Ed (float): Design shear force.
        A_sw (float): Area of shear reinforcement.
        s_w (float): Stirrup spacing.
        z_v (float): Effective shear depth.
        f_ywd (float): Design yield strength of the shear reinforcement.
        cot_theta (float): Cotangent of theta which describes the inclination
        of the compressive stress field
        depending on the LoA.
        a (float): Distance between axis of the support and the concentrated
        load, see Figure 30.1-9.
        x (float): Distance between axis of the support and the investigated
        cross-section load, see Figure 30.1-9.

    Returns:
        float: Increment moment due to shear in case loads are carried
        directly to a support through strut or arch
        action.
    """
    return (V_Ed - A_sw / s_w * z_v * f_ywd * cot_theta) * (a / 2 - x)


def cot_theta_direct_LoA1(
    cot_beta: float,
    cot_theta_min,
) -> float:
    """Calculate the optimized cot_theta in case of loads carried directly to
    a support through strut or arch action in
    LoA I approach.

    fib Model Code 2020: Eq. (30.1-50).

    Args:
        cot_beta (float): Cotangent of beta as shown in Fig. 30.1-9
        cot_theta_min (float): Cotangent of theta_min, where theta_min is the
        minimum inclination of the compressive
        stress field which depends on the LoA.

    Returns:
        float: Optimized cot_theta in case of loads carried directly to a
        support through strut or arch action in
    LoA I approach.
    """
    return min(cot_theta_min, (cot_beta + (1 + (cot_beta**2) ** 0.5)))


def cot_theta_direct_LoA2(
    a: float,
    z_v: float,
    cot_theta_y_c: float,
) -> float:
    """Calculate the optimized cot_theta in case of loads carried directly to
    a support through strut or arch action in
    LoA II approach.

    fib Model Code 2020: Eq. (30.1-51).

    Args:
        a (float): Distance between axis of the support and the concentrated
        load, see Figure 30.1-9.
        z_v (float): Effective shear depth.
        cot_theta_y_c (float): Cotangent of theta_y_c, which is the angle
        where simulatenuous shear reinforcement yield
        and concrete crushing is expected and can be calculated with Eq. (
        30.1-38).

    Returns:
        float: Optimized cot_theta in case of loads carried directly to a
        support through strut or arch action in
    LoA II approach.
    """
    return min(cot_theta_y_c, 1.3 * a / z_v)


def sig_swd(E_s: float, cot_theta: float, eps_x: float, f_ywd: float) -> float:
    """Calculate the design reinforcing steel stress, which replaces f_ywd in
    Eq. (30.1-48) and Eq. (30.1-49) if
    cot_theta < 1 is assumed.

    fib Model Code 2020: Eq. (30.1-52).

    Args:
        E_s (float): Modulus of elasticity of reinforcing steel.
        cot_theta (float): Cotangent of theta which describes the inclination
        of the compressive stress field
        depending on the LoA.
        eps_x (float): Longitudinal strain at mid-depth of the effective
        shear depth core layer as shown in Fig. 30.1-8.
        f_ywd (float): Design yield strength of the shear reinforcement.

    Returns:
        float: The design reinforcing steel stress, which replaces f_ywd in
        Eq. (30.1-48) and Eq. (30.1-49) if
    cot_theta < 1 is assumed.
    """
    return min(f_ywd, (E_s * (cot_theta**2 * (eps_x + 0.001) - 0.001)))


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                           Section 30.1.3.4                          #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
def V_Rd_1(V_Rdc: float, V_Rd_FRC: float) -> float:
    """Calculate the design shear resistance of a FRC member without
        shear reinforcement.

    fib Model Code 2020: Eq. (30.1-54).

    Args:
        V_Rd_c (float): Contribution of the concrete in N.
        V_Rd_FRC (float): Contribution of the residual tensile strength
            related to fibres in N.

    Returns:
        float: The design shear resistance of a FRC member without
            shear reinforcement in N.
    """
    return V_Rdc + V_Rd_FRC


def V_Rd_FRC(
    k_s: float,
    f_FTU: float,
    gamma_F: float,
    cot_theta_F: float,
    z_v: float,
    b_w: float,
) -> float:
    """Calculate contribution of residual tensile strength related to fibres.

    fib Model Code 2020: Eq. (30.1-56).

    Args:
        k_s (float): Parameter k_s should be taken as 0.64.
        f_FTU (float): Characteristic value of the ultimate residual tensile
        strength for FRC in MPa.
        gamma_F (float): Partial factor for FRC (unitless).
        cot_theta_F (float): Cotangent of inclination SFRC compressive stress
        field.
        z_v (float): Effective shear depth.
        b_w (float): Width of cross-section at centroidal axis.

    Returns:
        float: The design shear resistance contribution of fibers in FRC member
         without shear reinforcement.
    """
    return k_s * f_FTU / gamma_F * cot_theta_F * z_v * b_w


def f_FTU_LOA_1(f_R3k: float) -> float:
    """Calculate the characteristic value of the ultimate residual
        tensile strength for FRC in MPa, assuming Level I
        approximation.

    fib Model Code 2020: Eq. (30.1-57).

    Args:
        f_R3k (float): Characteristic residual strength of fibre
            reinforced concrete significant for ultimate conditions.

    Returns:
        float: The characteristic value of the ultimate residual
            tensile strength for FRC in MPa.
    """
    return f_R3k / 3


def k_v_FRC_LOA_1(f_yd: float, E_s: float, zeta: float, z_v: float) -> float:
    """Calculate design shear resistance for fibre reinforced concrete,
        assuming Level I approximation.

    fib Model Code 2020: Eq. (30.1-58).

    Args:
        f_yd (float): Design value of the yield strength of the
            flexural tensile reinforcement in MPa.
        E_s (float): Young's modulus of reinforcement steel in MPa.
        zeta (float): Reduction coefficient of the concrete component
            due to the increased crack width that occurs in the
            presence of fibres at ULS.
        z_v (float): Effective shear depth in mm.

    Returns:
        float: Design shear resistance for fibre reinforced concrete in
            MPa.
    """
    return 0.4 / ((1 + 750 * f_yd / E_s) + zeta) * 1300 / (1000 + 1.25 * z_v)


def zeta_LOA_1(f_yd: float, E_s: float, f_FTU: float, f_ck: float) -> float:
    """Calculate reduction coefficient of the concrete component
            due to the increased crack width that occurs in the
            presence of fibres at ULS, assuming Level I approximation.

    fib Model Code 2020: Eq. (30.1-59).

    Args:
        f_yd (float): Design value of the yield strength of the
            flexural tensile reinforcement in MPa.
        E_s (float): Young's modulus of reinforcement steel in MPa.
        f_FTU(float): The characteristic value of the ultimate residual
            tensile strength for FRC in MPa.
        f_ck (float): Characteristic compressive strength of concrete
            in MPa.

    Returns:
        float: Reduction coefficient of the concrete component
            due to the increased crack width that occurs in the
            presence of fibres at ULS (unitless).
    """
    return max((20 - 7000 * f_yd / E_s) * f_FTU / f_ck, 0)


def theta_F_LOA_1(f_yd: float, E_s: float) -> float:
    """Calculate the inclination of the steel fibre reinforced
        concrete (SFRC) compressive stress field in degrees, assuming
        Level I approximation.

    fib Model Code 2020: Eq. (30.1-60).

    Args:
        f_yd (float): Design value of the yield strength of the
            flexural tensile reinforcement in MPa.
        E_s (float): Young's modulus of reinforcement steel in MPa.

    Returns:
        float: The inclination of the steel fibre reinforced
        concrete (SFRC) compressive stress field in degrees.
    """
    return 29 + 3500 * f_yd / E_s


def tau_RD_cF_LOA_1(
    eta: float,
    gamma_c: float,
    rho_l: float,
    f_ck: float,
    d_dg: float,
    d: float,
    f_Ftud: float,
    tau_Rdc_min: float,
) -> float:
    """Calculate the design shear strength of FRC not requiring shear
        reinforcement, and having ordinary longitudinal reinforcement,
        assuming Level I approximation.

    fib Model Code 2020: Eq. (30.1-61).

    Args:
        eta (float): Factor to account for strain-rate effects in shear
            strength.
        gamma_c (float): Partial factor for concrete (unitless).
        rho_l (float): Longitudinal reinforcement ratio (unitless).
        f_ck (float): Characteristic compressive strength of concrete
            in MPa
        d_dg (float): Parameter related to aggregate size d_g in
            mm.
        d (float): Effective depth in mm.
        f_Ftud (float): Design value of post-cracking strength for
            ultimate crack opening for fibre-reinforced concrete in
            MPa.
        tau_Rdc_min (float): Minimum design value of the shear
            resistance of the concrete in MPa.

    Returns:
        float: The design shear strength of FRC not requiring shear
        reinforcement, and having ordinary longitudinal reinforcement
        in MPa.
    """
    return max(
        eta * 0.6 / gamma_c * (100 * rho_l * f_ck * d_dg / d) ** (1 / 3)
        + f_Ftud,
        eta * tau_Rdc_min + f_Ftud,
    )


def eta_LOA_1(f_Ftuk: float) -> float:
    """Calculate factor to account for strain-rate effects in shear
        strength, assuming Level I approximation.

    fib Model Code 2020: Eq. (30.1-62).

    Args:
        f_Ftuk (float): Characteristic value of post-cracking strength
            for ultimate crack opening for fibre-reinforced concrete in
            MPa.

    Returns:
        float: Factor to account for strain-rate effects in shear
        strength.
    """
    eta = 1.2 - 0.5 * f_Ftuk
    if eta < 0.4:
        return 0.4
    elif eta > 1:
        return 1
    return eta


def f_FTU_alt(f_R3k: float, f_R1k: float) -> float:
    """Use an alternative approach to calculate The characteristic
        value of the ultimate residual tensile strength for FRC in MPa.
        Method might be applied for LOA 1 and 2.

    fib Model Code 2020: Eq. (30.1-63).

    Args:
        f_R3k (float): Characteristic residual strength of fibre
            reinforced concrete significant for ultimate conditions.
        f_R1k (float): Characteristic residual strength of fibre
            reinforced concrete significant for serviceability
            conditions.

    Returns
        float: The characteristic value of the ultimate residual
            tensile strength for FRC in MPa.
    """
    return max(0.46 * f_R3k - 0.13 * f_R1k, 0)


def k_v_FRC_LOA_2(eps_x: float, zeta: float, k_dg: float, z_v: float) -> float:
    """Calculate the design shear resistance for fibre reinforced
        concrete, assuming Level II approximation.

    fib Model Code 2020: Eq. (30.1-64).

    Args:
        eps_x (float): Longitudinal strain at mid-depth of the
            effective shear depth (unitless).
        zeta (float): Reduction coefficient of the concrete component
            due to the increased crack width that occurs in the
            presence of fibres at ULS.
        k_dg (float): Multiplication factor to account for roughness of
            critical shear crack (unitless).
        z_v (float): Effective shear depth in mm.

    Returns:
        float: Design shear resistance for fibre reinforced concrete in
            MPa.
    """
    return 0.4 / (1 + 1500 * eps_x + zeta) * 1300 / (1000 + k_dg * z_v)


def zeta_LOA_2(eps_x: float, f_FTU: float, f_ck: float) -> float:
    """Calculate reduction coefficient of the concrete component
            due to the increased crack width that occurs in the
            presence of fibres at ULS, assuming Level II approximation.

    fib Model Code 2020: Eq. (30.1-65).

    Args:
        eps_x (float): Longitudinal strain at mid-depth of the
            effective shear depth (unitless).
        f_FTU(float): The characteristic value of the ultimate residual
            tensile strength for FRC in MPa.
        f_ck (float): Characteristic compressive strength of concrete
            in MPa.

    Returns:
        float: Reduction coefficient of the concrete component
            due to the increased crack width that occurs in the
            presence of fibres at ULS.
    """
    return max((20 - 14000 * eps_x) * f_FTU / f_ck, 0)


def theta_F_LOA_2_3(eps_x: float) -> float:
    """Calculate the inclination of the steel fibre reinforced
        concrete (SFRC) compressive stress field in degrees, assuming
        Level II approximation.

    fib Model Code 2020: Eqs. (30.1-66), (30.1-68), (30.1-73) and (30.1-75).

    Args:
        eps_x (float): Longitudinal strain at mid-depth of the
            effective shear depth (unitless).

    Returns:
        float: The inclination of the steel fibre reinforced
        concrete (SFRC) compressive stress field in degrees.
    """
    return 29 + 7000 * eps_x


def k_v_FRC_LOA_3_1(eps_x: float, k_dg: float, z_v: float) -> float:
    """Calculate the design shear resistance for fibre reinforced
        concrete, assuming Level III approximation.

    fib Model Code 2020: Eq. (30.1-67).

    Args:
        eps_x (float): Longitudinal strain at mid-depth of the
            effective shear depth (unitless).
        k_dg (float): Multiplication factor to account for roughness of
            critical shear crack (unitless).
        z_v (float): Effective shear depth in mm.

    Returns:
        float: Design shear resistance for fibre reinforced concrete in
            MPa.
    """
    return 0.4 / (1 + 1500 * eps_x) * 1300 / (1000 + k_dg * z_v)


def w_u_LOA_3_1(
    eps_x: float, k_dg: float, d_v: float, theta_F: float
) -> float:
    """Calculate the crack width required for determination of f_FTU in
        mm, assuming Level III approximation.

    fib Model Code 2020: Eq. (30.1-69).

    Args:
        eps_x (float): Longitudinal strain at mid-depth of the
            effective shear depth (unitless).
        k_dg (float): Multiplication factor to account for roughness of
            critical shear crack (unitless).
        d_v (float): Effective depth in mm.
        theta_F (float): The inclination of the steel fibre reinforced
            concrete (SFRC) compressive stress field in degrees.

    Returns:
        float: The crack width required for determination of f_FTU in
            mm.
    """
    return max(
        (0.2 + 1000 * eps_x)
        * ((1000 + k_dg * d_v) / 1300)
        / math.cos(math.radians(theta_F)),
        0.5,
    )


def tau_Rd_sF(
    A_sv: float,
    s: float,
    z: float,
    b: float,
    d: float,
    f_yd: float,
    f_Ftud: float,
    eta_sw: float = 0.75,
    eta_F: float = 1,
) -> float:
    """Calculate the design shear strength for members with shear
        reinforcement and classical longitudinal reinforcement.

    fib Model Code 2020: Eq. (30.1-70).

    Args:
        A_sv (float): Area of shear reinforcement in mm2.
        s (float): Spacing of bars in mm.
        z (float): Internal lever arm in mm.
        b (float): Width of the element in mm.
        d (float): Effective depth in mm.
        f_yd (float): Design value of the yield strength of the
            flexural tensile reinforcement in MPa.
        f_Ftud (float): Design value of post-cracking strength for
            ultimate crack opening for fibre-reinforced concrete in
            MPa.
        eta_sw (float): Parameter expressing the shear contribution of
            ordinary steel in the design shear strength of FRC. Default
            value = 0.75.
        eta_F (float): Parameter espressing the shear contribution of
            steel fibres in the design shear strength. Default value =
            1.0.

    Returns:
        float: The design shear strength for members with shear
            reinforcement and classical longitudinal reinforcement.
    """
    return eta_sw * A_sv / s * z / (b * d) * f_yd + eta_F * f_Ftud


def V_Rd_2(
    V_Rdc: float, V_Rd_FRC: float, V_Rds: float, V_RD_max: float
) -> float:
    """Calculate the design shear resistance for FRC members with shear
        reinforcement in N.

    fib Model Code 2020: Eq. (30.1-71).

    Args:
        V_Rd_c (float): Design shear resistance attributed to the
            concrete as defined by Eq. (30.1-55), depends on LoA.
        V_Rd_FRC (float): Design shear resistance provided by (steel)
            fibres as defined in Eq. (30.1-56).
        V_Rds (float): Design shear resistance provided by shear
            reinforcement as defined in Eq. (30.1-24).
        V_Rd_max (float): Maximum shear resistance related to crushing
            of concrete carrying the compression field, as defined in
            Eq. (30.1-25).

    Returns:
        float: Design shear resistance for FRC members with shear
            reinforcement.
    """
    return min(V_Rdc + V_Rd_FRC + V_Rds, V_RD_max)


def f_FTU_LOA_2(f_R3k: float, f_R1k: float) -> float:
    """Calculate the characteristic value of the ultimate residual
        tensile strength for FRC in MPa, assuming to Level II
        Approximation.

    fib Model Code 2020: Eq. (30.1-72).

    Args:
        f_R3k (float): Characteristic residual strength of fibre
            reinforced concrete significant for ultimate conditions.
        f_R1k (float): Characteristic residual strength of fibre
            reinforced concrete significant for serviceability
            conditions.

    Returns:
        float: The characteristic value of the ultimate residual
            tensile strength for FRC in MPa, assuming to Level II
            Approximation.
    """
    return max(0.46 * f_R3k - 0.13 * f_R1k, 0)


def k_v_FRC_LOA_3_2(eps_x: float) -> float:
    """Calculate the design shear resistance for fibre reinforced
        concrete, assuming Level III approximation.

    fib Model Code 2020: Eq. (30.1-74).

    Args:
        eps_x (float): Longitudinal strain at mid-depth of the
            effective shear depth (unitless).

    Returns:
        float: The design shear resistance for fibre reinforced
        concrete, assuming Level III approximation.
    """
    return 0.4 / (1 + 1500 * eps_x)


def w_u_LOA_3_2(eps_x: float, theta_F: float) -> float:
    """Calculate the crack width required for determination of f_FTU in
        mm, assuming Level III approximation.

    fib Model Code 2020: Eq. (30.1-76).

    Args:
        eps_x (float): Longitudinal strain at mid-depth of the
            effective shear depth (unitless).
        theta_F (float): The inclination of the steel fibre reinforced
            concrete (SFRC) compressive stress field in degrees.

    Returns:
        float: The crack width required for determination of f_FTU in
        mm, assuming Level III approximation.
    """
    return max((0.2 + 1000 * eps_x) / math.cos(math.radians(theta_F)), 0.5)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                           Section 30.1.3.5                          #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
def tau_Rd(
    tau_Rdc: float, rho_w: float, f_fw_Rd: float, nu: float, fcd: float
) -> float:
    """Calculate the design shear stress resistance of concrete with FRP
    reinforcement.

    fib Model Code 2020: Eq. (30.1-77).

    Args:
        tau_Rdc (float): Design value of shear stress resistance of concrete.
        rho_w (float): Reinforcement ratio provided by the FRP stirrups.
        f_fw_Rd (float): Design value of the effective capacity in FRP stirrups.
        nu (float): Strength reduction factor.
        fcd (float): Design value of concrete strength.

    Returns:
        float: Design shear stress resistance of concrete with FRP
        reinforcement, tau_Rd.
    """
    tau_rd = tau_Rdc + rho_w * f_fw_Rd
    tau_rd_max = 0.5 * nu * fcd

    return min(tau_rd, tau_rd_max)


def f_fw_Rd(f_fd: float, E_f: float) -> float:
    """Calculate the design value of the effective capacity in FRP stirrups.

    fib Model Code 2020: Eq. (30.1-78).

    Args:
        f_fd (float): Design tensile capacity of the FRP.
        E_f (float): Modulus of elasticity of the FRP reinforcement.

    Returns:
        float: Design value of the effective capacity in FRP stirrups, f_fw_Rd.
    """
    return min(f_fd, 0.004 * E_f)


def eta_bend(r_f: float, d_f: float) -> float:
    """Calculate the reduction factor that accounts for the effective tensile
    capacity by the bended part of the FRP
    stirrup.

    fib Model Code 2020: Eq. (30.1-79).

    Args:
        r_f (float): radius of the bended part of FRP.
        d_f (float): diameter of the FRP bar.

    Returns:
        float: Reduction factor that accounts for the effective tensile
        capacity by the bended part of the FRP stirrup,
        eta_bend.
    """
    return min(0.03 + (0.05 * r_f / d_f), 1.0)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                           Section 30.1.3.6                          #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
def alpha_l(l_x: float, l_bpt95: float) -> float:
    """Calculate the ratio between distance point of failure and transmission
    length.

    fib Model Code 2020: part of Eq. (30.1-80).

    Args:
        l_x (float): Distance point of failure, see Fig. 30.1-10.
        l_bpt95 (float): Transmission length.

    Returns:
        float: Ratio between distance point of failure and transmission
        length, alpha_l.
    """
    return l_x / l_bpt95


def V_Rd_ct_loa1(
    I_c: float,
    b_w: float,
    S_c: float,
    f_ctd: float,
    alpha_l: float,
    sigma_cp: float,
) -> float:
    """Calculate the design shear resistance for a prestressed hollow core
    slab using the level I approximation.

    fib Model Code 2020: Eq. (30.1-80).

    Args:
        I_c (float): Second moment of area.
        b_w (float): Width of cross-section at centroidal axis.
        S_c (float): First moment of area above and about the centroidal axis.
        f_ctd (float): Design value of concrete tensile strength.
        alpha_l (float): Ratio distance point of failure and transmission
        length.
        sigma_cp (float): Concrete compressive stress at centroidal axis due
        to prestressing in the area where the
        prestressing force is fully introduced.

    Returns:
        float: Design shear resistance for a prestressed hollow core slab,
        V_Rd_ct_loa1.
    """
    V_Rd_ct = (
        0.8
        * I_c
        * b_w
        / S_c
        * (f_ctd**2.0 + (alpha_l * sigma_cp * f_ctd)) ** 0.5
    )

    return V_Rd_ct


def V_Rd_ct_loa2(
    I_c: float,
    b_w_y: float,
    S_c_y: float,
    f_ctd: float,
    alpha_l: float,
    sigma_cp_y: float,
    tau_cp_y: float,
) -> float:
    """Calculate the design shear resistance for a prestressed hollow core
    slab using the level II approximation.

    fib Model Code 2020: Eq. (30.1-81).

    Args:
        I_c (float): Second moment of area.
        b_w_y (float): Width of cross-section at height y, where y is the
        critical point at the line of failure.
        S_c_y (float): First moment of area above height y and about the
        centroidal axis.
        f_ctd (float): Design value of concrete tensile strength.
        alpha_l (float): Ratio distance point of failure and transmission
        length.
        sigma_cp_y (float): Concrete compressive stress at height y and
        distance l_x
        tau_cp_y (float): Shear stress in concrete due to transmission of
        prestress at height y and distance l_x.

    Returns:
        float: Design shear resistance for a prestressed hollow core slab,
        V_Rd_ct_loa2.
    """
    V_Rd_ct = (
        0.8
        * I_c
        * b_w_y
        / S_c_y
        * ((f_ctd**2.0 + (alpha_l * sigma_cp_y * f_ctd)) ** 0.5 - tau_cp_y)
    )

    return V_Rd_ct


def sigma_cp_y(
    A_c: float,
    y_c: float,
    y: float,
    y_pt: float,
    I: float,
    F_p_lx: float,
    M_Ed: float,
) -> float:
    """Calculate the concrete compressive stress at height y and distance l_x.

    fib Model Code 2020: Eq. (30.1-82) --> ADJUSTED! Reason: the equation
    provided is dimensionally inconsistent, and
    probably incorrect. Instead, we took the eq. from EN 1168:2005 in
    4.3.3.2.2.2.

    Args:
        A_c (float): Cross-sectional area of concrete.
        y_c (float): Height of concrete centroidal axis.
        y (float): Height.
        y_pt (float): Height of centroidal axis of prestressing steel.
        I (float): Second moment of area.
        F_p_lx (float): Prestressing force at distance l_x.
        M_Ed (float): Design value of the applied bending moment.

    Returns:
        float: Concrete compressive stress at height y and distance l_x,
        sigma_cp_y.
    """
    # sigma_cp_y = ((1.0 / A_c) + ((y_c - y) / I)) * F_p_lx  # Eq. (30.1-82)
    # MC2020
    sigma_cp_y = (((1.0 / A_c) + ((y_c - y) * (y_c - y_pt) / I)) * F_p_lx) - (
        M_Ed * (y_c - y) / I
    )

    return sigma_cp_y


def tau_cp_y(
    b_w_y: float,
    A_c_y: float,
    A_c: float,
    S_c_y: float,
    y_c: float,
    y_pt: float,
    I: float,
    dFplx_dx: float,
) -> float:
    """Calculate the shear stress in concrete due to transmission of
    prestress at height y and distance l_x.

    fib Model Code 2020: Eq. (30.1-83).

    Args:
        b_w_y (float): Width of cross-section at height y, where y is the
        critical point at the line of failure.
        A_c_y (float): Concrete area above height y.
        A_c (float): Cross-sectional area of concrete.
        S_c_y (float): First moment of area above height y and about the
        centroidal axis.
        y_c (float): Height of concrete centroidal axis.
        y_pt (float): Height of centroidal axis of prestressing steel.
        I (float): Second moment of area.
        dFplx_dx (float): Gradient of prestressing force at distance l_x.

    Returns:
        float: Shear stress in concrete due to transmission of prestress at
        height y and distance l_x., tau_cp_y.
    """
    tau_cp_y = (
        (1.0 / b_w_y) * ((A_c_y / A_c) - (S_c_y * (y_c - y_pt) / I)) * dFplx_dx
    )

    return tau_cp_y


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                           Section 30.1.3.8                          #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
def unity_check_tau(tau_Rdi: float, tau_Edi: float) -> float:
    """Perform a unity check on shear stress at the interface between
    concrete cast at different times.

    fib Model Code 2020: Eqs. (30.1-84).

    Args:
        tau_Rdi (float): Design shear stress resistance at the interface.
        tau_Edi (float): Design shear stress at the interface.

    Returns:
        float: The unity check (unitless).
    """
    return tau_Edi / tau_Rdi


def tau_Edi(beta: float, V_Ed: float, z: float, bi: float) -> float:
    """Calculate the shear stress at the interface between concrete cast at
    different times.

    fib Model Code 2020: Eqs. (30.1-85).

    Args:
        beta (float): Ratio of the longitudinal force in the new concrete to
        the total longitudinal force.
        V_Ed (float): Shear force on the composite section.
        z (float): Inner lever arm of the composite section.
        bi (float): Width of the interface.

    Returns:
        float: Design value of shear stress at the interface.
    """
    return beta * V_Ed / (z * bi)


def tau_Rdi(
    c_a: float, f_ctd: float, mu: float, sigma_n: float, f_cd: float
) -> float:
    """Calculate the design shear stress resistance at the interface for
    scenarios without reinforcement.

    fib Model Code 2020: Eqs. (30.1-86).

    Args:
        c_a (float): Coefficient for the adhesive bond from Table 30.1-4.
        f_ctd (float): Design value of axial tensile strength of concrete.
        mu (float): Friction coefficient from Table 30.1-5.
        sigma_n (float): Lowest expected compressive stress resulting from an
        eventual normal force acting on the
        interface.
        f_cd (float): Design value of cylinder compressive strength of concrete.

    Returns:
        float: Design value of shear stress resistance at the interface.
    """
    tau_Rdi = c_a * f_ctd + mu * sigma_n
    tau_Rdi_lim = 0.25 * f_cd

    return min(tau_Rdi, tau_Rdi_lim)


def tau_Rdi_dowel(
    c_r: float,
    f_cd: float,
    mu: float,
    sigma_n: float,
    kappa_1: float,
    rho: float,
    f_yd: float,
    alpha_deg: float,
    kappa_2: float,
) -> float:
    """Calculate the design shear stress resistance at the interface for
    scenarios with dowels or reinforcement.

    fib Model Code 2020: Eqs. (30.1-87).

    Args:
        c_r (float): Coefficient for aggregate interlock effects from Table
        30.1-5.
        f_cd (float): Design value of cylinder compressive strength of concrete.
        mu (float): Friction coefficient from Table 30.1-5.
        sigma_n (float): Lowest expected compressive stress resulting from an
        eventual normal force acting on the
        interface.
        kappa_1 (float): interaction coefficient for tensile force activated
        in the dowels or reinforcement.
        rho: reinforcement ratio of reinforcing steel crossing the interface.
        f_yd (float): Design yield strength of the reinforcement.
        alpha_deg (float): Inclination of the reinforcement crossing the
        interface in degrees, according to fig.
        30.1-12.
        kappa_2 (float): interaction coefficient for flexural resistance.

    Returns:
        float: Design value of shear stress resistance at the interface.
    """
    alpha_rad = math.radians(alpha_deg)

    tau_Rdi = (
        c_r * f_cd ** (1 / 3)
        + mu * sigma_n
        + kappa_1
        * rho
        * f_yd
        * (mu * math.sin(alpha_rad) + math.cos(alpha_rad))
        + kappa_2 * rho * math.sqrt(f_yd * f_cd)
    )
    tau_Rdi_lim = 0.25 * f_cd

    return min(tau_Rdi, tau_Rdi_lim)


def R_f(V: float, D: float) -> float:
    """Calculate the interface surface roughness.

    fib Model Code 2020: Eqs. (30.1-88).

    Args:
        V (float): Volume.
        D (float): Diameter.

    Returns:
        float: Surface roughness in mm.
    """
    return 40 * V / (math.pi * D**2)
