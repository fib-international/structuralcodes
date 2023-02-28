"""A collection of shear formulas for concrete"""
import warnings
from typing import Optional
from math import pi, tan, sin


def epsilon_x(
    E_s: float,
    As: float,
    Med: float,
    Ved: float,
    Ned: float,
    z: float,
    delta_e: float,
) -> float:
    """The maximum allowed shear resistance

    fib Model Code 2010, eq. (7.3-16)

    Args:
        E (float): The E-modulus to the material in MPa
        As (Float): The cross-section area of reinforcement in mm^2
        Med (Float): The positive moment working on the material in Nmm
        Ved (float): The positive shear force working on the material in N
        Ned: (float): The normal force working on the material in N with
        positive sign for tension and negative sign for compression
        z: (float): The effective shear depth in mm
        delta_E (float): The eccentricity of the axial load due to
        imperfection in the construction with distance in mm as a positive
        value
    Returns:
        float: The longitudinal strain"""
    if Ned >= 0:
        return max(
            (
                (1 / (2 * E_s * As))
                * ((abs(Med) / z) + abs(Ved) + Ned * ((1 / 2) + (delta_e / z)))
            ),
            0,
        )

    return max(
        (
            (1 / (2 * E_s * As))
            * ((abs(Med) / z) + abs(Ved) + Ned * ((1 / 2) - (delta_e / z)))
        ),
        0,
    )


def eta_fc(fck: float):
    """returns N_fc that is used to determin the strength reduction factor"""
    return min((30 / fck) ** (1 / 3), 1)


def v_rd(
    approx_lvl: int,
    with_shear_reinforcment: bool,
    fck: float,
    z: float,
    bw: float,
    dg: float,
    E_s: float,
    As: float,
    Med: float,
    Ved: float,
    Ned: float,
    delta_e: float,
    asw: Optional[float],
    sw: Optional[float],
    f_ywd: Optional[float],
    theta: Optional[float],
    alfa: float = 90.0,
    gamma_c: float = 1.5,
) -> float:
    """Compute the shear resistance of a web or slab.

    fib Model Code 2010, Eq. (7.3-11)

    Args:
        approx_lvl_c (int): Approximation level for concrete
        approx_lvl_s (int): Approximation level for reinforcment
        reinforcment (bool): Shear reinforced concrete or no shear
        reinforcement
        fck (float): Characteristic strength in MPa
        gamma_c (float): Concrete safety factor
        z: (float): distances between the centerline of the
        compressive chord and the reinforcement in mm
        bw: (float): Thickness of web in cross section in mm
        dg: (float): Maximum size of aggregate in mm
        E_s: (float): The E_s-modulus to the materialb in MPa
        As: (float): The cross-section area of reinforcement in mm^2
        Med (Float): The positive moment working on the material in Nmm
        Ved (float): The positive shear force working on the material in N
        Ned: (float): The normal force working on the material in N with
        positive sign for tension and negative sign for compression
        delta_E (float): The eccentricity of the axial load due to
        imperfection in the construction with distance in mm as a positive
        value
        alfa (float): Inclination of the stirrups in degrees
        asw (float): Area of shear reinforcement in mm^2
        sw (float): Senter distance between the shear reinforcement in mm
        f_ywd (float): The design yield strength of the shear reinforcement
        theta (float): Inclitantion of the compression stressfield in degrees

    Returns:
        float: Design shear resistance
    """

    if not with_shear_reinforcment:
        if not approx_lvl == (1 or 2):
            warnings.warn(
                "Choosen approximation is not suited without reinforcment"
            )
        return v_rdc(
            approx_lvl,
            fck,
            z,
            bw,
            dg,
            E_s,
            As,
            Med,
            Ved,
            Ned,
            delta_e,
            alfa,
            gamma_c,
        )

    if approx_lvl == (3 or 4):
        return min(
            v_rds(asw, sw, z, f_ywd, theta, alfa),
            v_rd_max(
                approx_lvl,
                fck,
                bw,
                theta,
                z,
                E_s,
                As,
                Med,
                Ved,
                Ned,
                delta_e,
                alfa,
                gamma_c,
            ),
        )

    if approx_lvl == 5:
        V_rdc = v_rdc(
            approx_lvl,
            fck,
            z,
            bw,
            dg,
            E_s,
            As,
            Med,
            Ved,
            Ned,
            delta_e,
            alfa,
            gamma_c,
        )
        V_rd_max = v_rd_max(
            approx_lvl,
            fck,
            bw,
            theta,
            z,
            E_s,
            As,
            Med,
            Ved,
            Ned,
            delta_e,
            alfa,
            gamma_c,
        )
        if (V_rdc + v_rds(asw, sw, z, f_ywd, theta, alfa)) < V_rd_max:
            return V_rdc + v_rds(asw, sw, z, f_ywd, theta, alfa)
        v_rd(
            4,
            with_shear_reinforcment,
            fck,
            z,
            bw,
            dg,
            E_s,
            As,
            Med,
            Ved,
            Ned,
            delta_e,
            asw,
            sw,
            f_ywd,
            theta,
            alfa,
            gamma_c,
        )
    raise ValueError("invalid approx level")


def v_rdc(
    approx_lvl,
    fck: float,
    z: float,
    bw: float,
    dg: float,
    E_s: float,
    As: float,
    Med: float,
    Ved: float,
    Ned: float,
    delta_e: float,
    alfa: float = 90.0,
    gamma_c: float = 1.5,
) -> float:
    """The design shear resistance of a web or a slab without
    shear reinforcement.

    fib Model Code 2010, Eq. (7.3-17)

    Args:
        approx_lvl_c (int): Approximation level for concrete
        approx_lvl_s (int): Approximation level for reinforcement
        with_shear_reinforcment (bool): Shear reinforced concrete or no shear
        reinforcement
        fck (float): Characteristic strength in MPa
        gamma_c (float): Concrete safety factor
        z: (float): The length to the areasenter of cross-section in mm
        bw: (float): Thickness of web in cross section in mm
        dg: (float): Maximum size of aggregate
        E_s: (float): The E_s-modulus to the materialb in MPa
        As: (float): The cross-section area in mm^2
        Med (Float): The positive moment working on the material in Nmm
        Ved (float): The positive shear force working on the material in N
        Ned: (float): The normal force working on the material in N with
        positive sign for tension and negative sign for compression
        delta_E (float): The eccentricity of the axial load due to
        imperfection in the construction with distance in mm as a positive
        value
        alfa (float): Inclination of the stirrups in degrees

    Returns:
        float: The design shear resistance attributed to the concrete
    """

    if approx_lvl == 1:
        return v_rdc_approx1(fck, z, bw, gamma_c)

    if approx_lvl == 2:
        return v_rdc_approx2(
            fck, z, bw, dg, E_s, As, Med, Ved, Ned, delta_e, gamma_c
        )
    if approx_lvl == 5:
        return v_rdc_approx3(
            approx_lvl,
            fck,
            z,
            bw,
            E_s,
            As,
            Med,
            Ved,
            Ned,
            delta_e,
            alfa,
            gamma_c,
        )
    raise ValueError("Invalid approx level")


def v_rdc_approx1(
    fck: float,
    z: float,
    bw: float,
    gamma_c: float = 1.5,
) -> float:
    """Gives the shear resistance for concrete with approx level 1

    For members with no segnificant axal load, with fyk <= 600 Mpa,
    fck <= 70 Mpa and with maximum aggrigate size of not less then 10mm.

    fib Model Code 2010, Eq. (7.3-17) and (7.3-19)

    Args:
        fck (float): The characteristic compressive strength in MPa.
        z (float): The length to the areasenter of cross-section in mm
        gamma_c: Safety factor for concrete
        bw: Thickness of web in cross section

    Returns:
        float: Design shear resistance without shear reinforcement
    """
    fsqr = min(fck**0.5, 8)
    kv = 180 / (1000 + 1.25 * z)
    return (kv * fsqr * z * bw) / gamma_c


def v_rdc_approx2(
    fck: float,
    z: float,
    bw: float,
    dg: float,
    E_s: float,
    As: float,
    Med: float,
    Ved: float,
    Ned: float,
    delta_e: float,
    gamma_c: float = 1.5,
) -> float:
    """Gives the shear resistance for concrete with approx level 2

    In higher strength concrete and light-weight aggregate concretes,
    the fracture surface may go through the aggregate particles,
    rather then around, reducing the crack roughness

    fib Model Code 2010, Eq. (7.3-17), (7.3-20) and (7.3-21)

    Args:
        fck (float): Characteristic strength in MPa
        gamma_c (float): Concrete safety factor
        z (float): The length to the areasenter of cross-section in mm
        bw (float): Thickness of web in cross section
        dg (float): Maximum size of aggregate
        E_s (float): The E_s-modulus to the materialb in MPa
        As (float): The cross-section area in mm^2
        Med (Float): The positive moment working on the material in Nmm
        Ved (float): The positive shear force working on the material in N
        Ned (float): The normal force working on the material in N
        delta_E (float): The eccentricity of the axial load due to
        imperfection in the construction with distance in mm as a positive
        value

    Returns:
        float: Design shear resistance without shear reinforcement
    """
    fsqr = min(fck**0.5, 8)
    epsilonx = epsilon_x(E_s, As, Med, Ved, Ned, z, delta_e)
    k_dg = max(32 / (16 + dg), 0.75)
    kv = (0.4 / (1 + 1500 * epsilonx)) * (1300 / (1000 + k_dg * z))
    return (kv * fsqr * z * bw) / gamma_c


def v_rdc_approx3(
    approx_lvl: float,
    fck: float,
    z: float,
    bw: float,
    E_s: float,
    As: float,
    Med: float,
    Ved: float,
    Ned: float,
    delta_e: float,
    alfa: float = 90.0,
    gamma_c: float = 1.5,
) -> float:
    """Gives the shear resistance for concrete with approx level 3

    The design shear resistance of a web or a slab without
    shear reinforcement.

    fib Model Code 2010, Eq. (7.3-17), (7.3-39) and (7.3-43)

    Args:
        fck (float): Characteristic strength in MPa
        gamma_c (float): Concrete safety factor
        z: (float): The length to the areasenter of cross-section in mm
        bw: (float): Thickness of web in cross section
        dg: (float): Maximum size of aggregate
        E_s: (float): The E_s-modulus to the materialb in MPa
        As: (float): The cross-section area in mm^2
        Med (Float): The positive moment working on the material in Nmm
        Ved (float): The positive shear force working on the material in N
        Ned: (float): The normal force working on the material in N with
        positive sign for tension and negative sign for compression
        delta_E (float): The eccentricity of the axial load due to
        imperfection in the construction with distance in mm as a positive
        value
        alfa (float): Inclination of the stirrups in degrees

    Returns:
        float: Design shear resistance without shear reinforcement
    """
    fsqr = min(fck**0.5, 8)
    theta_min = 20 + 10000 * epsilon_x(E_s, As, Med, Ved, Ned, z, delta_e)
    V_rd_max = v_rd_max(
        approx_lvl,
        fck,
        bw,
        theta_min,
        z,
        E_s,
        As,
        Med,
        Ved,
        Ned,
        delta_e,
        alfa,
        gamma_c,
    )
    epsilonx = epsilon_x(E_s, As, Med, Ved, Ned, z, delta_e)
    kv = max(
        (0.4 / (1 + 1500 * epsilonx)) * (1 - Ved / V_rd_max),
        0,
    )

    return (kv * fsqr * z * bw) / gamma_c


def v_rds(
    asw: float,
    sw: float,
    z: float,
    f_ywd: float,
    theta: float,
    alpha: float,
) -> float:
    """The shear resistans that shear reinforcement gives

    fib Model Code 2010, Eq. (7.3-25) and (7.3-29)
    Args:
        asw (float): Area of shear reinforcement in mm
        sw (float): Senter distance between the shear reinforcement in mm
        z: (float): The length to the areasenter of cross-section in mm
        f_ywd (float): Design yield strength of the shear reinforcement in MPa
        theta (float): Inclitaniton of the compression stressfield in degrees
        alfa (float): Inclination of the stirrups

    Returns:
        The design shear resistance provided by shear reinforcement
    """
    if 45 < theta or theta < 20:
        warnings.warn("Too high or too low compression field angel")

    result = (
        (asw / sw)
        * z
        * f_ywd
        * ((1 / tan(theta * pi / 180)) + (1 / tan(alpha * pi / 180)))
        * sin(alpha * pi / 180)
    )
    return result


def v_rd_max(
    approx_lvl: int,
    fck: float,
    bw: float,
    theta: float,
    z: float,
    E_s: float,
    As: float,
    Med: float,
    Ved: float,
    Ned: float,
    delta_e: float,
    alfa: float = 90,
    gamma_c: float = 1.5,
) -> float:
    """The maximum allowed shear resistance, when there is shear reinforcment

    fib Model Code 2010, eq. (7.3-24) and (7.3-26)

    Args:
        approx_lvl_s (int): Approximation level for steel
        fck (float): Characteristic strength in MPa
        gamma_c (float): Concrete safety factor
        z: (float): The length to the areasenter of cross-section in mm
        bw: (float): Thickness of web in cross section
        dg: (float): Maximum size of aggregate
        E_s: (float): The E_s-modulus to the materialb in MPa
        As: (float): The cross-section area in mm^2
        Med (Float): The positive moment working on the material in Nmm
        Ved (float): The positive shear force working on the material in N
        Ned: (float): The normal force working on the material in N with
        positive sign for tension and negative sign for compression
        delta_E (float): The eccentricity of the axial load due to
        imperfection in the construction with distance in mm as a positive
        value
        alfa (float): Inclination of the stirrups in degrees
    Returns:
        float: The maximum allowed shear resistance regardless of
        approximation level"""
    if approx_lvl == 3:
        return v_rd_max_approx1(fck, bw, theta, z, alfa, gamma_c)

    if approx_lvl == 4:
        return v_rd_max_approx2(
            fck, bw, theta, z, E_s, As, Med, Ved, Ned, delta_e, alfa, gamma_c
        )

    if approx_lvl == 5:
        return v_rd_max_approx3(
            fck, bw, z, E_s, As, Med, Ved, Ned, delta_e, alfa, gamma_c
        )
    raise ValueError("invalid approx level")


def v_rd_max_approx1(
    fck: float,
    bw: float,
    theta: float,
    z: float,
    alfa: float = 90.0,
    gamma_c: float = 1.5,
) -> float:
    """The maximum allowed shear resistance, with level 1 approximation.

    fib Model Code 2010, eq. (7.3-37)

    Args:
        fck (float): Characteristic strength in MPa
        gamma_c (float): Concrete safety factor
        z: (float): The length to the areasenter of cross-section in mm
        bw: (float): Thickness of web in cross section
        delta_E (float): The eccentricity of the axial load due to
        imperfection in the construction with distance in mm as a positive
        value
        alfa (float): Inclination of the stirrups in degrees

    Returns:
        float: The maximum allowed shear resistance regardless of
        approximation level"""

    return (
        0.55
        * eta_fc(fck)
        * (fck / gamma_c)
        * bw
        * z
        * (
            ((1 / tan(theta * pi / 180)) + (1 / tan(alfa * pi / 180)))
            / (1 + (1 / tan(theta * pi / 180)) ** 2)
        )
    )


def v_rd_max_approx2(
    fck: float,
    bw: float,
    theta: float,
    z: float,
    E_s: float,
    As: float,
    Med: float,
    Ved: float,
    Ned: float,
    delta_e: float,
    alfa: float = 90.0,
    gamma_c: float = 1.5,
) -> float:
    """The maximum allowed shear resistance, with level 2 approximation

    fib Model Code 2010, eq. (7.3-24), (7.3-40) and (7.3-41)

    Args:
        fck (float): Characteristic strength in MPa
        gamma_c (float): Concrete safety factor
        z: (float): The length to the areasenter of cross-section in mm
        bw: (float): Thickness of web in cross section
        dg: (float): Maximum size of aggregate
        E_s: (float): The E_s-modulus to the materialb in MPa
        As: (float): The cross-section area of reinforcement in mm^2
        Med (Float): The positive moment working on the material in Nmm
        Ved (float): The positive shear force working on the material in N
        Ned: (float): The normal force working on the material in N with
        positive sign for tension and negative sign for compression
        delta_E (float): The eccentricity of the axial load due to
        imperfection in the construction with distance in mm as a positive
        value
        alfa (float): Inclination of the stirrups in degrees

    Returns:
        float: The maximum allowed shear resistance regardless of
        approximation level"""
    epsilonx = epsilon_x(E_s, As, Med, Ved, Ned, z, delta_e)
    epsilon_1 = epsilonx + (epsilonx + 0.002) * (
        (1 / tan(theta * pi / 180)) ** 2
    )
    k_epsilon = min(1 / (1.2 + 55 * epsilon_1), 0.65)

    return (
        k_epsilon
        * eta_fc(fck)
        * (fck / gamma_c)
        * bw
        * z
        * (
            ((1 / tan(theta * pi / 180)) + (1 / tan(alfa * pi / 180)))
            / (1 + (1 / tan(theta * pi / 180)) ** 2)
        )
    )


def v_rd_max_approx3(
    fck: float,
    bw: float,
    z: float,
    E_s: float,
    As: float,
    Med: float,
    Ved: float,
    Ned: float,
    delta_e: float,
    alfa: float = 90.0,
    gamma_c: float = 1.5,
) -> float:
    """The maximum allowed shear resistance, with level 3 approximation

    fib Model Code 2010, eq. (7.3-24), (7.3-40) and (7.3-41)

    Args:
        fck (float): Characteristic strength in MPa
        gamma_c (float): Concrete safety factor
        z: (float): The length to the areasenter of cross-section in mm
        bw: (float): Thickness of web in cross section
        dg: (float): Maximum size of aggregate
        E_s: (float): The E_s-modulus to the materialb in MPa
        As: (float): The cross-section area of reinforcement in mm^2
        Med (Float): The positive moment working on the material in Nmm
        Ved (float): The positive shear force working on the material in N
        Ned: (float): The normal force working on the material in N with
        positive sign for tension and negative sign for compression
        delta_E (float): The eccentricity of the axial load due to
        imperfection in the construction with distance in mm as a positive
        value
        alfa (float): Inclination of the stirrups in degrees

    Returns:
        float: The maximum allowed shear resistance regardless of
        approximation level"""
    epsilonx = epsilon_x(E_s, As, Med, Ved, Ned, z, delta_e)
    theta_min = 20 + 10000 * epsilonx

    return v_rd_max_approx2(
        fck, bw, theta_min, z, E_s, As, Med, Ved, Ned, delta_e, alfa, gamma_c
    )


def v_rd_ct(
    approx_lvl_h: int,
    f_ctd: float,
    i_c: float,
    s_c: float,
    b_w: float,
    sigma_cp: float,
    l_x: float,
    l_bd0: float,
    S_cy: float,
    b_wy: float,
    y: float,
    y_c: float,
    A_c: float,
    A_cy: float,
    y_pt: float,
    f_p_lx: float,
    f_p_lx_dx: float,
) -> float:
    """The shear resistance for a hollow core slab.

    fib Model Code 2010, eq. (7.3-44) and (7.3-45)

    Args:
        approx_lvl_h (int): What approximation level we want for hollow core
        f_ctd (float): The design value of concrete axial tensile strength
        i_c (float): Second moment of area in mm^4
        s_c (float): First moment of area, abouve and about the
        centriodal axis in mm^3
        b_w (float): The width of the cross-section at the centroidal axis
        sigma_cp (float): The compressive stress at centroidal axis
        due to prestress
        l_x (float): The distance between edge and point of
        failure (Figure: 7.3-12) in mm
        l_bd0 (float): follows 7.13-5 in mm
        S_cy (float): The first moment of area above y in mm^3
        b_wy (float): The width at hight y in mm
        y (float): Hight at of the critical point at the line of failure in mm
        y_c (float): The hight of the concrete centroidal axis in mm
        A_c (float): The area of concrete cross-section in mm^2
        A_cy (float): The area of concrete cross-section above y in mm^2
        y_pt (float): The hight of centroidal axis of prestressed steel in mm
        f_p_lx (float): The prestressing force at the distance l_x in N
        f_p_lx_dx (float): The change of prestressing force at the distance l_x

    Return:
        The maximum allowed shear force in a hollow core. Regardless of the
        approximation level"""
    if approx_lvl_h == 1:
        return v_rd_ct_approx1(f_ctd, i_c, s_c, b_w, sigma_cp, l_x, l_bd0)
    if approx_lvl_h == 2:
        return v_rd_ct_approx2(
            f_ctd,
            i_c,
            l_x,
            l_bd0,
            S_cy,
            b_wy,
            y,
            y_c,
            A_c,
            A_cy,
            y_pt,
            f_p_lx,
            f_p_lx_dx,
        )
    raise ValueError("Invalid approx level")


def v_rd_ct_approx1(
    f_ctd: float,
    i_c: float,
    s_c: float,
    b_w: float,
    sigma_cp: float,
    l_x: float,
    l_bd0: float,
) -> float:
    """Calculating level 1 approximation for shear
    resistance in hollow core slabs.

    fib Model Code 2010, eq. (7.3-44)

    Args:
        f_ctd (float): The design value of concrete axial tensile strength
        i_c (float): The second moment of area
        s_c (float): The first moment of area, above and about the
        centriodal axis in mm^3
        b_w (float): The width of the cross-section at the centroidal axis
        sigma_cp (float): The compressive stress at centroidal
        axis due to prestress in Mpa
        l_x (float): Distance from the edge to point of
        failure (Figure: 7.3-12)
        l_bd0 (float): follows 7.13-5
    return:
        Vrd Appoximation 1 for hollow slabs"""
    alfa_l = l_x / (1.2 * l_bd0)
    return (
        0.8
        * ((i_c * b_w) / s_c)
        * ((f_ctd**2) + alfa_l * sigma_cp * f_ctd) ** 0.5
    )


def v_rd_ct_approx2(
    f_ctd: float,
    i_c: float,
    l_x: float,
    l_bd0: float,
    S_cy: float,
    b_wy: float,
    y: float,
    y_c: float,
    A_c: float,
    A_cy: float,
    y_pt: float,
    f_p_lx: float,
    f_p_lx_dx: float,
) -> float:
    """Calculates the maximum shear force for level 2 approximation
       in hollow core slabs.

    fib Model Code 2010, eq. (7.3-45), (7.3-46) and (7.3-47)

       Args:
           f_ctd (float): Design value of concrete axial tensile strength in
           MPa
           i_c (float): The second moment of area in mm^4
           l_x (float): Distance from the edge to point of
           failure (Figure: 7.3-12)
           l_bd0 (float): follows 7.13-5
           S_cy (float): The first moment of area above y in mm^3
           b_wy (float): The width at hight y in mm
           y (float): The highh of the critical point at the line of failure
           in mm
           y_c (float): The hight of the concrete centroidal axis in mm
           A_c (float): The area of concrete cross-section in mm^2
           A_cy (float): The area of concrete cross-section above y in mm^2
           y_pt (float): The hight of centroidal axis of prestressed steel in
           mm
           f_p_lx (float): The prestressing force at the distance l_x in N
           f_p_lx (float): The derivative of prestressing force at the
           distance l_x with respect to dx

       sigma_cpy: The compressiv stress in the concrete at hight y and l_x
       tau_cpy: The shear stress due to prestress at hight y and l_x

       return:
           The maximum shear force for level 2 approximation"""

    alfa_l = l_x / (1.2 * l_bd0)
    sigma_cpy = ((1 / A_c) * ((y_c - y) / i_c)) * f_p_lx
    tau_cpy = (
        (1 / b_wy) * ((A_cy / A_c) - (S_cy * (y_c - y_pt)) / i_c) * f_p_lx_dx
    )
    return (i_c * b_wy / S_cy) * (
        ((f_ctd**2) + alfa_l * sigma_cpy * f_ctd) ** 0.5 - tau_cpy
    )
