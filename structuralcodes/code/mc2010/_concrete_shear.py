"""A collection of shear formulas for concrete"""
import typing as t
import warnings
from math import pi, tan, sin, cos


def epsilon_x(
    E: float,
    As: float,
    Med: float,
    Ved: float,
    Ned: float,
    z: float,
    delta_e: float,
) -> float:
    """The maximum allowed shear resistance

    fib Model Code 2010, eq. (7.3-26) and (7.3-24)

    Args:
        E (float): The E-modulus to the materialb in MPa
        AS (Float): The cross-section area in mm^2
        Med (Float): The moment working on the material in Nmm
        Ved (float): The shear working on the material in N
        Ned: (float): The normal force working on the material in N
        z: (float): The length to the areasenter of cross-section in mm
        delta_E (float): The exentricity of the load in mm
    Returns:
        float: The longitudinal strain"""
    return max(
        (1 / (2 * E * As)) * (
            (abs(Med) / z) + abs(Ved) + abs(Ned) * ((1 / 2) + (delta_e / z))),
        0
    )


def v_rd(
    approx_lvl_c: int,
    approx_lvl_s: int,
    reinforcment: bool,
    fck: float,
    z: float,
    bw: float,
    dg: float,
    E: float,
    As: float,
    Med: float,
    Ved: float,
    Ned: float,
    delta_e: float,
    alfa: float,
    gamma_c: float,
    asw: float,
    sw: float,
    fywd: float,
    theta: float,
) -> float:
    """Compute the shear resistance of a web or slab.

    fib Model Code 2010, Eq. (7.3-11)

    Args:
        approx_lvl_c (int): Approximation level for concrete
        approx_lvl_s (int): Approximation level for steel
        reinforcment (bool): Shear reinforced concrete or no shear
        reinforcement
        fck (float): Characteristic strength in MPa
        z: (float): The length to the areasenter of cross-section in mm
        bw: (float): Thickness of web in cross section
        dg: (float): Maximum size of aggregate
        E: (float): The E-modulus to the materialb in MPa
        As: (float): The cross-section area in mm^2
        Med: (float): The moment working on the material in Nmm
        Ved: (float): The shear working on the material in N
        Ned: (float): The normal force working on the material in N
        delta_e (float): The exentricity of the load in mm
        alfa (float): Inclination of the stirrups
        gamma_c (float): Concrete safety factor
        asw (float): Area of shear reinforcement
        sw (float): Senter distance between the shear reinforcement
        fywd (float): The design yield strength of the shear reinforcement
        theta (float): Inclitaniton of the compression stressfield

    Returns:
        float: Design shear resistance
    """
    if not (approx_lvl_c == 1 or 2) and not approx_lvl_s == 3:
        warnings.warn("Not a valid approximation level")

    if not reinforcment:
        return v_rdc(
            approx_lvl_c,
            fck,
            z,
            bw,
            dg,
            E,
            As,
            Med,
            Ved,
            Ned,
            delta_e,
            alfa,
            gamma_c,
            )
    if reinforcment and approx_lvl_s == 3:
        return min(v_rdc(
            approx_lvl_c,
            fck,
            z,
            bw,
            dg,
            E,
            As,
            Med,
            Ved,
            Ned,
            delta_e,
            alfa,
            gamma_c) +
            v_rds(asw, sw, z, fywd, theta, alfa),
            v_rd_max(
            approx_lvl_s,
            fck,
            bw,
            theta,
            z,
            E,
            As,
            Med,
            Ved,
            Ned,
            delta_e,
            alfa,
            gamma_c,
            )
            )
    elif reinforcment and approx_lvl_s == 1 or 2:
        return min(v_rds(asw, sw, z, fywd, theta, alfa), v_rd_max(
            approx_lvl_s,
            fck,
            bw,
            theta,
            z,
            E,
            As,
            Med,
            Ved,
            Ned,
            delta_e,
            alfa,
            gamma_c,
            )
            )


def v_rdc(
    approx_lvl_c: int,
    approx_lvl_s: int,
    fck: float,
    z: float,
    bw: float,
    dg: float,
    E: float,
    As: float,
    Med: float,
    Ved: float,
    Ned: float,
    delta_e: float,
    alfa: float,
    gamma_c: float = 1.5,
) -> float:

    """The design shear resistance of a web or a slab without
    shear reinforcement.

    fib Model Code 2010, Eq. (7.3-17)

    Args:
        approx_lvl_c (int): Approximation level for concrete
        approx_lvl_s (int): Approximation level for steel
        reinforcment (bool): Shear reinforced concrete or no shear
        reinforcement
        fck (float): Characteristic strength in MPa
        z: (float): The length to the areasenter of cross-section in mm
        bw: (float): Thickness of web in cross section
        dg: (float): Maximum size of aggregate
        E: (float): The E-modulus to the materialb in MPa
        As: (float): The cross-section area in mm^2
        Med: (float): The moment working on the material in Nmm
        Ved: (float): The shear working on the material in N
        Ned: (float): The normal force working on the material in N
        delta_e (float): The exentricity of the load in mm
        alfa (float): Inclination of the stirrups
        gamma_c (float): Safety factor

    Returns:
        float: The design shear resistance attributed to the concrete
    """

    if approx_lvl_c == 1:
        return v_rdc_approx1(fck, z, bw, gamma_c)

    elif approx_lvl_c == 2:
        return v_rdc_approx2(
            fck, z, bw, dg, E, As, Med, Ved, Ned, delta_e, gamma_c
        )

    elif approx_lvl_s == 3:
        return v_rds_approx3(
            approx_lvl_s, fck, z, bw, E, As, Med, Ved,
            Ned, delta_e, alfa, gamma_c
        )


def v_rdc_approx1(
    fck: float,
    z: float,
    bw: float,
    gamma_c: float = 1.5,
) -> float:
    """The design shear resistance of a web or a slab without
    shear reinforcement.

    fib Model Code 2010, Eq. (7.3-17)

    Args:
        fck (float): The characteristic compressive strength in MPa.
        z (float): The length to the areasenter of cross-section in mm
        gamma_c: Safety factor
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
    E: float,
    As: float,
    Med: float,
    Ved: float,
    Ned: float,
    delta_e: float,
    gamma_c: float = 1.5,
) -> float:
    """The design shear resistance of a web or a slab without
    shear reinforcement.

    fib Model Code 2010, Eq. (7.3-17)

    Args:
        fck (float): Characteristic strength in MPa
        z: (float): The length to the areasenter of cross-section in mm
        bw: (float): Thickness of web in cross section
        dg: (float): Maximum size of aggregate
        E: (float): The E-modulus to the materialb in MPa
        As: (float): The cross-section area in mm^2
        Med: (float): The moment working on the material in Nmm
        Ved: (float): The shear working on the material in N
        Ned: (float): The normal force working on the material in N
        delta_e (float): The exentricity of the load in mm
        gamma_c (float): Safety factor

    Returns:
        float: Design shear resistance without shear reinforcement
    """
    fsqr = min(fck**0.5, 8)

    kdg = max(32 / (16 + dg), 0.75)
    kv = (0.4 / (1 + 1500 * epsilon_x(E, As, Med, Ved, Ned, z, delta_e))) * (
        1300 / (1000 + kdg * z)
    )
    return (kv * fsqr * z * bw) / gamma_c


def v_rds_approx3(  # tror dette egentlig er vrds3
    approx_lvl_s: float,
    fck: float,
    z: float,
    bw: float,
    E: float,
    As: float,
    Med: float,
    Ved: float,
    Ned: float,
    delta_e: float,
    alfa: float,
    gamma_c: float = 1.5,
) -> float:

    """The design shear resistance of a web or a slab without
    shear reinforcement.

    fib Model Code 2010, Eq. (7.3-17)

    Args:
        fck (float): Characteristic strength in MPa
        z: (float): The length to the areasenter of cross-section in mm
        bw: (float): Thickness of web in cross section
        dg: (float): Maximum size of aggregate
        E: (float): The E-modulus to the materialb in MPa
        As: (float): The cross-section area in mm^2
        Med: (float): The moment working on the material in Nmm
        Ved: (float): The shear working on the material in N
        Ned: (float): The normal force working on the material in N
        delta_e (float): The exentricity of the load in mm
        alfa (float): Inclination of the stirrups
        gamma_c (float): Safety factor

    Returns:
        float: Design shear resistance without shear reinforcement
    """
    fsqr = min(fck**0.5, 8)
    theta_min = (20 + 10000 * epsilon_x(E, As, Med, Ved, Ned, z, delta_e))
    kv = max((0.4 / (1 + 1500 * epsilon_x(E, As, Med, Ved, Ned, z, delta_e)))*(
        1 - Ved / v_rd_max(
            approx_lvl_s, fck, bw, theta_min, z, E, As, Med, Ved,
            Ned, delta_e, alfa, gamma_c)),
            0)

    return (kv * fsqr * z * bw) / gamma_c


def v_rds(
    asw: float,
    sw: float,
    z: float,
    fywd: float,
    theta: float,
    alpha: t.Optional[float] = pi / 2,
) -> float:
    """fib Model Code 2010, Eq. (7.3-29)
    Args:
        asw (float): Area of shear reinforcement
        sw (float): Senter distance between the shear reinforcement
        z: (float): The length to the areasenter of cross-section in mm
        fywd (float): The design yield strength of the shear reinforcement
        theta (float): Inclitaniton of the compression stressfield
        alfa (float): Inclination of the stirrups

    Returns:
        The design shear resistance provided by shear reinforcement
        """
    if 45 < theta < 20:
        warnings.warn("Too high or too low compression field angel")
    return (
        (asw / sw) *
        z * fywd * ((1 / tan(theta*pi/180)) + (1 / tan(alpha*pi/180))) *
        sin(alpha*pi/180)
    )


def v_rd_max(
    approx_lvl_s: int,
    fck: float,
    bw: float,
    theta: float,
    z: float,
    E: float,
    As: float,
    Med: float,
    Ved: float,
    Ned: float,
    delta_e: float,
    alfa: float = 0,
    gamma_c: float = 1.5,
) -> float:
    """The maximum allowed shear resistance, when there is shear reinforcment

    fib Model Code 2010, eq. (7.3-26) and (7.3-24)

    Args:
        approx_lvl_s (int): Approximation level for steel
        fck (float): Characteristic strength in MPa
        z: (float): The length to the areasenter of cross-section in mm
        bw: (float): Thickness of web in cross section
        dg: (float): Maximum size of aggregate
        E: (float): The E-modulus to the materialb in MPa
        As: (float): The cross-section area in mm^2
        Med: (float): The moment working on the material in Nmm
        Ved: (float): The shear working on the material in N
        Ned: (float): The normal force working on the material in N
        delta_e (float): The exentricity of the load in mm
        alfa (float): Inclination of the stirrups
        gamma_c (float): Safety factor
    Returns:
        float: The maximum allowed shear resistance regardless of
        approximation level"""
    if approx_lvl_s == 1:
        return v_rd_max_approx1(fck, bw, theta, z, alfa, gamma_c)

    elif approx_lvl_s == 2:
        return v_rd_max_approx2(
            fck, bw, theta, z, E, As, Med, Ved, Ned, delta_e, alfa, gamma_c
        )

    elif approx_lvl_s == 3:
        return v_rd_max_approx3(
            fck, bw, theta, z, E, As, Med, Ved, Ned, delta_e, alfa, gamma_c
        )


def v_rd_max_approx1(
    fck: float,
    bw: float,
    theta: float,
    z: float,
    alfa: float = 0,
    gamma_c: float = 1.5,
) -> float:
    """The maximum allowed shear resistance, when there is shear reinforcment

    fib Model Code 2010, eq. (7.3-26) and (7.3-24)

    Args:
        fck (float): Characteristic strength in MPa
        z: (float): The length to the areasenter of cross-section in mm
        bw: (float): Thickness of web in cross section
        delta_e (float): The exentricity of the load in mm
        alfa (float): Inclination of the stirrups
        gamma_c (float): Safety factor

    Returns:
        float: The maximum allowed shear resistance regardless of
        approximation level"""
    nfc = min((30 / fck) ** (1 / 3), 1)

    return (
        0.55
        * nfc
        * (fck / gamma_c)
        * bw
        * z
        * (
            ((1 / tan(theta*pi/180)) + (1 / tan(alfa*pi/180)))
            / (1 + (1 / tan(theta*pi/180)) ** 2)
        )
    )


def v_rd_max_approx2(
    fck: float,
    bw: float,
    theta: float,
    z: float,
    E: float,
    As: float,
    Med: float,
    Ved: float,
    Ned: float,
    delta_e: float,
    alfa: float = 0,
    gamma_c: float = 1.5,
) -> float:
    """The maximum allowed shear resistance, when there is shear reinforcment

    fib Model Code 2010, eq. (7.3-26) and (7.3-24)

    Args:
        fck (float): Characteristic strength in MPa
        z: (float): The length to the areasenter of cross-section in mm
        bw: (float): Thickness of web in cross section
        dg: (float): Maximum size of aggregate
        E: (float): The E-modulus to the materialb in MPa
        As: (float): The cross-section area in mm^2
        Med: (float): The moment working on the material in Nmm
        Ved: (float): The shear working on the material in N
        Ned: (float): The normal force working on the material in N
        delta_e (float): The exentricity of the load in mm
        alfa (float): Inclination of the stirrups
        gamma_c (float): Concrete safety factor

    Returns:
        float: The maximum allowed shear resistance regardless of
        approximation level"""
    nfc = min((30 / fck) ** (1 / 3), 1)

    epsilon_1 = epsilon_x(E, As, Med, Ved, Ned, z, delta_e) + (
        epsilon_x(E, As, Med, Ved, Ned, z, delta_e) + 0.002
    ) * ((1 / tan(theta*pi/180)) ** 2)
    k_epsilon = min(1 / (1.2 + 55 * epsilon_1), 0.65)

    return (
        k_epsilon * nfc * (fck / gamma_c) * bw * z *(
            ((1 / tan(theta*pi/180)) + (1 / tan(alfa*pi/180)))
            / (1 + (1 / tan(theta*pi/180)) ** 2)
            )
        )


def v_rd_max_approx3(
    fck: float,
    bw: float,
    theta_min: float,
    z: float,
    E: float,
    As: float,
    Med: float,
    Ved: float,
    Ned: float,
    delta_e: float,
    alfa: float = 0,
    gamma_c: float = 1.5,
) -> float:
    """The maximum allowed shear resistance, when there is shear reinforcment

    fib Model Code 2010, eq. (7.3-26) and (7.3-24)

    Args:
        fck (float): Characteristic strength in MPa
        z: (float): The length to the areasenter of cross-section in mm
        bw: (float): Thickness of web in cross section
        dg: (float): Maximum size of aggregate
        E: (float): The E-modulus to the materialb in MPa
        As: (float): The cross-section area in mm^2
        Med: (float): The moment working on the material in Nmm
        Ved: (float): The shear working on the material in N
        Ned: (float): The normal force working on the material in N
        delta_e (float): The exentricity of the load in mm
        alfa (float): Inclination of the stirrups
        gamma_c (float): Concrete safety factor

    Returns:
        float: The maximum allowed shear resistance regardless of
        approximation level"""
    nfc = min((30 / fck) ** (1 / 3), 1)

    epsilon_1 = epsilon_x(E, As, Med, Ved, Ned, z, delta_e) + (
        epsilon_x(E, As, Med, Ved, Ned, z, delta_e) + 0.002
    ) * ((1 / tan(theta_min*pi/180)) ** 2)
    k_epsilon = min(1 / (1.2 + 55 * epsilon_1), 0.65)

    return (
        k_epsilon
        * nfc
        * (fck / gamma_c)
        * bw
        * z
        * (
            ((1 / tan(theta_min*pi/180)) + (1 / tan(alfa*pi/180)))
            / (1 + (1 / tan(theta_min*pi/180)) ** 2)
        )
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
        approx_lvl_h: What approximation level we want for hollow core
        f_ctd: The design value of concrete axial tensile strength
        i_c: The second moment of area
        s_c: The first moment of area, abouve and about the centriodal axis
        b_w: The width of the cross-section at the centroidal axis
        sigma_cp: The compressive stress at centroidal axis due to prestress
        alfa_l: l_x/(1.2*l_bd0)
        l_x: distance between edge and point of failure (Figure: 7.3-12)
        l_bd0: follows 7.13-5
        S_cy: The first moment of area above y
        b_wy: The width at hight y
        y: The hight at of the critical point at the line of failure
        y_c: The hight of the concrete centroidal axis
        A_c: The area of concrete cross-section
        A_cy: The area of concrete cross-section above y
        y_pt: The hight of centroidal axis of prestressed steel
        f_p_lx: The prestressing force at the distance l_x
        f_p_lx_dx: The derivative of prestressing force at the distance l_x

    Return:
        The maximum allowed shear force in a hollow core. Regardless of the
        approximation level"""
    if approx_lvl_h == 1:
        return v_rd_ct_approx1(f_ctd, i_c, s_c, b_w, sigma_cp, l_x, l_bd0)
    elif approx_lvl_h == 2:
        return v_rd_ct_approx2(
            f_ctd, i_c, l_x, l_bd0, S_cy, b_wy, y, y_c,
            A_c, A_cy, y_pt, f_p_lx, f_p_lx_dx
        )


def v_rd_ct_approx1(
    f_ctd: float,
    i_c: float,
    s_c: float,
    b_w: float,
    sigma_cp: float,
    l_x: float,
    l_bd0: float
) -> float:
    """Calculating level 1 approximation.

    Args:
        f_ctd: The design value of concrete axial tensile strength
        i_c: The second moment of area
        s_c: The first moment of area, abouve and about the centriodal axis
        b_w: The width of the cross-section at the centroidal axis
        sigma_cp: The compressive stress at centroidal axis due to prestress
        l_x: distance between edge and point of failure (Figure: 7.3-12)
        l_bd0: follows 7.13-5
    return:
        Vrd Appoximation 1 for hollow slabs"""
    alfa_l = l_x/(1.2*l_bd0)
    return 0.8*((i_c*b_w)/s_c)*((f_ctd**2)*+alfa_l*sigma_cp*f_ctd)**0.5


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
    Args:
        f_ctd: The design value of concrete axial tensile strength
        i_c: The second moment of area
        l_x: distance between edge and point of failure (Figure: 7.3-12)
        l_bd0: follows 7.13-5
        S_cy: The first moment of area above y
        b_wy: The width at hight y
        y: The hight at of the critical point at the line of failure
        y_c: The hight of the concrete centroidal axis
        A_c: The area of concrete cross-section
        A_cy: The area of concrete cross-section above y
        y_pt: The hight of centroidal axis of prestressed steel
        f_p_lx: The prestressing force at the distance l_x
        f_p_lx: The derivative of prestressing force at the distance l_x

    sigma_cpy: The compressiv stress in the concrete at hight y and l_x
    tau_cpy: The shear stress due to prestress at hight y and l_x

    return:
        The maximum shear force for level 2 approximation"""

    alfa_l = l_x/(1.2*l_bd0)
    sigma_cpy = ((1/A_c)*((y_c-y)/i_c))*f_p_lx
    tau_cpy = (1/b_wy)*((A_cy/A_c)-(S_cy*(y_c-y_pt))/i_c)*f_p_lx_dx
    return (i_c*b_wy/S_cy)*(((f_ctd**2)+alfa_l*sigma_cpy*f_ctd)**0.5-tau_cpy)


def tau_edi(beta: float, v_ed: float, z: float, b_i: float):
    """Shear at the interface between cocrete cast at different times
    fib Model Code 2010, eq. (7.3-49)
    Args:
        beta: The ratio of longitudinal force in the new concrete and
        the longitudinal force in either compression or tension zone
        z: The inner lever arm of composed section
        b_i: The width of the inerface
        v_ed: The shear force at the interface

    return:
        The shear force that should be used at the intersection"""
    return (beta*v_ed)/(z*b_i)


def tau_rdi_without_reinforceent(
    c_a: float,
    f_ctd: float,
    mu: float,
    sigma_n: float,
    f_ck: float,
    f_cd: float
):
    """Shear resistance without reinforcement at the intesection with
    different casting time

    fib Model Code 2010, eq. (7.3-50)

    Args:
        c_a: The coefficient for adhesive bond (tabel 7.3-1)
        mu: The friction coefficient
        sigma_n: The loweat expected compressiv stress from normal forces
        f_ck: Characteristic strength in MPa
        f_cd: The design value of cylinder compressive strength concrete

    return:
        Shear resistance without reinforcement at the intesection with
    different casting time"""

    v = min(0.55*(30/f_ck)**(1/3), 0.55)
    return min((c_a*f_ctd) + (mu * sigma_n), 0.5*v*f_cd)


def tau_rdi_with_reinforceent(
    c_r: float,
    k1: float,
    k2: float,
    mu: float,
    ro: float,
    sigma_n: float,
    alfa: float,
    beta_c: float,
    f_ck: float,
    f_yd: float,
    f_cd: float
):
    """Shear resistance with reinforcement or dowels at the intesection with
    different casting time

    fib Model Code 2010, eq. (7.3-51)

    Args:
        c_r: The coefficient for aggregate interlock effects (tabel 7.3-2)
        k1: The interction coefficient for tensile
        force activated in reinforcment (tabel 7.3-2)
        k2: The interction coeffiction for flexural resistance (tabel 7.3-2)
        mu: The friction coefficient (tabel 7.3-2)
        ro: reinforcement ratio of reinforcment crossing the interface
        sigma_n: The loweat expected compressiv stress from normal forces
        alfa: The inclination of reinforcement crossing the
        interface (tabel 7.3-14)
        beta_c: The coefficient for strength of
        compresstion strut (tabel 7.3-2)
        f_ck: Characteristic strength in MPa
        f_yd: design strength of reinforment steel
        f_cd: The design value of cylinder compressive strength concrete

    return:
        Shear resistance with reinforcement at intesection with
        different casting time"""
    v = min(0.55*(30/f_ck)**(1/3), 0.55)
    return min(
        (c_r*f_ck**(1/3))+(mu*sigma_n)+k1*ro * f_yd *
        (ro*sin(alfa)+cos(alfa)+k2*ro*(f_yd*f_cd)**0.5),
        beta_c*v*f_cd
        )


def v_ed_ti(t_ed: float, a_k: float, z_i: float):
    """Shear force from an torsion
    fib Model Code 2010, eq. (7.3-53)

    Args:
        t_ed: The acting torsyion force in the cross section
        z_i: Can be found in figure 7.3-18
        a_k: Can be found in figure 7.3-18
    
    Returns:
        The shear force that will ocurre due to torsion force."""
    
    return t_ed*z_i/(2*a_k)


def t_rd_max(k_c, f_ck, gamma_c, d_k, a_k, theta):
    """The maximum allowed torsion allowed
    fib Model Code 2010, eq. (7.3-56)
    
    args:
        f_ck: Characteristic strength in MPa
        gamma_c: Concrete safety factor
        d_k: Is the diameter in the smalest circel in the cross section
        a_k: Can be found in figure 7.3-18
        theta: Inclitaniton of the compression stressfield

    return:
    """
    nfc = min((30 / fck) ** (1 / 3), 1)

    if approx_lvl_s == 1:
        0.55
    elif approx_lvl_s == 2:
        epsilon_1 = epsilon_x(E, As, Med, Ved, Ned, z, delta_e) + (
        epsilon_x(E, As, Med, Ved, Ned, z, delta_e) + 0.002
        ) * ((1 / tan(theta*pi/180)) ** 2)
        k_epsilon = 1 / (1.2 + 55 * epsilon_1)
    elif approx_lvl_s == 3:
        epsilon_1 = epsilon_x(E, As, Med, Ved, Ned, z, delta_e) + (
            epsilon_x(E, As, Med, Ved, Ned, z, delta_e) + 0.002
            ) * ((1 / tan(theta_min*pi/180)) ** 2)
        k_epsilon = min(1 / (1.2 + 55 * epsilon_1), 0.65)


def t_rd(
    t_ed: float,
    v_ed: float,
    t_rd_max: float,
    )

    """hei
    
    Args:
        t
        
    return:
        e"""
    v_rd_max(approx_lvl_s, fck, bw, theta, z, E, As, Med,Ved, Ned, delta_e, alfa),
