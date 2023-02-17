"""A collection of shear formulas for concrete"""
import warnings
from math import pi, tan, sin, cos


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
        Med (Float): The moment working on the material in Nmm
        Ved (float): The shear force working on the material in N
        Ned: (float): The normal force working on the material in N
        z: (float): The effective shear depth in mm
        delta_E (float): The exentricity of the load in mm
    Returns:
        float: The longitudinal strain"""
    return max((1 / (2 * E_s * As)) * (
        (abs(Med) / z) + abs(Ved) + abs(Ned) * ((1 / 2) + (delta_e / z)))
    ), 0


def v_rd(
    approx_lvl_c: int,
    approx_lvl_s: int,
    reinforcment: bool,
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
    alfa: float,
    gamma_c: float,
    asw: float,
    sw: float,
    f_ywd: float,
    theta: float,
) -> float:
    """Compute the shear resistance of a web or slab.

    fib Model Code 2010, Eq. (7.3-11)

    Args:
        approx_lvl_c (int): Approximation level for concrete
        approx_lvl_s (int): Approximation level for reinforcment
        reinforcment (bool): Shear reinforced concrete or no shear
        reinforcement
        fck (float): Characteristic strength in MPa
        z: (float): distances between the centerline of the
        compressive chord and the reinforcement in mm
        bw: (float): Thickness of web in cross section in mm
        dg: (float): Maximum size of aggregate in mm
        E_s: (float): The E_s-modulus to the materialb in MPa
        As: (float): The cross-section area of reinforcement in mm^2
        Med: (float): The moment working on the material in Nmm
        Ved: (float): The shear working on the material in N
        Ned: (float): The normal force working on the material in N
        delta_e (float): The exentricity of the load in mm
        alfa (float): Inclination of the stirrups in degrees
        gamma_c (float): Concrete safety factor
        asw (float): Area of shear reinforcement in mm^2
        sw (float): Senter distance between the shear reinforcement in mm
        f_ywd (float): The design yield strength of the shear reinforcement
        theta (float): Inclitantion of the compression stressfield in degrees

    Returns:
        float: Design shear resistance
    """

    if not (approx_lvl_c == 1 or 2) and not approx_lvl_s == 3:
        warnings.warn("Not a valid approximation level")

    if not reinforcment:
        return v_rdc(
            approx_lvl_c, fck, z, bw, dg, E_s, As,
            Med, Ved, Ned, delta_e, alfa, gamma_c,
            )
    if reinforcment and approx_lvl_s == 3:
        return min(v_rdc(
            approx_lvl_c, fck, z, bw, dg, E_s, As, Med,
            Ved, Ned, delta_e, alfa, gamma_c) +
            v_rds(asw, sw, z, f_ywd, theta, alfa),
            v_rd_max(
            approx_lvl_s, fck, bw, theta, z, E_s, As, Med,
            Ved, Ned, delta_e, alfa, gamma_c,
            )
        )
    elif reinforcment and approx_lvl_s == 1 or 2:
        return min(v_rds(asw, sw, z, f_ywd, theta, alfa), v_rd_max(
            approx_lvl_s, fck, bw, theta, z, E_s, As,
            Med, Ved, Ned, delta_e, alfa, gamma_c,
            )
        )


def v_rdc(
    approx_lvl_c: int,
    approx_lvl_s: int,
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
    alfa: float,
    gamma_c: float = 1.5,
) -> float:

    """The design shear resistance of a web or a slab without
    shear reinforcement.

    fib Model Code 2010, Eq. (7.3-17)

    Args:
        approx_lvl_c (int): Approximation level for concrete
        approx_lvl_s (int): Approximation level for reinforcement
        reinforcment (bool): Shear reinforced concrete or no shear
        reinforcement
        fck (float): Characteristic strength in MPa
        z: (float): The length to the areasenter of cross-section in mm
        bw: (float): Thickness of web in cross section in mm
        dg: (float): Maximum size of aggregate
        E_s: (float): The E_s-modulus to the materialb in MPa
        As: (float): The cross-section area in mm^2
        Med: (float): The moment working on the material in Nmm
        Ved: (float): The shear working on the material in N
        Ned: (float): The normal force working on the material in N
        delta_e (float): The exentricity of the load in mm
        alfa (float): Inclination of the stirrups in degrees
        gamma_c (float): Safety factor

    Returns:
        float: The design shear resistance attributed to the concrete
    """

    if approx_lvl_c == 1:
        return v_rdc_approx1(fck, z, bw, gamma_c)

    elif approx_lvl_c == 2:
        return v_rdc_approx2(
            fck, z, bw, dg, E_s, As, Med, Ved, Ned, delta_e, gamma_c
        )

    elif approx_lvl_s == 3:
        return v_rds_approx3(
            approx_lvl_s, fck, z, bw, E_s, As, Med, Ved,
            Ned, delta_e, alfa, gamma_c
        )


def v_rdc_approx1(
    fck: float,
    z: float,
    bw: float,
    gamma_c: float = 1.5,
) -> float:
    """For members with no segnificant axal load, with fyk <= 600 Mpa,
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
    """In higher strength concrete and light-weight aggregate concretes,
    the fracture surface may go through the aggregate particles,
    rather then around, reducing the crack roughness

    fib Model Code 2010, Eq. (7.3-17), (7.3-20) and (7.3-21)

    Args:
        fck (float): Characteristic strength in MPa
        z (float): The length to the areasenter of cross-section in mm
        bw (float): Thickness of web in cross section
        dg (float): Maximum size of aggregate
        E_s (float): The E_s-modulus to the materialb in MPa
        As (float): The cross-section area in mm^2
        Med (float): The moment working on the material in Nmm
        Ved (float): The shear working on the material in N
        Ned (float): The normal force working on the material in N
        delta_e (float): The exentricity of the load in mm
        gamma_c (float): Safety factor

    Returns:
        float: Design shear resistance without shear reinforcement
    """
    fsqr = min(fck**0.5, 8)

    k_dg = max(32 / (16 + dg), 0.75)
    kv = (0.4 / (1 + 1500 * epsilon_x(E_s, As, Med, Ved, Ned, z, delta_e))) * (
        1300 / (1000 + k_dg * z)
    )
    return (kv * fsqr * z * bw) / gamma_c


def v_rds_approx3(  # tror dette egentlig er vrds3
    approx_lvl_s: float,
    fck: float,
    z: float,
    bw: float,
    E_s: float,
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

    fib Model Code 2010, Eq. (7.3-17), (7.3-39) and (7.3-43)

    Args:
        fck (float): Characteristic strength in MPa
        z: (float): The length to the areasenter of cross-section in mm
        bw: (float): Thickness of web in cross section
        dg: (float): Maximum size of aggregate
        E_s: (float): The E_s-modulus to the materialb in MPa
        As: (float): The cross-section area in mm^2
        Med: (float): The moment working on the material in Nmm
        Ved: (float): The shear working on the material in N
        Ned: (float): The normal force working on the material in N
        delta_e (float): The exentricity of the load in mm
        alfa (float): Inclination of the stirrups in degrees
        gamma_c (float): Safety factor

    Returns:
        float: Design shear resistance without shear reinforcement
    """
    fsqr = min(fck**0.5, 8)
    theta_min = (20 + 10000 * epsilon_x(E_s, As, Med, Ved, Ned, z, delta_e))
    kv = max((0.4 / (1 + 1500 * epsilon_x(E_s, As, Med, Ved, Ned, z, delta_e)))*(
        1 - Ved / v_rd_max(
            approx_lvl_s, fck, bw, theta_min, z, E_s, As, Med, Ved,
            Ned, delta_e, alfa, gamma_c)),
            0)

    return (kv * fsqr * z * bw) / gamma_c


def v_rds(
    asw: float,
    sw: float,
    z: float,
    f_ywd: float,
    theta: float,
    alpha: float,
) -> float:
    """fib Model Code 2010, Eq. (7.3-25) and (7.3-29)
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
    if 45 < theta < 20:
        warnings.warn("Too high or too low compression field angel")
    return (
        (asw / sw) *
        z * f_ywd * ((1 / tan(theta*pi/180)) + (1 / tan(alpha*pi/180))) *
        sin(alpha*pi/180)
    )


def v_rd_max(
    approx_lvl_s: int,
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
    alfa: float = 0,
    gamma_c: float = 1.5,
) -> float:
    """The maximum allowed shear resistance, when there is shear reinforcment

    fib Model Code 2010, eq. (7.3-24) and (7.3-26)

    Args:
        approx_lvl_s (int): Approximation level for steel
        fck (float): Characteristic strength in MPa
        z: (float): The length to the areasenter of cross-section in mm
        bw: (float): Thickness of web in cross section
        dg: (float): Maximum size of aggregate
        E_s: (float): The E_s-modulus to the materialb in MPa
        As: (float): The cross-section area in mm^2
        Med: (float): The moment working on the material in Nmm
        Ved: (float): The shear working on the material in N
        Ned: (float): The normal force working on the material in N
        delta_e (float): The exentricity of the load in mm
        alfa (float): Inclination of the stirrups in degrees
        gamma_c (float): Safety factor
    Returns:
        float: The maximum allowed shear resistance regardless of
        approximation level"""
    if approx_lvl_s == 1:
        return v_rd_max_approx1(fck, bw, theta, z, alfa, gamma_c)

    elif approx_lvl_s == 2:
        return v_rd_max_approx2(
            fck, bw, theta, z, E_s, As, Med, Ved, Ned, delta_e, alfa, gamma_c
        )

    elif approx_lvl_s == 3:
        return v_rd_max_approx3(
            fck, bw, z, E_s, As, Med, Ved, Ned, delta_e, alfa, gamma_c
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

    fib Model Code 2010, eq. (7.3-37)

    Args:
        fck (float): Characteristic strength in MPa
        z: (float): The length to the areasenter of cross-section in mm
        bw: (float): Thickness of web in cross section
        delta_e (float): The exentricity of the load in mm
        alfa (float): Inclination of the stirrups in degrees
        gamma_c (float): Safety factor

    Returns:
        float: The maximum allowed shear resistance regardless of
        approximation level"""
    nfc = min((30 / fck) ** (1 / 3), 1)

    return (
        0.55 * nfc * (fck / gamma_c) * bw * z * (
            ((1 / tan(theta*pi/180)) + (1 / tan(alfa*pi/180)))
            / (1 + (1 / tan(theta*pi/180)) ** 2)
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
    alfa: float = 0,
    gamma_c: float = 1.5,
) -> float:
    """The maximum allowed shear resistance, when there is shear reinforcment

    fib Model Code 2010, eq. (7.3-24), (7.3-40) and (7.3-41)

    Args:
        fck (float): Characteristic strength in MPa
        z: (float): The length to the areasenter of cross-section in mm
        bw: (float): Thickness of web in cross section
        dg: (float): Maximum size of aggregate
        E_s: (float): The E_s-modulus to the materialb in MPa
        As: (float): The cross-section area of reinforcement in mm^2
        Med: (float): The moment working on the material in Nmm
        Ved: (float): The shear working on the material in N
        Ned: (float): The normal force working on the material in N
        delta_e (float): The exentricity of the load in mm
        alfa (float): Inclination of the stirrups in degrees
        gamma_c (float): Concrete safety factor

    Returns:
        float: The maximum allowed shear resistance regardless of
        approximation level"""
    nfc = min((30 / fck) ** (1 / 3), 1)

    epsilon_1 = epsilon_x(E_s, As, Med, Ved, Ned, z, delta_e) + (
        epsilon_x(E_s, As, Med, Ved, Ned, z, delta_e) + 0.002
    ) * ((1 / tan(theta*pi/180)) ** 2)
    k_epsilon = min(1 / (1.2 + 55 * epsilon_1), 0.65)

    return (
        k_epsilon * nfc * (fck / gamma_c) * bw * z * (
            ((1 / tan(theta*pi/180)) + (1 / tan(alfa*pi/180)))
            / (1 + (1 / tan(theta*pi/180)) ** 2)
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
    alfa: float = 0,
    gamma_c: float = 1.5,
) -> float:
    """The maximum allowed shear resistance, when there is shear reinforcment

    fib Model Code 2010, eq. (7.3-24), (7.3-40) and (7.3-41)

    Args:
        fck (float): Characteristic strength in MPa
        z: (float): The length to the areasenter of cross-section in mm
        bw: (float): Thickness of web in cross section
        dg: (float): Maximum size of aggregate
        E_s: (float): The E_s-modulus to the materialb in MPa
        As: (float): The cross-section area of reinforcement in mm^2
        Med: (float): The moment working on the material in Nmm
        Ved: (float): The shear working on the material in N
        Ned: (float): The normal force working on the material in N
        delta_e (float): The exentricity of the load in mm
        alfa (float): Inclination of the stirrups in degrees
        gamma_c (float): Concrete safety factor

    Returns:
        float: The maximum allowed shear resistance regardless of
        approximation level"""
    nfc = min((30 / fck) ** (1 / 3), 1)
    theta_min = 20 + 10000 * epsilon_x(E_s, As, Med, Ved, Ned, z, delta_e)

    epsilon_1 = epsilon_x(E_s, As, Med, Ved, Ned, z, delta_e) + (
        epsilon_x(E_s, As, Med, Ved, Ned, z, delta_e) + 0.002
    ) * ((1 / tan(theta_min*pi/180)) ** 2)
    k_epsilon = min(1 / (1.2 + 55 * epsilon_1), 0.65)

    return (
        k_epsilon * nfc * (fck / gamma_c) * bw * z
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
    alfa_l = l_x/(1.2*l_bd0)
    return 0.8*((i_c*b_w)/s_c)*((f_ctd**2)+alfa_l*sigma_cp*f_ctd)**0.5


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
        f_ctd (float): Design value of concrete axial tensile strength in MPa
        i_c (float): The second moment of area in mm^4
        l_x (float): Distance from the edge to point of
        failure (Figure: 7.3-12)
        l_bd0 (float): follows 7.13-5
        S_cy (float): The first moment of area above y in mm^3
        b_wy (float): The width at hight y in mm
        y (float): The highh of the critical point at the line of failure in mm
        y_c (float): The hight of the concrete centroidal axis in mm
        A_c (float): The area of concrete cross-section in mm^2
        A_cy (float): The area of concrete cross-section above y in mm^2
        y_pt (float): The hight of centroidal axis of prestressed steel in mm
        f_p_lx (float): The prestressing force at the distance l_x in N
        f_p_lx (float): The derivative of prestressing force at the
        distance l_x with respect to dx

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
        beta (float): The ratio of longitudinal force in the new concrete and
        the longitudinal force in either compression or tension zone
        z (float): The inner lever arm of the composed section in mm
        b_i (float): The width of the inerface in mm
        v_ed (float): The shear force at the interface in N

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
        c_a (float): The coefficient for adhesive bond (tabel 7.3-1)
        mu (float): The friction coefficient
        sigma_n (float): The loweat expected compressiv stress from
        normal forces in MPa
        f_ck (float): Characteristic strength in MPa
        f_cd (float): The design value of cylinder compressive
        strength concrete in MPa

    return:
        The shear resistance without reinforcement at the intesection with
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
        c_r (float): Coefficient for aggregate interlock effects (tabel 7.3-2)
        k1 (float): The interction coefficient for tensile
        force activated in reinforcment (tabel 7.3-2)
        k2 (float): The interction coeffiction for flexural
        resistance (tabel 7.3-2)
        mu (float): The friction coefficient (tabel 7.3-2)
        ro (float): The reinforcement ratio of reinforing steel
        crossing the interface
        sigma_n (float): The loweat expected compressiv stress resulting
        from normal forces acting on the interface in MPa
        alfa (float): The inclination of reinforcement crossing the
        interface (tabel 7.3-14)
        beta_c (float): The coefficient for strength of
        compresstion strut (tabel 7.3-2)
        f_ck (float): Characteristic strength in MPa
        f_yd (float): design strength of reinforment steel in MPa
        f_cd (float): The design value of cylinder compressive
        strength concrete

    return:
        Shear resistance with reinforcement at intesection with
        different casting time"""
    v = min(0.55*(30/f_ck)**(1/3), 0.55)
    return min(
        (c_r*f_ck**(1/3))+(mu*sigma_n)+k1*ro * f_yd *
        (ro*sin(alfa*pi/180)+cos(alfa*pi/180)+k2*ro*(f_yd*f_cd)**0.5),
        beta_c*v*f_cd
        )


def v_ed_ti(t_ed: float, a_k: float, z_i: float):
    """Shear force due to torsion

    fib Model Code 2010, eq. (7.3-53)

    Args:
        t_ed: The acting torsion force in the cross section in Nmm
        z_i: Can be found in figure 7.3-18
        a_k: Can be found in figure 7.3-18

    Returns:
        The shear force that will ocurre due to torsion force."""
    return t_ed*z_i/(2*a_k)


def t_rd_max(
        f_ck: float, gamma_c: float, d_k: float, a_k: float,
        theta: float, approx_lvl_s: int, E_s: float, As: float,
        Med: float, Ved: float, Ned: float, z: float, delta_e: float
) -> float:
    """The maximum allowed torsion allowed
    fib Model Code 2010, eq. (7.3-56)
    args:
        f_ck (float): Characteristic strength in MPa
        gamma_c (float): Concrete safety factor
        d_k (float): Is the diameter in the smalest circel in the cross section
        a_k: Can be found in figure 7.3-18
        theta (float): Inclitaniton of the compression stressfield in degrees
        approx_lvl_s (int): Approximation method for cocrete with reinforcement
        E_s (float): The E_s-modulus to the materialb in MPa
        As (float): The cross-section reinforcement in mm^2
        Med (float): The moment working on the material in Nmm
        Ved (float): The shear working on the material in N
        Ned (float): The normal force working on the material in N
        z (float): distances between the centerline of the
        compressive chord and the reinforcement in mm
        delta_e (float): the exentrisity of the load in mm

    return:
        The maximum allowed torsion allowed
    """
    t_ef = d_k/8
    nfc = min((30 / f_ck) ** (1 / 3), 1)

    if approx_lvl_s == 1:
        k_epsilon = 0.55
    elif approx_lvl_s == 2:
        epsilon_1 = epsilon_x(E_s, As, Med, Ved, Ned, z, delta_e) + (
            epsilon_x(E_s, As, Med, Ved, Ned, z, delta_e) + 0.002
        ) * ((1 / tan(theta*pi/180)) ** 2)
        k_epsilon = 1 / (1.2 + 55 * epsilon_1)
    elif approx_lvl_s == 3:
        epsilon_1 = epsilon_x(E_s, As, Med, Ved, Ned, z, delta_e) + (
            epsilon_x(E_s, As, Med, Ved, Ned, z, delta_e) + 0.002
            ) * ((1 / tan(theta*pi/180)) ** 2)
        k_epsilon = min(1 / (1.2 + 55 * epsilon_1), 0.65)
    k_c = nfc*k_epsilon

    return k_c*f_ck*t_ef*2*a_k*sin(theta*pi/180)*cos(theta*pi/180)/gamma_c


def t_rd(
    t_ed: float,
    v_ed: float,
    approx_lvl_s: int,
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
    alfa: float,
    f_ck: float,
    d_k: float,
    a_k: float,
    gamma_c: float = 1.5,
) -> bool:
    """Checks if the combination of torstion ans shear is ok

    fib Model Code 2010, eq. (7.3-56)

    Args:
        approx_lvl_s (int): Approximation level for steel
        fck (float): Characteristic strength in MPa
        z: (float): The length to the areasenter of cross-section in mm
        bw: (float): Thickness of web in cross section
        dg: (float): Maximum size of aggregate
        E_s: (float): The E_s-modulus to the materialb in MPa
        As: (float): The cross-section area in mm^2
        Med: (float): The moment working on the material in Nmm
        Ved: (float): The shear working on the material in N
        Ned: (float): The normal force working on the material in N
        delta_e (float): The exentricity of the load in mm
        alfa (float): Inclination of the stirrups in degrees
        f_ck: Characteristic strength in MPa
        d_k: Is the diameter in the smalest circel in the cross section
        a_k: Can be found in figure 7.3-18
        gamma_c (float): Safety factor
    return:
        Returns a bool that is true if the criteria for torsion and
        shear is fulfilled"""
    if (
        (t_ed
         / t_rd_max(
            f_ck, gamma_c, d_k, a_k, theta, approx_lvl_s,
            E_s, As, Med, Ved, Ned, z, delta_e))**2 +
        (v_ed
         / v_rd_max(
            approx_lvl_s, fck, bw, theta, z, E_s, As, Med,
            Ved, Ned, delta_e, alfa, gamma_c))**2 <= 1):
        check = True

    else:
        check = False
    return check


def b_0(v_ed: float, v_prep_d_max: float) -> float:
    """Gives the general output for b_0, shear-resisting control perimeter.

    fib Model Code 2010, eq. (7.3-57)
    Args:
        V_ed (float): The acting shear force from the columns
        v_prep_d_max (float): The maximum shear force per unit length
        perpendiculerer to the basic control parameter (Figure 7.3-24)

    Return:
        The shear-resisting control perimeter, b_0"""

    return v_ed/v_prep_d_max


def m_ed(
    v_ed: float, e_u: float, r_sx: float, r_sy: float, l_min: float,
    inner: bool, edge_par: bool, edge_per: bool, corner: bool
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
    b_s = min(1.5 * (r_sx*r_sy)**0.5, l_min)
    if inner:
        return v_ed*((1/8)+e_u/(2*b_s))
    elif edge_par:
        return max(v_ed*((1/8)+e_u/(2*b_s)), v_ed/4)
    elif edge_per:
        return v_ed*((1/8)+e_u/(b_s))
    elif corner:
        return max(v_ed*((1/8)+e_u/(b_s)), v_ed/2)


def psi_punching(
    r_s: float, l_x: float, l_y: float, f_yd: float, d: float, e_s: float,
    approx_lvl_p: float, v_ed: float, e_u: float, r_sx: float, r_sy: float,
    l_min: float, inner: bool, edge_par: bool, edge_per: bool,
    corner: bool, m_rd: float, m_pd: float
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

    r_s = 0.22 * l_x

    if approx_lvl_p == 1:
        if not (0.5 < l_x/l_y < 2):
            warnings.warn("Reconsider maximum r_s value")
        psi = 1.5 * r_s * f_yd / (d*e_s)

    elif approx_lvl_p == 2 or 3:
        if 0.5 < l_x/l_y < 2:
            warnings.warn("Reconsider maximum r_s value")
        psi = (1.5 * r_s * f_yd / (d*e_s))*((m_ed(
            v_ed, e_u, r_sx, r_sy, l_min, inner, edge_par, edge_per, corner
        )-m_pd)/(m_rd-m_pd))**1.5

    return psi


def v_rdc_punching(
    r_s: float, l_x: float, l_y: float, f_yd: float, d: float, e_s: float,
    approx_lvl_p: float, dg: float, f_ck: float, d_v: float, v_ed: float,
    e_u: float, r_sx: float, r_sy: float, l_min: float, inner: bool,
    edge_par: bool, edge_per: bool, corner: bool, m_rd: float, m_pd: float,
    v_prep_d_max: float, gamma_c: float = 1.5
) -> float:

    """Punching resistance from the concrete
    fib Model Code 2010, eq. (7.3-61), (7.3-62) and (7.3-63)
    args:
        r_s (float): Denotes the position where the radial bending moment is
        zero with respect to support axis
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

    r_s = 0.22 * l_x
    k_dg = max(32/(16+dg), 0.75)
    k_psi = min(
        1/(1.5+0.9*k_dg*d*psi_punching(
            r_s, l_x, l_y, f_yd, d, e_s,
            approx_lvl_p, v_ed, e_u, r_sx, r_sy, l_min, inner, edge_par,
            edge_per, corner, m_rd, m_pd
        )), 0.6)

    return k_psi*b_0(v_ed, v_prep_d_max)*d_v * (f_ck**0.5)/gamma_c


def v_rds_punching(
    e_u: float, b_u: float, r_s: float, l_x: float, l_y: float, f_yd: float,
    d: float, e_s: float, approx_lvl_p: float, v_ed: float, r_sx: float,
    r_sy: float, l_min: float, inner: bool, edge_par: bool, edge_per: bool,
    corner: bool, m_rd: float, m_pd: float, alfa: float,
    f_bd: float, f_ywd: float, kam_w: float, a_sw: float
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
    k_e = 1 / (1+e_u/b_u)
    sigma_swd = min(
        (e_s*psi_punching(
            r_s, l_x, l_y, f_yd, d, e_s,
            approx_lvl_p, v_ed, e_u, r_sx, r_sy, l_min, inner, edge_par,
            edge_per, corner, m_rd, m_pd
        )/6)*(sin(alfa*pi/180)+cos(alfa*pi/180)) *
        (sin(alfa*pi/180)+f_bd*d/(f_ywd*kam_w)), f_ywd)

    if (a_sw*k_e*sigma_swd*sin(alfa*pi/180)) < 0.5 * v_ed:
        warnings.warn("""In order to ensure sufficent deformation capacity,
                      the shear resistance in punching most increase""")

    return a_sw*k_e*sigma_swd*sin(alfa*pi/180)


def v_rd_max_punching(
    r_s: float, l_x: float, l_y: float, f_yd: float, d: float, e_s: float,
    approx_lvl_p: float, v_ed: float, e_u: float, r_sx: float, r_sy: float,
    l_min: float, inner: bool, edge_par: bool, edge_per: bool, dg: float,
    corner: bool, m_rd: float, m_pd: float, v_prep_d_max: float,
    d_v: float, f_ck: float, d_head: bool, stirrups_compression: bool,
    gamma_c: float = 1.5
):
    """Finds the maximum value you can have for what the shear
    reinforcement can take.
     
    Args:
    
    Return
    """

    if d_head:
        k_sys = 2.8
    elif stirrups_compression:
        k_sys = 2.4
    else:
        k_sys = 2

    k_dg = max(32/(16+dg), 0.75)
    k_psi = min(
        1/(1.5+0.9*k_dg*d*psi_punching(
            r_s, l_x, l_y, f_yd, d, e_s,
            approx_lvl_p, v_ed, e_u, r_sx, r_sy, l_min, inner, edge_par,
            edge_per, corner, m_rd, m_pd
        )), 0.6)
    
    return min((k_sys*k_psi*b_0(v_ed, v_prep_d_max)*d_v*f_ck**0.5) /
               gamma_c, (b_0(v_ed, v_prep_d_max)*d_v*f_ck**0.5)/gamma_c)


def v_rd_punching(
    e_u, b_u, r_s: float, l_x: float, l_y: float, f_yd: float, d: float,
    e_s: float, approx_lvl_p: float, v_ed: float, r_sx: float, r_sy: float,
    l_min: float, inner: bool, edge_par: bool, edge_per: bool, corner: bool,
    m_rd: float, m_pd: float, alfa: float, f_bd: float, f_ywd: float,
    kam_w: float, a_sw: float, dg: float, f_ck: float, d_v: float,
    v_prep_d_max: float, d_head: bool, stirrups_compression: bool,
    gamma_c: float = 1.5
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

    return: The punching resistance for punching"""

    return v_rdc_punching(
        r_s, l_x, l_y, f_yd, d, e_s, approx_lvl_p, dg, f_ck, d_v, v_ed, e_u,
        r_sx, r_sy, l_min, inner, edge_par, edge_per, corner, m_rd,
        v_prep_d_max, gamma_c
        ) + min(v_rds_punching(
            e_u, b_u, r_s, l_x, l_y, f_yd, d, e_s, approx_lvl_p, v_ed, r_sx,
            r_sy, l_min, inner, edge_par, edge_per, corner, m_rd, m_pd, alfa,
            f_bd, f_ywd, kam_w, a_sw
        ), v_rd_max_punching(
            r_s, l_x, l_y, f_yd, d, e_s,
            approx_lvl_p, v_ed, e_u, r_sx, r_sy, l_min, inner, edge_par,
            edge_per, dg, corner, m_rd, m_pd, v_prep_d_max,
            d_v, f_ck, d_head, stirrups_compression, gamma_c)
)
