"Covers torsion in Model code 2010, 7.3.4"

from math import pi, tan, sin, cos
from ._concrete_shear import epsilon_x, v_rd_max, eta_fc


def v_ed_ti(t_ed: float, a_k: float, z_i: float):
    """Shear force due to torsion

    fib Model Code 2010, eq. (7.3-53)

    Args:
        t_ed: The acting torsion force in the cross section in Nmm
        z_i: Can be found in figure 7.3-18
        a_k: Can be found in figure 7.3-18

    Returns:
        The shear force that will ocurre due to torsion force."""
    return t_ed * z_i / (2 * a_k)


def t_rd_max(
    f_ck: float,
    d_k: float,
    a_k: float,
    theta: float,
    approx_lvl: int,
    z: float,
    E_s: float,
    As: float,
    loads: dict,
    gamma_c: float = 1.5,
) -> float:
    """The maximum allowed torsion allowed

    fib Model Code 2010, eq. (7.3-56)

    Args:
        f_ck (float): Characteristic strength in MPa
        gamma_c (float): Concrete safety factor
        d_k (float): Is the diameter in the smalest circel in the cross section
        a_k: Can be found in figure 7.3-18
        theta (float): Inclitaniton of the compression stressfield in degrees
        approx_lvl (int): Approximation method for cocrete with reinforcement
        E_s (float): The E_s-modulus to the materialb in MPa
        As (float): The cross-section reinforcement in mm^2
        Med (float): The moment working on the material in Nmm
        Ved (float): The shear working on the material in N
        Ned: (float): The normal force working on the material in N with
        positive sign for tension and negative sign for compression
        z (float): distances between the centerline of the
        compressive chord and the reinforcement in mm
        delta_E (float): The eccentricity of the axial load due to
        imperfection in the construction with distance in mm as a positive
        value

    Return:
        The maximum allowed torsion allowed
    """
    t_ef = d_k / 8

    if approx_lvl == 3:
        k_epsilon = 0.55
    elif approx_lvl == 4:
        epsilon_1 = epsilon_x(E_s, As, z, loads) + (
            epsilon_x(E_s, As, z, loads) + 0.002
        ) * ((1 / tan(theta * pi / 180)) ** 2)
        k_epsilon = min(1 / (1.2 + 55 * epsilon_1), 0.65)
    elif approx_lvl == 5:
        theta_min = 20 + 10000 * epsilon_x(E_s, As, z, loads)
        epsilon_1 = epsilon_x(E_s, As, z, loads) + (
            epsilon_x(E_s, As, z, loads) + 0.002
        ) * ((1 / tan(theta_min * pi / 180)) ** 2)
        k_epsilon = min(1 / (1.2 + 55 * epsilon_1), 0.65)
    k_c = eta_fc(f_ck) * k_epsilon
    result = (
        k_c
        * f_ck
        * t_ef
        * 2
        * a_k
        * sin(theta * pi / 180)
        * cos(theta * pi / 180)
        / gamma_c
    )

    return result


def t_rd(  # pylint: disable=r0801
    t_ed: float,
    approx_lvl: int,
    fck: float,
    bw: float,
    theta: float,
    z: float,
    E_s: float,
    As: float,
    loads: dict,
    d_k: float,
    a_k: float,
    alpha: float = 90.0,
    gamma_c: float = 1.5,
) -> bool:
    """Checks if the combination of torstion ans shear is ok

    fib Model Code 2010, eq. (7.3-56)

    Args:
        t_ed (float): The torsion working on the material in Nmm
        approx_lvl (int): Approximation level choosen for shear resistance
        fck (float): Characteristic strength in MPa
        bw: (float): Thickness of web in cross section
        z: (float): The length to the areasenter of cross-section in mm
        E_s: (float): The E_s-modulus to the materialb in MPa
        As: (float): The cross-section area in mm^2
        Med: (float): The moment working on the material in Nmm
        Ved: (float): The shear working on the material in N
        Ned: (float): The normal force working on the material in N with
        positive sign for tension and negative sign for compression
        delta_E (float): The eccentricity of the axial load due to
        imperfection in the construction with distance in mm as a positive
        value
        alpha (float): Inclination of the stirrups in degrees
        d_k (float): Is the diameter in the smalest circel in the cross section
        a_k (float): Can be found in figure 7.3-18
        gamma_c (float): Safety factor

    Return:
        Returns a bool that is true if the criteria for torsion and
        shear is fulfilled"""
    check = bool(
        (
            t_ed
            / t_rd_max(
                fck,
                d_k,
                a_k,
                theta,
                approx_lvl,
                z,
                E_s,
                As,
                loads,
                gamma_c,
            )
        )
        ** 2
        + (
            loads.get('Ved')
            / v_rd_max(
                approx_lvl,
                fck,
                bw,
                theta,
                z,
                E_s,
                As,
                loads,
                alpha,
                gamma_c,
            )
        )
        ** 2
        <= 1
    )
    return check
