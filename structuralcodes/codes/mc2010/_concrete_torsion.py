"""Covers torsion in Model code 2010, 7.3.4."""

from math import cos, pi, sin, tan

from ._concrete_shear import epsilon_x, eta_fc, v_rd_max


def v_ed_ti(t_ed: float, a_k: float, z_i: float):
    """Shear force due to torsion.

    fib Model Code 2010, eq. (7.3-53).

    Args:
        t_ed: The acting torsion force in the cross section in Nmm.
        z_i: Can be found in figure 7.3-18.
        a_k: Can be found in figure 7.3-18.

    Returns:
        float: The shear force that will occur due to torsion moment.
    """
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
    """The maximum allowed torsion.

    fib Model Code 2010, eq. (7.3-56).

    Args:
        f_ck (float): Characteristic strength in MPa.
        gamma_c (float): Concrete safety factor.
        d_k (float): Is the diameter in the smallest circel in the cross
            section.
        a_k: Can be found in figure 7.3-18.
        theta (float): inclination of the compression stressfield in degrees.
        approx_lvl (int): Approximation method for concrete with reinforcement.
        z (float): distances between the centerline of the compressive chord
            and the reinforcement in mm.
        E_s (float): The E_modulus to the material in MPa.
        As (float): The cross-section reinforcement in mm^2.
        loads (dict): The given loads in a dictionary. See create_load_dict.
        gamma_c (float): Concrete safety factor.

    Returns:
        float: The maximum allowed torsion.
    """
    t_ef = d_k / 8

    if approx_lvl == 1:
        k_epsilon = 0.55

    elif approx_lvl == 2:
        epsilonx = epsilon_x(E_s, As, z, loads)
        epsilon_1 = epsilonx + (epsilonx + 0.002) * (
            (1 / tan(theta * pi / 180)) ** 2
        )
        k_epsilon = min(1 / (1.2 + 55 * epsilon_1), 0.65)

    elif approx_lvl == 3:
        epsilonx = epsilon_x(E_s, As, z, loads)
        theta_min = 20 + 10000 * epsilonx
        epsilon_1 = epsilonx + (epsilonx + 0.002) * (
            (1 / tan(theta_min * pi / 180)) ** 2
        )
        k_epsilon = min(1 / (1.2 + 55 * epsilon_1), 0.65)

    k_c = eta_fc(f_ck) * k_epsilon
    return (
        k_c
        * f_ck
        * t_ef
        * 2
        * a_k
        * sin(theta * pi / 180)
        * cos(theta * pi / 180)
        / gamma_c
    )


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
    """Checks if the combination of torsion and shear is ok.

    fib Model Code 2010, eq. (7.3-56).

    Args:
        t_ed (float): The torsion working on the material in Nmm.
        approx_lvl (int): Approximation level chosen for shear resistance.
        fck (float): Characteristic strength in MPa.
        bw: (float): Thickness of web in cross section.
        theta (float): inclination of the compression stressfield in degrees.
        z: (float): The length to the areasenter of cross-section in mm.
        E_s: (float): The E_modulus to the material in MPa.
        As: (float): The cross-section area in mm^2.
        loads (dict): The given loads in a dictionary. See create_load_dict.
        d_k (float): Is the diameter in the smallest circel in the cross
            section.
        a_k (float): Can be found in figure 7.3-18.
        alpha (float): Inclination of the stirrups in degrees.
        gamma_c (float): Concrete safety factor.

    Returns:
        bool: Returns a bool that is true if the criterion for torsion and
        shear is fulfilled.
    """
    return bool(
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
