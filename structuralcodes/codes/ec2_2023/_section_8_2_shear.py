import math


def tao_Ed(VEd: float, bw: float, d: float) -> float:
    """Calculate the average shear stress over the cross-section
        for linear members.

    EN1992-1-1:2923 Eq. (8.18)

    Parameters:
        VEd (float): Design shear force at the control section in
            linear members in kN.
        bw (float): Width of the cross-section of linear members in mm.
        d (float): Effective depth of the cross-section in mm.

    Returns:
        float: Average shear stress over the cross-section
            for linear members in Mpa.
    """
    if bw < 0:
        raise ValueError(f'bw must not be negative. Got {bw}')
    if d < 0:
        raise ValueError(f'd must not be negative. Got {d}')

    z = 0.9 * d
    return VEd * 1000 / (bw * z)


def tao_Ed_planar(vEd: float, d: float) -> float:
    """Calculate the average shear stress over the cross-section
        for planar members.

    EN1992-1-1:2923 Eq. (8.19)

    Parameters:
        vEd (float): Design shear force per unit width in planar members
            in kN/m.
        d (float): Effective depth of the cross-section in mm.

    Returns:
        float: Average shear stress over the cross-section for
            planar members in MPa.
    """
    if d < 0:
        raise ValueError(f'd must not be negative. Got {d}')

    z = 0.9 * d
    return vEd / z


def tau_rdc_min(
    gamma_v: float,
    f_ck: float,
    f_yd: float,
    d: float,
    d_lower: float,
) -> float:
    """Calculate the minimum shear stress resistance.

    EN1992-1-1:2023 Eq. (8.20)

    Args:
        gamma_v (float): Partial factor for shear design.
        f_ck (float): Characteristic compressive strength of
            concrete in MPa (must be positive).
        f_yd (float): Design value of the yield strength
            in MPa (must be positive).
        d (float): Effective depth of the flexural
            reinforcement in mm (must be positive).
        d_lower (float): Smallest value of the upper sieve size D in an
            aggregate for the coarsest fraction of aggregates
            in mm (must be positive).

    Returns:
        float: Minimum shear stress resistance in MPa.

    Raises:
        ValueError: If any input value is non-positive.
    """
    if gamma_v < 0:
        raise ValueError(f'gamma_v must not be negative. Got {gamma_v}')
    if f_ck < 0:
        raise ValueError(f'f_ck must not be negative. Got {f_ck}')
    if f_yd < 0:
        raise ValueError(f'f_yd must not be negative. Got {f_yd}')
    if d < 0:
        raise ValueError(f'd must not be negative. Got {d}')
    if d_lower < 0:
        raise ValueError(f'd_lower must not be negative. Got {d_lower}')

    if f_ck <= 60:
        d_dg = min(16 + d_lower, 40)
    else:
        d_dg = min(16 + d_lower * (60 / f_ck) ** 2, 40)

    return 11 / gamma_v * math.sqrt(f_ck / f_yd * d_dg / d)


def v_Ed(vEd_x: float, vEd_y: float) -> float:
    """Calculate the design shear force per unit width (vEd).

    EN1992-1-1:2023 Eq. (8.21)

    Args:
        vEd_x (float): Shear force in x-direction in kN/m.
        vEd_y (float): Shear force in y-direction in kN/m.

    Returns:
        float: Design shear force per unit width (vEd) kN/m.

    Raises:
        ValueError: If vEd_x or vEd_y is negative.
    """
    if vEd_x < 0:
        raise ValueError(f'vEd_x must not be negative. Got {vEd_x}')
    if vEd_y < 0:
        raise ValueError(f'vEd_y must not be negative. Got {vEd_y}')

    return math.sqrt(vEd_x**2 + vEd_y**2)


def d_eff(dx: float, dy: float, vEd_x: float, vEd_y: float) -> float:
    """Calculate the effective depth (d) based on the ratio of shear forces.

    EN1992-1-1:2023 Eq. (8.22), (8.23), (8.24)

    Args:
        dx (float): Effective depth in x-direction in mm.
        dy (float): Effective depth in y-direction in mm.
        vEd_x (float): Shear force in x-direction in kN/m.
        vEd_y (float): Shear force in y-direction in kN/m.

    Returns:
        float: Effective depth (d) in mm.

    Raises:
        ValueError: If dx, dy, vEd_x, or vEd_y is negative.
    """
    if dx < 0:
        raise ValueError(f'dx must not be negative. Got {dx}')
    if dy < 0:
        raise ValueError(f'dy must not be negative. Got {dy}')

    vEd_x = abs(vEd_x)
    vEd_y = abs(vEd_y)

    ratio = vEd_y / vEd_x if vEd_x != 0 else float('inf')

    if ratio <= 0.5:
        return dx
    if 0.5 < ratio < 2:
        return 0.5 * (dx + dy)
    return dy


def d_eff_with_angle(
    dx: float, dy: float, vEd_x: float, vEd_y: float
) -> float:
    """Calculate the effective depth (d) based on the angle αv.

    EN1992-1-1:2023 Eq. (8.25), (8.26)

    Args:
        dx (float): Effective depth in x-direction in mm.
        dy (float): Effective depth in y-direction in mm.
        vEd_x (float): Shear force in x-direction in kN/m.
        vEd_y (float): Shear force in y-direction kN/m.

    Returns:
        float: Effective depth (d) mm.

    Raises:
        ValueError: If dx, dy, vEd_x, or vEd_y is negative.
    """
    if dx < 0:
        raise ValueError(f'dx must not be negative. Got {dx}')
    if dy < 0:
        raise ValueError(f'dy must not be negative. Got {dy}')

    vEd_x = abs(vEd_x)
    vEd_y = abs(vEd_y)

    alpha_v = math.atan2(vEd_y, vEd_x)

    return dx * math.cos(alpha_v) ** 2 + dy * math.sin(alpha_v) ** 2


def tau_Rdc(
    gamma_v: float,
    rho_l: float,
    f_ck: float,
    d: float,
    d_g: float,
    tau_rdc_min: float,
) -> float:
    """Calculate the design value of the shear stress resistance.

    EN1992-1-1:2023 Eq. (8.27)

    Args:
        gamma_v (float): Partial factor for shear (unitless).
        rho_l (float): Reinforcement ratio (unitless).
        f_ck (float): Characteristic compressive strength of concrete in MPa.
        d (float): Effective depth in mm.
        d_g (float): Maximum aggregate size in mm.
        tau_rdc_min (float): Minimum resistance neede in MPa.

    Returns:
        float: The design value of the shear stress resistance MPa.

    Raises:
        ValueError: If any of the input values are negative.
    """
    if gamma_v <= 0 or rho_l <= 0 or f_ck <= 0 or d <= 0 or d_g <= 0:
        raise ValueError(
            f'All input values must be positive. Got gamma_v={gamma_v}, '
            + f'rho_l={rho_l}, f_ck={f_ck}, d={d}, d_g={d_g}'
        )

    return max(
        0.66 / gamma_v * (100 * rho_l * f_ck * d_g / d) ** (1 / 3),
        abs(tau_rdc_min),
    )


def rho_l(A_sl: float, b_w: float, d: float) -> float:
    """Calculate the reinforcement ratio.

    EN1992-1-1:2023 Eq. (8.28)

    Args:
        A_sl (float): Effective area of tensile reinforcement in mm2.
        b_w (float): Width of the cross-section in mm.
        d (float): Effective depth in mm.

    Returns:
        float: The reinforcement ratio (unitless).

    Raises:
        ValueError: If any of the input values are negative.
    """
    if A_sl <= 0 or b_w <= 0 or d <= 0:
        raise ValueError(
            'All input values must be positive.'
            + f'Got A_sl={A_sl}, b_w={b_w}, d={d}'
        )

    return A_sl / (b_w * d)


def a_v(a_cs: float, d: float) -> float:
    """Calculate the effective shear span.

    EN1992-1-1:2023 Eq. (8.30)

    Args:
        M_Ed (float): Bending moment kN·m.
        V_Ed (float): Shear force kN.
        d (float): Effective depth mm.

    Returns:
        float: The effective shear span in mm.

    Raises:
        ValueError: If any of the input values are negative.
    """
    if a_cs <= 0 or d <= 0:
        raise ValueError(
            'All input values must be positive.' + f'Got a_cs={a_cs}, d={d}'
        )
    return math.sqrt(a_cs / 4 * d)


def a_cs(M_Ed: float, V_Ed: float, d: float) -> float:
    """Calculate the effective shear span.

    EN1992-1-1:2023 Eq. (8.30)

    Args:
        M_Ed (float): Bending moment kN·m.
        V_Ed (float): Shear force kN.
        d (float): Effective depth mm.

    Returns:
        float: The effective shear span in mm.

    Raises:
        ValueError: If any of the input values are negative.
    """
    if d <= 0:
        raise ValueError('All input values must be positive.' + f'Got d={d}')

    return max(abs(M_Ed * 1000 / V_Ed), d)


def k_vp(N_Ed: float, V_Ed: float, d: float, a_cs: float) -> float:
    """Calculate the coefficient k_vp.

    EN1992-1-1:2023 Eq. (8.31)

    Args:
        N_Ed (float): Axial force in kN.
        V_Ed (float): Shear force in kN.
        d (float): Effective depth in mm.
        a_cs (float): Effective shear span in mm.

    Returns:
        float: The coefficient k_vp (unitless).

    Raises:
        ValueError: If any of the input values are negative.

    """
    if d <= 0 or a_cs <= 0:
        raise ValueError(
            'All input values must be positive. '
            + f'Got N_Ed={N_Ed}, V_Ed={V_Ed}, d={d}, a_cs={a_cs}'
        )

    k_vp = 1 + (N_Ed / abs(V_Ed)) * (d / (3 * a_cs))
    return max(k_vp, 0.1)


def tau_Rdc_0(
    gamma_v: float, rho_l: float, f_ck: float, d: float, d_g: float
) -> float:
    """Calculate the design value of the shear stress
        resistance without axial force effects.

    EN1992-1-1:2023 Eq. (8.33)

    Args:
        gamma_v (float): Partial factor for shear (unitless).
        rho_l (float): Reinforcement ratio (unitless).
        f_ck (float): Characteristic compressive strength of concrete in MPa.
        d (float): Effective depth in mm.
        d_g (float): Maximum aggregate size in mm.

    Returns:
        float: The design value of the shear stress resistance in MPa.

    Raises:
        ValueError: If any of the input values are negative.
    """
    if gamma_v <= 0 or rho_l <= 0 or f_ck <= 0 or d <= 0 or d_g <= 0:
        raise ValueError(
            'All input values must be positive. '
            + f'Got gamma_v={gamma_v}, rho_l={rho_l}, '
            + f'f_ck={f_ck}, d={d}, d_g={d_g}'
        )

    return 0.66 / gamma_v * (100 * rho_l * f_ck * d_g / d) ** (1 / 3)


def tau_Rdc_comp(
    tau_Rdc_0: float, k1: float, sigma_cp: float, tau_Rdc_max: float
) -> float:
    """Calculate the design value of the shear stress
        resistance considering compressive normal forces.

    EN1992-1-1:2023 Eq. (8.32)

    Args:
        tau_Rdc_0 (float): Design value of the shear stress resistance
            without axial force effects in MPa.
        k1 (float): Factor considering the effect of compressive
            normal forces (unitless).
        sigma_cp (float): Compressive stress due to axial force in MPa.
        tau_Rdc_max (float): Maximum design value of the
            shear stress resistance in MPa.

    Returns:
        float: The design value of the shear stress resistance in MPa.

    Raises:
        ValueError: If any of the input values are negative.
    """
    if tau_Rdc_0 <= 0 or k1 <= 0 or sigma_cp < 0 or tau_Rdc_max <= 0:
        raise ValueError(
            'All input values must be positive. '
            + f'Got tau_Rdc_0={tau_Rdc_0}, k1={k1}, '
            + f'sigma_cp={sigma_cp}, tau_Rdc_max={tau_Rdc_max}'
        )

    tau_Rdc = tau_Rdc_0 - k1 * sigma_cp
    return min(tau_Rdc, tau_Rdc_max)


def k1(
    a_cs_0: float,
    e_p: float,
    A_c: float,
    b_w: float,
    z: float,
    d: float,
) -> float:
    """Calculate the factor k1 considering the effect of
        compressive normal forces.

    EN1992-1-1:2023 Eq. (8.34)

    Args:
        a_cs_0 (float): Effective shear span without
            considering prestressing effects in mm.
        e_p (float): Eccentricity of the prestressing
            force or external load in mm.
        A_c (float): Area of concrete cross-section in mm2.
        b_w (float): Width of the cross-section in mm.
        z (float): Lever arm in mm.
        d (float): Effective depth in mm.

    Returns:
        float: The factor k1 (unitless).

    Raises:
        ValueError: If any of the input values are negative.
    """
    if a_cs_0 <= 0 or e_p < 0 or A_c <= 0 or b_w <= 0 or z <= 0:
        raise ValueError(
            'All input values must be positive.'
            + f' Got a_cs_0={a_cs_0}, e_p={e_p}, A_c={A_c}, b_w={b_w}, z={z}'
        )

    k1 = 0.5 * a_cs_0 / (e_p + d / 3) * (A_c / (b_w * z))
    return min(k1, A_c * 0.18 / (b_w * z))


def tau_Rdc_max(tau_Rdc_0: float, a_cs_0: float, d: float) -> float:
    """Calculate the maximum design value of the shear stress resistance.

    EN1992-1-1:2023 Eq. (8.35)

    Args:
        tau_Rdc_0 (float): Design value of the shear stress
            resistance without axial force effects in MPa.
        a_cs_0 (float): Effective shear span without
            considering prestressing effects in mm.
        d (float): Effective depth in mm.

    Returns:
        float: The maximum design value of the shear stress resistance in MPa.

    Raises:
        ValueError: If any of the input values are negative.
    """
    if tau_Rdc_0 <= 0 or a_cs_0 <= 0 or d <= 0:
        raise ValueError(
            'All input values must be positive. '
            + f'Got tau_Rdc_0={tau_Rdc_0}, a_cs_0={a_cs_0}, d={d}'
        )

    tau_Rdc_max = 2.15 * tau_Rdc_0 * (a_cs_0 / d) ** (1 / 6)
    return min(tau_Rdc_max, 2.7 * tau_Rdc_0)


def d_eff_p(ds: float, As: float, dp: float, Ap: float) -> float:
    """Calculate the effective depth for prestressed
        members with bonded tendons.

    EN1992-1-1:2023 Eq. (8.36)
    Parameters:
        ds (float): Depth of the tension reinforcement in mm
        As (float): Area of the tension reinforcement in mm2
        dp (float): Depth of the prestressed reinforcement in mm
        Ap (float): Area of the prestressed reinforcement in mm2

    Returns:
        float: Effective depth in mm

    Raises:
        ValueError: If any of the input values are negative
    """
    if ds < 0:
        raise ValueError(f'ds must not be negative. Got {ds}')
    if As < 0:
        raise ValueError(f'As must not be negative. Got {As}')
    if dp < 0:
        raise ValueError(f'dp must not be negative. Got {dp}')
    if Ap < 0:
        raise ValueError(f'Ap must not be negative. Got {Ap}')

    return (ds**2 * As + dp**2 * Ap) / (ds * As + dp * Ap)


def rho_l_p(
    ds: float, As: float, dp: float, Ap: float, bw: float, d: float
) -> float:
    """Calculate the reinforcement ratio for prestressed
        members with bonded tendons.

    EN1992-1-1:2023 Eq. (8.37)
    Parameters:
        ds (float): Depth of the tension reinforcement in mm
        As (float): Area of the tension reinforcement in mm2
        dp (float): Depth of the prestressed reinforcement in mm
        Ap (float): Area of the prestressed reinforcement in mm2
        bw (float): Width of the member in mm
        d (float): Effective depth in mm

    Returns:
    float: Reinforcement ratio

    Raises:
        ValueError: If any of the input values are negative
    """
    if ds < 0:
        raise ValueError(f'ds must not be negative. Got {ds}')
    if As < 0:
        raise ValueError(f'As must not be negative. Got {As}')
    if dp < 0:
        raise ValueError(f'dp must not be negative. Got {dp}')
    if Ap < 0:
        raise ValueError(f'Ap must not be negative. Got {Ap}')
    if bw < 0:
        raise ValueError(f'bw must not be negative. Got {bw}')
    if d < 0:
        raise ValueError(f'd must not be negative. Got {d}')

    return (ds * As + dp * Ap) / (bw * d**2)


def rho_l_planar(
    vEd_y: float, vEd_x: float, rho_l_x: float, rho_l_y: float
) -> float:
    """Calculate the reinforcement ratio for planar members with
        different reinforcement ratios in both directions.

    EN1992-1-1:2023 Eq. (8.38), (8.39), (8.40)

    Parameters:
        vEd_y (float): Shear force in y-direction (kN)
        vEd_x (float): Shear force in x-direction (kN)
        rho_l_x (float): Reinforcement ratio in x-direction
        rho_l_y (float): Reinforcement ratio in y-direction

    Returns:
        float: Reinforcement ratio

    Raises:
    ValueError: If any of the input values are negative or if the
        shear force ratio is not in the valid range
    """
    if vEd_y < 0:
        raise ValueError(f'vEd_y must not be negative. Got {vEd_y}')
    if vEd_x < 0:
        raise ValueError(f'vEd_x must not be negative. Got {vEd_x}')
    if rho_l_x < 0:
        raise ValueError(f'rho_l_x must not be negative. Got {rho_l_x}')
    if rho_l_y < 0:
        raise ValueError(f'rho_l_y must not be negative. Got {rho_l_y}')

    ratio = vEd_y / vEd_x

    if ratio <= 0.5:
        rho_l = rho_l_x
    elif ratio >= 2:
        rho_l = rho_l_y
    else:
        alpha_v = math.atan(vEd_y / vEd_x)
        rho_l = (
            rho_l_x * math.cos(alpha_v) ** 4 + rho_l_y * math.sin(alpha_v) ** 4
        )

    return rho_l
