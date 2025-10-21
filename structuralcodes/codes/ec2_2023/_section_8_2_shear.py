"""Functions from Section 8.2 of EN 1992-1-1:2023."""

import math
from typing import List, Literal


def tau_Ed(VEd: float, bw: float, d: float) -> float:
    """Calculate the average shear stress over the cross-section for linear
    members.

    EN1992-1-1:2923 Eq. (8.18).

    Args:
        VEd (float): Design shear force at the control section in linear
            members in kN.
        bw (float): Width of the cross-section of linear members in mm.
        d (float): Effective depth of the cross-section in mm.

    Returns:
        float: Average shear stress over the cross-section for linear members
        in Mpa.
    """
    if bw < 0:
        raise ValueError(f'bw must not be negative. Got {bw}')
    if d < 0:
        raise ValueError(f'd must not be negative. Got {d}')

    z = 0.9 * d
    return VEd * 1000 / (bw * z)


def tau_Ed_planar(vEd: float, d: float) -> float:
    """Calculate the average shear stress over the cross-section for planar
    members.

    EN1992-1-1:2923 Eq. (8.19).

    Args:
        vEd (float): Design shear force per unit width in planar members in
            kN/m.
        d (float): Effective depth of the cross-section in mm.

    Returns:
        float: Average shear stress over the cross-section for planar members
        in MPa.
    """
    if d < 0:
        raise ValueError(f'd must not be negative. Got {d}')

    z = 0.9 * d
    return vEd / z


def d_dg(f_ck: float, d_lower: float) -> float:
    """Calculate the size parameter describing the failure zone roughness.

    EN1992-1-1:2023 Note 2 for Eq. (8.20).

    Args:
        f_ck (float): Characteristic compressive strength of concrete in MPa
            (must be positive).
        d_lower (float): Smallest value of the upper sieve size D in an
            aggregate for the coarsest fraction of aggregates in mm (must be
            positive).

    Returns:
        float: Size parameter ddg in mm.

    Raises:
        ValueError: If any input value is non-positive.
    """
    if f_ck < 0:
        raise ValueError(f'f_ck must not be negative. Got {f_ck}')
    if d_lower < 0:
        raise ValueError(f'd_lower must not be negative. Got {d_lower}')

    if f_ck <= 60:
        return min(16 + d_lower, 40)

    return min(16 + d_lower * (60 / f_ck) ** 2, 40)


def tau_rdc_min(
    gamma_v: float,
    f_ck: float,
    f_yd: float,
    d: float,
    d_dg: float,
) -> float:
    """Calculate the minimum shear stress resistance.

    EN1992-1-1:2023 Eq. (8.20).

    Args:
        gamma_v (float): Partial factor for shear design.
        f_ck (float): Characteristic compressive strength of concrete in MPa
            (must be positive).
        f_yd (float): Design value of the yield strength in MPa (must be
            positive).
        d (float): Effective depth of the flexural reinforcement in mm (must be
            positive).
        d_dg (float): Size parameter describing the failure zone roughness in
            mm (must be positive).

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
    if d_dg < 0:
        raise ValueError(f'd_dg must not be negative. Got {d_dg}')

    return 11 / gamma_v * math.sqrt(f_ck / f_yd * d_dg / d)


def v_Ed(vEd_x: float, vEd_y: float) -> float:
    """Calculate the design shear force per unit width (vEd).

    EN1992-1-1:2023 Eq. (8.21).

    Args:
        vEd_x (float): Shear force in x-direction in kN/m.
        vEd_y (float): Shear force in y-direction in kN/m.

    Returns:
        float: Design shear force per unit width (vEd) kN/m.
    """
    return math.sqrt(vEd_x**2 + vEd_y**2)


def d_eff(dx: float, dy: float, vEd_x: float, vEd_y: float) -> float:
    """Calculate the effective depth (d) based on the ratio of shear forces.

    EN1992-1-1:2023 Eq. (8.22), (8.23), (8.24).

    Args:
        dx (float): Effective depth in x-direction in mm (must be positive).
        dy (float): Effective depth in y-direction in mm (must be positive).
        vEd_x (float): Shear force in x-direction in kN/m.
        vEd_y (float): Shear force in y-direction in kN/m.

    Returns:
        float: Effective depth (d) in mm.

    Raises:
        ValueError: If dx or dy is negative.
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


def d_eff_angle(dx: float, dy: float, vEd_x: float, vEd_y: float) -> float:
    """Calculate the effective depth (d) based on the angle alpha_v.

    EN1992-1-1:2023 Eq. (8.25), (8.26).

    Args:
        dx (float): Effective depth in x-direction in mm.
        dy (float): Effective depth in y-direction in mm.
        vEd_x (float): Shear force in x-direction in kN/m.
        vEd_y (float): Shear force in y-direction kN/m.

    Returns:
        float: Effective depth (d) mm.

    Raises:
        ValueError: If dx or dy is negative.
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
    d_dg: float,
    tau_rdc_min: float,
) -> float:
    """Calculate the design value of the shear stress resistance.

    EN1992-1-1:2023 Eq. (8.27).

    Args:
        gamma_v (float): Partial factor for shear (unitless).
        rho_l (float): Reinforcement ratio (unitless).
        f_ck (float): Characteristic compressive strength of concrete in MPa.
        d (float): Effective depth in mm.
        d_dg (float): Size parameter describing the failure zone roughness in
            mm.
        tau_rdc_min (float): Minimum resistance neede in MPa.

    Returns:
        float: The design value of the shear stress resistance MPa.

    Raises:
        ValueError: If any input values are negative
            or if gamma_v or d are non-positive.
    """
    if gamma_v <= 0 or rho_l < 0 or f_ck < 0 or d <= 0 or d_dg < 0:
        raise ValueError(
            'gamma_v and d must be positive other values must be non-negative.'
            + f'Got gamma_v={gamma_v}, rho_l={rho_l}, '
            + f' f_ck={f_ck}, d={d}, d_dg={d_dg}'
        )

    return max(
        0.66 / gamma_v * (100 * rho_l * f_ck * d_dg / d) ** (1 / 3),
        abs(tau_rdc_min),
    )


def rho_l(A_sl: float, b_w: float, d: float) -> float:
    """Calculate the reinforcement ratio.

    EN1992-1-1:2023 Eq. (8.28).

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
    """Calculate the mechanical shear span.

    EN1992-1-1:2023 Eq. (8.29).

    Args:
        a_cs (float): Effective shear span in mm.
        d (float): Effective depth in mm.

    Returns:
        float: The mechanical shear span av in mm.

    Raises:
        ValueError: If any of the input values are non-positive.
    """
    if a_cs <= 0 or d <= 0:
        raise ValueError(
            'All input values must be positive.' + f'Got a_cs={a_cs}, d={d}'
        )
    return math.sqrt(a_cs / 4 * d)


def a_cs(M_Ed: float, V_Ed: float, d: float) -> float:
    """Calculate the effective shear span.

    EN1992-1-1:2023 Eq. (8.30).

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

    EN1992-1-1:2023 Eq. (8.31).

    Args:
        N_Ed (float): Axial force in kN
            (compression negative, tension positive).
        V_Ed (float): Shear force in kN.
        d (float): Effective depth in mm.
        a_cs (float): Effective shear span in mm.

    Returns:
        float: The coefficient k_vp (unitless).

    Raises:
        ValueError: If d or a_cs are non-positive.

    """
    if d <= 0 or a_cs <= 0:
        raise ValueError(
            'd and a_cs must be positive. ' + f'Got d={d}, a_cs={a_cs}'
        )

    k_vp = 1 + (N_Ed / abs(V_Ed)) * (d / (3 * a_cs))
    return max(k_vp, 0.1)


def tau_Rdc_0(
    gamma_v: float, rho_l: float, f_ck: float, d: float, d_dg: float
) -> float:
    """Calculate the design value of the shear stress resistance without axial
    force effects.

    EN1992-1-1:2023 Eq. (8.33).

    Args:
        gamma_v (float): Partial factor for shear (unitless).
        rho_l (float): Reinforcement ratio (unitless).
        f_ck (float): Characteristic compressive strength of concrete in MPa.
        d (float): Effective depth in mm.
        d_dg (float): Size parameter describing the failure zone roughness in
            mm.

    Returns:
        float: The design value of the shear stress resistance in MPa.

    Raises:
        ValueError: If any of the input values are negative.
    """
    if gamma_v <= 0 or rho_l <= 0 or f_ck <= 0 or d <= 0 or d_dg <= 0:
        raise ValueError(
            'All input values must be positive. '
            + f'Got gamma_v={gamma_v}, rho_l={rho_l}, '
            + f'f_ck={f_ck}, d={d}, d_dg={d_dg}'
        )

    return 0.66 / gamma_v * (100 * rho_l * f_ck * d_dg / d) ** (1 / 3)


def tau_Rdc_comp(
    tau_Rdc_0: float,
    k1: float,
    sigma_cp: float,
    tau_Rdc_max: float,
    tau_rdc_min: float,
) -> float:
    """Calculate the design value of the shear stress resistance considering
    compressive normal forces.

    EN1992-1-1:2023 Eq. (8.32).

    Args:
        tau_Rdc_0 (float): Design value of the shear stress resistance without
            axial force effects in MPa.
        k1 (float): Factor considering the effect of compressive normal forces
            (unitless).
        sigma_cp (float): Compressive stress due to axial force in MPa.
        tau_Rdc_max (float): Maximum design value of the shear stress
            resistance in MPa.
        tau_rdc_min (float): Minimum design value of the shear stress
            resistance in MPa.

    Returns:
        float: The design value of the shear stress resistance in MPa.

    Raises:
        ValueError: If any of the input values are negative.
    """
    if (
        tau_Rdc_0 <= 0
        or k1 <= 0
        or sigma_cp < 0
        or tau_Rdc_max <= 0
        or tau_rdc_min <= 0
    ):
        raise ValueError(
            'All input values must be positive. '
            + f'Got tau_Rdc_0={tau_Rdc_0}, k1={k1}, '
            + f'sigma_cp={sigma_cp}, tau_Rdc_max={tau_Rdc_max}, '
            + f'tau_rdc_min={tau_rdc_min}'
        )

    tau_Rdc = tau_Rdc_0 - k1 * sigma_cp
    return max(min(tau_Rdc, tau_Rdc_max), tau_rdc_min)


def k1(
    a_cs_0: float,
    e_p: float,
    A_c: float,
    b_w: float,
    z: float,
    d: float,
) -> float:
    """Calculate the factor k1 considering the effect of compressive normal
    forces.

    EN1992-1-1:2023 Eq. (8.34).

    Args:
        a_cs_0 (float): Effective shear span without considering prestressing
            effects in mm.
        e_p (float): Eccentricity of the prestressing force or external load in
            mm.
        A_c (float): Area of concrete cross-section in mm2.
        b_w (float): Width of the cross-section in mm.
        z (float): Lever arm in mm.
        d (float): Effective depth in mm.

    Returns:
        float: The factor k1 (unitless).

    Raises:
        ValueError: If any of the input values are negative.
    """
    if a_cs_0 <= 0 or e_p < 0 or A_c <= 0 or b_w <= 0 or z <= 0 or d <= 0:
        raise ValueError(
            'All input values must be positive.'
            + f' Got a_cs_0={a_cs_0}, e_p={e_p}, A_c={A_c}, '
            + f'b_w={b_w}, z={z}, d={d}'
        )

    k1 = 0.5 / a_cs_0 / (e_p + d / 3) * (A_c / (b_w * z))
    return min(k1, A_c * 0.18 / (b_w * z))


def tau_Rdc_max(tau_Rdc_0: float, a_cs_0: float, d: float) -> float:
    """Calculate the maximum design value of the shear stress resistance.

    EN1992-1-1:2023 Eq. (8.35).

    Args:
        tau_Rdc_0 (float): Design value of the shear stress resistance without
            axial force effects in MPa.
        a_cs_0 (float): Effective shear span without considering prestressing
            effects in mm.
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
    """Calculate the effective depth for prestressed members with bonded
    tendons.

    EN1992-1-1:2023 Eq. (8.36).

    Args:
        ds (float): Depth of the tension reinforcement in mm.
        As (float): Area of the tension reinforcement in mm2.
        dp (float): Depth of the prestressed reinforcement in mm.
        Ap (float): Area of the prestressed reinforcement in mm2.

    Returns:
        float: Effective depth in mm.

    Raises:
        ValueError: If any of the input values are negative or if
            ds*As + dp*Ap equals zero (division by zero).
    """
    if ds < 0:
        raise ValueError(f'ds must not be negative. Got {ds}')
    if As < 0:
        raise ValueError(f'As must not be negative. Got {As}')
    if dp < 0:
        raise ValueError(f'dp must not be negative. Got {dp}')
    if Ap < 0:
        raise ValueError(f'Ap must not be negative. Got {Ap}')

    denominator = ds * As + dp * Ap
    if denominator == 0:
        raise ValueError(
            'Division by zero: ds*As + dp*Ap cannot be zero. '
            + f'Got ds={ds}, As={As}, dp={dp}, Ap={Ap}'
        )

    return (ds**2 * As + dp**2 * Ap) / denominator


def rho_l_p(
    ds: float, As: float, dp: float, Ap: float, bw: float, d: float
) -> float:
    """Calculate the reinforcement ratio for prestressed members with bonded
    tendons.

    EN1992-1-1:2023 Eq. (8.37).

    Args:
        ds (float): Depth of the tension reinforcement in mm.
        As (float): Area of the tension reinforcement in mm2.
        dp (float): Depth of the prestressed reinforcement in mm.
        Ap (float): Area of the prestressed reinforcement in mm2.
        bw (float): Width of the member in mm.
        d (float): Effective depth in mm.

    Returns:
        float: Reinforcement ratio.

    Raises:
        ValueError: If any of the input values are negative or if bw is not
            positive.
    """
    if ds < 0:
        raise ValueError(f'ds must not be negative. Got {ds}')
    if As < 0:
        raise ValueError(f'As must not be negative. Got {As}')
    if dp < 0:
        raise ValueError(f'dp must not be negative. Got {dp}')
    if Ap < 0:
        raise ValueError(f'Ap must not be negative. Got {Ap}')
    if bw <= 0:
        raise ValueError(f'bw must be positive. Got {bw}')
    if d < 0:
        raise ValueError(f'd must not be negative. Got {d}')

    return (ds * As + dp * Ap) / (bw * d**2)


def rho_l_planar(
    vEd_y: float, vEd_x: float, rho_l_x: float, rho_l_y: float
) -> float:
    """Calculate the reinforcement ratio for planar members with different
    reinforcement ratios in both directions.

    EN1992-1-1:2023 Eq. (8.38), (8.39), (8.40).

    Args:
        vEd_y (float): Shear force in y-direction (kN).
        vEd_x (float): Shear force in x-direction (kN) (cannot be zero).
        rho_l_x (float): Reinforcement ratio in x-direction.
        rho_l_y (float): Reinforcement ratio in y-direction.

    Returns:
        float: Reinforcement ratio.

    Raises:
        ValueError: If any of the input values are negative or if vEd_x is zero
            (division by zero).
    """
    if vEd_y < 0:
        raise ValueError(f'vEd_y must not be negative. Got {vEd_y}')
    if vEd_x < 0:
        raise ValueError(f'vEd_x must not be negative. Got {vEd_x}')
    if rho_l_x < 0:
        raise ValueError(f'rho_l_x must not be negative. Got {rho_l_x}')
    if rho_l_y < 0:
        raise ValueError(f'rho_l_y must not be negative. Got {rho_l_y}')

    if vEd_x == 0:
        raise ValueError(
            'Division by zero: vEd_x cannot be zero for ratio calculation. '
            + f'Got vEd_x={vEd_x}'
        )

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


def cot_theta_min(
    NEd: float,
    VEd: float,
    x: float,
    d: float,
) -> float:
    """Calculate the minimum cotangent of the compression field inclination
    angle, thetamin, according to the conditions provided.

    EN1992-1-1:2023 Eq. (8.41).

    Args:
        NEd (float): Axial force in the member in kN.
        VEd (float): Shear force in the member in kN.
        x (float): Depth of the compression chord in mm.
        d (float): Effective depth of the member in mm.

    Returns:
        float: Minimum cotangent of the compression field inclination angle.

    Raises:
        ValueError: If any of the dimensions or forces are negative, or if d is
            zero.
    """
    if d <= 0 or x < 0:
        raise ValueError(
            'Dimensions and forces must be positive, and d must not be zero.'
        )

    if NEd > 0:
        cot_theta_min_value = 2.5 - 0.1 * NEd / abs(VEd)
        return max(cot_theta_min_value, 1.0)
    if NEd < 0:
        return 3.0

    return 2.5  # NEd == 0


def tau_Rd_sy(rho_w: float, fywd: float, cot_theta: float) -> float:
    """Calculate the shear stress resistance of yielding shear reinforcement.

    EN1992-1-1:2023 Eq. (8.42).

    Args:
        rho_w (float): Shear reinforcement ratio (unitless).
        fywd (float): Design yield strength of the shear reinforcement in MPa.
        cot_theta (float): Cotangent of the angle of the compression field.

    Returns:
        float: Shear stress resistance in MPa.

    Raises:
        ValueError: If rho_w or fywd is negative.
    """
    if rho_w < 0 or fywd < 0:
        raise ValueError(
            'Shear reinforcement ratio and yield strength must not be negative'
        )

    return rho_w * fywd * cot_theta


def rho_w(Asw: float, bw: float, s: float) -> float:
    """Calculate the shear reinforcement ratio rho_w.

    EN1992-1-1:2023 Eq. (8.43).

    Args:
        Asw (float): Area of shear reinforcement in mm2.
        bw (float): Width of the web in mm.
        s (float): Spacing of the shear reinforcement in mm.

    Returns:
        float: Shear reinforcement ratio, unitless.

    Raises:
        ValueError: If Asw, bw, or s is negative or zero.

    """
    if Asw <= 0 or bw <= 0 or s <= 0:
        raise ValueError('Asw, bw, and s must be positive and non-zero.')

    return Asw / (bw * s)


def sigma_cd(
    tau_Ed: float, cot_theta: float, tan_theta: float, nu: float, f_cd: float
) -> float:
    """Calculate the stress in the compression field sigma_cd and verify it.

    EN1992-1-1:2023 Eq. (8.44).

    Args:
        tau_Ed (float): Design value of the shear stress in MPa.
        cot_theta (float): Cotangent of the angle of the compression field.
        tan_theta (float): Tangent of the angle of the compression field.
        nu (float): Coefficient (usually 0.5 as per the note).
        f_cd (float): Design value of the concrete compressive strength in MPa.

    Returns:
        float: Stress in the compression field sigma_cd in MPa.

    Raises:
        ValueError: If any of the parameters are negative.
    """
    if tau_Ed < 0 or cot_theta < 0 or tan_theta < 0 or nu < 0 or f_cd < 0:
        raise ValueError('All parameters must be positive.')

    sigma_cd_value = tau_Ed * (cot_theta + tan_theta)
    return min(sigma_cd_value, nu * f_cd)


def tau_Rd(
    rho_w: float, fywd: float, cot_theta: float, nu: float, f_cd: float
) -> float:
    """Calculate the shear stress resistance tau_Rd considering the
    simultaneous yielding of the shear reinforcement and failure of the
    compression field.

    EN1992-1-1:2023 Eq. (8.42) and (8.44).

    Args:
        rho_w (float): Shear reinforcement ratio, unitless.
        fywd (float): Design yield strength of the shear reinforcement in MPa.
        cot_theta (float): Cotangent of the angle of the compression field.
        nu (float): Coefficient (usually 0.5 as per the note).
        f_cd (float): Design value of the concrete compressive strength in MPa.

    Returns:
        float: Shear stress resistance tau_Rd in MPa.

    Raises:
        ValueError: If any of the parameters are negative.
    """
    if rho_w < 0 or fywd < 0 or cot_theta < 0 or nu < 0 or f_cd < 0:
        raise ValueError('All parameters must be positive.')

    tau_Rd_value = rho_w * fywd * cot_theta
    return min(tau_Rd_value, nu * f_cd / 2)


def cot_theta(
    nu: float, f_cd: float, rho_w: float, fywd: float, cot_theta_min: float
) -> float:
    """Calculate the cotangent of the angle of the compression field
    considering the simultaneous yielding of the shear reinforcement and
    failure of the compression field.

    EN1992-1-1:2023 Eq. (8.44).

    Args:
        nu (float): Coefficient (usually 0.5 as per the note).
        f_cd (float): Design value of the concrete compressive strength in MPa.
        rho_w (float): Shear reinforcement ratio, unitless.
        fywd (float): Design yield strength of the shear reinforcement in MPa.
        cot_theta_min (float): Value of cot_theta_min.

    Returns:
        float: Cotangent of the angle of the compression field.

    Raises:
        ValueError: If any of the parameters are negative.
    """
    if nu < 0 or f_cd < 0 or rho_w < 0 or fywd < 0:
        raise ValueError('All parameters must be positive.')

    cot_theta_value = (nu * f_cd) / (rho_w * fywd) - 1
    return min(max(cot_theta_value, 1.0), cot_theta_min)


def epsilon_xt(Ftd: float, Est: float, Ast: float) -> float:
    """Calculate epsilon_xt.

    EN1992-1-1:2023 Eq. (8.47).

    Args:
        Ftd (float): Tensile force in the flexural tension chord in kN.
        Est (float): Modulus of elasticity of the steel in the tension chord in
            MPa.
        Ast (float): Area of the longitudinal reinforcement in the flexural
            tension chord in mm2.

    Returns:
        float: epsilon_xt (strain).
    """
    if Est < 0:
        raise ValueError(f'Est must not be negative. Got {Est}')
    if Ast < 0:
        raise ValueError(f'Ast must not be negative. Got {Ast}')
    return abs(Ftd) * 1000 / (Est * Ast)


def epsilon_xc_comp(Fcd: float, Ecc: float, Acc: float) -> float:
    """Calculate epsilon_xc for compression.

    EN1992-1-1:2023 Eq. (8.48).

    Args:
        Fcd (float): Compressive force in the flexural compression chord in kN.
        Ecc (float): Modulus of elasticity of the concrete in the compression
            chord in MPa.
        Acc (float): Area of the flexural compression chord in mm2.

    Returns:
        float: epsilon_xc (strain).
    """
    if Acc < 0:
        raise ValueError(f'Acc must not be negative. Got {Acc}')
    if Ecc < 0:
        raise ValueError(f'Ecc must not be negative. Got {Ecc}')
    return abs(Fcd) * 1000 / (Ecc * Acc)


def epsilon_xc_tens(Fcd: float, Esc: float, Asc: float) -> float:
    """Calculate epsilon_xc.

    EN1992-1-1:2023 Eq. (8.49).

    Args:
        Fcd (float): Tensile force in the flexural compression chord kN.
        Esc (float): Modulus of elasticity of the steel in the compression
            chord in MPa.
        Asc (float): Area of the longitudinal reinforcement in the flexural
            compression chord mm2.

    Returns:
        float: epsilon_xc (strain).
    """
    if Asc < 0:
        raise ValueError(f'Acc must not be negative. Got {Asc}')
    if Esc < 0:
        raise ValueError(f'Esc must not be negative. Got {Esc}')
    return abs(Fcd) * 1000 / (Esc * Asc)


def epsilon_x(epsilon_xt: float, epsilon_xc: float) -> float:
    """Calculate epsilon_x.

    EN1992-1-1:2023 Eq. (8.46).

    Args:
        epsilon_xt (float): Strain in the flexural tension chord.
        epsilon_xc (float): Strain in the flexural compression chord.

    Returns:
        float: epsilon_x (strain).
    """
    epsilon_x_value = (epsilon_xt + epsilon_xc) / 2
    return max(epsilon_x_value, 0)


def nu(epsilon_x: float, cot_theta: float) -> float:
    """Calculate nu.

    EN1992-1-1:2023 Eq. (8.45).

    Args:
        epsilon_x (float): Average strain of the bottom and top chords.
        cot_theta (float): cotan of the compression field inclination to the
            member axis.

    Returns:
        float: nu (dimensionless factor).
    """
    if epsilon_x < 0:
        raise ValueError(f'epsilon_x must not be negative. Got {epsilon_x}')

    nu_value = 1 / (
        1.0 + 110 * (epsilon_x + (epsilon_x + 0.001) * cot_theta**2)
    )
    return min(nu_value, 1.0)


def Nvd(VEd: float, cot_theta: float) -> float:
    """Calculate the additional tensile axial force NVd due to shear VEd.

    EN1992-1-1:2023 Eq. (8.50).

    Args:
        VEd (float): Shear force in kN.
        cot_theta (float): Cotangent of the angle.

    Returns:
        float: Additional tensile axial force NVd in kN.
    """
    return abs(VEd) * cot_theta


def Ftd(
    MEd: float,
    z: float,
    NVd: float,
    NE: float,
) -> float:
    """Calculate the chord force Ftd.

    EN1992-1-1:2023 Eq. (8.51).

    Args:
        MEd (float): Moment in kNm.
        z (float): Lever arm in mm.
        NVd (float): Additional tensile axial force in kN.
        NE (float): Axial force in kN.

    Returns:
        float: Chord force Ftd in kN.

    Raises:
        ValueError: If any input is negative.

    """
    return MEd * 1000 / z + (NVd + NE) / 2


def Fcd(
    MEd: float,
    z: float,
    NVd: float,
    NE: float,
) -> float:
    """Calculate the chord force Fcd.

    EN1992-1-1:2023 Eq. (8.52).

    Args:
        MEd (float): Moment in kNm.
        z (float): Lever arm in mm.
        NVd (float): Additional tensile axial force in kN.
        NE (float): Axial force in kN.

    Returns:
        float: Chord force Fcd in kN.
    """
    return MEd * 1000 / z - (NVd + NE) / 2


def k_duct(
    duct_material: Literal['steel', 'plastic'],
    is_grouted: bool,
    wall_thickness: float,
    duct_diameter: float,
) -> float:
    """Calculate the k_duct coefficient based on duct material, filling, and
    wall thickness.

    EN1992-1-1:2023 guidelines for k_duct.

    Args:
        duct_material (str): Material of the duct ('steel' or 'plastic').
        is_grouted (bool): True if the duct is grouted, False otherwise.
        wall_thickness (float): Wall thickness of the duct in mm.
        duct_diameter (float): Outer diameter of the duct in mm.

    Returns:
        float: Coefficient k_duct.

    Raises:
        ValueError: If wall_thickness or duct_diameter is negative.
        ValueError: If duct_material is not 'steel' or 'plastic'.

    """
    if wall_thickness < 0:
        raise ValueError(
            f'Wall thickness must not be negative. Got {wall_thickness}'
        )
    if duct_diameter < 0:
        raise ValueError(
            f'Duct diameter must not be negative. Got {duct_diameter}'
        )

    max_thickness = max(0.035 * duct_diameter, 2.0)

    if duct_material == 'steel' and is_grouted:
        return 0.5
    if duct_material == 'plastic':
        if is_grouted:
            if wall_thickness <= max_thickness:
                return 0.8
            return 1.2
        return 1.2
    raise ValueError(
        'Invalid duct material. Expected "steel" or "plastic". '
        + f'Got {duct_material}'
    )


def bw_nom(
    bw: float,
    duct_diameters: List[float],
    k_duct: float,
) -> float:
    """Calculate the nominal web width considering the presence of ducts.

    EN1992-1-1:2023 Eq. (8.54).

    Args:
        bw (float): Actual web width in mm.
        duct_diameters (List[float]): List of duct diameters in mm (each must
            be non-negative).
        k_duct (float): Coefficient depending on the material and filling of
            the duct.

    Returns:
        float: Nominal web width in mm.

    Raises:
        ValueError: If bw is negative or if any duct diameter is negative.
        ValueError: If the sum of duct diameters exceeds bw/8.
    """
    if bw < 0:
        raise ValueError(f'bw must not be negative. Got {bw}')
    for d in duct_diameters:
        if d < 0:
            raise ValueError(f'Duct diameters must not be negative. Got {d}')

    sum_phi_duct = sum(duct_diameters)
    if sum_phi_duct > bw / 8:
        raise ValueError(
            'Sum of duct diameters exceeds bw/8.'
            + f'Got sum {sum_phi_duct}, limit {bw/8}'
        )

    return bw - k_duct * sum_phi_duct


def tau_rd(
    nu: float,
    f_cd: float,
    cot_theta: float,
    cot_beta_incl: float,
    rho_w: float,
    f_ywd: float,
) -> float:
    """Calculate the enhanced shear stress resistance tau_Rd.

    EN1992-1-1:2023 Eq. (8.55).

    Args:
        nu (float): The factor nu (unitless).
        f_cd (float): Design value of concrete compressive strength in MPa.
        cot_theta (float): Cotangent of the inclination of compression field.
        cot_beta_incl (float): Cotangent of the inclination of load (cot
            beta_incl).
        rho_w (float): Reinforcement ratio rho_w.
        f_ywd (float): Design yield strength of shear reinforcement in MPa.

    Returns:
        float: Enhanced shear stress resistance tau_Rd in MPa.

    Raises:
        ValueError: If any input parameter is negative where applicable.
    """
    if f_cd < 0:
        raise ValueError(f'f_cd must not be negative. Got {f_cd}')
    if f_ywd < 0:
        raise ValueError(f'f_ywd must not be negative. Got {f_ywd}')

    tau_rd_value = (
        nu * f_cd * (cot_theta - cot_beta_incl) / (1 + cot_theta**2)
        + rho_w * f_ywd * cot_beta_incl
    )
    tau_rd_max = nu * f_cd * cot_theta / (1 + cot_theta**2)

    return min(tau_rd_value, tau_rd_max)


def sigma_swd(
    Es: float, eps_x: float, f_ywd: float, cot_theta: float
) -> float:
    """Calculate the stress sigma_swd in the shear reinforcement.

    EN1992-1-1:2023 Eq. (8.56).

    Args:
        Es (float): Modulus of elasticity of steel Es in MPa.
        eps_x (float): Longitudinal strain epsilon_x.
        f_ywd (float): Design yield strength of shear reinforcement in MPa.
        cot_theta (float): Cotangent of the inclination of compression field.

    Returns:
        float: Stress sigma_swd in the shear reinforcement in MPa.

    Raises:
        ValueError: If any input parameter is negative where applicable.
    """
    if Es < 0:
        raise ValueError(f'e_s must not be negative. Got {Es}')
    if f_ywd < 0:
        raise ValueError(f'f_ywd must not be negative. Got {f_ywd}')

    sigma_swd_value = Es * (cot_theta**2 * (eps_x + 0.001) - 0.001)

    return min(sigma_swd_value, f_ywd)


def delta_MEd(
    tau_ed: float,
    rho_w: float,
    f_ywd: float,
    cot_theta: float,
    z: float,
    b_w: float,
    a: float,
    x: float,
) -> float:
    """Calculate the additional moment delta MEd.

    EN1992-1-1:2023 Eq. (8.57).

    Args:
        tau_ed (float): Shear stress tau_Ed in MPa (should be positive).
        rho_w (float): Reinforcement ratio rho_w.
        f_ywd (float): Design yield strength of shear reinforcement in MPa.
        cot_theta (float): Cotangent of the inclination of compression field.
        z (float): Lever arm z in mm (should be positive).
        b_w (float): Width of the web bw in mm (should be positive).
        a (float): Distance between the axis of the support and the
            concentrated force in mm.
        x (float): Distance between the support and the investigated
            cross-section in mm.

    Returns:
        float: Additional moment delta MEd in kNm.

    Raises:
        ValueError: If any input parameter is negative where applicable.
    """
    if tau_ed < 0:
        raise ValueError(f'tau_ed must not be negative. Got {tau_ed}')
    if f_ywd < 0:
        raise ValueError(f'f_ywd must not be negative. Got {f_ywd}')
    if z < 0:
        raise ValueError(f'z must not be negative. Got {z}')
    if b_w < 0:
        raise ValueError(f'b_w must not be negative. Got {b_w}')
    if a < 0:
        raise ValueError(f'a must not be negative. Got {a}')
    if x < 0:
        raise ValueError(f'x must not be negative. Got {x}')

    return (
        (tau_ed - rho_w * f_ywd * cot_theta) * z * b_w * (a / 2 - x)
    ) / 1e6  # Convert to kNm


def tau_rd_sy(
    rho_w: float,
    f_ywd: float,
    cot_theta: float,
    alpha_w: float,
    cot_theta_min: float,
) -> float:
    """Calculate the shear stress resistance tau_Rd,sy for inclined shear
    reinforcement.

    EN1992-1-1:2023 Eq. (8.58), (8.59).

    Args:
        rho_w (float): Reinforcement ratio rho_w.
        f_ywd (float): Design yield strength of shear reinforcement in MPa.
        cot_theta (float): Cotangent of the inclination of compression field in
            radians (cottheta).
        alpha_w (float): Angle of inclined shear reinforcement alpha_w in
            degrees (45 ≤ alpha_w < 90).
        cot_theta_min (float): max value for cot_theta.

    Returns:
        float: Shear stress resistance tau_Rd,sy in MPa.

    Raises:
        ValueError: If any input parameter is negative where applicable.
    """
    if f_ywd < 0:
        raise ValueError(f'f_ywd must not be negative. Got {f_ywd}')
    if alpha_w < 45 or alpha_w >= 90:
        raise ValueError(
            f'alpha_w must be between 45 degrees and 90 degrees. Got {alpha_w}'
        )

    alpha_w_rad = math.radians(alpha_w)
    cot_theta = min(max(math.tan(alpha_w_rad / 2), cot_theta), cot_theta_min)
    cot_alpha_w = 1 / math.tan(alpha_w_rad)

    return rho_w * f_ywd * (cot_theta + cot_alpha_w) * math.sin(alpha_w_rad)


def sigma_cd_s(
    tau_ed: float,
    cot_theta: float,
    alpha_w: float,
    nu: float,
    f_cd: float,
    cot_theta_min: float,
) -> float:
    """Calculate the compression stress sigma_cd.

    EN1992-1-1:2023 Eq. (8.58), (8.60).

    Args:
        tau_ed (float): Shear stress tau_Ed in MPa.
        cot_theta (float): Cotangent of the inclination of compression field in
            radians (cottheta).
        alpha_w (float): Angle of inclined shear reinforcement alpha_w in
            degrees (45 ≤ alpha_w < 90).
        nu (float): The factor nu.
        f_cd (float): Design value of concrete compressive strength in MPa.
        cot_theta_min (float): max value for cot_theta.

    Returns:
        float: Compression stress sigma_cd in MPa.

    Raises:
        ValueError: If any input parameter is negative where applicable.
    """
    if tau_ed < 0:
        raise ValueError(f'tau_ed must not be negative. Got {tau_ed}')
    if f_cd < 0:
        raise ValueError(f'f_cd must not be negative. Got {f_cd}')
    if alpha_w < 45 or alpha_w >= 90:
        raise ValueError(
            f'alpha_w must be between 45 degrees and 90 degrees. Got {alpha_w}'
        )

    alpha_w_rad = math.radians(alpha_w)
    cot_theta = min(max(math.tan(alpha_w_rad / 2), cot_theta), cot_theta_min)
    cot_alpha_w = 1 / math.tan(alpha_w_rad)

    sigma_cd_value = tau_ed * (1 + cot_theta**2) / (cot_theta + cot_alpha_w)

    return min(sigma_cd_value, nu * f_cd)


def NVds(
    VEd: float, cot_theta: float, alpha_w: float, cot_theta_min: float
) -> float:
    """Calculate the axial tensile force NVd.

    EN1992-1-1:2023 Eq. (8.58), (8.61).

    Args:
        VEd (float): Design shear force VEd in kN (should be positive).
        cot_theta (float): Cotangent of the nclination of compression field in
            radians (cottheta).
        alpha_w (float): Angle of inclined shear reinforcement alpha_w in
            degrees (45 ≤ alpha_w < 90).
        cot_theta_min (float): max value for cot_theta.

    Returns:
        float: Axial tensile force NVd in kN.

    Raises:
        ValueError: If v_ed is negative.
    """
    if VEd < 0:
        raise ValueError(f'VEd must not be negative. Got {VEd}')
    if alpha_w < 45 or alpha_w >= 90:
        raise ValueError(
            f'alpha_w must be between 45 degrees and 90 degrees. Got {alpha_w}'
        )

    alpha_w_rad = math.radians(alpha_w)
    cot_alpha_w = 1 / math.tan(alpha_w_rad)
    cot_theta = min(max(math.tan(alpha_w_rad / 2), cot_theta), cot_theta_min)

    return abs(VEd) * (cot_theta - cot_alpha_w)


def tau_rd_incl(
    nu: float,
    f_cd: float,
    cot_theta: float,
    cot_beta_incl: float,
    rho_w: float,
    f_ywd: float,
    alpha_w: float,
    cot_theta_min: float,
) -> float:
    """Calculate the shear stress resistance tau_Rd for members with inclined
    shear reinforcement.

    EN1992-1-1:2023 Eq. (8.62).

    Args:
        nu (float): The factor nu (unitless).
        f_cd (float): Design value of concrete compressive strength in MPa.
        cot_theta (float): Cotangent of the inclination of compression field in
            radians.
        cot_beta_incl (float): Cotangent of the inclination of load in radians
            (cot beta_incl).
        rho_w (float): Reinforcement ratio rho_w (unitless).
        f_ywd (float): Design yield strength of shear reinforcement in MPa.
        alpha_w (float): Angle of inclined shear reinforcement alpha_w in
            degrees (45 ≤ alpha_w < 90).
        cot_theta_min (float): max value for cot_theta.

    Returns:
        float: Shear stress resistance tau_Rd in MPa.

    Raises:
        ValueError: If any input parameter is negative where applicable.
    """
    if f_cd < 0:
        raise ValueError(f'f_cd must not be negative. Got {f_cd}')
    if f_ywd < 0:
        raise ValueError(f'f_ywd must not be negative. Got {f_ywd}')
    if alpha_w < 45 or alpha_w >= 90:
        raise ValueError(
            f'alpha_w must be between 45 degrees and 90 degrees. Got {alpha_w}'
        )

    alpha_w_rad = math.radians(alpha_w)
    cot_theta = min(max(math.tan(alpha_w_rad / 2), cot_theta), cot_theta_min)
    cot_alpha_w = 1 / math.tan(alpha_w_rad)

    tau_rd_value = nu * f_cd * (cot_theta - cot_beta_incl) / (
        1 + cot_theta**2
    ) + rho_w * f_ywd * (cot_beta_incl + cot_alpha_w) * math.sin(alpha_w_rad)

    tau_rd_max = nu * f_cd * (cot_theta + cot_alpha_w) / (1 + cot_theta**2)

    return min(tau_rd_value, tau_rd_max)


def sigma_swd_v2(
    Es: float, eps_x: float, cot_theta: float, alpha_w: float, f_ywd: float
) -> float:
    """Calculate the stress sigma_swd in the shear reinforcement for
    compression field inclinations.

    EN1992-1-1:2023 Eq. (8.63).

    Args:
        Es (float): Modulus of elasticity of steel Es in MPa.
        eps_x (float): Longitudinal strain epsilon_x (unitless).
        cot_theta (float): Cotangent inclination of compression field in
            radians.
        alpha_w (float): Angle of inclined shear reinforcement alpha_w in
            degrees (45 ≤ alpha_w < 90).
        f_ywd (float): Design yield strength of shear reinforcement in MPa.

    Returns:
        float: Stress sigma_swd in the shear reinforcement in MPa.

    Raises:
        ValueError: If any input parameter is negative where applicable.
    """
    if Es < 0:
        raise ValueError(f'Es must not be negative. Got {Es}')
    if f_ywd < 0:
        raise ValueError(f'f_ywd must not be negative. Got {f_ywd}')

    alpha_w_rad = math.radians(alpha_w)
    cot_alpha_w = 1 / math.tan(alpha_w_rad)

    sigma_swd_value = Es * (
        (eps_x + 0.001)
        * ((cot_theta + cot_alpha_w) ** 2 / (1 + cot_alpha_w**2))
        - 0.001
    )

    return min(sigma_swd_value, f_ywd)


def tao_Rd_m(tau_rd: float, m_ed: float, m_rd: float) -> float:
    """Calculate the shear stress resistance reduced by the influence of
    transverse bending.

    EN1992-1-1:2023 Eq. (8.64).

    Args:
        tau_rd (float): Shear resistance tau_Rd in MPa.
        m_ed (float): Applied transverse bending moment mEd in kNm.
        m_rd (float): Bending resistance without interaction with shear mRd in
            kNm.

    Returns:
        float: Reduced shear stress resistance tau_Rdm in MPa.

    Raises:
        ValueError: If any of the input parameters are negative.
    """
    if tau_rd < 0:
        raise ValueError(f'tau_rd must not be negative. Got {tau_rd}')
    if m_ed < 0:
        raise ValueError(f'm_ed must not be negative. Got {m_ed}')
    if m_rd < 0:
        raise ValueError(f'm_rd must not be negative. Got {m_rd}')

    return tau_rd * (1 - (m_ed / m_rd))


def tao_Ed_flang(
    delta_Fd: float,
    hf: float,
    delta_x: float,
) -> float:
    """Calculate the longitudinal shear stress at the junction between a flange
    and web.

    EN1992-1-1:2023 Eq. (8.65).

    Args:
        delta_Fd (float): Change of axial force in the flange over length delta
            x in kN.
        hf (float): Thickness of the flange at the junction in mm.
        delta_x (float): Length under consideration delta x in mm.

    Returns:
        float: Longitudinal shear stress tau_Ed in MPa.

    Raises:
        ValueError: If any of the input parameters are negative.
    """
    if delta_Fd < 0:
        raise ValueError(f'delta_fd must not be negative. Got {delta_Fd}')
    if hf <= 0:
        raise ValueError(f'hf must be positive. Got {hf}')
    if delta_x <= 0:
        raise ValueError(f'delta_x must be positive. Got {delta_x}')

    return delta_Fd * 1000 / (hf * delta_x)


def Ast_min_flang(tau_ed: float, sf: float, hf: float, fyd: float) -> float:
    """Calculate the minimum transverse reinforcement in web/flanges.

    EN1992-1-1:2023 Eq. (8.66)

    Args:
        tau_ed (float): Longitudinal shear stress tau_Ed MPa.
        sf (float): Spacing of reinforcement sf mm.
        hf (float): Thickness of the flange at the junction mm.
        fyd (float): Design yield strength of reinforcement fyd MPa.

    Returns:
        float: Minimum required transverse reinforcement area Asf in mm2.

    Raises:
        ValueError: If any of the input parameters are negative or
            if cot_theta_f is out of the specified range.
    """
    if tau_ed < 0:
        raise ValueError(f'tau_ed must not be negative. Got {tau_ed}')
    if sf <= 0:
        raise ValueError(f'sf must be positive. Got {sf}')
    if hf <= 0:
        raise ValueError(f'hf must be positive. Got {hf}')
    if fyd < 0:
        raise ValueError(f'fyd must not be negative. Got {fyd}')

    return (tau_ed * sf * hf) / (fyd)


def Asf_flang(
    tau_ed: float, sf: float, hf: float, fyd: float, cot_theta_f: float
) -> float:
    """Calculate the transverse reinforcement.

    EN1992-1-1:2023 Eq. (8.69).

    Args:
        tau_ed (float): Longitudinal shear stress tau_Ed MPa.
        sf (float): Spacing of reinforcement sf mm.
        hf (float): Thickness of the flange at the junction mm.
        fyd (float): Design yield strength of reinforcement fyd MPa.
        cot_theta_f (float): Cotangent of the inclination angle of the
            compression field in the flange thetaf.

    Returns:
        float: Required transverse reinforcement area Asf in mm2.

    Raises:
        ValueError: If any of the input parameters are negative or if
            cot_theta_f is out of the specified range.
    """
    if tau_ed < 0:
        raise ValueError(f'tau_ed must not be negative. Got {tau_ed}')
    if sf <= 0:
        raise ValueError(f'sf must be positive. Got {sf}')
    if hf <= 0:
        raise ValueError(f'hf must be positive. Got {hf}')
    if fyd < 0:
        raise ValueError(f'fyd must not be negative. Got {fyd}')

    return (tau_ed * sf * hf) / (fyd * cot_theta_f)


def sigma_cd_flang(
    tau_ed: float,
    theta_f: float,
    fcd: float,
    nu: float = 0.5,
) -> float:
    """Computes the compression field stress in the flange.

    EN1992-1-1:2023 Eq. (8.70), (8.71).

    Args:
        tau_ed (float or int): Longitudinal shear stress tau_Ed in MPa.
        theta_f (float or int): Inclination angle of the compression field in
            the flange thetaf (degrees).
        fcd (float or int): Design compressive strength of concrete fcd in MPa.
        nu (float, optional): Strength reduction factor, default is 0.5.

    Returns:
        bool: The compression field stress in MPa.

    Raises:
        ValueError: If any of the input parameters are
            negative or if theta_f is out of the specified range.
    """
    if tau_ed < 0:
        raise ValueError(f'tau_ed must not be negative. Got {tau_ed}')
    if theta_f < 0:
        raise ValueError(f'theta_f must not be negative. Got {theta_f}')
    if fcd < 0:
        raise ValueError(f'fcd must not be negative. Got {fcd}')
    if nu <= 0:
        raise ValueError(f'nu must be positive. Got {nu}')

    theta_rad = math.radians(theta_f)
    sigma_cd = tau_ed * (1 / math.tan(theta_rad) + math.tan(theta_rad))

    return min(sigma_cd, nu * fcd)


def eps_x_flang(Ftd: float, Ast: float, Es: float) -> float:
    """Calculate the longitudinal strain in the tensile flange based on the
    thickness of the tensile flange and the distance from the neutral axis.

    EN1992-1-1:2023 Eq. (8.72)

    Args:
        Ftd (float): force in the tension chord in kN.
        As (float): Longitudinal reinforcement area in mm2.
        Es (float): Steel modulus of elasticity in MPa.

    Returns:
        float: The longitudinal strain in the tensile flange (unitless).

    Raises:
        ValueError: If td is non-positive or if sst is non-positive.
    """
    if Ftd <= 0:
        raise ValueError(f'Ftd must be positive. Got {Ftd}')
    if Ast <= 0:
        raise ValueError(f'Ast must be positive. Got {Ast}')
    if Es <= 0:
        raise ValueError(f'Ast must be positive. Got {Es}')

    return Ftd * 1000 / (Ast * Es)


def tau_Edi(VEdi: float, Ai: float) -> float:
    """Calculate the design value of the shear stress at an interface.

    EN1992-1-1:2023 Eq. (8.74).

    Args:
        VEdi (float): Shear force acting parallel to the interface in kN.
        Ai (float): Area of the interface in mm2.

    Returns:
        float: Shear stress at the interface in MPa.

    Raises:
        ValueError: If any input value is negative.
    """
    if VEdi < 0:
        raise ValueError(f'VEdi must not be negative. Got {VEdi}')
    if Ai < 0:
        raise ValueError(f'Ai must not be negative. Got {Ai}')

    return VEdi * 1000 / Ai


def tau_Edi_composite(
    beta_new: float, VEd: float, z: float, bi: float
) -> float:
    """Calculate the longitudinal shear stress between concrete interfaces due
    to composite action.

    EN1992-1-1:2023 Eq. (8.75).

    Args:
        beta_new (float): Ratio of the longitudinal force in the new concrete
            to the total longitudinal force, dimensionless.
        VEd (float): Shear force acting perpendicular to the interface in kN.
        z (float): Lever arm of the composite section in mm.
        bi (float): Width of the interface in mm.

    Returns:
        float: Longitudinal shear stress at the interface in MPa.

    Raises:
        ValueError: If any input value is negative or zero when it should not
            be.
    """
    if beta_new < 0:
        raise ValueError(f'beta_new must not be negative. Got {beta_new}')
    if VEd < 0:
        raise ValueError(f'Ved must not be negative. Got {VEd}')
    if z <= 0:
        raise ValueError(f'z must be positive. Got {z}')
    if bi < 0:
        raise ValueError(f'bi must not be negative. Got {bi}')

    return beta_new * VEd * 1000 / (z * bi)


def tau_Rdi(
    fck: float,
    sigma_n: float,
    Ai: float,
    Asi: float,
    fyd: float,
    alpha_deg: float,
    cv1: float,
    mu_v: float,
    gamma_c: float,
) -> float:
    """Calculate the design shear stress resistance at the interface for
    scenarios without reinforcement or where reinforcement is sufficiently
    anchored.

    EN1992-1-1:2023 Eq. (8.76).

    Args:
        fck (float): Lowest compressive strength of the concretes at the
            interface in MPa.
        sigma_n (float): Compressive or tensile stress over the interface area
            in MPa.
        Ai (float): Area of the interface in mm2.
        Asi (float): Cross-sectional area of bonded reinforcement crossing the
            interface in mm2.
        fyd (float): Design yield strength of the reinforcement in MPa.
        alpha_deg (float): Angle of reinforcement crossing the interface in
            degrees.
        cv1 (float): Coefficient depending on the roughness of the interface,
            dimensionless.
        mu_v (float): Friction coefficient depending on the roughness of the
            interface, dimensionless.
        gamma_c (float): safety factory for concrete.

    Returns:
        float: Shear stress resistance at the interface in MPa.

    Raises:
        ValueError: If any input value is negative or outside expected ranges.
    """
    if fck < 0 or Ai < 0 or Asi < 0 or fyd < 0 or gamma_c < 0:
        raise ValueError('Input values must not be negative.')
    if not (35 <= alpha_deg <= 135):
        raise ValueError('Alpha must be between 35 and 135 degrees.')

    sigma_n = max(0, min(sigma_n, 0.6 * fck / gamma_c))

    alpha_rad = math.radians(alpha_deg)
    rho_i = Asi / Ai

    tau_rdi = (
        cv1 * math.sqrt(fck) / gamma_c
        + mu_v * sigma_n
        + rho_i * fyd * (mu_v * math.sin(alpha_rad) + math.cos(alpha_rad))
    )

    # Limiting tau_rdi according to the specification
    return min(tau_rdi, 0.30 * fck + rho_i * fyd * math.cos(alpha_rad))


def cv1(
    surface_roughness: Literal[
        'very smooth', 'smooth', 'rough', 'very rough', 'keyed'
    ],
    tensile_stress: bool = False,
) -> float:
    """Get the cv1 coefficient based on the surface roughness and tensile
    stress condition.

    EC1992-1-1:2023 Table (8.2).

    Args:
        surface_roughness (str): Description of the surface roughness.
        tensile_stress (bool): True if tensile stresses are present.

    Returns:
        float: The cv1 coefficient.

    Raises:
        ValueError: If an unknown surface roughness is provided.
    """
    if tensile_stress:
        return 0

    coefficients = {
        'very smooth': 0.01,
        'smooth': 0.08,
        'rough': 0.15,
        'very rough': 0.19,
        'keyed': 0.37,
    }
    return coefficients[surface_roughness]


def mu_v(
    surface_roughness: Literal[
        'very smooth', 'smooth', 'rough', 'very rough', 'keyed'
    ],
) -> float:
    """Get the mu_v coefficient based on the surface roughness.

    EC1992-1-1:2023 Table (8.2).

    Args:
        surface_roughness (str): Description of the surface roughness.

    Returns:
        float: The mu_v coefficient.

    Raises:
        ValueError: If an unknown surface roughness is provided.
    """
    coefficients = {
        'very smooth': 0.5,
        'smooth': 0.6,
        'rough': 0.7,
        'very rough': 0.9,
        'keyed': 0.9,
    }
    return coefficients[surface_roughness]


def cv2(
    surface_roughness: Literal['very smooth', 'smooth', 'rough', 'very rough'],
    tensile_stress: bool = False,
) -> float:
    """Get the cv2 coefficient based on the surface roughness and tensile
    stress condition.

    EC1992-1-1:2023 Table (8.2).

    Args:
        surface_roughness (str): Description of the surface roughness.
        tensile_stress (bool): True if tensile stresses are present.

    Returns:
        float: The cv2 coefficient, or None for keyed surfaces.

    Raises:
        ValueError: If an unknown surface roughness is provided.
    """
    if tensile_stress:
        return 0

    coefficients = {
        'very smooth': 0,
        'smooth': 0,
        'rough': 0.08,
        'very rough': 0.15,
    }
    return coefficients[surface_roughness]


def kv(
    surface_roughness: Literal['very smooth', 'smooth', 'rough', 'very rough'],
) -> float:
    """Get the kv coefficient based on the surface roughness.

    EC1992-1-1:2023 Table (8.2).

    Args:
        surface_roughness (str): Description of the surface roughness.

    Returns:
        float: The kv coefficient, or None for keyed surfaces.

    Raises:
        ValueError: If an unknown surface roughness is provided.
    """
    coefficients = {
        'very smooth': 0,
        'smooth': 0.5,
        'rough': 0.5,
        'very rough': 0.5,
    }
    return coefficients[surface_roughness]


def kdowel(
    surface_roughness: Literal['very smooth', 'smooth', 'rough', 'very rough'],
) -> float:
    """Get the kdowel coefficient based on the surface roughness.

    EC1992-1-1:2023 Table (8.2).

    Args:
        surface_roughness (str): Description of the surface roughness.

    Returns:
        float: The kdowel coefficient, or None for keyed surfaces.

    Raises:
        ValueError: If an unknown surface roughness is provided.
    """
    coefficients = {
        'very smooth': 1.5,
        'smooth': 1.1,
        'rough': 0.9,
        'very rough': 0.9,
    }
    return coefficients[surface_roughness]


def tau_Rdi_ny(
    cv2: float,
    fck: float,
    gamma_c: float,
    mu_v: float,
    sigma_n: float,
    kv: float,
    rho_i: float,
    fyd: float,
    kdowel: float,
) -> float:
    """Calculate the shear stress resistance at the interface when yielding is
    not ensured at the interface.

    EN 1992-1-1:2022 Eq. (8.77).

    Args:
        cv2 (float): Coefficient depending on the roughness of the interface
            (unitless).
        fck (float): Concrete compressive resistance in MPa.
        gamma_c (float): Partial safety factor for concrete (unitless).
        mu_v (float): Coefficient mu_v from the Eurocode (unitless).
        sigma_n (float): Normal stress in the interface MPa.
        kv (float): Coefficient kv from the Eurocode (unitless).
        rho_i (float): Reinforcement ratio at the interface (unitless).
        fyd (float): Design yield strength of reinforcement MPa.
        kdowel (float): Coefficient for dowel action of reinforcement
            (unitless).

    Returns:
        float: Shear stress resistance tau_Rdi MPa.

    Raises:
        ValueError: If any of the dimensions or resistances are negative.
    """
    if any(x < 0 for x in [fck, gamma_c, kv, rho_i, fyd, kdowel]):
        raise ValueError('Dimensions and resistances must not be negative.')

    sigma_n = max(0, sigma_n)
    tau_rdi = (
        cv2 * math.sqrt(fck) / gamma_c
        + mu_v * sigma_n
        + kv * rho_i * fyd * mu_v
        + kdowel * rho_i * math.sqrt(fyd * fck / gamma_c)
    )
    return min(
        tau_rdi, 0.25 * fck / gamma_c
    )  # Cap tau_Rdi according to the formula


def as_min(tmin: float, fctm: float, fyk: float) -> float:
    """Calculate the minimum interface reinforcement per unit length along the
    edge of composite slabs.

    EN 1992-1-1:2022 Eq. (8.78).

    Args:
        tmin (float): Smaller value of the thickness of new and old concrete
            layers in mm.
        fctm (float): Mean tensile strength of the respective concrete layer in
            MPa.
        fyk (float): Characteristic yield strength of the reinforcement in MPa.

    Returns:
        float: Minimum interface reinforcement per unit length in mm/mm.

    Raises:
        ValueError: If any input is negative or zero where it shouldn't be.
    """
    if tmin <= 0:
        raise ValueError(f'tmin must be positive. Got {tmin}')
    if fctm <= 0:
        raise ValueError(f'fctm must be positive. Got {fctm}')
    if fyk <= 0:
        raise ValueError(f'fyk must be positive. Got {fyk}')

    return tmin * fctm / fyk
