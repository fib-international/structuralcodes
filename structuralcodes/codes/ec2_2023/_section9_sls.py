"""Functions from Section 9 of EN 1992-1-1:2023."""

import math
from typing import Tuple

from ._annexB_time_dependent import alpha_c
from ._section5_materials import fcm, fctm


def Ec_eff(fcm: float, phi: float, kE: float = 9500) -> float:
    """Returns de effective modulus of elasticity from fcm and phi.

    EN 1992-1-1:2023, Eq. (9.1).

    Args:
        fcm (float): The mean compressive strength in MPa.
        phi (float): The creep coefficient.

    Keyword Args:
        kE (float): Constant to account for the type of aggregate.

    Returns:
        float: The effective modulus of elastiticy in MPa.
    """
    Ecm = kE * fcm ** (1 / 3)
    return alpha_c(fcm) * Ecm / (1 + phi)


def As_min_y(
    NEd: float, b: float, h: float, fct_eff: float, fyk: float
) -> Tuple[float, float]:
    """Returns the minimum reinforcement to avoid yielding of steel. Box or T
    sections are to be divided into rectangles.

    EN 1992-1-1:2023, Eq. (9.4)

    Eq. (9.2) and (9.3) are particular cases of the general equation

    Eq. (9.2) is valid for pure bending, hence NEd=0

    Eq. (9.3) is valid for pure tension. The general expression has an upper
    limit that equals the values of Eq. (9.3)

    Args:
        NEd (float): SLS axial force applied on the section or rectangle
            (compressions are negative) in kN.
        b (float): The width of the section or rectangle in meters.
        h (float): The height of the section or rectange in meters.
        fct_eff (float): Effective tension strength of concrete (can normally
            be taken as the mean tensile strength) in MPa.
        fyk (float): Characteristic yield strength of steel in MPa.

    Returns:
        tuple(float, float): The minimum tensile reinforcement to avoid
        yielding of steel on the most tensioned fibre of the rectangle
        (As_min_y1) in cm2, and the minimum tensile reinforcement to avoid
        yielding of steel on the most tensioned fibre of the rectangle
        (As_min_y2) in cm2.
    """
    As_min_y1 = (
        min(
            max(
                (0.3 * NEd / 1000 + 0.2 * kh(b, h) * fct_eff * b * h) / fyk, 0
            ),
            0.5 * kh(b, h) * fct_eff * b * h / fyk,
        )
        * 1e4
    )
    return As_min_y1, min(max(NEd / fyk * 10 - As_min_y1, 0), As_min_y1)


def kh(b: float, h: float) -> float:
    """Returns factor kh, which reduces the tensile strength of concrete to
    account for imposed restrained deformations due to shrinkage.

    EN 1992-1-1:2023, Eq. (9.5).

    Args:
        b (float): Width of the rectangle in meters.
        h (float): Height of the rectangle in meters.

    Returns:
        float: Factor kh which reduces the tensile strength of concrete to
        account for imposed restrained deformations due to shrinkage.
    """
    return min(max(0.8 - 0.6 * (min(b, h) - 0.3), 0.5), 0.8)


def wk_cal2(
    kw: float, k_1_r: float, srm_cal: float, epssm_epscm: float
) -> float:
    """Returns the calculated characteristic crack width.

    EN 1992-1-1:2023, Eq. (9.8).

    Args:
        kw (float): Factor that converts the mean crack spacing to a
            characteristic value.
        k_1_r (float): Factor accounting for the effect of curvature on crack
            width - can be determined using the function k_1_r.
        srm_cal (float): Mean crack spacing - can be determined using the
            function srm_cal.
        epssm_epscm (float): Mean diference of strain between steel and
            concrete - can be determined using the function epssm_epscm.

    Returns:
        float: The calculated characteristic crack width in in units consistent
        with srm_cal.
    """
    return kw * k_1_r * srm_cal * epssm_epscm


def k_1_r(h: float, x: float, ay: float) -> float:
    """Returns k1/r factor to account for increase in crack width due to
    curvature of the section in bending.

    EN 1992-1-1:2023, Eq. (9.9).

    Args:
        h (float): Height of the section in consistent units (e.g. meters).
        x (float): Distance from most compressed fibre to neutra axis in
            consistent units (e.g. meters).
        ay (float): Cover to centre of tensioned reinforcement closest to most
            tensioned face in consistent units (e.g. meters).

    Returns:
        float: Factor k1/r (non-dimensional) which accounts for the increase in
        crack width due to curvature of the section in bending.
    """
    return (h - x) / (h - ay - x)


def epssm_epscm(
    sigma_s: float,
    kt: float,
    fct_eff: float,
    rho_eff: float,
    alphae: float,
    Es: float,
) -> float:
    """Returns the mean strain difference between steel and concrete along 2
    transfer lengths.

    EN 1992-1-1:2023, Eq. (9.11).

    Args:
        sigma_s (float): The stress in steel at the section of the crack.
        kt (float): An integration factor to account for the variation in
            strain in steel and concrete it is to be taken as 0.6 for short
            term loading or instantaneous loading and equal to 0.4 for long
            term or repeated loading.
        fct_eff (float): The effective cracking stress, which can be taken
            equal to the mean tensile strength of concrete, fctm.
        rho_eff (float): The effective reinforcement ratio in the tension zone.
        alphae (float): The equivalence factor equal to Es/Ecm.
        Es (float): The modulus of elasticity of steel, normally taken as 200
            GPa.

    Returns:
        float: The mean strain difference bewteen steel and concrete along 2
        transfer lengths.
    """
    return max(
        (sigma_s - kt * fct_eff / rho_eff * (1 + alphae * rho_eff)) / Es,
        (1 - kt) * sigma_s / Es,
    )


def kfl(h: float, xg: float, hceff: float) -> float:
    """Returns factor kfl which accounts for the distribution of stresses
    before cracking.

    EN 1992-1-1:2023, Eq. (9.17).

    Args:
        h (float): Height of the cross section.
        xg (float): Distance from the compressed fibre to the centroid of the
            uncracked section.
        hceff (float): Height of the effective tension area.

    Returns:
        float: Returns factor kfl which accounts for the distribution of
        stresses before cracking.
    """
    return max(0.5 * (1 + (h - xg - hceff) / (h - xg)), 0.5)


def srm_cal(
    c: float,
    kfl_: float,
    kb: float,
    phi: float,
    rho_eff: float,
    kw: float,
    h,
    x: float,
) -> float:
    """Returns the mean crack spacing.

    EN 1992-1-1:2023, Eq. (9.15).

    Args:
        c (float): Concrete cover of reinforcement to bar surface. Larger value
            of lateral and vertical cover should be applied.
        kfl (float): Factor accounting for distribution of stresses prior to
            cracking.
        kb (float): Factor accounting for bond conditions.
        phi (float): Bar diameter.
        rho_eff(float): Effective reinforcement ratio in the tension zone.
        kw (float): Factor converting the mean crack spacing into a
            characteristic crack spacing, with a reocmmended value of 1.3
            (NDP).
        h (float): Height of the cross section.
        x (float): Depth of the neutral axis measured form the most compressed
            fibre.

    Returns:
        float: The mean crack spacing in units consistent with c and phi.
    """
    return min(1.5 * c + kfl_ * kb / 7.2 * phi / rho_eff, 1.3 / kw * (h - x))


def wk_cal(
    kw: float,
    h: float,
    xg: float,
    hc_eff: float,
    c: float,
    kb: float,
    phi: float,
    rho_eff: float,
    x: float,
    sigma_s: float,
    kt: float,
    fct_eff: float,
    alphae: float,
    Es: float,
) -> Tuple[float, float, float, float]:
    """Returns the characteristic crack width, wk,cal, as well as auxiliary
    variables, 1/r, srm,cal and epssm-epscm.

    EN1992-1-1:2023 Eq. (9.8), complemented with Eq. (9.11), Eq. (9.15), Eq.
    (9.17).

    Args:
        kw (float): Factor that converts the mean crack spacing to a
            characteristic value.
        h (float): Height of cross section.
        xg (float): Depth of centroid of section measured from compressed
            fibre.
        hc_eff (float): Height of the effective tensioned concrete area.
        c (float): Concrete cover of reinforcement to bar surface. Larger
            value of lateral and vertical cover should be applied.
        kb (float): Factor account for bond conditions of bar.
        phi (float): Diameter of tensioned bars (for different bar diameters,
            equivalent diameter according to Eq. (9.19).
        rho_eff (float): Effective tension reinforcement ratio.
        x (float): Depth of the neutral axis of the cracked section measured
            from compressed fibre.
        sigma_s (float): Tension in most tensioned bar according to fully
            cracked analysis.
        kt (float): Factor accounting for tension stiffening.
        fct_eff (float): Effective tensile strength of concrete.
        alphae (float): Modular ratio Es/Ecm.
        Es (float): Modulus of elasticity of steel bars (normally Es=200 MPa).

    Returns:
        Tuple[float, float, float, float]: The characteristic crack width,
        wk,cal, in consistent units, as well as auxiliary variables, 1/r,
        srm,cal and epssm-epscm.
    """
    k_1_r_ = k_1_r(h, x, c + phi / 2)
    srm_cal_ = srm_cal(c, kfl(h, xg, hc_eff), kb, phi, rho_eff, kw, h, x)
    epssm_epscm_ = epssm_epscm(sigma_s, kt, fct_eff, rho_eff, alphae, Es)
    wk_cal_ = kw * k_1_r_ * srm_cal_ * epssm_epscm_
    return wk_cal_, k_1_r_, srm_cal_, epssm_epscm_


def delta_simpl(
    delta_loads: float,
    delta_shr: float,
    fck1: float,
    phi1: float,
    b1: float,
    h: float,
    d: float,
    As1: float,
    Mk: float,
) -> float:
    """Simplified calculation of the deflection for rectangular sections.

    EN1992-1-1:2023, Eq. (9.23).

    Args:
        delta_loads (float): Linear elastic deflection due to loads.
        delta_shr (float): Linear elastic deflection due to shrinkage.
        fck1 (float): Characteristic concrete strength in MPa.
        phi1 (float): Weighted mean value of the creep coefficient.
        b1 (float): Width of rectangular cross-section in m.
        h (float): Height of rectanguar cross-section in m.
        d (float): Effective height of cross-section in m.
        As1 (float): Tension reinforcement at centre span for continuous in cm2
            beams or at the embedment for a cantilever.
        Mk (float): Characteristic moment at centre span for continuous
            beams or at the embedment for a cantilever.

    Returns:
        float: The deflection of the beam in units consistent with delta_loads
        and delta_shr.
    """
    Mcr = fctm(fck1) * 1000 * b1 * h ** (2) / 6
    zeta = 1 - 0.5 * (Mcr / Mk) ** 2
    rho_l = As1 / b1 / d * 1e-4
    alpha_e_eff = 200000 / Ec_eff(fcm(fck1), phi1)
    Ig_Icr = 1 / (
        2.7 * math.pow(alpha_e_eff * rho_l, 0.6) * math.pow(d / h, 3)
    )

    if Mk < Mcr:
        kS = 1.00
        kI = 1.00
    else:
        kS = 455 * rho_l**2 - 35 * rho_l + 1.6
        kI = zeta * Ig_Icr + (1 - zeta)
    return kI * (delta_loads + kS * delta_shr)


def rho_p_eff(As: float, xi1: float, Ap: float, Ac_eff: float) -> float:
    """Effective bond ratio between areas.

    EN 1992-1-1:2023, Eq. (9.12).

    Args:
        As (float): Steel area in mm2.
        xi1 (float): The adjusted ratio of bond according to expression (9.6).
        Ap (float): The area in mm2 of post-tensioned tendons in ac_eff.
        Ac_eff (float): Effective area of concrete in tension surrounding the
            reinforcement or prestressing tendons of depth hc_eff.

    Returns:
        float: With the retio between areas.

    Raises:
        ValueError: If any of As, xi1, Ap or Ac_eff is less than 0.
    """
    if As < 0:
        raise ValueError(f'As={As} cannot be less than 0')
    if xi1 < 0:
        raise ValueError(f'xi1={xi1} cannot be less than 0')
    if Ap < 0:
        raise ValueError(f'Ap={Ap} cannot be less than 0')
    if Ac_eff < 0:
        raise ValueError(f'Ac_eff={Ac_eff} cannot be less than 0')

    return (As + xi1**2 * Ap) / Ac_eff


def xi1(xi: float, phi_p: float, phi_s: float) -> float:
    """Computes the adjusted ratio of bond strength taking into account
    the different diameters of prestressing and reinforcing steel.

    EN 1992-1-1:2023, Eq. (9.6).

    Args:
        xi (float): ratio of bond strength of prestressing and reinforcing
            steel, according to Table 10.1 in 10.3(2)
        phi_p (float): largest bar diameter in mm of reinforcing steel.
            Equal to 0 if only reinforced steel is used in control cracking.
        phi_s (float): equivalent diameter in mm of tendon acoording
            to 9.19.

    Returns:
        float: With the value of the ratio.

    Raises:
        ValueError: If diameters phi_s or phi_p are lower than 0. If ratio of
            bond strength xi is less than 0.15 or larger than 0.8.
    """
    if phi_p == 0:
        return 0
    if phi_p < 0:
        raise ValueError(f'phi_p={phi_p} cannot be less than 0')
    if phi_s < 0:
        raise ValueError(f'phi_s={phi_s} cannot be less than 0')
    if xi < 0.15:
        raise ValueError(f'The minimum value for xi={xi} is 0.15')
    if xi > 0.8:
        raise ValueError(f'The maximum value for xi={xi} is 0.8')

    return ((xi * phi_s / phi_p) ** 0.5) if phi_s > 0 else xi**0.5


def _lower_circular_segment_area(d, x):
    """Calculates the area of the lower circular segment cut by a horizontal line at a height `x`
    from the top of the circle.

    Args:
        d (float): Diameter of the circle. Must be greater than 0.
        x (float): Vertical distance from the top edge of the circle to the horizontal chord.
                   Must satisfy 0 < x < d.

    Returns:
        float: Area of the circular segment below the chord.

    Raises:
        ValueError: If `x` is not in the range (0, d).
    """
    if x < 0 or x > d:
        raise ValueError(
            'x must be between 0 and d (the diameter of the circle).'
        )

    r = d / 2
    y = r - x  # vertical distance from the center to the chord
    theta = math.acos(y / r)  # central angle in radians
    upper_area = r**2 * theta - y * math.sqrt(r**2 - y**2)
    return r**2 * math.pi - upper_area  # lower area


def Ac_eff(
    x: float,
    ay,
    phi,
    h: float = None,
    b: float = None,
    diameter: float = None,
    n: int = 1,
    sy: float = None,
    loading_type: str = 'bending',
    section_type: str = 'rectangular',
    bar_spacing: float = None,
    ax: float = None,
) -> tuple[float, float]:
    """Returns the effective area of concrete in tension surrounding the
    reinforcement or prestressing tendons.

    EN 1992-1-1:2023, Figure 9.3 (a to f).

    Args:
        x (float): distance in mm to the zero tensile stress line.
        ay (float): distance from extreme fibre of concrete to centroid of reinforcement bars in mm.
        phi (float): diameter of the tensioned bars in mm.
        h (float): total heigth of the rectangular section in mm.
        b (float): total width of the rectangular section in mm.
        diameter (float): diameter of the circular section in mm.
        n (int): Number of layers. Default is 1.
        sy (float): Spacing of the layers in mm. Default is None.
        loading_type (str): Type of loading. Default is 'bending'. Can be 'bending' or 'tension'.
        section_type (str): Type of section. Default is 'rectangular'. Can be 'rectangular' or 'circular'.
        bar_spacing (float): Spacing of the bars in mm. Default is less than 10 times phi.
        ax (float): distance from extreme fibre of concrete to centroid of reinforcement bars in mm (Figure 9.3 f).

    Returns:
        tuple: (Ac_eff, hc_eff)
            Ac_eff (float): The effective area in mm2.
            hc_eff (float): The effective height in mm.

    Raises:
        ValueError: If any of h, ay or x is lower than zero.
        ValueError: If x is greater than h.
        ValueError: If ay is greater than h.
        ValueError: If phi is lower than zero.
        ValueError: If n is lower than 1.
        ValueError: If sy is None and n > 1.
        ValueError: If loading_type is not 'bending' or 'tension'.
        ValueError: If section_type is not 'rectangular' or 'circular'.
        ValueError: If missing a parameter gor rectangular or circular section.
    """
    if section_type not in ['rectangular', 'circular']:
        raise ValueError(
            f'section_type={section_type} not implemented. It should be "rectangular" or "circular for this method"'
        )
    if section_type == 'rectangular':
        if b is None:
            raise ValueError('b must be provided for rectangular sections')
        if h is None:
            raise ValueError('h must be provided for rectangular sections')
        if h < 0:
            raise ValueError(f'h={h} cannot be less than 0')
        if b < 0:
            raise ValueError(f'b={b} cannot be less than 0')
        if ay > h:
            raise ValueError(f'ay={ay} cannot be larger than h={h}')
        if x > h:
            raise ValueError(f'x={x} cannot be larger than h={h}')
        bc_eff = (
            b if bar_spacing is None else min(b, 10 * phi / bar_spacing * b)
        )
    if section_type == 'circular':
        if diameter is None:
            raise ValueError('diameter must be provided for circular sections')
        if diameter < 0:
            raise ValueError(f'diameter={diameter} cannot be less than 0')
        if ay > diameter:
            raise ValueError(f'ay={ay} cannot be larger than D={diameter}')
        if x > diameter:
            raise ValueError(f'x={x} cannot be larger than D={diameter}')

    if ay < 0:
        raise ValueError(f'ay={ay} cannot be less than 0')
    if x < 0:
        raise ValueError(f'x={x} cannot be less than zero')
    if phi < 0:
        raise ValueError(f'phi={phi} cannot be less than 0')
    if n < 1:
        raise ValueError(f'n={n} cannot be less than 1')
    if sy is None and n > 1:
        raise ValueError(f'sy cannot be None if n={n} > 1')
    if loading_type not in ['bending', 'tension']:
        raise ValueError(
            f'loading_type={loading_type} is not valid. It should be "bending" or "tension"'
        )

    if section_type == 'rectangular' and loading_type == 'bending':
        if n == 1:
            hc_eff = min(ay + 5 * phi, 10 * phi, 3.5 * ay, h - x, h / 2)
        else:
            hc_eff = min(
                min(ay + 5 * phi, 10 * phi, 3.5 * ay) + (n - 1) * sy,
                h - x,
                h / 2,
            )
        Ac_eff = bc_eff * hc_eff
    if section_type == 'circular' and loading_type == 'bending':
        hc_eff = min(ay + 5 * phi, 10 * phi, 3.5 * ay)
        Ac_eff = _lower_circular_segment_area(
            diameter, x
        ) - _lower_circular_segment_area(
            diameter - hc_eff * 2, max(0, x - hc_eff)
        )
    if section_type == 'rectangular' and loading_type == 'tension':  # d,e,f
        hc_eff = min(ay + 5 * phi, 10 * phi, 3.5 * ay, h / 2)
        if ax is not None:  # f
            bc_eff = min(ax + 5 * phi, 10 * phi, 3.5 * ax, b / 2)
            Ac_eff = h * b - (b - 2 * bc_eff) * (h - 2 * hc_eff)
        else:  # d,e
            Ac_eff = bc_eff * hc_eff

    return (Ac_eff, hc_eff)
