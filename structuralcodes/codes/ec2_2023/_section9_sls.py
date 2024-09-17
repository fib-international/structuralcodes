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
