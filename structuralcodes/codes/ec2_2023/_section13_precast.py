"""Functions from Section 13 of EN 1992-1-1:2023."""

import math
from typing import List, Literal


def fcm_t_precast(
    fcmp: float, fcm_tref: float, t: float, tref: float, tp: float
) -> float:
    """Calculate the mean compressive strength of concrete at an age `t` less
        than `tref` after heat curing.

    EN1992-1-1:2023 Eq. (13.1)

    Args:
        fcmp (float): Mean compressive strength after heat curing in MPa.
        fcm_tref (float): Mean compressive strength
            at the reference age `tref` in MPa.
        t (float): Age of the concrete in days.
        tref (float): Reference age of the concrete in days.
        tp (float): Age at which the heat curing is completed in days.

    Returns:
        float: Mean compressive strength of concrete at age `t` in days.

    Raises:
        ValueError: If any input value is negative.
    """
    if fcmp < 0 or fcm_tref < 0 or t < 0 or tref < 0 or tp < 0:
        raise ValueError(
            f'All input values must be non-negative. Got fcmp={fcmp}, '
            + f' fcm_tref={fcm_tref}, t={t}, tref={tref}, tp={tp}'
        )

    if t <= tp:
        raise ValueError(
            'The age t must be greater than tp (t > tp).'
            + f' Got t={t}, tp={tp}'
        )

    numerator = math.log10(t - tp + 1)
    denominator = math.log10(tref - tp + 1)

    return fcmp + (fcm_tref - fcmp) * (numerator / denominator)


def t_eq_precast(
    temperatures: List[float],
    time_intervals: List[float],
    tmax: float,
) -> float:
    """Calculate the equivalent time (teq) to account for the
    effect of heat treatment on relaxation losses in precast concrete.

    EN1992-1-1:2023 Eq. (13.2)

    Args:
        teq_max (float): Maximum temperature during the heat treatment in °C.
        temperatures (List[float]): List of temperatures during each
            time interval in °C.
        time_intervals (List[float]): List of corresponding time
            intervals in hours.
        tmax (float): Maximum temperature during heat treatment in °C.

    Returns:
        float: Equivalent time t_eq_precast in hours.

    Raises:
        ValueError: If any temperature or time interval is
            negative, if lists are of unequal length or Tmax
            is less or equal to 20 °C.
    """
    if len(temperatures) != len(time_intervals):
        raise ValueError(
            'The lengths of temperatures and'
            + 'time_intervals lists must be the same.'
        )

    if any(t < 0 for t in temperatures) or any(
        ti < 0 for ti in time_intervals
    ):
        raise ValueError(
            'Temperatures and time intervals must be non-negative.'
        )
    if tmax <= 20:
        raise ValueError('Tmax should be larger than 20 °C')

    term1 = math.pow(1.14, tmax - 20) / (tmax - 20)
    _sum = 0
    for t, delta_t in zip(temperatures, time_intervals):
        _sum += delta_t * (t - 20)

    return term1 * _sum


def delta_signma_th_precast(
    Ep: float, alpha_c_th: float, tmax: float, t0: float
) -> float:
    """Calculate the specific thermal loss Δσθ induced by heat treatment
        in precast concrete.

    EN1992-1-1:2023 Eq. (13.3)

    Args:
        Ep (float): Modulus of elasticity of the prestressing tendons in MPa.
        alpha_c_th (float): Coefficient of thermal expansion of
            concrete in 1/°C.
        tmax (float): Maximum concrete temperature near
            the tendons during heat treatment in °C.
        T0 (float): Initial concrete temperature
            near the tendons in °C.

    Returns:
        float: Specific thermal loss Δσθ in MPa.

    Raises:
        ValueError: If any temperature is negative, or if Tmax < T0.
    """
    if tmax < t0:
        raise ValueError(
            'Tmax must be greater than or equal to T0. '
            + f'Got Tmax={tmax}, T0={t0}'
        )
    if Ep < 0 or alpha_c_th < 0:
        raise ValueError(
            'Modulus of elasticity and coefficient of thermal '
            + 'expansion must be non-negative. '
            + f'Got Ep={Ep}, alpha_c_th={alpha_c_th}'
        )

    return 0.5 * Ep * alpha_c_th * (tmax - t0)


def min_csx_precast(csx: float, phi_p: float, Dupper: float) -> float:
    """Calculate the minimum clear spacing in horizontal direction
        for bond for pre-tensioning tendons in precast concrete.

    EN1992-1-1:2023 Section 13.5.1(2)

    Args:
        csx (float): Clear spacing in the x-direction in mm.
        phi_p (float): Diameter of the pre-tensioning tendon in mm.
        Dupper (float): Diameter of the upper reinforcement in mm.

    Returns:
        float: Minimum clear spacing in x direction in mm.

    Raises:
        ValueError: If any input value is negative.
    """
    if csx < 0 or phi_p < 0 or Dupper < 0:
        raise ValueError(
            f'All input values must be non-negative. Got csx={csx}, '
            + f' phi_p={phi_p}, Dupper={Dupper}'
        )

    return max(2 * phi_p, Dupper + 5, 20)


def min_csy_precast(csy: float, phi_p: float, Dupper: float) -> float:
    """Calculate the minimum clear spacing in vertical direction
        for bond for pre-tensioning tendons in precast concrete.

    EN1992-1-1:2023 Section 13.5.1(2)

    Args:
        csy (float): Clear spacing in the y-direction in mm.
        phi_p (float): Diameter of the pre-tensioning tendon in mm.
        Dupper (float): Diameter of the upper reinforcement in mm.

    Returns:
        float: Minimum clear spacing in y direction mm.

    Raises:
        ValueError: If any input value is negative.
    """
    if csy < 0 or phi_p < 0 or Dupper < 0:
        raise ValueError(
            f'All input values must be non-negative. Got csy={csy}, '
            + f' phi_p={phi_p}, Dupper={Dupper}'
        )

    return max(2 * phi_p, Dupper + 5, 20)


def cmin_b_precast(
    s: float, phi_p: float, tendon_type: Literal['strand', 'indented_wire']
) -> float:
    """Calculate the minimum concrete cover cmin,b for pre-tensioning tendons
        in precast concrete.

    EN1992-1-1:2023 Table 13.1

    Args:
        s (float): Clear spacing between tendons in mm.
        phi_p (float): Diameter of the pre-tensioning tendon in mm.
        tendon_type (str): Type of tendon ('strand' or 'indented_wire').

    Returns:
        float: Minimum concrete cover cmin,b in mm.

    Raises:
        ValueError: If any input value is negative or
            if tendon_type is not recognized.
    """
    if s < 0 or phi_p < 0 or phi_p < 0:
        raise ValueError(
            f'All input values must be non-negative. Got s={s}, '
            + f' phi_p={phi_p}'
        )

    if tendon_type == 'strand':
        if s == 2 * phi_p:
            return 3.0 * phi_p
        if s > 2 * phi_p:
            return 2.5 * phi_p
        raise ValueError(f's should be larger or equal than {2*phi_p}')

    # If indented_wire
    if s == 2 * phi_p:
        return 4.5 * phi_p
    if s > 2 * phi_p:
        return 4.0 * phi_p
    raise ValueError(f's should be larger or eequal than {2*phi_p}')


def lpt_precast(
    gamma_c: float,
    release: Literal['gradual', 'sudden'],
    tendon: Literal['indented_wires', 'stands'],
    sigma_pm0: float,
    position: Literal['favorable', 'unfavorable'],
    fck_t: float,
    phi_p: float,
) -> float:
    """Calculate the basic value of the transmission
        length lpt for prestressing tendons in precast concrete.

    EN1992-1-1:2023 Eq. (13.4)

    Args:
        gamma_c (float): Partial safety factor for concrete.
        release (float): Gradual or sudden release.
        tendon (str): type of tendon,
        sigma_pm0 (float): Tendon stress just after release in MPa.
        position: favorable or unfavorable.
        fck_t (float): Concrete compressive
            strength at the time of release in MPa.
        phi_p (float): Nominal diameter of the tendon in mm.

    Returns:
        float: Transmission length lpt in mm.

    Raises:
        ValueError: If any input value is negative.
    """
    if gamma_c < 0 or sigma_pm0 < 0 or fck_t < 0 or phi_p < 0:
        raise ValueError(
            f'Input values must be non-negative. Got gamma_c={gamma_c},'
            + f' sigma_pm0={sigma_pm0},'
            + f' fck_t={fck_t}, phi_p={phi_p}'
        )

    alpha1 = 1.0 if release == 'gradual' else 1.25
    alpha2 = 0.4 if tendon == 'indented_wires' else 0.26
    eta1 = 1.0 if position == 'favorable' else 0.7

    return (
        gamma_c
        / 1.5
        * (alpha1 * alpha2 * sigma_pm0)
        / (eta1 * math.sqrt(fck_t))
        * phi_p
    )


def lpt1_precast(lpt: float) -> float:
    """Calculate the design transmission length lpt1
    for the verification of local stresses at release.

    EN1992-1-1:2023 Eq. (13.6)

    Args:
        lpt (float): Basic transmission length in mm.

    Returns:
        float: Design transmission length lpt1 in mm.
    """
    if lpt < 0:
        raise ValueError(f'lpt must be non-negative. Got lpt={lpt}')

    return 0.8 * lpt


def lpt2_precast(lpt: float) -> float:
    """Calculate the design transmission
        length lpt2 for ultimate limit states.

    EN1992-1-1:2023 Eq. (13.7)

    Args:
        lpt (float): Basic transmission length in mm.

    Returns:
        float: Design transmission length lpt2 in mm.
    """
    if lpt < 0:
        raise ValueError(f'lpt must be non-negative. Got lpt={lpt}')

    return 1.2 * lpt


def ldisp_precast(lpt: float, d: float) -> float:
    """Calculate the dispersion length ldisp where concrete
        stresses have a linear distribution.

        EN1992-1-1:2023 Eq. (13.8)

    Args:
        lpt (float): Basic transmission length in mm.
        d (float): Depth of the member cross-section in mm.

    Returns:
        float: Dispersion length ldisp in mm.

    Raises:
        ValueError: If any input value is negative.
    """
    if lpt < 0 or d < 0:
        raise ValueError(
            f'Input values must be non-negative. Got lpt={lpt}, d={d}'
        )

    return math.sqrt(lpt**2 + d**2)


def lbpd_precast(
    lpt2: float,
    gamma_c: float,
    tendon: Literal['indented_wires', 'stands'],
    fatigue: bool,
    position: Literal['favorable', 'unfavorable'],
    sigma_pd: float,
    sigma_pm_inf: float,
    fck: float,
    phi_p: float,
) -> float:
    """Calculate the total anchorage length lbpd for
        anchoring a tendon at ultimate limit state in precast concrete.

    EN1992-1-1:2023 Eq. (13.9)

    Args:
        lpt2 (float): Upper design value of transmission length in mm.
        gamma_c (float): Partial safety factor for concrete.
        tendon (float): Gradual or sudden release.
        fatigue (bool): True if requires fatigue verification. False otherwise.
        sigma_pm0 (float): Tendon stress just after release in MPa.
        position: favorable or unfavorable.
        sigma_pd (float): Tendon stress at ultimate limit state in MPa.
        sigma_pm_inf (float): Tendon stress after all losses in MPa.
        fck (float): Concrete compressive strength in MPa.
        phi_p (float): Nominal diameter of the tendon in mm.

    Returns:
        float: Total anchorage length lbpd in mm.

    Raises:
        ValueError: If any input value is negative.
    """
    if any(
        param < 0
        for param in [
            lpt2,
            gamma_c,
            sigma_pd,
            sigma_pm_inf,
            fck,
            phi_p,
        ]
    ):
        raise ValueError('All input values must be non-negative.')

    alpha2 = 0.4 if tendon == 'indented_wires' else 0.26
    alpha3 = 1.5 if fatigue else 1.0
    eta1 = 1.0 if position == 'favorable' else 0.7

    part1 = lpt2**2
    part2 = (
        gamma_c
        / 1.5
        * 2
        * alpha2
        * alpha3
        * (sigma_pd - sigma_pm_inf)
        / (eta1 * math.sqrt(fck))
        * phi_p
    )
    return part1 + part2


def sigma_1_Ed_precast(sigma_xEd_y: float, tau_Ed_y: float) -> float:
    """Calculate the principal tensile stress σ1Ed
        at a given distance y from the centroidal axis.

    EN1992-1-1:2023 Eq. (13.11)

    Args:
        sigma_xEd_y (float): Normal stress in the
            longitudinal direction at distance y in MPa.
        tau_Ed_y (float): Shear stress at distance y in MPa.

    Returns:
        float: Principal tensile stress σ1Ed in MPa.

    Raises:
        ValueError: If any input value is negative.

    """
    if sigma_xEd_y < 0 or tau_Ed_y < 0:
        raise ValueError(
            'Stress values must be non-negative. '
            + f' Got sigma_xEd_y={sigma_xEd_y}, tau_Ed_y={tau_Ed_y}'
        )

    return sigma_xEd_y / 2 + math.sqrt((sigma_xEd_y / 2) ** 2 + tau_Ed_y**2)


def tau_Ed_y_precast(VEd: float, Sy: float, In: float, by: float) -> float:
    """Calculate the shear stress τEd in a fiber at distance
        y from the centroidal axis.

    EN1992-1-1:2023 Eq. (13.12)

    Args:
        VEd (float): Shear force in kN.
        Sy (float): First moment of area of the cross-section
            above the fiber at distance y in mm3.
        In (float): Second moment of area of the
            concrete cross-section in mm4.
        by (float): Width of the concrete cross-section
            at distance y from the centroidal axis in mm.

    Returns:
        float: Shear stress τEd MPa.

    Raises:
        ValueError: If any input value is negative.
    """
    if VEd < 0 or Sy < 0 or In < 0 or by < 0:
        raise ValueError(
            'All input values must be non-negative. '
            + f' Got VEd={VEd}, Sy={Sy}, In={In}, by={by}'
        )

    tau_Ed = (VEd * Sy) / (In * by)
    return tau_Ed * 1000  # Convert to MPa


def sigma_1_Rd_precast(f_ctk_005: float, gamma_c: float) -> float:
    """Get the resistant principal stress in precast concrete.

    EN1992-1-1:2023 Eq. (13.10)

    Args:
        sigma_1Ed (float): Principal tensile stress σ1Ed in MPa.
        f_ctk_005 (float): Characteristic tensile
            strength of concrete fctk,0.05 in MPa.
        gamma_c (float): Partial safety factor for concrete.

    Returns:
        bool: True if σ1Ed ≤ fctk,0.05/γC, False otherwise.
    """
    if f_ctk_005 < 0 or gamma_c <= 0:
        raise ValueError(
            'Values must be non-negative, and gamma_c must be positive. '
            + f'Got f_ctk_005={f_ctk_005}, gamma_c={gamma_c}'
        )

    return f_ctk_005 / gamma_c


def vEd_joints_precast(qEd: float, be: float) -> float:
    """Calculate the shear force vEd acting in joints
        between adjacent precast elements per unit length
        of the longitudinal joint.

    EN1992-1-1:2023 Eq. (13.13)

    Args:
        qEd (str): Design value of variable load in N/mm2.
        be (str): Width of the element in mm.

    Returns:
        float: Shear force vEd in N/mm.

    Raises:
        ValueError: If any input value is negative.
    """
    if qEd < 0 or be < 0:
        raise ValueError(
            f'Input values must be non-negative. Got qEd={qEd}, be={be}'
        )

    return qEd * be / 3


def sT_transverse_ribs(
    sL: float,
    lL: float,
    h: float,
    load_type: Literal['residential_snow', 'other'],
) -> float:
    """Calculate the maximum spacing of transverse
        ribs sT based on type of load and geometry.

    EN1992-1-1:2023 Table 13.2

    Args:
        sL (float): Spacing of longitudinal ribs in mm.
        lL (float): Span of longitudinal ribs in mm.
        h (float): Depth of the floor in mm.
        load_type (str): Type of load ('residential_snow' or 'other').

    Returns:
        float: Maximum spacing of transverse ribs sT in mm.

    Raises:
        ValueError: If any input value is negative or
            if load_type is not recognized.
    """
    if sL < 0 or lL < 0 or h < 0:
        raise ValueError(
            f'Input values must be non-negative. Got sL={sL}, lL={lL}, h={h}'
        )

    if load_type == 'residential_snow':
        if sL <= lL / 8:
            return float('inf')  # No requirement for transverse ribs
        return 12 * h

    # If other
    if sL <= lL / 8:
        return 10 * h
    return 8 * h


def As_conn_precast(FEd: float, fyd: float, t: float, h: float) -> float:
    """Calculates the required reinforcement area
        As for transverse tensile stresses near
        connections transmitting compression.

    EN1992-1-1:2023 Eq. (13.14)

    Args:
        FEd (float): Design force in kN.
        fyd (float): Design yield strength of the reinforcement in MPa.
        t (float): support thickness in mm.
        h (float): Thickness of the wall in mm.

    Returns:
        float: Required reinforcement area As in mm2.

    Raises:
        ValueError: If any input is negative.
    """
    if FEd < 0:
        raise ValueError(f'FEd must not be negative. Got {FEd}')
    if fyd < 0:
        raise ValueError(f'fyd must not be negative. Got {fyd}')
    if h < 0:
        raise ValueError(f'h must not be negative. Got {h}')

    return 0.25 * t / h * FEd * 1000 / fyd


def FRd_conn_precast(
    FEd: float,
    h: float,
    fcd: float,
    reinforced: bool = True,
) -> bool:
    """Computes the vertical load resistance for a wall connection.

    EN1992-1-1:2023 Eq. (13.15) and (13.16)

    Args:
        FEd (float or int): Design force in kN.
        h (float or int): Wall thickness in mm.
        fcd (float or int): Design value of concrete
            compressive strength in MPa.
        reinforced (bool): True if the wall is reinforced, False otherwise.

    Returns:
        float: the resistance force in N/mm.

    Raises:
        ValueError: If any input is negative.
    """
    if FEd < 0:
        raise ValueError(f'FEd must not be negative. Got {FEd}')
    if h < 0:
        raise ValueError(f'h must not be negative. Got {h}')
    if fcd < 0:
        raise ValueError(f'fcd must not be negative. Got {fcd}')

    limit_factor = 0.6 if reinforced else 0.5
    return limit_factor * h * fcd


def l_emb_min_precast(hcol: float, MEd: float, NEd: float) -> float:
    """Calculates the minimum embedded length l
        for pocket foundations with smooth or rough surfaces.

    EN1992-1-1:2023 Eq. (13.17) and (13.18)

    Args:
        hcol (float): Largest side of the column's section in mm.
        MEd (float): Design bending moment in kNm.
        NEd (float): Design axial force in kN.

    Returns:
        float: Minimum embedded length l in mm.

    Raises:
        ValueError: If any input is negative.
    """
    if hcol < 0:
        raise ValueError(f'hcol must not be negative. Got {hcol}')
    if MEd < 0:
        raise ValueError(f'MEd must not be negative. Got {MEd}')
    if NEd < 0:
        raise ValueError(f'NEd must not be negative. Got {NEd}')

    ratio = MEd * 1000 / NEd  # Convert to mm

    if ratio <= 0.15 * hcol:
        l_emb = 1.2 * hcol
    elif ratio >= 2.0 * hcol:
        l_emb = 2.0 * hcol
    else:
        # Linear interpolation between the two conditions
        l_emb = 1.2 * hcol + (
            (ratio - 0.15 * hcol) / (2.0 * hcol - 0.15 * hcol)
        ) * (2.0 * hcol - 1.2 * hcol)

    return l_emb
