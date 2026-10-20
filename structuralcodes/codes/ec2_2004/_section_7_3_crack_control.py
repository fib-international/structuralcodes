"""Collection of functions from EUROCODE 1992-1-1:2004
Chapter 7.3 - Crack control.
"""

import math
import typing as t

import numpy as np
import scipy.interpolate


def w_max(
    exposure_class: t.Literal[
        'X0', 'XC1', 'XC2', 'XC3', 'XC4', 'XD1', 'XD2', 'XS1', 'XS2', 'XS3'
    ],
    load_combination: t.Literal['f', 'qp'],
) -> float:
    """Computes the recommended value of the maximum crack width.

    EUROCODE 2 1992-1-1:2004, Table (7.1N).

    Args:
        exposure_class (str): The exposure class. Possible values: X0, XC1,
            XC2, XC3, XC4, XD1, XD2, XS1, XS2, XS3.
        load_combination (str): The characteristic of the load combination.
            Frequent (f), or quasi-permanent (qp).

    Returns:
        float: The maximum recommended value for the crack width wmax in mm.

    Raises:
        ValueError: if not valid exposure_class or load_combination values.
    """
    _load_combination = load_combination.lower().strip()
    _exposure_class = exposure_class.upper().strip()
    if _load_combination == 'f':
        if _exposure_class in ('X0', 'XC1'):
            return 0.2
        if _exposure_class in ('XC2', 'XC3', 'XC4'):
            return 0.2
    if _load_combination == 'qp':
        if _exposure_class in ('X0', 'XC1'):
            return 0.4
        if _exposure_class in (
            'XC2',
            'XC3',
            'XC4',
            'XD1',
            'XD2',
            'XS1',
            'XS2',
            'XS3',
        ):
            return 0.3
        raise ValueError(
            f'{exposure_class} is not a valid value for exposure_class.'
            + ' Please enter one of the following: X0, XC1, XC2, XC3, XC4, XD1'
            + ',XD2, XS1, XS2, XS3'
        )
    raise ValueError(
        f'{load_combination} is not a valid value for load_combination.'
        + 'Please enter "f" for frequent load combination or "qp" for'
        + 'quasi-permanent load combination.'
    )


def As_min(
    A_ct: float, sigma_s: float, fct_eff: float, k: float, kc: float
) -> float:
    """Computes the minimum area of reinforcing steel within the tensile zone
    for control of cracking areas.

    EUROCODE 2 1992-1-1:2004, Eq. (7.1).

    Args:
        A_ct (float): is the area of concrete within the tensile zone in mm2.
            The tensile zone is that part of the section which is calculated
            to be in tension just before the formation of the first crack.
        sigma_s (float): is the absolute value of the maximum stress in MPa
            permitted in the reinforcement immediately after the formation
            of the crack. This may be taken as theyield strength of the
            reinforcement, fyk. A lower value may, however, be needed to
            satisfy the crack width limits according to the maximum
            bar size of spacing (see 7.3.3 (2)).
        fct_eff (float): is the mean value of the tensile strength in MPa of
            the concrete effective at the time when the cracks may first be
            expected to occur: fct,eff=fct or lower (fct(t)), is cracking
            is expected earlier than 28 days.
        k (float): is the coefficient which allow for the effect of
            non-uniform self-equilibrating stresses, which lead to a
            reduction of restraint forces.
            k=1 for webs w<=300mm or flanges widths less than 300mm
            k=0.65 for webs w>=800mm or flanges with widths greater than 800mm
            Intermediate values may be interpolated.
        kc (float): is a coefficient which takes account of the stress
            distribution within the section immediately prior to cracking and
            the change of the lever arm.

    Returns:
        float: The minimum area of reinforcing steel within the tensile zone in
        mm2.

    Raises:
        ValueError: if k value is not between 0.65 and 1 or kc is not
            larger than 0 and lower than 1.
    """
    fct_eff = abs(fct_eff)

    if A_ct <= 0:
        raise ValueError(f'A_ct={A_ct} must be larger than 0')
    if sigma_s < 0:
        raise ValueError(f'sigma_s={sigma_s} must be equal or larger than 0')
    if k < 0.65 or k > 1.0:
        raise ValueError(f'k={k} must be between 0.65 and 1')
    if kc > 1 or kc < 0:
        raise ValueError(f'kc={kc} must be lower than 1 and larger than 0')

    return kc * k * fct_eff * A_ct / sigma_s


def k(h: float) -> float:
    """Is the coefficient which allow for the effect of
    non-uniform self-equilibrating stresses, which lead to a
    reduction of restraint forces.
    k=1 for webs w<=300mm or flanges widths less than 300mm
    k=0.65 for webs w>=800mm or flanges with widths greater than 800mm.

    EUROCODE 2 1992-1-1:2004, Eq. (7.1).

    Args:
        h (float): flange length or flange width in mm

    Returns:
        float: k coefficient value.

    Raises:
        ValueError: if h is less than 0
    """
    if h < 0:
        raise ValueError(f'h={h} cannot be less than 0mm')
    if h <= 300:
        return 1
    if h < 800:
        interpol = scipy.interpolate.interp1d((300, 800), (1, 0.65))
        return interpol(h)
    return 0.65


def kc_tension() -> float:
    """Computes the coefficient which takes account of the stress
    distribution within the section immediately prior to cracking and
    the change of the lever arm in pure tension.

    EUROCODE 2 1992-1-1:2004, Eq. (7.1).

    Returns:
        float: Value of the kc coefficient in pure tension.
    """
    return 1


def kc_rect_area(h: float, b: float, fct_eff: float, N_ed: float) -> float:
    """Computes the coefficient which takes account of the stress distribution
    within the section immediately prior to cracking and the change of the
    lever arm for bending+axial combination in rectangular sections and webs of
    box sections and T-sections.

    EUROCODE 2 1992-1-1:2004, Eq. (7.2).

    Args:
        h (float): heigth of the element in mm
        b (float): width of the element in mm
        fct_eff (float): is the mean value of the tensile strength in MPa of
            the concrete effective at the time when the cracks may first be
            expected to occur: fct,eff=fct or lower (fct(t)), is cracking is
            expected earlier than 28 days.
        N_ed (str): axial force at the serviceability limit state acting on the
            part of the cross-section under consideration (compressive force
            positive). n_ed should be determined considering the characteristic
            values of prestress and axial forces under the relevant combination
            of actions.

    Returns:
        float: Value of the kc coefficient.

    Raises:
        ValueError: If h or b are less than 0.
    """
    if h < 0:
        raise ValueError(f'h={h} should be larger than 0mm')
    if b < 0:
        raise ValueError(f'b={b} should be larger than 0mm')

    h_s = min(h, 1000)
    k1 = 1.5 if N_ed >= 0 else 2 * h_s / 3 / h
    s_concrete = N_ed * 1000 / b / h
    h_ratio = h / h_s
    return min(max(0.4 * (1 - s_concrete / k1 / h_ratio / fct_eff), 0), 1)


def kc_flanges_area(f_cr: float, A_ct: float, fct_eff: float) -> float:
    """Computes the coefficient which takes account of the stress
    distribution within the section immediately prior to cracking and
    the change of the lever arm for bending+axial combination
    in rectangular sections for flanges of box sections and T-sections.

    EUROCODE 2 1992-1-1:2004, Eq. (7.3).

    Args:
        f_cr: is the absolute value in kN of the tensile force within the
            flange immediately prior to cracking due to cracking moment
            calculated with fct,eff.
        A_ct (float): is the area of concrete within the tensile zone in mm2.
            The tensile zone is that part of the section which is calculated
            to be in tension just before the formation of the first crack.
        fct_eff (float): is the mean value of the tensile strength in MPa of
            the concrete effective at the time when the cracks may first be
            expected to occur: fct,eff=fct or lower (fct(t)), is cracking
            is expected earlier than 28 days.

    Returns:
        float: Value of the kc coefficient.

    Raises:
        ValueError: If A_ct is less than 0mm2.
    """
    f_cr = abs(f_cr)
    return max(0.9 * f_cr * 1000 / A_ct / fct_eff, 0.5)


def xi1(xi: float, phi_p: float, phi_s: float) -> float:
    """Computes the adjusted ratio of bond strength taking into account
    the different diameters of prestressing and reinforcing steel.

    EUROCODE 2 1992-1-1:2004, Eq. (7.5).

    Args:
        xi (float): ratio of bond strength of prestressing and reinforcing
            steel, according to Table 6.2 in 6.8.2.
        phi_p (float): largest bar diameter in mm of reinforcing steel.
            Equal to 0 if only prestressing is used in control cracking.
        phi_s (float): equivalent diameter in mm of tendon acoording
            to 6.8.2.

    Returns:
        float: With the value of the ratio.

    Raises:
        ValueError: If diameters phi_s or phi_p are lower than 0. If ratio of
            bond strength xi is less than 0.15 or larger than 0.8.
    """
    if phi_p <= 0:
        raise ValueError(f'phi_p={phi_p} cannot be less than 0')
    if phi_s < 0:
        raise ValueError(f'phi_s={phi_s} cannot be less than 0')
    if xi < 0.15:
        raise ValueError(f'The minimum value for xi={xi} is 0.15')
    if xi > 0.8:
        raise ValueError(f'The maximum value for xi={xi} is 0.8')

    return ((xi * phi_s / phi_p) ** 0.5) if phi_s > 0 else xi**0.5


def hc_eff(h: float, d: float, x: float) -> float:
    """Returns the effective height of concrete in tension surrounding the
    reinforcement or prestressing tendons.

    EUROCODE 2 1992-1-1:2004, Section (7.3.2-3).

    Args:
        h (float): total depth of the element in mm.
        d (float): distance in mm to the level of the steel centroid.
        x (float): distance in mm to the zero tensile stress line.

    Returns:
        float: The effective height in mm.

    Raises:
        ValueError: If any of h, d or x is lower than zero.
        ValueError: If d is greater than h.
        ValueError: If x is greater than h.
    """
    if h < 0:
        raise ValueError(f'h={h} cannot be less than 0')
    if d < 0:
        raise ValueError(f'd={d} cannot be less than 0')
    if x < 0:
        raise ValueError(f'x={x} cannot be less than zero')
    if d > h:
        raise ValueError(f'd={d} cannot be larger than h={h}')
    if x > h:
        raise ValueError(f'x={x} cannot be larger than h={h}')

    return min(2.5 * (h - d), (h - x) / 3, h / 2)


def As_min_p(
    A_ct: float,
    sigma_s: float,
    fct_eff: float,
    k: float,
    kc: float,
    Ap: float,
    phi_s: float,
    phi_p: float,
    xi: float,
    delta_s: float,
) -> float:
    """Computes the minimum area of reinforcing steel within the tensile zone
    for control of cracking areas in addition with bonded tendons.

    EUROCODE 2 1992-1-1:2004, Eq. (7.1).

    Args:
        A_ct (float): is the area of concrete within the tensile zone in mm2.
            The tensile zone is that part of the section which is calculated to
            be in tension just before the formation of the first crack.
        sigma_s (float): is the absolute value of the maximum stress in MPa
            permitted in the reinforcement immediately after the formation of
            the crack. This may be taken as the yield strength of the
            reinforcement, fyk. A lower value may, however, be needed to
            satisfy the crack width limits according to the maximum bar size of
            spacing (see 7.3.3 (2)).
        fct_eff (float): is the mean value of the tensile strength in MPa of
            the concrete effective at the time when the cracks may first be
            expected to occur: fct,eff=fct or lower (fct(t)), is cracking is
            expected earlier than 28 days.
        k (float): is the coefficient which allow for the effect of non-
            uniform self-equilibrating stresses, which lead to a reduction of
            restraint forces. k=1 for webs w<=300mm or flanges widths less than
            300mm. k=0.65 for webs w>=800mm or flanges with widths greater than
            800mm. Intermediate values may be interpolated.
        kc (float): is a coefficient which takes account of the stress
            distribution within the section immediately prior to cracking and
            the change of the lever arm.
        Ap (float): is the area in mm2 of pre or post-tensioned tendons within
            ac_eff.
        phi_s (float): largest bar diameter in mm of reinforcing steel. Equal
            to 0 if only prestressing is used in control cracking.
        phi_p (float): equivalent diameter in mm of tendon according to 6.8.2.
        xi (float): ratio of bond strength of prestressing and reinforcing
            steel, according to Table 6.2 in 6.8.2.
        delta_s (float): stress variation in MPa in prestressing tendons from
            the state of zero strain of the concrete at the same level.

    Returns:
        float: The minimm area of reinforcing steel within the tensile zone in
        mm2.

    Raises:
        ValueError: If k value is not between 0.65 and 1 or kc is not larger
        than 0 and lower than 1. If diameters phi_s or phi_p are lower than 0.
        If ratio of bond xi strength is less than 0.15 or larger than 0.8. If
        stress variation incr_stress is less than 0.
    """
    fct_eff = abs(fct_eff)

    if Ap < 0:
        raise ValueError(f'Ap={Ap} cannot be less than 0')
    if delta_s < 0:
        raise ValueError(f'delta_s={delta_s} cannot be less than 0')
    if A_ct <= 0:
        raise ValueError(f'A_ct={A_ct} must be larger than 0')
    if sigma_s < 0:
        raise ValueError(f'sigma_s={sigma_s} must be equal or larger than 0')
    if k < 0.65 or k > 1.0:
        raise ValueError(f'k={k} must be between 0.65 and 1')
    if kc > 1 or kc < 0:
        raise ValueError(f'kc={kc} must be lower than 1 and larger than 0')

    a1 = kc * k * fct_eff * A_ct
    e1 = xi1(xi, phi_p, phi_s)
    a2 = e1 * Ap * delta_s
    a = a1 - a2

    return a / sigma_s


def As_min_2(
    wk: float,
    sigma_s: float,
    fct_eff: float,
    h_cr: float,
    h: float,
    d: float,
    delta_s: float = 0,
    kc: t.Optional[float] = None,
) -> t.Tuple[float, float]:
    """Computes the minimum area of reinforcing steel within the tensile zone
    for control of cracking areas.

    EUROCODE 2 1992-1-1:2004, Table (7.2N), Table (7.3N).

    Args:
        wk (float): the characteristic crack width value in mm.
        sigma_s (float): the steel stress value in MPa under the relevant
            combination of actions.
        fct_eff (float): is the mean value of the tensile strength in MPa of
            the concrete effective at the time when the cracks may first be
            expected to occur: fct,eff=fct or lower (fct(t)), is cracking is
            expected earlier than 28 days.
        h_cr (float): is the depth of the tensile zone immediately prior to
            cracking, considering the characteristic values of prestress and
            axial forces under the quasi-permanent combination of actions.
        h (float): the overall depth of the section in mm.
        d (float): is the effective depth to the centroid of the outer layer of
            the reinforcement.

    Keyword Args:
        delta_s (float, optional): value of prestressed stress in MPa if
            applicable.
        kc (float, optional): is a coefficient which takes account of the
            stress distribution within the section immediately prior to
            cracking and the change of the lever arm in a bending section. None
            for pure tensile uniform axial section.

    Returns:
        tuple(float, float): With the value of the maximum bar diameters in mm
        in the first position and the maximum bar spacing in mm in the second
        position.

    Raises:
        ValueError: If wk, fct_eff, h_cr, h or d are less than 0.
        ValueError: If kc is not between 0 and 1.
        ValueError: If combination of wk and stress values are out of scope.
    """
    if wk < 0:
        raise ValueError(f'wk={wk} cannot be less than 0')
    if fct_eff < 0:
        raise ValueError(f'fct_eff={fct_eff} is less than 0')
    if h_cr < 0:
        raise ValueError(f'h_cr={h_cr} is less than 0')
    if h < 0:
        raise ValueError(f'h={h} is less than 0')
    if d < 0:
        raise ValueError(f'd={d} is less than 0')
    if kc is not None and (kc < 0 or kc > 1):
        raise ValueError(f'kc={kc} is not between 0 and 1')

    s = sigma_s - delta_s
    if s <= 0:
        return (0, 0)

    x = (0.4, 0.3, 0.2)
    y_phi = (160, 200, 240, 280, 320, 360, 400, 450)
    y_spa = (160, 200, 240, 280, 320, 360)
    phi_s_v = (
        40,
        32,
        25,
        32,
        25,
        16,
        20,
        16,
        12,
        16,
        12,
        8,
        12,
        10,
        6,
        10,
        8,
        5,
        8,
        6,
        4,
        6,
        5,
        None,
    )
    spa_v = (
        300,
        300,
        200,
        300,
        250,
        150,
        250,
        200,
        100,
        200,
        150,
        50,
        150,
        100,
        None,
        100,
        50,
        None,
    )

    points_phi = np.array(np.meshgrid(y_phi, x)).T.reshape(-1, 2)
    points_spa = np.array(np.meshgrid(y_spa, x)).T.reshape(-1, 2)
    xi = (s, wk)

    phi_star = float(
        scipy.interpolate.griddata(points_phi, phi_s_v, xi, method='linear')
    )
    if kc is not None:
        phi = phi_star * (fct_eff / 2.9) * kc * h_cr / (2 * (h - d))
    else:
        phi = phi_star * (fct_eff / 2.9) * h_cr / (8 * (h - d))

    spa = float(
        scipy.interpolate.griddata(points_spa, spa_v, xi, method='linear')
    )

    if math.isnan(phi) or math.isnan(spa):
        raise ValueError('Combination of wk or stress values out of scope')

    return phi, spa


def alpha_e(Es: float, Ecm: float) -> float:
    """Compute the ratio between the steel and mean concrete elastic modules.

    EUROCODE 2 1992-1-1:2004, Section 7.3.4-2

    Args:
        Es (float): Steel elastic modulus in MPa.
        Ecm (float): Concrete mean elastic modulus in MPa.

    Returns:
        float: Ratio between modules.

    Raises:
        ValueError: If any of es or ecm is lower than 0.
    """
    if Es < 0:
        raise ValueError(f'Es={Es} cannot be less than 0')
    if Ecm < 0:
        raise ValueError(f'Ecm={Ecm} cannot be less than 0')

    return Es / Ecm


def rho_p_eff(As: float, xi1: float, Ap: float, Ac_eff: float) -> float:
    """Effective bond ratio between areas.

    EUROCODE 2 1992-1-1:2004, Eq. (7.10).

    Args:
        As (float): Steel area in mm2.
        xi1 (float): The adjusted ratio of bond according to expression (7.5).
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


def kt(load_type: t.Literal['short', 'long']) -> float:
    """Returns the kt factor dependent on the load duration for the crack width
    calculation.

    Args:
        load_type (str): The load type, 'short' for term loading, 'long' for
            long term loading.

    Returns:
        float: With the kt factor.

    Raises:
        ValueError: If load_type is not 'short' and not 'long'.
    """
    if not isinstance(load_type, str):
        raise TypeError

    load_type = load_type.lower().strip()
    if load_type not in ('short', 'long'):
        raise ValueError(
            f'load_type={load_type} can only have "short" or "long" as a value'
        )

    return 0.6 if load_type == 'short' else 0.4


def eps_sm_eps_cm(
    sigma_s: float,
    alpha_e: float,
    rho_p_eff: float,
    kt: float,
    fct_eff: float,
    Es: float,
) -> float:
    """Returns the strain difference (epsilon_sm - epsilon_cm) needed to
    compute the crack width. esm is the mean strain in the reinforcement under
    the relevant combination of loads of imposed deformations and taking into
    account the effects of tension stiffening. Only the additional tensile
    strain beyond the state of zero strain of the concrete is considered.
    epsilon_cm is the mean strain in the concrete between the cracks.

    EUROCODE 2 1992-1-1:2004, Eq. (7.9).

    Args:
        sigma_s (float): Is the stress in MPa in the tension reinforcement
            assuming a cracked section. For pretensioned members, s_steel may
            be replaced by increment of s_steel stress variation in
            prestressing tendons from the state of zero strain of the concrete
            at the same level.
        alpha_e (float): Is the ratio Es/Ecm.
        rho_p_eff (float): Effective bond ratio between areas given by Eq.
            (7.10).
        kt (float): Is a factor dependent on the load duration.
        fct_eff (float): Is the mean value of the tensile strength in MPa of
            the concrete effective at the time when the cracks may first be
            expected to occur: fct_eff=fctm or fctm(t) if crack is expected
            earlier than 28 days.
        Es (float): Steel elastic modulus in MPa.

    Returns:
        float: The strain difference between concrete and steel.

    Raises:
        ValueError: If any sigma_s, alpha_e, rho_p_eff, fct_eff or Es is less
            than 0.
        ValueError: if kt is not 0.6 and not 0.4.
    """
    if sigma_s < 0:
        raise ValueError(f'sigma_s={sigma_s} cannot be less than 0')
    if alpha_e < 0:
        raise ValueError(f'alpha_e={alpha_e} cannot be less than 0')
    if rho_p_eff < 0:
        raise ValueError(f'rho_p_eff={rho_p_eff} cannot be less than 0')
    if fct_eff < 0:
        raise ValueError(f'fct_eff={fct_eff} cannot be less than 0')
    if Es < 0:
        raise ValueError(f'Es={Es} cannot be less than 0')
    if kt not in (0.6, 0.4):
        raise ValueError(f'kt={kt} can only take as values 0.4 and 0.6')

    min_val = 0.6 * sigma_s / Es

    a = 1 + alpha_e * rho_p_eff
    b = kt * fct_eff / rho_p_eff * a
    c = (sigma_s - b) / Es

    return max(c, min_val)


def w_spacing(c: float, phi: float) -> float:
    """Computes the distance threshold from which the maximum crack spacing is
    constant.

    EUROCODE 2 1992-1-1:2004, Sect. (7.3.4-3).

    Args:
        c (float): Cover of the longitudinal reinforcement in mm.
        phi (float): Is the bar diameter in mm. Where mixed bar diameters used,
            then it should be replaced for an equivalent bar diameter.

    Returns:
        float: Threshold distance in mm.

    Raises:
        ValueError: If any of c or phi is less than 0.
    """
    if c < 0:
        raise ValueError(f'c={c} cannot be less than 0')
    if phi < 0:
        raise ValueError(f'phi={phi} cannot be less than 0')

    return 5 * (c + phi / 2)


def phi_eq(n1: int, n2: int, phi1: float, phi2: float) -> float:
    """Computes the equivalent diameter. For a section with n1 bars of diameter
    phi1 and n2 bars of diameter phi2.

    EUROCODE 2 1992-1-1:2004, Sect. (7.12).

    Args:
        n1 (int): Number of bars with diameter phi1.
        n2 (int): Number of bars with diameter phi2.
        phi1 (float): Diameter of n1 bars in mm.
        phi2 (float): Diamater of n2 bars in mm.

    Returns:
        float: The equivalent diameter in mm.

    Raises:
        ValueError: If any of n1 or n2 is less than 0.
        ValueError: If any of phi1 or phi2 is less than 0.
        TypeError: If any of n1 or n2 is not an integer.
    """
    if n1 < 0:
        raise ValueError(f'n1={n1} cannot be less than 0')
    if not isinstance(n1, int):
        raise TypeError(f'n1={n1} needs to be an integer value')
    if n2 < 0:
        raise ValueError(f'n2={n2} cannot be less than 0')
    if not isinstance(n2, int):
        raise TypeError(f'n2={n2} needs to be an integer value')
    if phi1 < 0:
        raise ValueError(f'phi1={phi1} cannot be less than 0')
    if phi2 < 0:
        raise ValueError(f'phi2={phi2} cannot be less than 0')

    a = n1 * phi1**2 + n2 * phi2**2
    b = n1 * phi1 + n2 * phi2
    return a / b


def k1(bond_type: t.Literal['bond', 'plain']) -> float:
    """Get the k1 coefficient which takes account of the bond properties of the
    bounded reinforcement.

    EUROCODE 2 1992-1-1:2004, Eq. (7.11-k1).

    Args:
        bond_type (str): The bond property of the reinforcement. High bond bars
            (bond), or bars with an effectively plain surface (plain).

    Returns:
        float: Value of the k1 coefficient.

    Raises:
        ValueError: If bond_type is neither 'bond' nor 'plain'.
        TypeError: If bond_type is not an str.
    """
    if not isinstance(bond_type, str):
        raise TypeError(f'bond_type={bond_type} is not an str')

    bond_type = bond_type.lower().strip()
    if bond_type not in ('bond', 'plain'):
        raise ValueError(
            f'bond_type={bond_type} can only have "bond" or "plain" as values'
        )

    return 0.8 if bond_type == 'bond' else 1.6


def k2(eps_r: float) -> float:
    """Computes a coeff. which takes into account the distribution of strain.

    EUROCODE 2 1992-1-1:2004, Eq. (7.13).

    Args:
        eps_r (float): ratio epsilon_2/epsilon_1 where epsilon_1 is the greater
            and epsilon_2 is the lesser strain at the boundaries of the section
            considered, assessed on the basis of a cracked section. epsilon_r=0
            for bending and epsilon_r=1 for pure tension.

    Returns:
        float: The k2 coefficient value.

    Raises:
        ValueError: If eps_r is not between 0 and 1.
    """
    if eps_r < 0 or eps_r > 1:
        raise ValueError(f'eps_r={eps_r} must be between 0 and 1')

    return (1 + eps_r) / 2


def k3():
    """Returns the k3 coefficient for computing sr_max.

    Returns:
        float: Value for the coefficient.
    """
    return 3.4


def k4():
    """Returns the k4 coefficient for computing sr_max.

    Returns:
        float: Value for the coefficient.
    """
    return 0.425


K_COEFFICIENTS = {'k1': k1, 'k2': k2, 'k3': k3, 'k4': k4}


def sr_max_close(
    c: float,
    phi: float,
    rho_p_eff: float,
    k1: float,
    k2: float,
    k3: t.Optional[float] = None,
    k4: t.Optional[float] = None,
) -> float:
    """Computes the maximum crack spacing in cases where bonded reinforcement
    is fixed at reasonably close centres within the tension zone
    (w_spacing<=5(c+phi/2)).

    EUROCODE 2 1992-1-1:2004, Eq. (7.11).

    Args:
        c (float): Is the cover in mm of the longitudinal reinforcement.
        phi (float): Is the bar diameter in mm. Where mixed bar diameters used,
            then it should be replaced for an equivalent bar diameter.
        rho_p_eff (float): Effective bond ratio between areas given by Eq.
            (7.10).
        k1 (float): Coefficient that takes into account the bound properties
            of the bonded reinforcement.
        k2 (float): Coefficient that takes into account the distribution of of
            the strain.

    Keyword Args:
        k3 (float, optional): Coefficient from the National Annex. If not
            specified then k3=3.4.
        k4 (float): Coefficient from the National Annex. If not specified then
            k4=0.425.

    Returns:
        float: The maximum crack spaing in mm.

    Raises:
        ValueError: If one or more of c, phi, rho_p_eff, k3 or k4
            is lower than zero.
        ValueError: If k1 is not 0.8 or 1.6.
        ValueError: If k2 is not between 0.5 and 1.0.
    """
    if k3 is None:
        k3 = K_COEFFICIENTS['k3']()
    if k4 is None:
        k4 = K_COEFFICIENTS['k4']()

    if c < 0:
        raise ValueError(f'c={c} cannot be less than zero')
    if phi < 0:
        raise ValueError(f'phi={phi} cannot be less than zero')
    if rho_p_eff < 0:
        raise ValueError(f'rho_p_eff={rho_p_eff} cannot be less than zero')
    if k3 < 0:
        raise ValueError(f'k3={k3} cannot be less than zero')
    if k4 < 0:
        raise ValueError(f'k4={k4} cannot be less than zero')
    if k1 not in (0.8, 1.6):
        raise ValueError(f'k1={k1} can only take as values 0.8 and 1.6')
    if k2 < 0.5 or k2 > 1:
        raise ValueError(f'k2={k2} is not between 0.5 and 1.0')

    return k3 * c + k1 * k2 * k4 * phi / rho_p_eff


def sr_max_far(h: float, x: float) -> float:
    """Computes the maximum crack spacing in cases where bonded reinforcement
    exceeds (w_spacing>5(c+phi/2)) or where there is no bonded reinforcement at
    all.

    EUROCODE 2 1992-1-1:2004, Eq. (7.14).

    Args:
        h (float): Total depth of the beam in mm.
        x (float): Distance to non tension area of the element mm.

    Returns:
        float: Maximum crack spacing in mm.

    Raises:
        ValueError: If one of h or x is less than zero.
        ValueError: If x is greater than h.
    """
    if x < 0:
        raise ValueError(f'x={x} cannot be less than zero')
    if h < 0:
        raise ValueError(f'h={h} cannot be less than zero')
    if x > h:
        raise ValueError(f'x={x} cannot be larger than h={h}')

    return 1.3 * (h - x)


def sr_max_theta(sr_max_y: float, sr_max_z: float, theta: float) -> float:
    """Computes the crack spacing sr_max when there is an angle between the
    angle of principal stress and the direction of the reinforcement, for
    members in two orthogonal directions, that is significant (> 15 degrees).

    EUROCODE 2 1992-1-1:2004, Eq. (7.15).

    Args:
        sr_max_y (float): Crack spacing in mm in the y-direction.
        sr_max_z (float): Crack spacing in mm in the z-direction.
        theta (float): Angle in radians between the reinforcement in the
            y-direction and the direction of the principal tensile stress.

    Returns:
        float: The crack spacing in mm.

    Raises:
        ValueError: If sr_max_y or sr_max_z is negative.
        ValueError: If theta is not between 0 and pi/2.
    """
    if sr_max_y < 0:
        raise ValueError(f'sr_max_y={sr_max_y} cannot be less than zero')
    if sr_max_z < 0:
        raise ValueError(f'sr_max_z={sr_max_z} cannot be less than zero')

    a = math.cos(theta) / sr_max_y
    b = math.sin(theta) / sr_max_z
    return 1 / (a + b)


def wk(sr_max: float, eps_sm_eps_cm: float) -> float:
    """Computes the crack width.

    EUROCODE 2 1992-1-1:2004, Eq. (7.8).

    Args:
        sr_max (float): The maximum crack length spacing in mm.
        eps_sm_eps_cm (float): the difference between the mean strain in the
            reinforcement under relevant combination of loads, including the
            effect of imposed deformations and taking into account tension
            stiffening and the mean strain in the concrete between cracks.

    Returns:
        float: Crack width in mm.

    Raises:
        ValueError: If any of sr_max or esm_ecm is less than zero.
    """
    if sr_max < 0:
        raise ValueError(f'sr_max={sr_max} cannot be less than zero')
    if eps_sm_eps_cm < 0:
        raise ValueError(
            f'eps_sm_eps_cm={eps_sm_eps_cm} cannot be less than zero'
        )

    return sr_max * eps_sm_eps_cm
