"""Collection of functions from EUROCODE 1992-1-1:2004
Chapter 7.3 - Crack control"""
import scipy.interpolate


def w_max(exposure_class: str, load_combination: str) -> float:
    """Computes the recomended value of the maximum crack width.

    EUROCODE 2 1992-1-1:2004, Table (7.1N)

    Args:
        exposure_class (str): The exposure class.
            Possible values: X0, XC1, XC2, XC3, XC4, XD1, XD2, XS1, XS2, XS3
        load_combination (str):
            - f: for frequent load combination
            - qp: for quasi-permanent load combination

    Returns:
        float: The maximum recommended value for the crack width wmax in mm.

    Raises:
        ValueError: if not valid exposure_class or load_combination values.
    """
    _load_combination = load_combination.lower()
    _exposure_class = exposure_class.upper()
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


def crack_min_steel_area(
    a_ct: float, s_steel: float, fct_eff: float, k: float, kc: float
) -> float:
    """Computes the minimum area of reinforcing steel within the tensile zone
    for control of cracking areas

    EUROCODE 2 1992-1-1:2004, Eq. (7.1)

    Args:
        a_ct (float): is the area of concrete within the tensile zone in mm2.
            The tensile zone is that parg of the section which is calculated
            to be in tension just before the formation of the first crack.
        s_steel (float): is the absolute value of the maximum stress in MPa
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
            reduction of restraint forces. Use 'k_crack_min_steel_area'
            to compute it
            k=1 for webs w<=300mm or flanges widths less than 300mm
            k=0.65 for webs w>=800mm or flanges with widths greater than 800mm
            Intermediate values may be interpolated.
        kc (float): is a coefficient which takes account of the stress
            distribution within the section immediately prior to cracking and
            the change of the lever arm.

    Returns:
        float: the minimm area of reinforcing steel within the tensile
            zone in mm2.

    Raises:
        ValueError: if k value is not between 0.65 and 1 or kc is not
            larger than 0 and lower than 1.
    """
    fct_eff = abs(fct_eff)

    if a_ct <= 0:
        raise ValueError(f'a_ct={a_ct} must be larger than 0')
    if s_steel < 0:
        raise ValueError(f's_steel={s_steel} must be equal or larger than 0')
    if k < 0.65 or k > 1.0:
        raise ValueError(f'k={k} must be between 0.65 and 1')
    if kc > 1 or kc < 0:
        raise ValueError(f'kc={kc} must be lower than 1 and larger than 0')

    return kc * k * fct_eff * a_ct / s_steel


def k_crack_min_steel_area(h: float) -> float:
    """Is the coefficient which allow for the effect of
    non-uniform self-equilibrating stresses, which lead to a
    reduction of restraint forces. Use 'k_crack_min_steel_area'
    to compute it
    k=1 for webs w<=300mm or flanges widths less than 300mm
    k=0.65 for webs w>=800mm or flanges with widths greater than 800mm

    EUROCODE 2 1992-1-1:2004, Eq. (7.1)

    Args:
        h (float): flange length or flange width in mm

    Returns:
        float: k coefficient value

    Raises:
        ValueError: if h is less than 0
    """
    if h < 0:
        raise ValueError(f'h={h} cannot be less than 0mm')
    if h <= 300:
        return 1
    if h < 800:
        interpol = scipy.interpolate.interp1d((300, 800), (1, 0.65))
        return (float)(interpol(h))
    return 0.65


def kc_crack_min_steel_area_pure_tension() -> float:
    """Computes the coefficient which takes account of the stress
    distribution within the section immediately prior to cracking and
    the change of the lever arm in pure dtension.

    EUROCODE 2 1992-1-1:2004, Eq. (7.1)

    Returns:
        float: value of the kc coefficient in pure tension
    """
    return 1


def kc_crack_min_steel_area_rectangular(
    h: float, b: float, fct_eff: float, n_ed: float
) -> float:
    """Computes the coefficient which takes account of the stress
    distribution within the section immediately prior to cracking and
    the change of the lever arm for bending+axial combination
    in rectangular sections and webs of box sections and T-sections.

    EUROCODE 2 1992-1-1:2004, Eq. (7.2)

    Args:
        h (float): heigth of the element in mm
        b (float): width of the element in mm
        fct_eff (float): is the mean value of the tensile strength in MPa of
            the concrete effective at the time when the cracks may first be
            expected to occur: fct,eff=fct or lower (fct(t)), is cracking
            is expected earlier than 28 days.
        n_ed (str): axial force at the serviceability limit state acting on
            the part of the cross-section under consideration (compressive
            force positive). n_ed should be determined considering the
            characteristic values of prestress and axial forces under the
            relevant combination of actions

    Returns:
        float: value of the kc coefficient

    Raises:
        ValueError: is h or b are less than 0
    """
    if h < 0:
        raise ValueError(f'h={h} should be larger than 0mm')
    if b < 0:
        raise ValueError(f'b={b} should be larger than 0mm')

    h_s = min(h, 1000)
    k1 = 1.5 if n_ed >= 0 else 2 * h_s / 3 / h
    s_concrete = n_ed * 1000 / b / h
    h_ratio = h / h_s
    return min(max(0.4 * (1 - s_concrete / k1 / h_ratio / fct_eff), 0), 1)


def kc_crack_min_steel_area_flanges(
    f_cr: float, a_ct: float, fct_eff: float
) -> float:
    """Computes the coefficient which takes account of the stress
    distribution within the section immediately prior to cracking and
    the change of the lever arm for bending+axial combination
    in rectangular sections for flanges of box sections and T-sections.

    EUROCODE 2 1992-1-1:2004, Eq. (7.3)

    Args:
        f_cr: is the absolute value in kN of the tensile force within the
            flange immediately prior to cracking due to cracking moment
            calculated with fct,eff
        a_ct (float): is the area of concrete within the tensile zone in mm2.
            The tensile zone is that part of the section which is calculated
            to be in tension just before the formation of the first crack.
        fct_eff (float): is the mean value of the tensile strength in MPa of
            the concrete effective at the time when the cracks may first be
            expected to occur: fct,eff=fct or lower (fct(t)), is cracking
            is expected earlier than 28 days.

    Returns:
        float: value of the kc coefficient

    Raises:
        ValueError: is a_ct is less than 0mm2
    """
    f_cr = abs(f_cr)
    return max(0.9 * f_cr * 1000 / a_ct / fct_eff, 0.5)


def crack_min_steel_area_with_prestresed_tendons(
    a_ct: float,
    s_steel: float,
    fct_eff: float,
    k: float,
    kc: float,
    ap: float,
    d_steel: float,
    d_press: float,
    e: float,
    incr_stress: float,
) -> float:
    """Computes the minimum area of reinforcing steel within the tensile zone
    for control of cracking areas in addition with bonded tendons

    EUROCODE 2 1992-1-1:2004, Eq. (7.1)

    Args:
        a_ct (float): is the area of concrete within the tensile zone in mm2.
            The tensile zone is that part of the section which is calculated
            to be in tension just before the formation of the first crack.
        s_steel (float): is the absolute value of the maximum stress in MPa
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
            reduction of restraint forces. Use 'k_crack_min_steel_area'
            to compute it
            k=1 for webs w<=300mm or flanges widths less than 300mm
            k=0.65 for webs w>=800mm or flanges with widths greater than 800mm
            Intermediate values may be interpolated.
        kc (float): is a coefficient which takes account of the stress
            distribution within the section immediately prior to cracking and
            the change of the lever arm.
        ac_eff (float): is the effective area in mm2 of concrete in tension
            surrounding or prestressing tendons if depth hc,ef
        ap (float): is the area in mm2 of pre or post-tensioned tendons
            within ac_eff
        d_steel (float): largest bar diameter in mm of reinforcing steel.
            Equal to zero if only prestressing is used in control cracking
        d_press (float): equivalent diameter in mm of tendon acoording
            to 6.8.2
        e (float): ratio of bond strength of prestressing and reinforcing
            steel, according to Table 6.2 in 6.8.2
        incr_stress (float): stress variation in MPa in prestressing tendons
            from the state of zero strain of the concrete at the same level

    Returns:
        float: the minimm area of reinforcing steel within the tensile
            zone in mm2.

    Raises:
        ValueError: if k value is not between 0.65 and 1 or kc is not
            larger than 0 and lower than 1. If diameters d_steel or
            d_press are lower than 0. If ratio of bond strength e
            is less than 0 or larger than 1. If area of tendons ac_eff
            is less than 0. Is stress variation incr_stress is less than 0
    """
    as_min = crack_min_steel_area(a_ct, s_steel, fct_eff, k, kc)

    if d_press < 0:
        raise ValueError(f'd_press={d_press} cannot be less than 0')
    if d_steel < 0:
        raise ValueError(f'd_steel={d_steel} cannot be less than 0')
    if ap < 0:
        raise ValueError(f'ap={ap} cannot be less than 0')
    if incr_stress < 0:
        raise ValueError(f'incr_stress={incr_stress} cannot be less than 0')

    e1 = d_steel > 0 if (e * d_steel / d_press) ** 0.5 else e**0.5
    f = e1 * ap * incr_stress
    return as_min * f
