"""Provides functions to determine parameters of the bond stress-slip
relations as shown in Figure 6.1-1 of the fib
ModelCode 2010.

- Eq. 6.1-5.
- Eq. 6.1-6.
- Table 6.1-1 s_1, s_2 and s_3
- Solve s for Eq. 6.1-1 with substitution of tau_bu_split for tau_b.
- Eq. 6.1-19.
"""

import warnings


def eta_2(bond: str) -> float:
    """eta_2 according to Eq. 6.1-5.

    Args:
        bond: Bond condition according to Subsection 6.1.3.2. Input must be
        'Good' or 'Other'.

    Returns:
        float: eta_2 value.
    """
    eta_2 = {'Good': 1.0, 'Other': 0.7}

    try:
        return eta_2[bond]
    except KeyError:
        raise ValueError("Invalid bond condition. Must be 'Good' or 'Other'.")


def tau_bu_split(
    eta_2: float,
    f_cm: float,
    phi: float,
    c_min: float,
    c_max: float,
    k_m: float,
    K_tr: float,
) -> float:
    """Maximum bond stress in case of splitting failure
    according to Eq.(6.1-5).

    Args:
        eta_2: Parameter based on bond condition. 1.0 for 'Good' bond
        condition and 0.7 for 'Other' conditions.
        f_cm: Mean cylinder concrete compressive strength in MPa.
        phi: Nominal bar diameter in mm.
        c_min: Parameter according to MC2010 Figure 6.1-2.
        c_max: Parameter according to MC2010 Figure 6.1-2.
        k_m: Parameter according to MC2010 Figure 6.1-3.
        K_tr: To be calculated with MC2010 Eq.(6.1-6).

    Returns:
        float: tau_bu_split in MPa.
    """
    return (
        eta_2
        * 6.5
        * (f_cm / 25) ** 0.25
        * (25 / phi) ** 0.2
        * ((c_min / phi) ** 0.33 * (c_max / c_min) ** 0.1 + k_m * K_tr)
    )


def K_tr(
    n_t: float,
    A_st: float,
    n_b: float,
    phi: float,
    s_t: float,
) -> float:
    """K_tr according to Eq. 6.1-6.

    Args:
        n_t: Number of legs of confining reinforcement crossing a potential
        splitting failure surface at a section.
        A_st: Cross-sectional area of one leg of a confining bar in mm^2.
        n_b: Number of anchored bars or pairs of lapped bars in potential
        splitting surface.
        phi: Nominal bar diameter in mm.
        s_t: Longitudnal spacing of confining reinforcement in mm.

    Returns:
        float: K_tr as value.
    """
    return min((n_t * A_st / (n_b * phi * s_t)), 0.05)


def tau_bmax(bond: str, f_cm: float) -> float:
    """tau_bmax according to Table 6.1-1.

    Args:
        bond: Bond condition according to Subsection 6.1.3.2. Must be
        'Good' or 'Other'.
        f_cm: Mean cylinder concrete compressive strength in MPa.

    Returns:
        float: tau_bmax in MPa.
    """
    tau_bmax = {'Good': 2.5 * f_cm**0.5, 'Other': 1.25 * f_cm**0.5}

    try:
        return tau_bmax[bond]
    except KeyError:
        raise ValueError("Invalid input bond. Must be 'Good' or 'Other'.")


def s_1(bond: str) -> float:
    """s_1 according to Table 6.1-1.

    Args:
        bond: Bond condition according to Subsection 6.1.3.2. Must be
        'Good' or 'Other'.

    Returns:
        float: s_1 in mm.
    """
    s_1_values = {'Good': 1.0, 'Other': 1.8}

    try:
        return s_1_values[bond]
    except KeyError:
        raise ValueError("Invalid input bond. Must be 'Good' or 'Other'.")


def s_2(bond: str) -> float:
    """s_2 according to Table 6.1-1.

    Args:
        bond: Bond condition according to Subsection 6.1.3.2. Must be
        'Good' or 'Other'.

    Returns:
        float: s_2 in mm.
    """
    s_2_values = {'Good': 2.0, 'Other': 3.6}

    try:
        return s_2_values[bond]
    except KeyError:
        raise ValueError("Invalid input bond. Must be 'Good' or 'Other'.")


def s_3(
    failmod: str,
    bond: str,
    confin: str,
    c_clear: float,
    s_1: float,
) -> float:
    """s_3 according to Table 6.1-1.

    Args:
        failmod: Failure mode. Must be "PO" for Pull-out of "SP" for Splitting.
        bond: Bond condition according to Subsection 6.1.3.2. Must be
        'Good' or 'Other'.
        confin: Confinement conditions. Must be "Unconfined" or "Stirrups"
        c_clear: clear distance between ribs in mm.
        s_1: s_1 according to Table 6.1-1 columns 1 and 2 for Pull-out failure.

    Returns:
        float: s_3 in mm.
    """
    valid_failmods = ['PO', 'SP']
    valid_bonds = ['Good', 'Other']
    valid_confins = ['Unconfined', 'Stirrups']

    if failmod not in valid_failmods:
        raise ValueError("Invalid failmod value. Must be 'PO' or 'SP'.")
    if bond not in valid_bonds:
        raise ValueError("Invalid bond value. Must be 'Good' or 'Other'.")
    if confin not in valid_confins:
        raise ValueError(
            "Invalid confin value. Must be 'Unconfined' or 'Stirrups'."
        )

    s_3_values = {
        'PO': {
            'Good': c_clear,
            'Other': c_clear,
        },
        'SP': {
            'Good': {'Unconfined': 1.2 * s_1, 'Stirrups': 0.5 * c_clear},
            'Other': {'Unconfined': 1.2 * s_1, 'Stirrups': 0.5 * c_clear},
        },
    }

    if failmod == 'PO':
        return s_3_values[failmod][bond]
    if failmod == 'SP':
        return s_3_values[failmod][bond][confin]
    return 0.0

def s_tau_bu_split(
    tau_bmax: float,
    tau_bu_split: float,
    alpha: float,
    s_1: float,
) -> float:
    """Calculates the slip at tau_bu_split by rewriting Eq. 6.1-1 to
    solve for s and substituting tau_bu_split for tau_b.

    Args:
        tau_bmax: Maximum bond stress according to Table 6.1-1.
        tau_bu_split: Concrete bond stress at splitting failure.
        alpha: Parameter in Eq. 6.1-1. Alpha is 0.4 in Table 6.1-1.
        s_1: s_1 according to Table 6.1-1.

    Returns:
        float: slip at tau_bu_split in mm.
    """
    return (tau_bu_split / tau_bmax) ** (1 / alpha) * s_1


def f_stm(
    f_cm: float,
    phi: float,
    l_b: float,
    c_min: float,
    c_max: float,
    k_m: float,
    K_tr: float,
) -> float:
    """f_stm according to Eq. 6.1-19.

    Args:
        f_cm: Mean cylinder concrete compressive strength in MPa.
        phi: Nominal bar diameter in mm.
        l_b: Bond length in mm.
        c_min: Parameter shown in Figure 6.1-2.
        c_max: Parameter shown in Figure 6.1-2.
        k_m: Parameter shown in Figure 6.1-3.
        K_tr: Parameter calculated with Eq.(6.1-6).

    Returns:
        f_stm in MPa.

    """
    if not 15 < f_cm < 110:
        warnings.warn(
            'Warning: Eq.(6.1-19) is valid for 15 MPa < f_cm < 110 MPa.',
            UserWarning,
        )
    if not 0.5 < c_min / phi < 3.5:
        warnings.warn(
            'Warning: Eq.(6.1-19) is valid for 0.5 < c_min / phi < 3.5.',
            UserWarning,
        )
    if not 1.0 < c_max / c_min < 5.0:
        warnings.warn(
            'Warning: Eq.(6.1-19) is valid for 1.0 < c_max / c_min < 5.0',
            UserWarning,
        )
    if not K_tr <= 0.05:
        warnings.warn(
            'Warning: Eq.(6.1-19) is valid for K_tr <= 0.05.',
            UserWarning,
        )

    return (
        54
        * (f_cm / 25) ** 0.25
        * (25 / phi) ** 0.2
        * (l_b / phi) ** 0.55
        * ((c_min / phi) ** 0.25 * (c_max / c_min) ** 0.1 + k_m * K_tr)
    )


def tau_yield(
    f_y: float,
    l_b: float,
    phi: float,
) -> float:
    """Calculates the bond stress at yield of the reinforcement bar based
    on a uniform bond stress over l_b.

    The bond stress at yield is a limit for:
    - tau_bmax, given with eps_s < eps_s,y in Table 6.1-1.
    - tau_bu_split, explicitly stated for Eq. 6.1-19.

    Args:
        f_y: Reinforcement bar yield stress in MPa.
        l_b: Bond length with uniform bond stress assumed in mm.
        phi: Nominal bar diameter in mm.

    Returns:
        Bond stress at yield of rebar in MPa.

    """
    return f_y / (4 * l_b / phi)
