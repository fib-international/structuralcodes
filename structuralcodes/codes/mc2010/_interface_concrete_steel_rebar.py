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


def eta_2(bond: t.Literal['good', 'other']) -> float:
    """Bond coefficient eta_2.

    fib Model Code 2010, Eq. (6.1-5).

    Args:
        bond (str): Bond condition according to Subsection 6.1.3.2. Input must
            be 'good' or 'other'.

    Returns:
        float: eta_2 value.

    Raises:
        ValueError: If a proper value for 'bond' is not input.
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
    """Maximum bond stress in case of splitting failure.

    fib Model Code 2010, Eq. (6.1-5).

    Args:
        eta_2 (float): Parameter based on bond condition. 1.0 for 'good' bond
            condition and 0.7 for 'other' conditions.
        f_cm (float): Mean cylinder concrete compressive strength in MPa.
        phi (float): Nominal bar diameter in mm.
        c_min (float): Parameter according to MC2010 Figure 6.1-2.
        c_max (float): Parameter according to MC2010 Figure 6.1-2.
        k_m (float): Parameter according to MC2010 Figure 6.1-3.
        K_tr (float): To be calculated with MC2010 Eq.(6.1-6).

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
    """Coefficient accounting for transverse stresses, K_tr.

    fib Model Code 2010, Eq. (6.1-6).

    Args:
        n_t (float): Number of legs of confining reinforcement crossing a
            potential splitting failure surface at a section.
        A_st (float): Cross-sectional area of one leg of a confining bar in
            mm^2.
        n_b (float): Number of anchored bars or pairs of lapped bars in
            potential splitting surface.
        phi (float): Nominal bar diameter in mm.
        s_t (float): Longitudnal spacing of confining reinforcement in mm.

    Returns:
        float: K_tr as value.
    """
    return min((n_t * A_st / (n_b * phi * s_t)), 0.05)


def tau_bmax(bond: t.Literal['good', 'other'], f_cm: float) -> float:
    """Maximum bond stress tau_bmax.

    fib Model Code 2010, Table 6.1-1.

    Args:
        bond (str): Bond condition according to Subsection 6.1.3.2. Must be
            'good' or 'other'.
        f_cm (float): Mean cylinder concrete compressive strength in MPa.

    Returns:
        float: tau_bmax in MPa.

    Raises:
        ValueError: If a proper value for 'bond' is not input.
    """
    tau_bmax = {'Good': 2.5 * f_cm**0.5, 'Other': 1.25 * f_cm**0.5}

    try:
        return tau_bmax[bond]
    except KeyError:
        raise ValueError("Invalid input bond. Must be 'Good' or 'Other'.")


def s_1(bond: t.Literal['good', 'other']) -> float:
    """Slip at maximum bond stress, s_1.

    fib Model Code 2010, Table 6.1-1.

    Args:
        bond (str): Bond condition according to Subsection 6.1.3.2. Must be
            'good' or 'other'.

    Returns:
        float: s_1 in mm.

    Raises:
        ValueError: If a proper value for 'bond' is not input.
    """
    s_1_values = {'Good': 1.0, 'Other': 1.8}

    try:
        return s_1_values[bond]
    except KeyError:
        raise ValueError("Invalid input bond. Must be 'Good' or 'Other'.")


def s_2(bond: t.Literal['good', 'other']) -> float:
    """s_2 according to Table 6.1-1.

    Args:
        bond (str): Bond condition according to Subsection 6.1.3.2. Must be
            'good' or 'other'.

    Returns:
        float: s_2 in mm.

    Raises:
        ValueError: If a proper value for 'bond' is not input.
    """
    s_2_values = {'Good': 2.0, 'Other': 3.6}

    try:
        return s_2_values[bond]
    except KeyError:
        raise ValueError("Invalid input bond. Must be 'Good' or 'Other'.")


def s_3(
    failmod: t.Literal['PO', 'SP'],
    bond: t.Literal['good', 'other'],
    confinement: t.Literal['unconfined', 'stirrups'],
    c_clear: float,
    s_1: float,
) -> float:
    """s_3 according to Table 6.1-1.

    Args:
        failmod (str): Failure mode. Must be "PO" for Pull-out of "SP" for
            Splitting.
        bond (str): Bond condition according to Subsection 6.1.3.2. Must be
            'Good' or 'Other'.
        confinement (str): Confinement conditions. Must be "Unconfined" or
            "Stirrups"
        c_clear (float): Clear distance between ribs in mm.
        s_1 (float): s_1 according to Table 6.1-1 columns 1 and 2 for Pull-out
            failure.

    Returns:
        float: s_3 in mm.

    Raises:
        ValueError: If a proper value for 'failmod' is not input.
        ValueError: If a proper value for 'bond' is not input.
        ValueError: If a proper value for 'confinement' is not input.
    """
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
    """Calculates the slip at tau_bu_split by rewriting Eq. 6.1-1 to solve for
    s and substituting tau_bu_split for tau_b.

    Args:
        tau_bmax (float): Maximum bond stress according to Table 6.1-1.
        tau_bu_split (float): Concrete bond stress at splitting failure.
        alpha (float): Parameter in Eq. 6.1-1. Alpha is 0.4 in Table 6.1-1.
        s_1 (float): s_1 according to Table 6.1-1.

    Returns:
        float: Slip at tau_bu_split in mm.
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
    """Reinforcement stress, f_stm.

    fib Model Code 2010, Eq. (6.1-19).

    Args:
        f_cm (float): Mean cylinder concrete compressive strength in MPa.
        phi (float): Nominal bar diameter in mm.
        l_b (float): Bond length in mm.
        c_min (float): Parameter shown in Figure 6.1-2.
        c_max (float): Parameter shown in Figure 6.1-2.
        k_m (float): Parameter shown in Figure 6.1-3.
        K_tr (float): Parameter calculated with Eq.(6.1-6).

    Returns:
        float: f_stm in MPa.

    Raises:
        UserWarning: If not 15 MPa < f_cm < 110 MPa.
        UserWarning: If not 0.5 < c_min / phi < 3.5.
        UserWarning: If not 1.0 < c_max / c_min < 5.0.
        UserWarning: If not K_tr <= 0.05.
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
        f_y (float): Reinforcement bar yield stress in MPa.
        l_b (float): Bond length with uniform bond stress assumed in mm.
        phi (float): Nominal bar diameter in mm.

    Returns:
        float: Bond stress at yield of rebar in MPa.

    """
    return f_y / (4 * l_b / phi)
