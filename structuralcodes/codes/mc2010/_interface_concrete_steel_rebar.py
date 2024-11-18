"""
A collection of functions for determining the parameters of the analytical bond stress-slip relation, shown in Figure
6.1-1 of the fib ModelCode 2010.
"""

import warnings


def tau_bu_split(
        f_cm: float,
        phi: float,
        c_min: float,
        c_max: float,
        k_m: float,
        K_tr: float,
        bond: str,
        f_ym: float = None,
) -> dict:
    """
    Eq. (6.1-5) of MC2010. This equation is derived from Eq. (6.1-19). Unlike in the ModelCode, the limit values and
    input validity of Eq. (6.1-19) are also included in this function.

    Args:
        bond: "Good" or "Poor" bond conditions, described in section 6.1.3.2 of the fib ModelCode 2010;
        f_cm: Mean value of compressive strength of concrete in MPa;
        phi: Nominal bar diameter in mm;
        c_min: Figure 6.1-2 of the fib ModelCode 2010;
        c_max: Figure 6.1-2 of the fib ModelCode 2010;
        k_m: Figure 6.1-3 of the fib ModelCode 2010;
        K_tr: Eq (6.1-6) of the fib ModelCode 2010;
        f_ym: Mean yield stress of the reinforcement in MPa;

    Returns:
        Dictionary with ultimate bond stress corresponding to splitting, yield and pull-out, and the failure mode and
        failure bond stress.
    """

    # Raise warnings
    if not 15 < f_cm < 110:
        warnings.warn("Warning: Eq. (6.1-19) is valid for 15 MPa < f_cm < 110 MPa")
    if not 0.5 < c_min / phi < 3.5:
        warnings.warn("Warning: Eq. (6.1-19) is valid for 0.5 < c_min / phi < 3.5")
    if not 1.0 < c_max / c_min < 5.0:
        warnings.warn("Warning: Eq. (6.1-19) is valid for 1.0 < c_max / c_min < 5.0")
    if not K_tr <= 0.05:
        warnings.warn("Warning: Eq. (6.1-19) is valid for K_tr <= 0.05")
    if f_ym is None:
        warnings.warn(
            "Warning: No input for f_ym.")


    bond_param = {
        'Good': {
            'tau_bmax': 2.5 * (f_cm) ** 0.5,
            'eta_2': 1.0
        },
        'Poor': {
            'tau_bmax': 1.25 * (f_cm) ** 0.5,
            'eta_2': 0.7
        }
    }

    tau_bu_split = (
            bond_param[bond]['eta_2']
            * 6.5
            * (f_cm / 25) ** 0.25
            * (25 / phi) ** 0.2
            * (
                    (c_min / phi) ** 0.25  # <-- Deviates from MC2010, but follows fib Bulletin 72, Appendix A.
                    * (c_max / c_min) ** 0.1
                    + k_m * K_tr
            )
    )

    tau_all = {
        'Pull-out': bond_param[bond]['tau_bmax'],
        'Splitting': tau_bu_split,
    }

    if f_ym is None:
        msg = 'Potential yielding of rebar prior to reaching tau_bu_split or tau_bmax was not checked.'
        if len(tau_all) != len(set(tau_all)):
            warnings.warn("Warning: The tau of two or more failure modes is identical")
        else:
           tau_fail = min(tau_all.values())
           fail_mode = f'{min(tau_all, key=tau_all.get)}. Warning: {msg}'
    else:
        tau_all['Yield'] = f_ym / 20,  # Based on the assumption that the bond stress is uniform over l_b = 5 * phi
        if len(tau_all) != len(set(tau_all)):
            warnings.warn("Warning: The tau of two or more failure modes is identical")
        else:
           tau_fail = min(tau_all.values())
           fail_mode = min(tau_all, key=tau_all.get)

    return {
        'tau_bu_split': tau_all['Splitting'],
        'tau_yield': tau_all['Yield'],
        'tau_bmax': tau_all['Pull-out'],
        'fail_mode': fail_mode,
        'tau_fail': tau_fail,
    }


def s_tau_bu_split(
        f_cm: float,
        bond: str,
        tau_bu_split: float,
) -> float:
    """
    Calculate s(tau_bu,split) of ModelCode 2010 Table 6.1-1.


    Args:
        f_cm: Mean value of compressive strength of concrete in MPa;
        bond: "Good" or "Poor" bond conditions, as defined in fib ModelCode 2010 section 6.1.3.2;
        tau_bu_split: Concrete bond stress at splitting failure.

    Returns:
        s(tau_bu,split)
    """

    alpha = 0.4

    if bond == 'Good':
        s1_PO = 1.0
        tau_bmax = 2.5 * f_cm ** 0.5
    elif bond == 'Poor':
        s1_PO = 2.0
        tau_bmax = 1.25 * f_cm ** 0.5
    else:
        warnings.ValueError("Invalid bond condition. Use 'Good' or 'Poor'")

    return s1_PO * (tau_bu_split / tau_bmax) ** (1 / alpha)
