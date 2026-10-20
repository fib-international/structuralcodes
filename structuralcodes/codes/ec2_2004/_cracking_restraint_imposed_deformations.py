"""Collection of functions from EUROCODE 1992-3:2006."""


def eps_sm_eps_cm_restraint_end(
    alpha_e: float,
    rho_p_eff: float,
    kc: float,
    k: float,
    fct_eff: float,
    Es: float,
) -> float:
    """Returns the strain difference (epsilon_sm - epsilon_cm) needed to
    compute the crack width for restraint member at its end.

    EN 1992-3:2006, Eq. (M.1).

    Args:
        alpha_e (float): Is the ratio Es/Ecm.
        rho_p_eff (float): Effective bond ratio between areas given by Eq.
            (7.10).
        kc (float): is a coefficient which takes account of the stress
            distribution within the section immediately prior to cracking and
            the change of the lever arm.
        k (float): is the coefficient which allow for the effect of
            non-uniform self-equilibrating stresses, which lead to a
            reduction of restraint forces.
            k=1 for webs w<=300mm or flanges widths less than 300mm
            k=0.65 for webs w>=800mm or flanges with widths greater than 800mm
            Intermediate values may be interpolated.
        fct_eff (float): Is the mean value of the tensile strength in MPa of
            the concrete effective at the time when the cracks may first be
            expected to occur: fct_eff=fctm or fctm(t) if crack is expected
            earlier than 28 days.
        Es (float): Steel elastic modulus in MPa.

    Returns:
        float: The calculated strain difference.
    """
    return (
        0.5 * alpha_e * kc * k * fct_eff * (1 + 1 / (alpha_e * rho_p_eff)) / Es
    )
