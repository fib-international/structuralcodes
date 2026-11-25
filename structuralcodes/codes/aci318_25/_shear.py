"""One-way shear strength functions according to ACI 318-25, Ch. 22.5.

Units: every function takes and returns SI units (MPa, mm, mm2, N), as the
rest of the library does. The expressions of ACI 318-25 are written in
inch-pound units, with empirical coefficients on sqrt(f'c) in psi. Each
function applies the Code expression in that form and converts with exact
factors, so that every coefficient below is the one printed in the Code.

Edition: the expressions were read in ACI 318-25 and compared with
ACI 318-19. Differences are recorded in the docstring of the function they
affect. The size effect factor lambda_s departs from ACI 318-25 as printed,
see :func:`lambda_s`.
"""

import math
import typing as t

# Exact conversion factors
_MPA_PER_PSI = 0.006894757293168361
_MM_PER_IN = 25.4

# Sec. 22.5.3.1: limit on sqrt(f'c) for one-way shear, in psi
_SQRT_FC_LIMIT_PSI = 100.0


def _sqrt_fc(fc: float, limit_psi: t.Optional[float] = None) -> float:
    """Return sqrt(f'c), evaluated in psi and expressed in MPa.

    Args:
        fc (float): Specified compressive strength of concrete in MPa.
        limit_psi (Optional[float]): Upper limit on sqrt(f'c) in psi. No
            limit is applied if None.

    Returns:
        float: sqrt(f'c in psi), converted from psi to MPa.

    Raises:
        ValueError: If fc is not positive.
    """
    if fc <= 0.0:
        raise ValueError(f'fc must be positive, got {fc}')
    root_psi = math.sqrt(fc / _MPA_PER_PSI)
    if limit_psi is not None:
        root_psi = min(root_psi, limit_psi)
    return root_psi * _MPA_PER_PSI


def _axial_stress(Nu: float, Ag: float, fc: float) -> float:
    """Return the term Nu / (6 Ag) of Table 22.5.5.1 in MPa.

    The term is limited to 0.05 f'c by Sec. 22.5.5.1.2.

    Args:
        Nu (float): Factored axial force in N, positive in compression.
        Ag (float): Gross area of the section in mm2. Required if Nu is not
            zero.
        fc (float): Specified compressive strength of concrete in MPa.

    Returns:
        float: The axial term in MPa.

    Raises:
        ValueError: If Nu is not zero and Ag is not positive.
    """
    if Nu == 0.0:
        return 0.0
    if Ag <= 0.0:
        raise ValueError(
            f'Ag must be positive when Nu is not zero, got Ag = {Ag}'
        )
    return min(Nu / (6.0 * Ag), 0.05 * fc)


# ---------------------------------------------------------------------------
# Size-effect factor
# ---------------------------------------------------------------------------


def lambda_s(d: float) -> float:
    """Size effect modification factor for shear.

    Evaluated in the form of ACI 318-19, Eq. (22.5.5.1.3), with d in
    inches: ``lambda_s = sqrt(2 / (1 + d / 10)) <= 1``.

    Note:
        ACI 318-25 prints Eq. (22.5.5.1.3) without the square root, as
        ``2 / (1 + d / 10) <= 1``. This function does not follow the
        equation as printed, for these reasons:

        - ACI 318-19 prints the equation under a square root.
        - The clause carries no change bar in ACI 318-25, and the
          commentary R22.5.5.1.3 is the same in the two editions. It
          describes the factor as consistent with fracture mechanics
          theory, whose size effect law has the square root form.
        - No erratum to ACI 318-25 on this equation had been published
          when this function was written (checked November 2026).

        The printed form would give a smaller factor for d greater than
        10 in, for example 0.667 against 0.816 at d = 20 in. Review this
        function if ACI confirms the printed form.

    Args:
        d (float): Effective depth of the member in mm.

    Returns:
        float: Size effect factor lambda_s, not greater than 1.0.

    Raises:
        ValueError: If d is not positive.
    """
    if d <= 0.0:
        raise ValueError(f'd must be positive, got {d}')
    d_in = d / _MM_PER_IN
    return min(math.sqrt(2.0 / (1.0 + d_in / 10.0)), 1.0)


# ---------------------------------------------------------------------------
# Concrete shear strength
# ---------------------------------------------------------------------------


def Vc_detailed(
    fc: float,
    bw: float,
    d: float,
    rho_w: float,
    Nu: float = 0.0,
    Ag: float = 0.0,
    lambda_concrete: float = 1.0,
    has_min_shear_reinforcement: bool = False,
    apply_lower_bound: bool = True,
) -> float:
    """Nominal one-way shear strength provided by concrete.

    ACI 318-25, Table 22.5.5.1, with f'c in psi and lengths in inches:

    - Eq. (b), for Av >= Av,min:
      ``Vc = [8 lambda rho_w^(1/3) sqrt(f'c) + Nu / (6 Ag)] bw d``
    - Eq. (c), for Av < Av,min:
      ``Vc = [8 lambda_s lambda rho_w^(1/3) sqrt(f'c) + Nu / (6 Ag)] bw d``

    Limits applied:

    - Sec. 22.5.3.1: sqrt(f'c) is not taken greater than 100 psi. The
      relaxation of Sec. 22.5.3.2 is not implemented.
    - Sec. 22.5.5.1.2: Nu / (6 Ag) is not taken greater than 0.05 f'c.
    - Sec. 22.5.5.1.1: Vc is not taken greater than
      ``5 lambda sqrt(f'c) bw d``, and need not be taken less than
      ``lambda sqrt(f'c) bw d`` unless the element is in net axial tension
      or Sec. 18.6.5.2 or 18.7.6.2.1 applies.
    - Table 22.5.5.1, Note 2: Vc is not taken less than zero.

    Note:
        The lower bound of Sec. 22.5.5.1.1 is not in ACI 318-19. Set
        ``apply_lower_bound=False`` to obtain the ACI 318-19 result, or
        where Sec. 18.6.5.2 or 18.7.6.2.1 applies.

    Args:
        fc (float): Specified compressive strength of concrete in MPa.
        bw (float): Web width in mm.
        d (float): Effective depth in mm.
        rho_w (float): Longitudinal reinforcement ratio As / (bw d).
        Nu (float): Factored axial force in N, positive in compression and
            negative in tension. Defaults to 0.
        Ag (float): Gross area of the section in mm2. Required if Nu is not
            zero.
        lambda_concrete (float): Lightweight concrete factor lambda.
            Defaults to 1.0.
        has_min_shear_reinforcement (bool): True if Av >= Av,min, which
            selects Eq. (b). Defaults to False, which selects Eq. (c).
        apply_lower_bound (bool): Apply the lower bound of Sec. 22.5.5.1.1.
            The bound is never applied when Nu is negative.

    Returns:
        float: Nominal concrete shear strength Vc in N.

    Raises:
        ValueError: If fc, bw, d or rho_w is not positive, or if Nu is not
            zero and Ag is not positive.
    """
    if bw <= 0.0 or d <= 0.0:
        raise ValueError(f'bw and d must be positive, got {bw}, {d}')
    if rho_w <= 0.0:
        raise ValueError(f'rho_w must be positive, got {rho_w}')

    sqrt_fc = _sqrt_fc(fc, _SQRT_FC_LIMIT_PSI)
    axial = _axial_stress(Nu, Ag, fc)
    size = 1.0 if has_min_shear_reinforcement else lambda_s(d)

    stress = (
        8.0 * size * lambda_concrete * rho_w ** (1.0 / 3.0) * sqrt_fc + axial
    )
    Vc = stress * bw * d

    # Sec. 22.5.5.1.1
    Vc = min(Vc, 5.0 * lambda_concrete * sqrt_fc * bw * d)
    if apply_lower_bound and Nu >= 0.0:
        Vc = max(Vc, lambda_concrete * sqrt_fc * bw * d)

    return max(Vc, 0.0)


def Vc_simplified(
    fc: float,
    bw: float,
    d: float,
    Nu: float = 0.0,
    Ag: float = 0.0,
    lambda_concrete: float = 1.0,
) -> float:
    """Nominal one-way shear strength provided by concrete, Eq. (a).

    ACI 318-25, Table 22.5.5.1 Eq. (a), with f'c in psi and lengths in
    inches. Permitted only where Av >= Av,min:
    ``Vc = [2 lambda sqrt(f'c) + Nu / (6 Ag)] bw d``.

    The limits of Sec. 22.5.3.1, 22.5.5.1.1 (upper limit), 22.5.5.1.2 and
    Table 22.5.5.1 Note 2 are applied as in :func:`Vc_detailed`. The same
    expression and limits are in ACI 318-19.

    Args:
        fc (float): Specified compressive strength of concrete in MPa.
        bw (float): Web width in mm.
        d (float): Effective depth in mm.
        Nu (float): Factored axial force in N, positive in compression and
            negative in tension. Defaults to 0.
        Ag (float): Gross area of the section in mm2. Required if Nu is not
            zero.
        lambda_concrete (float): Lightweight concrete factor lambda.
            Defaults to 1.0.

    Returns:
        float: Nominal concrete shear strength Vc in N.

    Raises:
        ValueError: If fc, bw or d is not positive, or if Nu is not zero
            and Ag is not positive.
    """
    if bw <= 0.0 or d <= 0.0:
        raise ValueError(f'bw and d must be positive, got {bw}, {d}')

    sqrt_fc = _sqrt_fc(fc, _SQRT_FC_LIMIT_PSI)
    axial = _axial_stress(Nu, Ag, fc)

    Vc = (2.0 * lambda_concrete * sqrt_fc + axial) * bw * d
    Vc = min(Vc, 5.0 * lambda_concrete * sqrt_fc * bw * d)
    return max(Vc, 0.0)


# ---------------------------------------------------------------------------
# Steel shear contribution and nominal shear strength
# ---------------------------------------------------------------------------


def Vs(Av: float, fyt: float, d: float, s: float) -> float:
    """Nominal shear strength provided by transverse reinforcement.

    ACI 318-19, Eq. (22.5.8.5.3), for shear reinforcement perpendicular to
    the axis of the member: ``Vs = Av fyt d / s``. The expression is
    dimensionally consistent and needs no unit conversion.

    Args:
        Av (float): Area of shear reinforcement within spacing s in mm2.
        fyt (float): Specified yield strength of transverse reinforcement in
            MPa.
        d (float): Effective depth in mm.
        s (float): Center-to-center spacing of transverse reinforcement in
            mm.

    Returns:
        float: Nominal shear strength Vs in N.
    """
    return Av * fyt * d / s


def Vn(Vc: float, Vs: float) -> float:  # noqa: N803
    """Nominal one-way shear strength of a section, Vn = Vc + Vs.

    Args:
        Vc (float): Nominal shear strength provided by concrete in N.
        Vs (float): Nominal shear strength provided by shear reinforcement
            in N.

    Returns:
        float: Nominal shear strength Vn in N.
    """
    return Vc + Vs


# ---------------------------------------------------------------------------
# Cross-section size check
# ---------------------------------------------------------------------------


def check_cross_section(
    Vu: float,
    phi: float,
    Vc: float,
    fc: float,
    bw: float,
    d: float,
) -> bool:
    """Check the limit on cross-sectional dimensions.

    ACI 318-25, Eq. (22.5.1.2), with f'c in psi and lengths in inches:
    ``Vu <= phi (Vc + 8 sqrt(f'c) bw d)``. The same expression is in
    ACI 318-19. The limit of Sec. 22.5.3.1 on sqrt(f'c) is written for the
    calculation of Vc and is not applied to this expression.

    Args:
        Vu (float): Factored shear force in N.
        phi (float): Strength reduction factor for shear.
        Vc (float): Nominal shear strength provided by concrete in N.
        fc (float): Specified compressive strength of concrete in MPa.
        bw (float): Web width in mm.
        d (float): Effective depth in mm.

    Returns:
        bool: True if the cross-section satisfies Eq. (22.5.1.2).
    """
    return Vu <= phi * (Vc + 8.0 * _sqrt_fc(fc) * bw * d)


# ---------------------------------------------------------------------------
# Minimum transverse reinforcement
# ---------------------------------------------------------------------------


def Av_min_per_s(fc: float, bw: float, fyt: float) -> float:
    """Minimum area of shear reinforcement per unit length of a beam.

    ACI 318-25, Table 9.6.3.4 (a) and (b), nonprestressed beams, with
    stresses in psi: ``Av,min / s`` is the greater of
    ``0.75 sqrt(f'c) bw / fyt`` and ``50 bw / fyt``. The same table is in
    ACI 318-19.

    Args:
        fc (float): Specified compressive strength of concrete in MPa.
        bw (float): Web width in mm.
        fyt (float): Specified yield strength of transverse reinforcement in
            MPa.

    Returns:
        float: Minimum Av / s in mm2/mm.
    """
    stress = max(0.75 * _sqrt_fc(fc), 50.0 * _MPA_PER_PSI)
    return stress * bw / fyt


def shear_reinforcement_required(Vu: float, phi_Vc: float) -> bool:
    """Determine whether Vu exceeds the design strength of the concrete.

    Returns True when ``Vu > phi Vc``. This is a comparison only. The Code
    thresholds at which minimum shear reinforcement is required in beams
    and slabs are separate provisions and are not evaluated here.

    Args:
        Vu (float): Factored shear force in N.
        phi_Vc (float): Design shear strength of the concrete, phi Vc, in N.

    Returns:
        bool: True if Vu exceeds phi Vc.
    """
    return Vu > phi_Vc


# ---------------------------------------------------------------------------
# Maximum stirrup spacing
# ---------------------------------------------------------------------------


def max_stirrup_spacing(
    d: float,
    Vs: float,
    fc: float,
    bw: float,
) -> float:
    """Maximum spacing of legs of shear reinforcement along a beam.

    ACI 318-25, Table 9.7.6.2.2, nonprestressed beam, along the length,
    with f'c in psi and lengths in inches. Where the required
    ``Vs <= 4 sqrt(f'c) bw d`` the maximum spacing is the lesser of d/2
    and 24 in. Above that it is the lesser of d/4 and 12 in. The same
    table is in ACI 318-19.

    Args:
        d (float): Effective depth in mm.
        Vs (float): Required nominal shear strength of the shear
            reinforcement in N.
        fc (float): Specified compressive strength of concrete in MPa.
        bw (float): Web width in mm.

    Returns:
        float: Maximum spacing in mm.
    """
    threshold = 4.0 * _sqrt_fc(fc) * bw * d
    if Vs <= threshold:
        return min(d / 2.0, 24.0 * _MM_PER_IN)
    return min(d / 4.0, 12.0 * _MM_PER_IN)
