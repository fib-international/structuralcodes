"""One-way shear strength functions according to ACI 318-25, Ch. 22.5."""

import math


# ---------------------------------------------------------------------------
# Size-effect factor
# ---------------------------------------------------------------------------


def lambda_s(d: float) -> float:
    """Size-effect modification factor for shear.

    ACI 318-25, Eq. 22.5.5.1.3.  Converts *d* from mm to inches and applies:
    ``lambda_s = min(2 / (1 + d_in / 10), 1.0)``.

    Args:
        d (float): Effective depth of the member in mm.

    Returns:
        float: Size-effect factor lambda_s (dimensionless, <= 1.0).
    """
    d_in = d / 25.4
    return min(2.0 / (1.0 + d_in / 10.0), 1.0)


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
    Av_provided: float = 0.0,
    Av_min: float = 0.0,
) -> float:
    """Nominal concrete shear strength — detailed method.

    ACI 318-25, Table 22.5.5.1.  The formula depends on whether minimum
    transverse reinforcement is provided:

    - **With** minimum reinforcement (``Av_provided >= Av_min``):
      ``Vc = (8*lambda*(rho_w)^(1/3)*sqrt_fc + axial_term) * bw * d``
    - **Without** minimum reinforcement (``Av_provided < Av_min``):
      ``Vc = (8*lambda_s(d)*lambda*(rho_w)^(1/3)*sqrt_fc + axial_term) * bw * d``

    Limits applied (Sec. 22.5.5.1.1):

    - Cap: ``Vc <= 5 * lambda * sqrt_fc * bw * d``
    - Floor: ``Vc >= lambda * sqrt_fc * bw * d`` (only when ``Nu >= 0``)
    - ``Vc >= 0`` always

    The axial term is ``min(Nu/(6*Ag), 0.05*fc)`` when ``Ag > 0``, else 0.
    *Nu* is positive for compression, negative for tension (Sec. 22.5.3.2).

    Args:
        fc (float): Specified compressive strength of concrete in MPa.
        bw (float): Web width in mm.
        d (float): Effective depth in mm.
        rho_w (float): Longitudinal reinforcement ratio ``As / (bw * d)``
            (dimensionless).
        Nu (float): Factored axial force in N; positive = compression,
            negative = tension. Defaults to 0.
        Ag (float): Gross cross-sectional area in mm². Required when *Nu* != 0.
            Defaults to 0.
        lambda_concrete (float): Concrete lightweight modification factor
            (dimensionless). Defaults to 1.0 (normal-weight).
        Av_provided (float): Provided transverse reinforcement area per unit
            length (mm²/mm), or just area (mm²) matched against *Av_min* on
            the same basis. Defaults to 0.
        Av_min (float): Minimum required transverse reinforcement area on the
            same basis as *Av_provided*. Defaults to 0.

    Returns:
        float: Nominal concrete shear strength *Vc* in N.
    """
    # Sec. 22.5.3.1 — limit on sqrt(fc)
    sqrt_fc = min(math.sqrt(fc), 8.3)

    # Axial-load term (Sec. 22.5.3.2)
    if Ag > 0.0:
        axial_term = min(Nu / (6.0 * Ag), 0.05 * fc)
    else:
        axial_term = 0.0

    # Reinforcement-ratio term
    rho_term = rho_w ** (1.0 / 3.0)

    if Av_provided >= Av_min:
        # Minimum reinforcement provided — no size-effect penalty
        base = 8.0 * lambda_concrete * rho_term * sqrt_fc + axial_term
    else:
        # Below minimum — apply size-effect factor
        ls = lambda_s(d)
        base = 8.0 * ls * lambda_concrete * rho_term * sqrt_fc + axial_term

    Vc = base * bw * d

    # Upper cap (Sec. 22.5.5.1.1)
    Vc_max = 5.0 * lambda_concrete * sqrt_fc * bw * d
    Vc = min(Vc, Vc_max)

    # Lower floor — only apply when there is no net tension (Nu >= 0)
    if Nu >= 0.0:
        Vc_floor = lambda_concrete * sqrt_fc * bw * d
        Vc = max(Vc, Vc_floor)

    # Absolute minimum
    return max(Vc, 0.0)


def Vc_simplified(
    fc: float,
    bw: float,
    d: float,
    Nu: float = 0.0,
    Ag: float = 0.0,
    lambda_concrete: float = 1.0,
) -> float:
    """Nominal concrete shear strength — simplified method.

    ACI 318-25, Table 22.5.5.1(a).  Valid only when ``Av >= Av_min``:
    ``Vc = (2 * lambda * sqrt(fc) + axial_term) * bw * d``.

    Args:
        fc (float): Specified compressive strength of concrete in MPa.
        bw (float): Web width in mm.
        d (float): Effective depth in mm.
        Nu (float): Factored axial force in N; positive = compression,
            negative = tension. Defaults to 0.
        Ag (float): Gross cross-sectional area in mm². Defaults to 0.
        lambda_concrete (float): Concrete lightweight modification factor
            (dimensionless). Defaults to 1.0.

    Returns:
        float: Nominal concrete shear strength *Vc* in N.
    """
    sqrt_fc = math.sqrt(fc)

    if Ag > 0.0:
        axial_term = min(Nu / (6.0 * Ag), 0.05 * fc)
    else:
        axial_term = 0.0

    return (2.0 * lambda_concrete * sqrt_fc + axial_term) * bw * d


# ---------------------------------------------------------------------------
# Steel shear contribution and nominal shear strength
# ---------------------------------------------------------------------------


def Vs(Av: float, fyt: float, d: float, s: float) -> float:
    """Nominal shear strength provided by transverse reinforcement.

    ACI 318-25, Eq. 22.5.8.5.3:
    ``Vs = Av * fyt * d / s``.

    Args:
        Av (float): Area of shear reinforcement within spacing *s* in mm².
        fyt (float): Specified yield strength of transverse reinforcement in
            MPa.
        d (float): Effective depth in mm.
        s (float): Center-to-center spacing of transverse reinforcement in mm.

    Returns:
        float: Steel shear contribution *Vs* in N.
    """
    return Av * fyt * d / s


def Vn(Vc: float, Vs: float) -> float:  # noqa: N803
    """Nominal shear strength of a section.

    ACI 318-25, Sec. 22.5.1.1:
    ``Vn = Vc + Vs``.

    Args:
        Vc (float): Nominal concrete shear strength in N.
        Vs (float): Nominal steel shear strength in N.

    Returns:
        float: Total nominal shear strength *Vn* in N.
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
    """Check that the cross-section is large enough.

    ACI 318-25, Eq. 22.5.1.2.  The factored shear must not exceed:
    ``Vu <= phi * (Vc + 8 * sqrt(fc) * bw * d)``.

    Args:
        Vu (float): Factored shear force in N.
        phi (float): Strength reduction factor for shear (dimensionless).
        Vc (float): Nominal concrete shear strength in N.
        fc (float): Specified compressive strength of concrete in MPa.
        bw (float): Web width in mm.
        d (float): Effective depth in mm.

    Returns:
        bool: ``True`` if the cross-section size is adequate; ``False``
        otherwise.
    """
    limit = phi * (Vc + 8.0 * math.sqrt(fc) * bw * d)
    return Vu <= limit


# ---------------------------------------------------------------------------
# Minimum transverse reinforcement
# ---------------------------------------------------------------------------


def Av_min_per_s(fc: float, bw: float, fyt: float) -> float:
    """Minimum area of shear reinforcement per unit length.

    ACI 318-25, Sec. 9.6.3.4:
    ``Av/s = max(0.062 * sqrt(fc), 0.35) * bw / fyt``.

    Args:
        fc (float): Specified compressive strength of concrete in MPa.
        bw (float): Web width in mm.
        fyt (float): Specified yield strength of transverse reinforcement in
            MPa.

    Returns:
        float: Minimum ``Av/s`` in mm²/mm.
    """
    return max(0.062 * math.sqrt(fc), 0.35) * bw / fyt


def shear_reinforcement_required(Vu: float, phi_Vc: float) -> bool:
    """Determine whether shear reinforcement is required.

    ACI 318-25, Sec. 9.6.3.1.  Reinforcement is required when:
    ``Vu > phi * Vc``.

    Args:
        Vu (float): Factored shear force in N.
        phi_Vc (float): Design concrete shear strength ``phi * Vc`` in N.

    Returns:
        bool: ``True`` if shear reinforcement is required; ``False``
        otherwise.
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
    """Maximum permitted stirrup spacing.

    ACI 318-25, Sec. 9.7.6.2.2.  When ``Vs <= 4*sqrt(fc)*bw*d`` the maximum
    spacing is ``min(d/2, 600 mm)``; when *Vs* exceeds that threshold the
    limit is halved to ``min(d/4, 300 mm)``.

    Args:
        d (float): Effective depth in mm.
        Vs (float): Nominal steel shear contribution in N.
        fc (float): Specified compressive strength of concrete in MPa.
        bw (float): Web width in mm.

    Returns:
        float: Maximum stirrup spacing in mm.
    """
    threshold = 4.0 * math.sqrt(fc) * bw * d
    if Vs <= threshold:
        return min(d / 2.0, 600.0)
    else:
        return min(d / 4.0, 300.0)
