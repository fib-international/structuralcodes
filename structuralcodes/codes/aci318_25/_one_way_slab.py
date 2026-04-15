"""One-way slab design rules according to ACI 318-25, Chapter 7."""

# ---------------------------------------------------------------------------
# Module-level constants
# ---------------------------------------------------------------------------

THICKNESS_RATIOS = {
    'simply_supported': 20,
    'one_end_continuous': 24,
    'both_ends_continuous': 28,
    'cantilever': 10,
}


# ---------------------------------------------------------------------------
# Minimum thickness
# ---------------------------------------------------------------------------


def min_thickness(
    span: float,
    support_condition: str,
    fy: float = 420.0,
    lightweight: bool = False,
    wc: float = 2320.0,
) -> float:
    """Minimum slab thickness for one-way slabs not supporting partitions.

    ACI 318-25, Table 7.3.1.1.  Returns the minimum overall thickness *h*
    (in mm) for the given span and support condition.

    Modifications applied when applicable:

    - **Sec. 7.3.1.1.1**: For reinforcement with *fy* other than 60 ksi,
      multiply the table value by ``(0.4 + fy_psi / 100 000)``.
    - **Sec. 7.3.1.1.2**: For lightweight concrete, multiply the table value
      by ``max(1.65 - 0.005 * wc_pcf, 1.09)``.

    Both factors are applied cumulatively when both conditions are present.

    Args:
        span (float): Clear span length in mm.
        support_condition (str): One of ``'simply_supported'``,
            ``'one_end_continuous'``, ``'both_ends_continuous'``, or
            ``'cantilever'``.
        fy (float): Specified yield strength of reinforcement in MPa.
            Defaults to 420 MPa (~60 ksi).
        lightweight (bool): ``True`` if lightweight concrete is used.
            Defaults to ``False``.
        wc (float): Unit weight of concrete in kg/m³.  Only used when
            *lightweight* is ``True``.  Defaults to 2320 kg/m³ (~145 pcf).

    Returns:
        float: Minimum slab thickness *h* in mm.

    Raises:
        ValueError: If *support_condition* is not one of the recognised keys
            in :data:`THICKNESS_RATIOS`.
    """
    if support_condition not in THICKNESS_RATIOS:
        raise ValueError(
            f"Unknown support condition '{support_condition}'. "
            f"Must be one of: {list(THICKNESS_RATIOS.keys())}."
        )

    ratio = THICKNESS_RATIOS[support_condition]
    h = span / ratio

    # Sec. 7.3.1.1.1 — fy adjustment (applies when fy != 60 ksi)
    fy_psi = fy * 145.038
    if not (59000.0 <= fy_psi <= 61000.0):
        h *= 0.4 + fy_psi / 100000.0

    # Sec. 7.3.1.1.2 — lightweight concrete adjustment
    if lightweight:
        wc_pcf = wc / 16.0185
        h *= max(1.65 - 0.005 * wc_pcf, 1.09)

    return h


# ---------------------------------------------------------------------------
# Shrinkage and temperature reinforcement
# ---------------------------------------------------------------------------


def As_shrinkage_temperature(
    fy: float,
    b: float,
    h: float,
) -> float:
    """Minimum shrinkage and temperature reinforcement area for a slab.

    ACI 318-25, Sec. 24.4.3.2.  The minimum steel ratio depends on *fy*
    in psi:

    - *fy* <= 50 000 psi : rho_min = 0.0020
    - *fy* <= 60 000 psi : rho_min = 0.0018
    - *fy*  > 60 000 psi : rho_min = max(0.0014, 0.0018 × 60 000 / fy_psi)

    Args:
        fy (float): Specified yield strength of reinforcement in MPa.
        b (float): Width of slab strip in mm.
        h (float): Overall slab thickness in mm.

    Returns:
        float: Minimum shrinkage and temperature steel area in mm².
    """
    fy_psi = fy * 145.038  # 1 MPa = 145.038 psi
    if fy_psi <= 50000.0:
        ratio = 0.0020
    elif fy_psi <= 60000.0:
        ratio = 0.0018
    else:
        ratio = max(0.0014, 0.0018 * 60000.0 / fy_psi)
    return ratio * b * h


# ---------------------------------------------------------------------------
# Bar spacing limits
# ---------------------------------------------------------------------------


def max_bar_spacing_flexure(h: float) -> float:
    """Maximum centre-to-centre spacing of flexural reinforcement in a slab.

    ACI 318-25, Sec. 7.7.2.3:  ``s_max = min(3h, 450 mm)``.

    Args:
        h (float): Overall slab thickness in mm.

    Returns:
        float: Maximum bar spacing in mm.
    """
    return min(3.0 * h, 450.0)


def max_bar_spacing_shrinkage(h: float) -> float:
    """Maximum centre-to-centre spacing of shrinkage and temperature reinforcement.

    ACI 318-25, Sec. 7.7.6.2.1:  ``s_max = min(5h, 450 mm)``.

    Args:
        h (float): Overall slab thickness in mm.

    Returns:
        float: Maximum bar spacing in mm.
    """
    return min(5.0 * h, 450.0)


# ---------------------------------------------------------------------------
# Shear critical section
# ---------------------------------------------------------------------------


def shear_critical_section_offset(d: float) -> float:
    """Distance from face of support to the shear critical section.

    ACI 318-25, Sec. 7.4.3.2.  For non-prestressed slabs the critical section
    for shear is located at a distance *d* from the face of the support.

    Args:
        d (float): Effective depth of the slab in mm.

    Returns:
        float: Critical section offset in mm (equal to *d*).
    """
    return d
