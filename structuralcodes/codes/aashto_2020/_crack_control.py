## AASHTO LRFD Functions for Control of Cracking ##


def beta_s(h: float, dc) -> float:
    """Determines the flexural strain ratio.

    AASHTO LRFD 2020 9th Edition, Eq (5.6.7-2)

    Args:
        h (float): height of the cross section in (in)
        d (float): the thickness of the concrete cover from the surface
        to the center of the reinforcement in (in)

    Returns:
        The flexural strain ratio

    Raises:
        ValueError: If h is less than 0
        ValueError: If dc is less than 0
    """
    if h < 0:
        raise ValueError(f'h={h} cannot be less than 0')
    if dc < 0:
        raise ValueError(f'd={dc} cannot be less than 0')

    return 1 + (dc / (0.7 * (h - dc)))


def s(fyk: float, beta_s: float, gamma_e: float, dc: float) -> float:
    """Determines the spacing of nonprestressed reinforcement.

    AASHTO LRFD 2020 9th Edition, Eq. (5.6.7-1)

    Args:
        fyk (float): The steel reinforcement yielding strength in ksi
        beta_s (float): The flexural strain ratio
        gamma_e (float): The exposure factor
        dc (float): The thickness of concrete cover from surface to center
        of the reinforcement in (in)

    Returns:
        The spacing of nonprestressed reinforcement in (in)

    Raises:
        ValueError: If fyk is less than 0
        ValueError: If gamma_e is less than 0
        ValueError: If dc is less than 0
    """
    if fyk < 0:
        raise ValueError(f'fyk={fyk} cannot be less than 0')
    if gamma_e < 0:
        raise ValueError(f'gamma_e={gamma_e} cannot be less than 0')
    if dc < 0:
        raise ValueError(f'dc={dc} cannot be less than 0')
