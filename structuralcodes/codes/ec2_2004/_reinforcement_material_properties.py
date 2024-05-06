"""Material properties for reinforcement steel."""


def fyd(fyk: float, gamma_s: float = 1.15) -> float:
    """Calculate the design value of the reinforcement yield strength.

    EUROCODE 2 1992-1-1:2004, Fig. 3.8

    Args:
        fyk (float): The characteristic yield strength in MPa.
        gamma_s (float): The partial factor.

    Returns:
        float: The design yield strength in MPa.

    Raises:
        ValueError: if fyk is less than 0
        ValueError: if gamma_s is less than or equal to 0
    """
    if fyk < 0:
        raise ValueError(f'fyk={fyk} cannot be less than 0')
    if gamma_s <= 0:
        raise ValueError(f'gamma_s={gamma_s} must be larger than 0')
    return fyk / gamma_s
