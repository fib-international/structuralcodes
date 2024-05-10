"""Material properties for reinforcement steel."""

import typing as t

DUCTILITY_CLASSES = {
    'A': {
        'epsuk': 2.5e-2,
        'k': 1.05,
    },
    'B': {
        'epsuk': 5.0e-2,
        'k': 1.08,
    },
    'C': {
        'epsuk': 7.5e-2,
        'k': 1.15,
    },
}


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


def reinforcement_duct_props(
    fyk: float,
    ductility_class: t.Literal['A', 'B', 'C'],
) -> t.Dict[str, float]:
    """Return a dict with the minimum characteristic ductility properties for
    reinforcement ductility class.

    EUROCODE 2 1992-1-1:2004, Tab. C.1

    Args:
        fyk (float): The characteristic yield strength.
        ductility_class (Literal['A', 'B', 'C']): The reinforcement ductility
            class designation.

    Returns:
        Dict[str, float]: A dict with the characteristik strain value at the
        ultimate stress level (epsuk), and the characteristic ultimate stress
        (ftk).
    """
    duct_props = DUCTILITY_CLASSES.get(ductility_class.upper(), None)
    if duct_props is None:
        raise ValueError(
            'The no properties was found for the provided ductility class '
            f'({ductility_class}).'
        )
    return {
        'epsuk': duct_props['epsuk'],
        'ftk': duct_props['k'] * fyk,
    }
