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
    'D': {
        'epsuk': 8.0e-2,
        'k': 1.25,
    },
}


def fyd(fyk: float, gamma_s: float = 1.15) -> float:
    """Calculate the design value of the reinforcement yield strength.

    fib Model Code 2010, Sec. 4.5.2.2.3

    Args:
        fyk (float): The characteristic yield strength in MPa.
        gamma_s (float): The partial factor. Default value 1.15.

    Returns:
        float: The design yield strength in MPa.

    Raises:
        ValueError: if fyk is less than 0
        ValueError: if gamma_s is less than 1
    """
    if fyk < 0:
        raise ValueError(f'fyk={fyk} cannot be less than 0')
    if gamma_s < 1:
        raise ValueError(f'gamma_s={gamma_s} must be larger or equal to 1')
    return fyk / gamma_s


def reinforcement_duct_props(
    fyk: float,
    ductility_class: t.Literal['A', 'B', 'C', 'D'],
) -> t.Dict[str, float]:
    """Return a dict with the minimum characteristic ductility properties for
    reinforcement ductility class.

    fib Model Code 2010, Sec. 5.2.5.4

    Args:
        fyk (float): The characteristic yield strength.
        ductility_class (Literal['A', 'B', 'C', 'D']): The reinforcement
            ductility class designation.

    Returns:
        Dict[str, float]: A dict with the characteristik strain value at the
        ultimate stress level (epsuk), and the characteristic ultimate stress
        (ftk).

    Raises:
        ValueError: when the ductility_class does not define a valid ductility
            class
    """
    duct_props = DUCTILITY_CLASSES.get(ductility_class.upper(), None)
    if duct_props is None:
        raise ValueError(
            'No properties were found for the provided ductility class '
            f'({ductility_class}).'
        )
    return {
        'epsuk': duct_props['epsuk'],
        'ftk': duct_props['k'] * fyk,
    }
