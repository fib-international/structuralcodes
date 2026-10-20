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


def fyd(fyk: float, gamma_s: float) -> float:
    """Calculate the design value of the reinforcement yield strength.

    EUROCODE 2 1992-1-1:2004, Fig. 3.8

    Args:
        fyk (float): The characteristic yield strength in MPa.
        gamma_s (float): The partial factor.

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


def epsud(epsuk: float, gamma_eps: float = 0.9) -> float:
    """Calculate the design value of the reinforcement ultimate strain.

    EUROCDE 2 1992-1-1:2004, Fig 3.8

    Args:
        epsuk (float): The characteristic ultimate strain

    Keyword Args:
        gamma_eps (float): The partial factor specified in NA.
            Default value 0.9.

    Returns:
        float: The design ultimate strain

    Raises:
        ValueError: if epsuk is less than 0
        ValueError: if gamma_eps is greater than 1
    """
    if epsuk < 0:
        raise ValueError(f'epsuk={epsuk} cannot be less than 0')
    if gamma_eps > 1:
        raise ValueError(
            f'gamma_eps={gamma_eps} must be smaller or equal to 1'
        )
    return epsuk * gamma_eps


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
