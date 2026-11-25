"""Two-way shear strength functions according to ACI 318-19, Sec. 22.6.

Units: every function takes and returns SI units (MPa, mm). The expressions
of ACI 318-19 are written in inch-pound units, with empirical coefficients
on sqrt(f'c) in psi. Each function applies the Code expression in that form
and converts with exact factors, so that every coefficient is the one
printed in the Code.
"""

import typing as t

from ._shear import _sqrt_fc, lambda_s

# Sec. 22.6.3.1: limit on sqrt(f'c) for two-way shear, in psi
_SQRT_FC_LIMIT_PSI = 100.0

# Sec. 22.6.5.3
_ALPHA_S = {'interior': 40.0, 'edge': 30.0, 'corner': 20.0}


def alpha_s(location: t.Literal['interior', 'edge', 'corner']) -> float:
    """Constant alpha_s used to calculate vc for two-way shear.

    ACI 318-19, Sec. 22.6.5.3: 40 for interior columns, 30 for edge columns
    and 20 for corner columns. The commentary to the provision relates the
    three terms to a critical section with the slab continuous on four,
    three and two sides.

    Args:
        location (str): One of 'interior', 'edge' or 'corner'.

    Returns:
        float: The value of alpha_s.

    Raises:
        ValueError: If location is not one of the three terms.
    """
    value = _ALPHA_S.get(location)
    if value is None:
        raise ValueError(
            f'location must be one of {list(_ALPHA_S)}, got {location!r}'
        )
    return value


def vc_two_way(
    fc: float,
    d: float,
    bo: float,
    beta: float,
    alpha_s: float,
    lambda_concrete: float = 1.0,
) -> float:
    """Two-way shear stress strength of a member without shear reinforcement.

    ACI 318-19, Table 22.6.5.2, with f'c in psi. vc is the least of:

    - (a) ``4 lambda_s lambda sqrt(f'c)``
    - (b) ``(2 + 4 / beta) lambda_s lambda sqrt(f'c)``
    - (c) ``(2 + alpha_s d / bo) lambda_s lambda sqrt(f'c)``

    sqrt(f'c) is not taken greater than 100 psi (Sec. 22.6.3.1). lambda_s
    is the size effect factor of Eq. (22.5.5.1.3).

    Note:
        Table 22.6.5.2 and Sec. 22.6.5.3 are the same in ACI 318-25. See
        ``structuralcodes.codes.aci318_25.lambda_s`` for a difference in
        how the two editions print Eq. (22.5.5.1.3).

    Args:
        fc (float): Specified compressive strength of concrete in MPa.
        d (float): Effective depth of the member in mm.
        bo (float): Perimeter of the critical section in mm.
        beta (float): Ratio of the long side to the short side of the
            column, concentrated load or reaction area. Not less than 1.
        alpha_s (float): Constant of Sec. 22.6.5.3.
        lambda_concrete (float): Lightweight concrete factor lambda.
            Defaults to 1.0.

    Returns:
        float: Stress strength vc in MPa, acting on the area bo d.

    Raises:
        ValueError: If fc, d or bo is not positive, or if beta is less
            than 1.
    """
    if bo <= 0.0:
        raise ValueError(f'bo must be positive, got {bo}')
    if beta < 1.0:
        raise ValueError(
            f'beta is the long side over the short side, got {beta}'
        )
    coefficient = min(4.0, 2.0 + 4.0 / beta, 2.0 + alpha_s * d / bo)
    return (
        coefficient
        * lambda_s(d)
        * lambda_concrete
        * _sqrt_fc(fc, _SQRT_FC_LIMIT_PSI)
    )
