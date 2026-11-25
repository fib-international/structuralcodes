"""One-way shear functions according to ACI 318-19, Sec. 22.5.

Units: every function takes and returns SI units (MPa, mm, mm2, N), as the
rest of the library does. The expressions of ACI 318-19 are written in
inch-pound units, with empirical coefficients on sqrt(f'c) in psi. Each
function applies the Code expression in that form and converts with exact
factors, so that every coefficient is the one printed in the Code.
"""

import math
import typing as t

# Exact conversion factors
_MPA_PER_PSI = 0.006894757293168361
_MM_PER_IN = 25.4


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


def lambda_s(d: float) -> float:
    """Size effect modification factor for shear.

    ACI 318-19, Eq. (22.5.5.1.3), with d in inches:
    ``lambda_s = sqrt(2 / (1 + d / 10)) <= 1``.

    Note:
        ACI 318-25 prints this equation without the square root. See
        ``structuralcodes.codes.aci318_25.lambda_s``.

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
