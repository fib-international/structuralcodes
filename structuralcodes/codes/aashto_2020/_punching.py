# Functions for AASHTO LRFD 2020 9th Edition Shear Punching Design
import math


def b0_edge(W: float, L: float, df: float) -> float:
    """Determines the critical perimeter for a bearing on the edge of the beam.

    AASHTO LRFD 2020 9th Edition, Eq. (5.8.4.4-4)

    Args:
        W (float) = width of the bearing plate or pad (in)
        L (float) = length of the bearing pad (in)
        df (float) = distance from the top of the ledge to the bottom
        longitudinal reinforcement (in)

    Returns:
        The critical perimeter for a bearing on the edge of the beam in (in)

    Raises:
        ValueError: If W is less than 0
        ValueError: If L is less than 0
        ValueError: If df is less than 0
    """
    if W < 0:
        raise ValueError(f'W={W} cannot be less than 0')
    if L < 0:
        raise ValueError(f'L={L} cannot be less than 0')
    if df < 0:
        raise ValueError(f'df={df} cannot be less than 0')

    return W + 2 * L + 2 * df


def b0_corner(W: float, L: float, df: float, c: float) -> float:
    """Determines the critical perimeter for a bearing in the corner of
    the beam.

    AASHTO LRFD 2020 9th Edition, Eq. (5.8.4.4-5)

    Args:
        W (float) = width of the bearing plate or pad (in)
        L (float) = length of the bearing pad (in)
        df (float) = distance from the top of the ledge to the bottom
        longitudinal reinforcement (in)
        c (float) = spacing from the centerline of the bearing to the end of
        beam ledge (in)

    Returns:
        The critical perimeter for a bearing in the corner of the beam in (in)

    Raises:
        ValueError: If W is less than 0
        ValueError: If L is less than 0
        ValueError: If df is less than 0
        ValueError: If c is less than 0
    """
    if W < 0:
        raise ValueError(f'W={W} cannot be less than 0')
    if L < 0:
        raise ValueError(f'L={L} cannot be less than 0')
    if df < 0:
        raise ValueError(f'df={df} cannot be less than 0')
    if c < 0:
        raise ValueError(f'c={c} cannot be less than 0')

    return 0.5 * W + L + df + c


def b0_interior(W: float, L: float, df: float) -> float:
    """Determines the critical perimeter for a bearing in the interior of
    the beam.

    Was interpolated based on the edge bearing formula

    Args:
        W (float) = width of the bearing plate or pad (in)
        L (float) = length of the bearing pad (in)
        df (float) = distance from the top of the ledge to the bottom
        longitudinal reinforcement (in)

    Returns:
        The critical perimeter for a bearing in the interior of the beam
        in (in)

    Raises:
        ValueError: If W is less than 0
        ValueError: If L is less than 0
        ValueError: If df is less than 0
    """
    if W < 0:
        raise ValueError(f'W={W} cannot be less than 0')
    if L < 0:
        raise ValueError(f'L={L} cannot be less than 0')
    if df < 0:
        raise ValueError(f'df={df} cannot be less than 0')

    return 2 * W + 2 * L + 4 * df


def Vn(fc_prime: float, b0: float, df: float) -> float:
    """Determines the nominal punching shear resistance.

    Args:
        fc_prime (float): compressive strength of concrete in ksi
        b0: the critical perimeter of the bearing in (in)
        df: distance from the top of the ledge to the bottom
        longitudinal reinforcement (in)

    Returns:
        The nominal punching shear resitance in kips

    Raises:
        ValueError: If fc_prime is less than 0
        ValueError: If b0 is less than 0
        ValueError: If df is less than 0
    """
    return 0.125 * math.sqrt(fc_prime) * b0 * df
