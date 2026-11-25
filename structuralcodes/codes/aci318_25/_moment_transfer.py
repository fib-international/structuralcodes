"""Two-way shear with moment transfer according to ACI 318-25, Sec. 8.4.4.2.

The expressions of this module are dimensionally consistent, so any
consistent set of units may be used. With the rest of the ACI 318-25
functions in mm, N and MPa, moments are in N mm and Jc is in mm4.

The critical section is the geometry returned by
``critical_section_rectangular``: straight parts in plan, each of depth d.
"""

import math
import typing as t

from shapely.geometry import MultiLineString

Direction = t.Literal['x', 'y']


def gamma_f(b1: float, b2: float) -> float:
    """Fraction of Msc transferred by flexure.

    ACI 318-25, Eq. (8.4.2.2.2):
    ``gamma_f = 1 / (1 + (2/3) sqrt(b1 / b2))``.

    Args:
        b1 (float): Dimension of the critical section measured in the
            direction of the span for which moments are determined.
        b2 (float): Dimension of the critical section measured
            perpendicular to b1.

    Returns:
        float: The fraction gamma_f.

    Raises:
        ValueError: If b1 or b2 is not positive.
    """
    if b1 <= 0.0 or b2 <= 0.0:
        raise ValueError(f'b1 and b2 must be positive, got {b1}, {b2}')
    return 1.0 / (1.0 + (2.0 / 3.0) * math.sqrt(b1 / b2))


def gamma_v(b1: float, b2: float) -> float:
    """Fraction of Msc transferred by eccentricity of shear.

    ACI 318-25, Eq. (8.4.4.2.2): ``gamma_v = 1 - gamma_f``, with gamma_f
    from Eq. (8.4.2.2.2). The modified values of gamma_f permitted by
    Sec. 8.4.2.2.4 are not considered here.

    Args:
        b1 (float): Dimension of the critical section measured in the
            direction of the span for which moments are determined.
        b2 (float): Dimension of the critical section measured
            perpendicular to b1.

    Returns:
        float: The fraction gamma_v.

    Raises:
        ValueError: If b1 or b2 is not positive.
    """
    return 1.0 - gamma_f(b1, b2)


def _segments(
    section: MultiLineString,
) -> t.List[t.Tuple[t.Tuple[float, float], t.Tuple[float, float]]]:
    """Return the straight segments of a critical section.

    Args:
        section (MultiLineString): The critical section.

    Returns:
        List: The end points of each straight segment.

    Raises:
        ValueError: If the section has no length.
    """
    segments = []
    for line in section.geoms:
        coords = list(line.coords)
        segments.extend(
            ((xa, ya), (xb, yb))
            for (xa, ya), (xb, yb) in zip(coords[:-1], coords[1:])
            if math.hypot(xb - xa, yb - ya) > 0.0
        )
    if not segments:
        raise ValueError('The critical section has no length.')
    return segments


def critical_section_centroid(
    section: MultiLineString,
) -> t.Tuple[float, float]:
    """Centroid of a critical section in plan.

    ACI 318-25, Sec. 8.4.4.2.2 and 8.4.4.2.3 refer the transferred moment
    and the variation of shear stress to the centroid of the critical
    section. The depth of the section is constant, so the centroid is that
    of the lines in plan.

    Args:
        section (MultiLineString): The critical section.

    Returns:
        Tuple[float, float]: The coordinates (x, y) of the centroid.

    Raises:
        ValueError: If the section has no length.
    """
    total = sx = sy = 0.0
    for (xa, ya), (xb, yb) in _segments(section):
        length = math.hypot(xb - xa, yb - ya)
        total += length
        sx += length * (xa + xb) / 2.0
        sy += length * (ya + yb) / 2.0
    return sx / total, sy / total


def Jc(section: MultiLineString, d: float, direction: Direction) -> float:
    """Property of the critical section analogous to polar moment of inertia.

    ACI 318-25, R8.4.4.2.3 gives Jc for an interior column and notes that
    similar expressions may be developed for other sections. This function
    evaluates Jc for any section made of straight parts of depth d, about
    the axis through the centroid of the section.

    With s the coordinate along ``direction`` measured from the centroid, a
    part of length L running from s_a to s_b contributes
    ``d L (s_a^2 + s_a s_b + s_b^2) / 3 + d^3 (s_b - s_a)^2 / (12 L)``.
    The first term is the second moment of the face of the part. The second
    is the term in d^3 that the commentary expression carries for the faces
    parallel to the direction. For an interior rectangular column the sum
    equals the commentary expression.

    Args:
        section (MultiLineString): The critical section.
        d (float): Effective depth of the slab.
        direction (str): 'x' or 'y', the direction along which the
            distances from the centroid are measured. This is the
            direction of the span for which the moment is determined.

    Returns:
        float: Jc, in length to the fourth power.

    Raises:
        ValueError: If d is not positive, the direction is not 'x' or 'y',
            or the section has no length.
    """
    if d <= 0.0:
        raise ValueError(f'd must be positive, got {d}')
    if direction not in ('x', 'y'):
        raise ValueError(f"direction must be 'x' or 'y', got {direction!r}")
    index = 0 if direction == 'x' else 1
    origin = critical_section_centroid(section)[index]
    total = 0.0
    for a, b in _segments(section):
        length = math.hypot(b[0] - a[0], b[1] - a[1])
        s_a, s_b = a[index] - origin, b[index] - origin
        total += d * length * (s_a**2 + s_a * s_b + s_b**2) / 3.0
        total += d**3 * (s_b - s_a) ** 2 / (12.0 * length)
    return total


def vu_max(
    Vu: float,
    section: MultiLineString,
    d: float,
    gamma_v_Msc_x: float = 0.0,
    gamma_v_Msc_y: float = 0.0,
) -> float:
    """Maximum factored two-way shear stress on a critical section.

    ACI 318-25, Sec. 8.4.4.2.1 to 8.4.4.2.3: the stress is the combination
    of ``vuv = Vu / (bo d)`` and the stress produced by gamma_v Msc, which
    is applied at the centroid of the critical section and varies linearly
    about it. At a point of the section at distance c from the centroid,
    measured along the direction of the moment arms, the stress is
    ``vuv + gamma_v Msc c / Jc`` (R8.4.4.2.3).

    Note:
        The moments are those at the centroid of the critical section,
        which is not the centroid of the column at an edge or corner
        column or where an opening has reduced the section.

        Sec. 8.4.4.2 is written for moment about one axis. Where moments
        are given for both directions this function adds the two linear
        distributions. That is an assumption of this function.

    Args:
        Vu (float): Factored shear force on the critical section.
        section (MultiLineString): The critical section.
        d (float): Effective depth of the slab.
        gamma_v_Msc_x (float): The moment gamma_v Msc for the span in the x
            direction, at the centroid of the critical section. A positive
            value increases the stress on the +x side of the centroid.
        gamma_v_Msc_y (float): The same for the span in the y direction. A
            positive value increases the stress on the +y side.

    Returns:
        float: The maximum stress vu over the section.

    Raises:
        ValueError: If d is not positive or the section has no length.
    """
    if d <= 0.0:
        raise ValueError(f'd must be positive, got {d}')
    segments = _segments(section)
    x0, y0 = critical_section_centroid(section)
    vuv = Vu / (section.length * d)
    slope_x = gamma_v_Msc_x / Jc(section, d, 'x') if gamma_v_Msc_x else 0.0
    slope_y = gamma_v_Msc_y / Jc(section, d, 'y') if gamma_v_Msc_y else 0.0
    # A linear distribution on straight parts has its maximum at an end.
    return max(
        vuv + slope_x * (x - x0) + slope_y * (y - y0)
        for segment in segments
        for x, y in segment
    )
