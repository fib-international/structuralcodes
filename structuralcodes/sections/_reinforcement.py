"""Functions related to reinforcement definition."""

import typing as t

import numpy as np
from shapely import Point

from structuralcodes.core.base import ConstitutiveLaw, Material
from structuralcodes.geometry import (
    CompoundGeometry,
    PointGeometry,
    SurfaceGeometry,
)


def add_reinforcement(
    geo: t.Union[SurfaceGeometry, CompoundGeometry],
    coords: t.Tuple[float, float],
    diameter: float,
    material: t.Union[Material, ConstitutiveLaw],
) -> CompoundGeometry:
    """Add a single bar given coordinate.

    Args:
        geo: a geometry (SurfaceGeometry or CompoundGeometry) to
            which add reinforcement
        coords: a tuple with cordinates of bar
        diameter: the diameter of the reinforcement
        material: a material or a constitutive law for the behavior
            of the reinforcement

    Returns:
        the resulting compound geometry
    """
    bar = PointGeometry(Point(coords), diameter, material)
    return CompoundGeometry([geo, bar])


def add_reinforcement_line(
    geo: t.Union[SurfaceGeometry, CompoundGeometry],
    coords_i: t.Tuple[float, float],
    coords_j: t.Tuple[float, float],
    diameter: float,
    material: t.Union[Material, ConstitutiveLaw],
    n: int = 0,
    s: float = 0.0,
    first: bool = True,
    last: bool = True,
) -> CompoundGeometry:
    """Adds a set of bars distributed in a line.

    Args:
        geo: the geometry used as input
        coords_i: coordinates of the initial point of line
        coords_j: coordinates of the final point of line
        diamter: the diameter of the bars
        material: a valid material or constitutive law
        n: the number of bars to be distributed inside the line (default = 0)
        s: the distance between the bars (default = 0)
        first: boolean indicating if placing the first bar (default = True)
        last: boolean indicating if placing the last bar (default = True)

    Note:
        at least n or s should be greater than zero

    Returns:
        compound geometry
    """
    from math import floor

    p1 = np.array(coords_i)
    p2 = np.array(coords_j)
    distance = np.linalg.norm(p2 - p1)
    v = (p2 - p1) / distance
    if n > 0 and s > 0:
        # Provided both the number of bars and spacing
        # 1. Check the is enough space for fitting the bars
        d = (n - 1) * s
        if d > distance:
            raise ValueError(
                f'There is not room to fit {n} bars with a spacing of {s} \
                in {distance}'
            )
        p1 = p1 + v * (distance - d) / 2.0
    elif n > 0:
        # Provided the number of bars
        s = distance / (n - 1)
    elif s > 0:
        # Provided the spacing
        # 1. Compute the number of bars
        n = floor(distance / s) + 1
        # 2. Distribute the bars centered in the segment
        d = (n - 1) * s
        p1 = p1 + v * (distance - d) / 2.0
    else:
        raise ValueError('At least n or s should be provided')
    # add the bars
    for i in range(n):
        if i == 0 and not first:
            continue
        if i == n - 1 and not last:
            continue
        coords = p1 + v * s * i
        geo = add_reinforcement(
            geo, (coords[0], coords[1]), diameter, material
        )
    return geo
