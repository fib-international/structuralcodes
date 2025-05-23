"""Functions related to reinforcement definition."""

import math
import typing as t

import numpy as np
from shapely import Point

from structuralcodes.core.base import ConstitutiveLaw, Material

from ._geometry import CompoundGeometry, PointGeometry, SurfaceGeometry


def add_reinforcement(
    geo: t.Union[SurfaceGeometry, CompoundGeometry],
    coords: t.Tuple[float, float],
    diameter: float,
    material: t.Union[Material, ConstitutiveLaw],
    group_label: t.Optional[str] = None,
) -> CompoundGeometry:
    """Add a single bar given coordinate.

    Arguments:
        geo (Union(SurfaceGeometry, CompoundGeometry)): A geometry to which add
            reinforcement.
        coords (Tuple(float, float)): A tuple with cordinates of bar.
        diameter (float): The diameter of the reinforcement.
        material (Union(Material, ConstitutiveLaw)): A material or a
            constitutive law for the behavior of the reinforcement.
        group_label (Optional(str)): A label for grouping several objects
            (default is None).

    Returns:
        CompoundGeometry: A compound geometry with the original geometry and
        the reinforcement.
    """
    bar = PointGeometry(
        Point(coords), diameter, material, group_label=group_label
    )
    return geo + bar


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
    group_label: t.Optional[str] = None,
) -> CompoundGeometry:
    """Adds a set of bars distributed in a line.

    Arguments:
        geo (Union(SurfaceGeometry, CompoundGeometry)): The geometry used as
            input.
        coords_i (Tuple(float, float)): Coordinates of the initial point of
            line.
        coords_j (Tuple(float, float)): Coordinates of the final point of line.
        diamter (float): The diameter of the bars.
        material (Union(Material, ConstitutiveLaw)): A valid material or
            constitutive law.
        n (int): The number of bars to be distributed inside the line (default
            = 0).
        s (float): The distance between the bars (default = 0).
        first (bool): Boolean indicating if placing the first bar (default =
            True).
        last (bool): Boolean indicating if placing the last bar (default =
            True).
        group_label (Optional(str)): A label for grouping several objects
            (default is None).

    Note:
        At least n or s should be greater than zero.

    Returns:
        CompoundGeometry: A compound geometry with the original geometry and
        the reinforcement.
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
            geo,
            (coords[0], coords[1]),
            diameter,
            material,
            group_label=group_label,
        )
    return geo


def add_reinforcement_circle(
    geo: t.Union[SurfaceGeometry, CompoundGeometry],
    center: t.Tuple[float, float],
    radius: float,
    diameter: float,
    material: t.Union[Material, ConstitutiveLaw],
    n: int = 0,
    s: float = 0.0,
    first: bool = True,
    last: bool = True,
    start_angle: float = 0.0,
    stop_angle: float = 2 * np.pi,
    group_label: t.Optional[str] = None,
) -> CompoundGeometry:
    """Adds a set of bars distributed in a circular arch line.
    By default the whole circle is considered. If one wants to specify a
    circular arch, the `start_angle` and `stop_angle` attributes need to be
    specified.

    Arguments:
        geo (Union(SurfaceGeometry, CompoundGeometry)): The geometry used as
            input.
        center (Tuple(float, float)): Coordinates of the center point of
            the circle line where reinforcement will be added.
        radius (float): Radius of the circle line where reinforcement will be
            added.
        diameter (float): The diameter of the bars.
        material (Union(Material, ConstitutiveLaw)): A valid material or
            constitutive law.
        n (int): The number of bars to be distributed inside the line (default
            = 0).
        s (float): The distance between the bars (default = 0).
        first (bool): Boolean indicating if placing the first bar (default =
            True).
        last (bool): Boolean indicating if placing the last bar (default =
            True).
        start_angle (float): Start angle (respect to X axis) for defining the
            arch where to add bars in radians (default = 0)
        stop_angle (float): Stop angle (respect to X axis) for defining the
            arch where to add bars in radians (default = 2pi)
        group_label (Optional(str)): A label for grouping several objects
            (default is None).

    Note:
        At least n or s should be greater than zero.
        Attribues start_angle and stop_angle by default are 0 and 2pi
        respectively, so that bars will be distributed along the whole circle.
        If only a portion of the circle must be used (i.e. an arch of
        circumference), then set start and stop angles correspondingly.
        stop_angle must always be larger than start_angle.

    Returns:
        CompoundGeometry: A compound geometry with the original geometry and
        the reinforcement.
    """
    # Check that difference between stop and start angle is
    # positive and less than 2pi
    if stop_angle - start_angle <= 0 or stop_angle - start_angle > 2 * np.pi:
        raise ValueError(
            'Stop angle should be larger than start angle and difference \
            them should be at most 2pi.'
        )
    # Calculate length from start_angle to stop_angle
    length = radius * (stop_angle - start_angle)

    # If the whole circle, than deal with the case that would add an extra bar
    whole = math.isclose(length - 2 * np.pi * radius, 0, abs_tol=1e-4)
    add_n = 0 if not whole else 1

    # delta_angle is used if we need to center the set of bars
    # in the curve.
    delta_angle = 0
    if n > 0 and s > 0:
        # Provided both the number of bars and spacing
        # Check there is enough space for fitting the bars
        n += add_n
        needed_length = (n - 1) * s
        if needed_length > length:
            raise ValueError(
                f'There is not room to fit {n} bars with a spacing of {s} \
                in {length}'
            )
        # Compute delta_angle to make bars centered in the curvilinear segment
        delta_angle = (length - needed_length) / 2.0 / radius
    elif n > 0:
        # Provided the number of bars
        s = length / (n - 1)
        # If we are distributing bars i the whole circle add a fictitious extra
        # bar (than later will be removed).
        n += add_n
    elif s > 0:
        # Provided the spacing
        # 1. Compute the number of bars
        n = math.floor(length / s) + 1
        # 2. Distribute the bars centered in the curvilinear segment
        needed_length = (n - 1) * s
        delta_angle = (length - needed_length) / 2.0 / radius
        # set whole to False bacause in this case we don't need to deal with
        # the special case
        whole = False
    else:
        raise ValueError('At least n or s should be provided')

    phi_rebars = np.linspace(
        start_angle + delta_angle, stop_angle - delta_angle, n
    )
    if whole:
        n -= 1
        phi_rebars = phi_rebars[:-1]

    x = center[0] + radius * np.cos(phi_rebars)
    y = center[1] + radius * np.sin(phi_rebars)

    # add the bars
    for i in range(n):
        if i == 0 and not first:
            continue
        if i == n - 1 and not last:
            continue
        geo = add_reinforcement(
            geo,
            (x[i], y[i]),
            diameter,
            material,
            group_label=group_label,
        )
    return geo
