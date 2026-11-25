"""Critical sections for two-way shear according to ACI 318-19, Sec. 22.6.4.

The functions of this module are geometric. Any consistent length unit may
be used. The rest of the ACI 318-19 functions take mm.
"""

import math
import typing as t

from shapely.geometry import LineString, MultiLineString, Point, Polygon, box

_FACES = ('+x', '-x', '+y', '-y')

# Sec. 22.6.4.3: reach of the opening provision as a multiple of h
_OPENING_REACH = 4.0

# Number of chords used to close the far side of a wedge
_WEDGE_STEPS = 8

Opening = t.Union[Polygon, t.Sequence[t.Tuple[float, float]]]


def _sides(
    c1: float,
    c2: float,
    d: float,
    edge_distances: t.Mapping[str, float],
) -> t.List[LineString]:
    """Return the straight sides of the section before openings are applied.

    Args:
        c1 (float): Column dimension along x.
        c2 (float): Column dimension along y.
        d (float): Effective depth of the slab.
        edge_distances (Mapping[str, float]): Distance from a column face to
            a free slab edge, by face name.

    Returns:
        List[LineString]: One line for each side that is not a free edge.
    """
    half = {'+x': c1 / 2.0, '-x': c1 / 2.0, '+y': c2 / 2.0, '-y': c2 / 2.0}
    # Distance from the column centroid to the section, or to the slab edge
    reach = {
        face: half[face] + edge_distances.get(face, d / 2.0) for face in _FACES
    }
    x0, x1 = -reach['-x'], reach['+x']
    y0, y1 = -reach['-y'], reach['+y']
    lines = {
        '+x': LineString([(x1, y0), (x1, y1)]),
        '+y': LineString([(x1, y1), (x0, y1)]),
        '-x': LineString([(x0, y1), (x0, y0)]),
        '-y': LineString([(x0, y0), (x1, y0)]),
    }
    return [line for face, line in lines.items() if face not in edge_distances]


def _wedge(opening: Polygon, radius: float) -> Polygon:
    """Return the region between the tangent lines from the origin.

    Args:
        opening (Polygon): The opening, clear of the origin.
        radius (float): Distance from the origin to which the wedge is
            drawn.

    Returns:
        Polygon: The wedge enclosed by the two lines from the origin tangent
        to the opening.

    Raises:
        ValueError: If the origin is inside or on the boundary of the
            convex hull of the opening.
    """
    hull = opening.convex_hull
    if hull.intersects(Point(0.0, 0.0)):
        raise ValueError(
            'The centroid of the column is inside an opening or its convex '
            'hull, so no tangent lines can be drawn.'
        )
    angles = sorted(
        math.atan2(y, x) % (2.0 * math.pi) for x, y in hull.exterior.coords
    )
    # The wedge is the complement of the largest angular gap between
    # vertices of the hull.
    gaps = [
        (angles[(i + 1) % len(angles)] - a) % (2.0 * math.pi)
        for i, a in enumerate(angles)
    ]
    i_gap = max(range(len(gaps)), key=gaps.__getitem__)
    start = angles[(i_gap + 1) % len(angles)]
    width = 2.0 * math.pi - gaps[i_gap]

    points = [(0.0, 0.0)]
    for step in range(_WEDGE_STEPS + 1):
        angle = start + width * step / _WEDGE_STEPS
        points.append((radius * math.cos(angle), radius * math.sin(angle)))
    return Polygon(points)


def critical_section_rectangular(
    c1: float,
    c2: float,
    d: float,
    edge_distances: t.Optional[t.Mapping[str, float]] = None,
    openings: t.Optional[t.Sequence[Opening]] = None,
    h: t.Optional[float] = None,
) -> MultiLineString:
    """Critical section for two-way shear at a rectangular column.

    ACI 318-19, Sec. 22.6.4.1 and 22.6.4.1.1: the section is located at d/2
    from the column faces and is defined with straight sides.

    Sec. 22.6.4.3: where an opening is closer than 4h to the periphery of
    the column, the portion of the section enclosed by the straight lines
    from the centroid of the column tangent to the opening is ineffective
    and is removed.

    The section is drawn in plan with the origin at the centroid of the
    column, c1 along x and c2 along y. The perimeter bo is the ``length``
    of the returned geometry.

    Note:
        Sec. 22.6.4.1 requires the section that makes bo a minimum. This
        function does not search for it. Whether a nearby slab edge, or an
        opening long enough to act as one, is treated as a free edge is
        for the caller to decide, by comparing the perimeters with and
        without it.

    Args:
        c1 (float): Column dimension along x.
        c2 (float): Column dimension along y.
        d (float): Effective depth of the slab.
        edge_distances (Optional[Mapping[str, float]]): Free slab edges.
            Each key is the column face the edge is beyond, one of '+x',
            '-x', '+y' or '-y'. Each value is the distance from that
            column face to the slab edge, 0 for a column flush with the
            edge. The section has no side on that face and the adjacent
            sides run to the slab edge.
        openings (Optional[Sequence]): Openings in the slab, each a shapely
            Polygon or a sequence of (x, y) points.
        h (Optional[float]): Overall thickness of the slab. Required if
            openings are given.

    Returns:
        MultiLineString: The effective straight parts of the section.

    Raises:
        ValueError: If a dimension is not positive, a face name is not
            recognised, an edge distance is negative, h is missing when
            openings are given, the column centroid is inside an opening,
            or no part of the section remains.
    """
    if c1 <= 0.0 or c2 <= 0.0 or d <= 0.0:
        raise ValueError(f'c1, c2 and d must be positive, got {c1}, {c2}, {d}')
    edge_distances = dict(edge_distances or {})
    unknown = [face for face in edge_distances if face not in _FACES]
    if unknown:
        raise ValueError(f'Faces must be named from {_FACES}, got {unknown}')
    if any(distance < 0.0 for distance in edge_distances.values()):
        raise ValueError('An edge distance cannot be negative.')
    openings = list(openings or [])
    if openings and (h is None or h <= 0.0):
        raise ValueError(
            'The slab thickness h is required to apply Sec. 22.6.4.3.'
        )

    sides = _sides(c1, c2, d, edge_distances)
    if not sides:
        raise ValueError('No side of the critical section remains.')

    column = box(-c1 / 2.0, -c2 / 2.0, c1 / 2.0, c2 / 2.0)
    extent = max(math.hypot(x, y) for side in sides for x, y in side.coords)
    for opening in openings:
        polygon = opening if isinstance(opening, Polygon) else Polygon(opening)
        if polygon.distance(column) >= _OPENING_REACH * h:
            continue
        wedge = _wedge(polygon, radius=100.0 * extent)
        sides = [side.difference(wedge) for side in sides]

    parts: t.List[LineString] = []
    for side in sides:
        pieces = getattr(side, 'geoms', [side])
        parts.extend(
            piece
            for piece in pieces
            if isinstance(piece, LineString) and piece.length > 0.0
        )
    if not parts:
        raise ValueError('No part of the critical section remains.')
    return MultiLineString(parts)
