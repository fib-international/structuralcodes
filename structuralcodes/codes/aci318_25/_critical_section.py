"""Critical sections for two-way shear according to ACI 318-25, Sec. 22.6.4.

Sec. 22.6.4.1, 22.6.4.1.1 and 22.6.4.3 of ACI 318-25 are worded as in
ACI 318-19, so the geometry is evaluated by the ACI 318-19 function.
"""

import typing as t

from shapely.geometry import MultiLineString

from ..aci318_19 import _critical_section as _aci318_19


def critical_section_rectangular(
    c1: float,
    c2: float,
    d: float,
    edge_distances: t.Optional[t.Mapping[str, float]] = None,
    openings: t.Optional[t.Sequence[_aci318_19.Opening]] = None,
    h: t.Optional[float] = None,
) -> MultiLineString:
    """Critical section for two-way shear at a rectangular column.

    ACI 318-25, Sec. 22.6.4.1 and 22.6.4.1.1: the section is located at d/2
    from the column faces and is defined with straight sides.

    Sec. 22.6.4.3: where an opening is closer than 4h to the periphery of
    the column, the portion of the section enclosed by the straight lines
    from the centroid of the column tangent to the opening is ineffective
    and is removed.

    The section is drawn in plan with the origin at the centroid of the
    column, c1 along x and c2 along y. The perimeter bo is the ``length``
    of the returned geometry. Any consistent length unit may be used.

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
    return _aci318_19.critical_section_rectangular(
        c1, c2, d, edge_distances=edge_distances, openings=openings, h=h
    )
