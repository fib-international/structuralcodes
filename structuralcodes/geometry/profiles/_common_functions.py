"""Common functions for creating profiles."""

import math

import numpy as np
from shapely import (
    LineString,
    Point,
    Polygon,
    get_geometry,
    polygonize,
    set_precision,
)
from shapely.affinity import scale, translate
from shapely.geometry.polygon import orient
from shapely.ops import linemerge, unary_union


def _create_I_section(h: float, b: float, tw: float, tf: float, r: float):
    """Returns a polygon for a I section."""
    # top flange
    top_flange = Polygon(
        [
            (-b / 2, -h / 2),
            (b / 2, -h / 2),
            (b / 2, -h / 2 + tf),
            (-b / 2, -h / 2 + tf),
        ]
    )
    # bottom flange
    bottom_flange = translate(top_flange, xoff=0, yoff=h - tf)
    web = Polygon(
        [
            (-tw / 2, -h / 2 + tf),
            (tw / 2, -h / 2 + tf),
            (tw / 2, h / 2 - tf),
            (-tw / 2, h / 2 - tf),
        ]
    )
    # fillets
    p = Point([tw / 2 + r, -h / 2 + tf + r]).buffer(r)
    s = Polygon(
        [
            (tw / 2, -h / 2 + tf),
            (tw / 2 + r, -h / 2 + tf),
            (tw / 2 + r, -h / 2 + tf + r),
            (tw / 2, -h / 2 + tf + r),
        ]
    )
    fillet = s.difference(p)
    p = Point([-tw / 2 - r, -h / 2 + tf + r]).buffer(r)
    s = Polygon(
        [
            (-tw / 2 - r, -h / 2 + tf),
            (-tw / 2, -h / 2 + tf),
            (-tw / 2, -h / 2 + tf + r),
            (-tw / 2 - r, -h / 2 + tf + r),
        ]
    )
    fillet = s.difference(p).union(fillet)
    fillet = translate(
        scale(fillet, 1, -1), xoff=0, yoff=h - 2 * tf - r
    ).union(fillet)

    # Estimate grid_size value
    # Tentative geometry (due to approximations can be a MultiPolygon)
    geom_trial = unary_union([fillet, top_flange, bottom_flange, web])
    # minx, miny, maxx, maxy
    bounds = geom_trial.bounds
    min_size = min(bounds[2] - bounds[0], bounds[3] - bounds[1])
    grid_size = 10 ** int(math.floor(math.log10(abs(min_size)))) * 1e-12
    # Create the geometry
    geometries = [
        set_precision(geometry, grid_size=grid_size)
        for geometry in [fillet, top_flange, bottom_flange, web]
    ]
    return orient(unary_union(geometries), 1)


def _create_taper_I_section(
    h: float,
    b: float,
    tw: float,
    tf: float,
    r1: float,
    r2: float,
    slope: float,
) -> Polygon:
    """Returns a shapely polygon representing a Taper Flange I Section."""
    # Create first part of line
    ls = [
        set_precision(
            LineString([[0, -h / 2], [b / 2, -h / 2]]), grid_size=1e-13
        )
    ]
    # Create first fillet
    xy = np.array(
        [
            [b / 2, -h / 2],
            [b / 2, -h / 2 + tf - (b / 4 * slope)],
            [b / 4, -h / 2 + tf],
        ]
    )
    ls.append(
        set_precision(
            LineString(xy).offset_curve(r2).offset_curve(-r2), grid_size=1e-13
        )
    )
    # Create second fillet
    xy = np.array(
        [
            [b / 4, -h / 2 + tf],
            [tw / 2, -h / 2 + tf + (b / 4 - tw / 2) * slope],
            [tw / 2, 0],
        ]
    )
    ls.append(
        set_precision(
            LineString(xy).offset_curve(-r1).offset_curve(r1), grid_size=1e-13
        )
    )

    # Merge filleted
    merged_ls = linemerge(ls)

    # mirror the parts
    merged_ls = linemerge(
        [
            merged_ls,
            translate(scale(merged_ls, -1, 1), -b / 2, 0),
        ]
    )
    merged_ls = linemerge(
        [
            merged_ls,
            translate(scale(merged_ls, 1, -1), 0, h / 2),
        ]
    )

    # Create a polygon
    poly = polygonize([merged_ls])

    # Return the first and only polygon of this collection
    return orient(get_geometry(poly, 0), 1)


def _create_taper_U_section(
    h: float,
    b: float,
    tw: float,
    tf: float,
    r1: float,
    r2: float,
    slope: float,
    u: float,
) -> Polygon:
    """Returns a shapely polygon representing a Taper Flange U Section."""
    # Create first part of line
    ls = [
        set_precision(
            LineString([[0, h / 2], [0, 0], [b, 0]]), grid_size=1e-13
        )
    ]
    # Create first fillet
    xy = np.array([[b, 0], [b, tf - slope * u], [b - u, tf]])
    ls.append(
        set_precision(
            LineString(xy).offset_curve(r2).offset_curve(-r2), grid_size=1e-13
        )
    )
    # Create second fillet
    xy = np.array(
        [
            [b - u, tf],
            [tw, tf + slope * (b - u - tw)],
            [tw, h / 2],
        ]
    )
    ls.append(
        set_precision(
            LineString(xy).offset_curve(-r1).offset_curve(r1), grid_size=1e-13
        )
    )

    # Merge filleted
    merged_ls = linemerge(ls)

    # mirror the parts
    merged_ls = linemerge(
        [
            merged_ls,
            translate(scale(merged_ls, 1, -1), 0, h / 2),
        ]
    )

    # Create a polygon
    poly = polygonize([merged_ls])

    # Get the first and only polygon of this collection
    poly = get_geometry(poly, 0)
    # Translate it to the centroid when returning
    return orient(
        translate(poly, xoff=-poly.centroid.x, yoff=-poly.centroid.y), 1
    )


def _create_parallel_U_section(
    h: float,
    b: float,
    tw: float,
    tf: float,
    r: float,
) -> Polygon:
    """Returns a shapely polygon representing a Parallel Flange U Section."""
    # top flange
    top_flange = Polygon(
        [
            (0, h / 2 - tf),
            (b, h / 2 - tf),
            (b, h / 2),
            (0, h / 2),
        ]
    )
    # bottom flange
    bottom_flange = translate(top_flange, xoff=0, yoff=-h + tf)
    web = Polygon(
        [
            (0, -h / 2 + tf),
            (tw, -h / 2 + tf),
            (tw, h / 2 - tf),
            (0, h / 2 - tf),
        ]
    )
    # fillets
    p = Point([tw + r, -h / 2 + tf + r]).buffer(r)
    s = Polygon(
        [
            (tw, -h / 2 + tf),
            (tw + r, -h / 2 + tf),
            (tw + r, -h / 2 + tf + r),
            (tw, -h / 2 + tf + r),
        ]
    )
    fillet = s.difference(p)
    fillet = translate(
        scale(fillet, 1, -1), xoff=0, yoff=h - 2 * tf - r
    ).union(fillet)
    # Estimate grid_size value
    # Tentative geometry (due to approximations can be a MultiPolygon)
    geom_trial = unary_union([fillet, top_flange, bottom_flange, web])
    # minx, miny, maxx, maxy
    bounds = geom_trial.bounds
    min_size = min(bounds[2] - bounds[0], bounds[3] - bounds[1])
    grid_size = 10 ** int(math.floor(math.log10(abs(min_size)))) * 1e-12
    # Create the geometry
    geometries = [
        set_precision(geometry, grid_size=grid_size)
        for geometry in [fillet, top_flange, bottom_flange, web]
    ]
    geometry = orient(unary_union(geometries), 1)
    # Return the geometry centered at the origin
    return translate(
        geometry, xoff=-geometry.centroid.x, yoff=-geometry.centroid.y
    )


def _create_L_section(
    h: float, b: float, t1: float, t2: float, r1: float, r2: float
) -> Polygon:
    """Returns a shapely polygon representing a L Section."""
    # Create first part of line
    ls = [set_precision(LineString([[0, h], [0, 0], [b, 0]]), grid_size=1e-13)]
    # Create fillet
    xy = np.array([[b, 0], [b, t1], [b / 2, t1]])
    ls.append(
        set_precision(
            LineString(xy).offset_curve(r2).offset_curve(-r2), grid_size=1e-13
        )
    )
    # Create second fillet
    xy = np.array([[b / 2, t1], [t2, t1], [t2, h / 2]])
    ls.append(
        set_precision(
            LineString(xy).offset_curve(-r1).offset_curve(r1), grid_size=1e-13
        )
    )
    # Last fillet
    xy = np.array([[t2, h / 2], [t2, h], [0, h]])
    ls.append(
        set_precision(
            LineString(xy).offset_curve(r2).offset_curve(-r2), grid_size=1e-13
        )
    )
    # Merge filleted
    merged_ls = linemerge(ls)

    merged_ls

    # Create a polygon
    poly = polygonize([merged_ls])

    # Return the first and only polygon of this collection
    poly = get_geometry(poly, 0)
    # Translate it to the centroid when returning
    return orient(
        translate(poly, xoff=-poly.centroid.x, yoff=-poly.centroid.y), 1
    )
