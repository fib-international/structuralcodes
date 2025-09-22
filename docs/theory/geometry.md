# Geometry

In Structuralcodes, a geometry object is a polygon or a point with a material assigend to it. Structuralcodes uses Shapely for creating polygons. Shapely works in screen **XY** coordinates, which are mapped to **yz** in the **GRS** system. Geometries are essential for defining structural sections and are categorized into the following classes:

## PointGeometry
Defines individual points in the plane. These typically represents reinforcing bars. A point geometry is represented by its coordinates combined with a material.

## SurfaceGeometry
Represents simple surface-based geometries such as rectangles, circles or generic polygons, also with holes. A surface geometry is represented as a shapely polygon combined with a material.

## CompoundGeometry
Combines multiple geometries into a single entity. Allows the creation of complex sections by grouping simpler geometries. A typical use is combining one or more surface geometries

:::{seealso}
For a description of the API for geometry creation and the different classes involved, refer to the [API reference](api-geometry-creation).
:::

**Figure Placeholder**: Example of geometries (SurfaceGeometry, PointGeometry, and CompoundGeometry) and their translations to the GRS coordinate system.