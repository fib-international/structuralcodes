"""Generic class section implemenetation."""

from __future__ import annotations  # To have clean hints of ArrayLike in docs

import typing as t
import warnings

import numpy as np
from numpy.typing import ArrayLike
from shapely import affinity
from shapely.geometry import (
    LinearRing,
    LineString,
    MultiLineString,
    MultiPolygon,
    Point,
    Polygon,
)
from shapely.ops import split

from structuralcodes.core.base import ConstitutiveLaw, Material
from structuralcodes.materials.concrete import Concrete
from structuralcodes.materials.constitutive_laws import Elastic

# Useful classes and functions: where to put?????? (core?
# utility folder in sections? here in this file?)
# Polygons, LineStrings, Points, MultiLyneStrings, MultiPolygons etc.

# to think: dataclass or class?
# Note that some things are already computed (like area) by shapely
# like: polygon.area, polygon.centroid, etc.

# For now dataclass, if we need convert to regular class,
# init commented for now


class Geometry:
    """Base class for a geometry object."""

    section_counter: t.ClassVar[int] = 0

    def __init__(
        self, name: t.Optional[str] = None, group_label: t.Optional[str] = None
    ) -> None:
        """Initializes a geometry object.

        The name and grouplabel serve for filtering in a compound object. By
        default it creates a new name each time.

        Arguments:
            name (Optional(str)): The name to be given to the object.
            group_label (Optional(str)): A label for grouping several objects.
        """
        if name is not None:
            self._name = name
        else:
            counter = Geometry.return_global_counter_and_increase()
            self._name = f'Geometry_{counter}'
        self._group_label = group_label

    @property
    def name(self):
        """Returns the name of the Geometry."""
        return self._name

    @property
    def group_label(self):
        """Returns the group_label fo the Geometry."""
        return self._group_label

    @classmethod
    def _increase_global_counter(cls):
        """Increases the global counter by one."""
        cls.section_counter += 1

    @classmethod
    def return_global_counter_and_increase(cls):
        """Returns the current counter and increases it by one."""
        counter = cls.section_counter
        cls._increase_global_counter()
        return counter

    @staticmethod
    def from_geometry(
        geo: Geometry,
        new_material: t.Optional[t.Union[Material, ConstitutiveLaw]] = None,
    ) -> Geometry:
        """Create a new geometry with a different material."""
        raise NotImplementedError(
            'This method should be implemented by subclasses'
        )


class PointGeometry(Geometry):
    """Class for a point geometry with material.

    Basically it is a wrapper for shapely Point including the material (and
    other parameters that may be needed).
    """

    def __init__(
        self,
        point: t.Union[Point, ArrayLike],
        diameter: float,
        material: t.Union[Material, ConstitutiveLaw],
        density: t.Optional[float] = None,
        name: t.Optional[str] = None,
        group_label: t.Optional[str] = None,
    ):
        """Initializes a PointGeometry object.

        The name and group_label serve for filtering in a compound object. By
        default it creates a new name each time.

        Arguments:
            point (Union(Point, ArrayLike)): A couple of coordinates or a
                shapely Point object.
            diameter (float): The diameter of the point.
            material (Union(Material, ConstitutiveLaw)): The material for the
                point (this can be a Material or a ConstitutiveLaw).
            density (Optional(float)): When a ConstitutiveLaw is passed as
                material, the density can be providen by this argument. When
                the material is a Material object the density is taken from the
                material.
            name (Optional(str)): The name to be given to the object.
            group_label (Optional(str)): A label for grouping several objects
                (default is None).
        """
        super().__init__(name, group_label)
        # I check if point is a shapely Point or an ArrayLike object
        if not isinstance(point, Point):
            # It is an ArrayLike object -> create the Point given the
            # coordinates x and y (coordinates can be a List, Tuple, np.array,
            # ...)
            coords = np.atleast_1d(point)
            num = len(coords)
            if num < 2:
                raise ValueError('Two coordinates are needed')
            if num > 2:
                warn_str = f'Two coordinates are needed. {num}'
                warn_str += ' coords provided. The extra entries will be'
                warn_str += ' discarded'
                warnings.warn(warn_str)
            point = Point(coords)
        if not isinstance(material, Material) and not isinstance(
            material, ConstitutiveLaw
        ):
            raise TypeError(
                f'mat should be a valid structuralcodes.base.Material \
                or structuralcodes.base.ConstitutiveLaw object. \
                {repr(material)}'
            )
        # Pass a constitutive law to the PointGeometry
        self._density = density
        if isinstance(material, Material):
            self._density = material.density
            material = material.constitutive_law

        self._point = point
        self._diameter = diameter
        self._material = material
        self._area = np.pi * diameter**2 / 4.0

    @property
    def diameter(self) -> float:
        """Returns the point diameter."""
        return self._diameter

    @property
    def area(self) -> float:
        """Returns the point area."""
        return self._area

    @property
    def material(self) -> Material:
        """Returns the point material."""
        return self._material

    @property
    def density(self) -> float:
        """Returns the density."""
        return self._density

    @property
    def x(self) -> float:
        """Returns the x coordinate of the point."""
        return self._point.x

    @property
    def y(self) -> float:
        """Returns the y coordinate of the point."""
        return self._point.y

    @property
    def point(self) -> Point:
        """Returns the shapely Point object."""
        return self._point

    def _repr_svg_(self) -> str:
        """Returns the svg representation."""
        return str(self._point._repr_svg_())

    def translate(self, dx: float = 0.0, dy: float = 0.0) -> PointGeometry:
        """Returns a new PointGeometry that is translated by dx, dy.

        Arguments:
            dx (float): Translation in x direction.
            dy (float): Translation in y direction.

        Returns:
            PointGeometry: A new, translated point.
        """
        return PointGeometry(
            point=affinity.translate(self._point, dx, dy),
            diameter=self._diameter,
            material=self._material,
            density=self._density,
            name=self._name,
            group_label=self._group_label,
        )

    def rotate(
        self,
        angle: float = 0.0,
        point: t.Union[t.Tuple[float, float], Point] = (0.0, 0.0),
        use_radians: bool = True,
    ) -> PointGeometry:
        """Returns a new PointGeometry that is rotated by angle.

        Arguments:
            angle (float): Rotation angle in radians (if use_radians = True),
                or degress (if use_radians = False).
            point (Union(Point, Tuple(float, float))): The origin of the
                rotation.
            use_radians (bool): True if angle is in radians, and False if angle
                is in degrees.
        """
        return PointGeometry(
            point=affinity.rotate(
                self._point, angle, origin=point, use_radians=use_radians
            ),
            diameter=self._diameter,
            material=self._material,
            density=self._density,
            name=self._name,
            group_label=self._group_label,
        )

    @staticmethod
    def from_geometry(
        geo: PointGeometry,
        new_material: t.Optional[t.Union[Material, ConstitutiveLaw]] = None,
    ) -> PointGeometry:
        """Create a new PointGeometry with a different material.

        Arguments:
            geo (PointGeometry): The geometry.
            new_material (Optional(Union(Material, ConstitutiveLaw))): A new
                material to be applied to the geometry. If new_material is
                None an Elastic material with same stiffness as the original
                material is created.

        Returns:
            PointGeometry: The new PointGeometry.

        Note:
            The polygon is not copied, but just referenced in the returned
            PointGeometry object.
        """
        if not isinstance(geo, PointGeometry):
            raise TypeError('geo should be a PointGeometry')
        if new_material is not None:
            # provided a new_material
            if not isinstance(new_material, Material) and not isinstance(
                new_material, ConstitutiveLaw
            ):
                raise TypeError(
                    f'new_material should be a valid structuralcodes.base.\
                    Material or structuralcodes.base.ConstitutiveLaw object. \
                    {repr(new_material)}'
                )
        else:
            # new_material not provided, assume elastic material with same
            # elastic modulus
            new_material = Elastic(E=geo.material.get_tangent(eps=0)[0])

        return PointGeometry(
            point=geo._point,
            diameter=geo._diameter,
            material=new_material,
            density=geo._density,
            name=geo._name,
            group_label=geo._group_label,
        )


def create_line_point_angle(
    point: t.Union[Point, t.Tuple[float, float]],
    theta: float,
    bbox: t.Tuple[float, float, float, float],
) -> LineString:
    """Creates a line from point and angle within the bounding box.

    Arguments:
        point (Union(Point, Tuple(float, float))): A Point or a coordinate the
            line should pass through.
        theta (float): The angle of the line in radians.
        bbox (Tuple(float, float, float, float)): Bounds for the created line.

    Returns:
        LineString: The created line.
    """
    # create a unit vector defining the line
    v = (np.cos(theta), np.sin(theta))

    # check if the line is vertical to avoid div by zero
    if abs(v[0]) > 1e-8:
        # it is a non vertical line
        tg = v[1] / v[0]
        x1 = bbox[0] - 1e-3
        x2 = bbox[2] + 1e-3
        y1 = point[1] + (x1 - point[0]) * tg
        y2 = point[1] + (x2 - point[0]) * tg
    else:
        # it is a near-vertical line
        # tg is almost zero
        ctg = v[0] / v[1]
        y1 = bbox[1] - 1e-3
        y2 = bbox[2] + 1e-3
        x1 = point[0] + (y1 - point[1]) * ctg
        x2 = point[0] + (y2 - point[1]) * ctg
    # create the line
    return LineString([(x1, y1), (x2, y2)])


class SurfaceGeometry:
    """Class for a surface geometry with material.

    Basically it is a wrapper for shapely polygon including the material (and
    other parameters needed). As a shapely polygon it can contain one or more
    holes.
    """

    def __init__(
        self,
        poly: Polygon,
        mat: t.Union[Material, ConstitutiveLaw],
        density: t.Optional[float] = None,
        concrete: bool = False,
    ) -> None:
        """Initializes a SurfaceGeometry object.

        Arguments:
            poly (shapely.Polygon): A Shapely polygon.
            mat (Union(Material, ConstitutiveLaw)): A Material or
                ConsitutiveLaw class applied to the geometry.
            density (Optional(float)): When a ConstitutiveLaw is passed as mat,
                the density can be provided by this argument. When mat is a
                Material object the density is taken from the material.
            concrete (bool): Flag to indicate if the geometry is concrete.
        """
        # Check if inputs are of the correct type, otherwise return error
        if not isinstance(poly, Polygon):
            raise TypeError(
                f'poly need to be a valid shapely.geometry.Polygon object. \
                {repr(poly)}'
            )
        if not isinstance(mat, Material) and not isinstance(
            mat, ConstitutiveLaw
        ):
            raise TypeError(
                f'mat should be a valid structuralcodes.base.Material \
                or structuralcodes.base.ConstitutiveLaw object. \
                {repr(mat)}'
            )
        self.polygon = poly
        # Pass a constitutive law to the SurfaceGeometry
        self._density = density
        if isinstance(mat, Material):
            self._density = mat.density
            if isinstance(mat, Concrete):
                concrete = True
            mat = mat.constitutive_law

        self.material = mat
        self.concrete = concrete

    @property
    def area(self) -> float:
        """Returns the area of the geometry.

        Returns:
            float: The area of the geometry.
        """
        return self.polygon.area

    @property
    def centroid(self) -> t.Tuple[float, float]:
        """Returns the centroid of the geometry.

        Returns:
            Tuple(float, float): x and y coordinates of the centroid.
        """
        return self.polygon.centroid.coords[0]

    @property
    def density(self) -> float:
        """Returns the density."""
        return self._density

    def calculate_extents(self) -> t.Tuple[float, float, float, float]:
        """Calculate extents of SurfaceGeometry.

        Calculates the minimum and maximum x and y values.

        Returns:
            Tuple(float, float, float, float): Minimum and maximum x and y
            values (``x_min``, ``x_max``, ``y_min``, ``y_max``).
        """
        min_x, min_y, max_x, max_y = self.polygon.bounds
        return min_x, max_x, min_y, max_y

    def split(
        self, line: t.Union[LineString, t.Tuple[t.Tuple[float, float], float]]
    ) -> t.Tuple[t.List[SurfaceGeometry], t.List[SurfaceGeometry]]:
        """Splits the geometry using a line.

        Arguments:
            line (Union(LineString, Tuple(Tuple(float, float), float))): A line
                either represented by a LineString shapely object, or a tuple
                (point, theta) where point is a coordinate pair and theta is
                the angle respect the horizontal axis in radians.

        Returns:
            Tuple(List(SurfaceGeometry), List(SurfaceGeometry)): The
            SurfaceGeometries above and below the line.
        """
        if not isinstance(line, LineString):
            point = line[0]
            theta = line[1]

            # get boundingbox of polygon
            bbox = self.polygon.bounds

            line = create_line_point_angle(point, theta, bbox)

        # split the geometry about the line
        above_polygons = []
        below_polygons = []
        if line.intersects(self.polygon):
            result = split(self.polygon, line)
            # divide polygons "above" and "below" line
            for geom in result.geoms:
                if LinearRing(
                    (line.coords[0], line.coords[1], geom.centroid.coords[0])
                ).is_ccw:
                    above_polygons.append(geom)
                else:
                    below_polygons.append(geom)
        else:
            # not intersecting, all the polygon is above or below the line
            geom = self.polygon
            if LinearRing(
                (line.coords[0], line.coords[1], geom.centroid.coords[0])
            ).is_ccw:
                above_polygons.append(geom)
            else:
                below_polygons.append(geom)

        return above_polygons, below_polygons

    def split_two_lines(
        self, lines: t.Union[t.Tuple[LineString, LineString], MultiLineString]
    ) -> t.Union[Polygon, MultiPolygon]:
        """Splits the geometry using two lines.

        Arguments:
            lines (Union(Tuple(LineString, Linestring)), MultiLineString): Two
                lines either represented by a tuple of two LineString shapely
                objects, or a MultiLineString shapely object.

        Returns:
            Union(Polygon, Multipolygon): The polygon or multipolygon obtained
            by splitting the SurfaceGeometry with the two lines.
        """
        if isinstance(lines, MultiLineString):
            multi_line = lines
        elif isinstance(lines, tuple):
            if len(lines) != 2:
                raise RuntimeError('Two lines must be input')
            multi_line = MultiLineString(lines)
        lines_polygon = multi_line.convex_hull

        # get the intersection
        return self.polygon.intersection(lines_polygon)

    def __add__(self, other: Geometry) -> CompoundGeometry:
        """Add operator "+" for geometries.

        Arguments:
            other (Geometry): The other geometry to add.

        Returns:
            CompoundGeometry: A new CompoundGeometry.
        """
        return CompoundGeometry([self, other])

    def __sub__(self, other: Geometry) -> SurfaceGeometry:
        """Add operator "-" for geometries.

        Arguments:
            other (Geometry): The other geometry to be subtracted.

        Returns:
            SurfaceGeometry: The resulting SurfaceGeometry.
        """
        material = self.material
        density = self._density

        # if we subtract a point from a surface we obtain the same surface
        sub_polygon = self.polygon
        if isinstance(other, SurfaceGeometry):
            # We are subtracting a surface from another surface
            sub_polygon = self.polygon - other.polygon
        if isinstance(other, CompoundGeometry):
            # We are subtracting a compound from a surface
            sub_polygon = self.polygon
            for g in other.geometries:
                sub_polygon = sub_polygon - g.polygon

        return SurfaceGeometry(poly=sub_polygon, mat=material, density=density)

    def _repr_svg_(self) -> str:
        """Returns the svg representation."""
        return str(self.polygon._repr_svg_())

    def translate(self, dx: float = 0.0, dy: float = 0.0) -> SurfaceGeometry:
        """Returns a new SurfaceGeometry that is translated by dx, dy.

        Arguments:
            dx (float): Translation in x direction.
            dy (float): Translation in y direction.

        Returns:
            SurfaceGeometry: The translated SurfaceGeometry.
        """
        return SurfaceGeometry(
            poly=affinity.translate(self.polygon, dx, dy),
            mat=self.material,
            density=self._density,
        )

    def rotate(
        self,
        angle: float = 0.0,
        point: t.Tuple[float, float] = (0.0, 0.0),
        use_radians: bool = True,
    ) -> SurfaceGeometry:
        """Returns a new SurfaceGeometry that is rotated by angle.

        Arguments:
            angle (float): Rotation angle in radians (if use_radians = True),
                or degress (if use_radians = False).
            point (Union(Point, Tuple(float, float))): The origin of the
                rotation.
            use_radians (bool): True if angle is in radians, and False if angle
                is in degrees.

        Returns:
            SurfaceGeometry: The rotated SurfaceGeometry.
        """
        return SurfaceGeometry(
            poly=affinity.rotate(
                self.polygon, angle, origin=point, use_radians=use_radians
            ),
            mat=self.material,
            density=self._density,
        )

    @staticmethod
    def from_geometry(
        geo: SurfaceGeometry,
        new_material: t.Optional[t.Union[Material, ConstitutiveLaw]] = None,
    ) -> SurfaceGeometry:
        """Create a new SurfaceGeometry with a different material.

        Arguments:
            geo (SurfaceGeometry): The geometry.
            new_material: (Optional(Union(Material, ConstitutiveLaw))): A new
                material to be applied to the geometry. If new_material is None
                an Elastic material with same stiffness of the original
                material is created.

        Returns:
            SurfaceGeometry: The new SurfaceGeometry.

        Note:
            The polygon is not copied, but just referenced in the returned
            SurfaceGeometry object.
        """
        if not isinstance(geo, SurfaceGeometry):
            raise TypeError('geo should be a SurfaceGeometry')
        if new_material is not None:
            # provided a new_material
            if not isinstance(new_material, Material) and not isinstance(
                new_material, ConstitutiveLaw
            ):
                raise TypeError(
                    f'new_material should be a valid structuralcodes.base.\
                    Material or structuralcodes.base.ConstitutiveLaw object. \
                    {repr(new_material)}'
                )
        else:
            # new_material not provided, assume elastic material with same
            # elastic modulus
            new_material = Elastic(E=geo.material.get_tangent(eps=0)[0])

        return SurfaceGeometry(
            poly=geo.polygon, mat=new_material, density=geo._density
        )

    # here we can also add static methods like:
    # from_points
    # from_points_and_facets
    # from_surface_geometry
    # from_dxf
    # from_ascii
    # ...
    # we could also add methods wrapping shapely function, like:
    # mirror, translation, rotation, etc.


def _process_geometries_multipolygon(
    geometries: MultiPolygon,
    materials: t.Optional[
        t.Union[t.List[Material], Material, ConstitutiveLaw]
    ],
) -> list:
    """Process geometries for initialization."""
    checked_geometries = []
    # a MultiPolygon is provided
    if isinstance(materials, (ConstitutiveLaw, Material)):
        for g in geometries.geoms:
            checked_geometries.append(SurfaceGeometry(poly=g, mat=materials))
    elif isinstance(materials, list):
        # the list of materials is provided, one for each polygon
        if len(geometries.geoms) != len(materials):
            raise ValueError(
                'geometries and materials should have the same length'
            )
        for g, m in zip(geometries.geoms, materials):
            checked_geometries.append(SurfaceGeometry(poly=g, mat=m))
    return checked_geometries


def _process_geometries_list(
    geometries: t.List[Geometry],
) -> t.Tuple[list, list]:
    """Process geometries for initialization."""
    # a list of SurfaceGeometry is provided
    checked_geometries = []
    checked_point_geometries = []
    for geo in geometries:
        if isinstance(geo, SurfaceGeometry):
            checked_geometries.append(geo)
        elif isinstance(geo, CompoundGeometry):
            for g in geo.geometries:
                checked_geometries.append(g)
            for pg in geo.point_geometries:
                checked_point_geometries.append(pg)
        elif isinstance(geo, PointGeometry):
            checked_point_geometries.append(geo)
    return (checked_geometries, checked_point_geometries)


class CompoundGeometry(Geometry):
    """Class for a compound geometry.

    It is basicaly a set of geometries, each one with its own materials and
    properties.
    """

    def __init__(
        self,
        geometries: t.Union[t.List[Geometry], MultiPolygon],
        materials: t.Optional[t.Union[t.List[Material], Material]] = None,
    ) -> None:
        """Creates a compound geometry.

        Arguments:
            geometries (Union(List(Geometry), MultiPolygon)): A list of
                SurfaceGeometry objects or a shapely MultiPolygon object
                (in this case also a list of materials should be given).
            materials (Optional(List(Material), Material)): A material (applied
                to all polygons) or a list of materials. In this case the
                number of polygons should match the number of materials.
        """
        if isinstance(geometries, MultiPolygon):
            # a MultiPolygon is provided
            self.geometries = _process_geometries_multipolygon(
                geometries, materials
            )
            self.point_geometries = []
            # useful for representation in svg
            geoms_representation = [g.polygon for g in self.geometries]
            self.geom = MultiPolygon(geoms_representation)
            return
        if isinstance(geometries, list):
            self.geometries, self.point_geometries = _process_geometries_list(
                geometries
            )
            # useful for representation in svg
            geoms_representation = [g.polygon for g in self.geometries]
            geoms_representation += [
                pg._point.buffer(pg._diameter / 2)
                for pg in self.point_geometries
            ]
            self.geom = MultiPolygon(geoms_representation)
        self._reinforced_concrete = None

    # we can add here static methods like
    # from_dxf
    # from_ascii
    # ...

    def _repr_svg_(self) -> str:
        """Returns the svg representation."""
        return str(self.geom._repr_svg_())

    @property
    def reinforced_concrete(self) -> bool:
        """Returns True if it is a Reinforced Concrete section."""
        if self._reinforced_concrete is None:
            self._reinforced_concrete = False
            for geo in self.geometries:
                if geo.concrete:
                    self._reinforced_concrete = True
                    break
        return self._reinforced_concrete

    @property
    def area(self) -> float:
        """Return the area of the compound geometry."""
        area = 0
        for geo in self.geometries:
            area += geo.area
        return area

    def calculate_extents(self) -> t.Tuple[float, float, float, float]:
        """Calculate extents of CompundGeometry.

        Calculates the minimum and maximum x and y-values.

        Returns:
            Tuple(float, float, float, float): Minimum and maximum x and y
            values (``x_min``, ``x_max``, ``y_min``, ``y_max``).

        Note:
            Considers only SurfaceGeometries and not PointGeometries!
        """
        min_x = 1e16
        max_x = -1e16
        min_y = 1e16
        max_y = -1e16
        for geo in self.geometries:
            xmin, xmax, ymin, ymax = geo.calculate_extents()
            min_x = min(min_x, xmin)
            min_y = min(min_y, ymin)
            max_x = max(max_x, xmax)
            max_y = max(max_y, ymax)
        return min_x, max_x, min_y, max_y

    def translate(self, dx: float = 0.0, dy: float = 0.0) -> CompoundGeometry:
        """Returns a new CompountGeometry that is translated by dx, dy.

        Arguments:
            dx (float): Translation in x direction.
            dy (float): Translation in y direction.

        Returns:
            CompoundGeometry: The translated CompoundGeometry.
        """
        processed_geoms = []
        for g in self.geometries:
            processed_geoms.append(g.translate(dx, dy))
        for pg in self.point_geometries:
            processed_geoms.append(pg.translate(dx, dy))
        return CompoundGeometry(geometries=processed_geoms)

    def rotate(
        self,
        angle: float = 0.0,
        point: t.Tuple[float, float] = (0.0, 0.0),
        use_radians: bool = True,
    ) -> CompoundGeometry:
        """Returns a new CompoundGeometry that is rotated by angle.

        Arguments:
            angle (float): Rotation angle in radians (if use_radians = True),
                or degress (if use_radians = False).
            point (Union(Point, Tuple(float, float))): The origin of the
                rotation.
            use_radians (bool): True if angle is in radians, and False if angle
                is in degrees.

        Returns:
            CompoundGeometry: The rotated CompoundGeometry.
        """
        processed_geoms = []
        for g in self.geometries:
            processed_geoms.append(g.rotate(angle, point, use_radians))
        for pg in self.point_geometries:
            processed_geoms.append(pg.rotate(angle, point, use_radians))
        return CompoundGeometry(geometries=processed_geoms)

    def __add__(self, other: Geometry) -> CompoundGeometry:
        """Add operator "+" for geometries.

        Arguments:
            other (Geometry): The other geometry to add.

        Returns:
            CompoundGeometry: A new CompoundGeometry.
        """
        return CompoundGeometry([self, other])

    def __sub__(self, other: Geometry) -> CompoundGeometry:
        """Add operator "-" for geometries.

        Arguments:
            other (Geometry): The other geometry to be subtracted.

        Returns:
            CompoundGeometry: A new CompoundGeometry.
        """
        # if we subtract a point from a surface we obtain the same surface
        if isinstance(other, PointGeometry):
            return self
        # Otherwise perform subtraction
        processed_geoms = []
        for g in self.geometries:
            processed_geoms.append(g - other)
        for pg in self.point_geometries:
            processed_geoms.append(pg)
        return CompoundGeometry(geometries=processed_geoms)

    @staticmethod
    def from_geometry(
        geo: CompoundGeometry,
        new_material: t.Optional[t.Union[Material, ConstitutiveLaw]] = None,
    ) -> CompoundGeometry:
        """Create a new CompoundGeometry with a different material.

        Arguments:
            geo (CompoundGeometry): The geometry.
            new_material (Optional(Union(Material, ConstitutiveLaw))): A new
                material to be applied to the geometry. If new_material is None
                an Elastic material with same stiffness of the original
                material is created.

        Returns:
            CompoundGeometry: The new CompoundGeometry.
        """
        if not isinstance(geo, CompoundGeometry):
            raise TypeError('geo should be a CompoundGeometry')
        processed_geoms = []
        for g in geo.geometries:
            processed_geoms.append(
                SurfaceGeometry.from_geometry(geo=g, new_material=new_material)
            )
        for pg in geo.point_geometries:
            processed_geoms.append(
                PointGeometry.from_geometry(geo=pg, new_material=new_material)
            )
        return CompoundGeometry(geometries=processed_geoms)
