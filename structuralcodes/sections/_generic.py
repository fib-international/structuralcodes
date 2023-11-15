""""Generic class section implemenetation"""
from __future__ import annotations

import typing as t
import warnings
from numpy.typing import ArrayLike
import numpy as np
from shapely.geometry import (Polygon,
                              Point,
                              LineString,
                              LinearRing,
                              MultiPolygon,
                              MultiLineString)
from shapely.ops import split
from structuralcodes.core.base import Material
from structuralcodes.core.base import Section, SectionCalculator
import structuralcodes.sections._section_results as s_res

# Useful classes and functions: where to put?????? (core?
# utility folder in sections? here in this file?)
# Polygons, LineStrings, Points, MultiLyneStrings, MultiPolygons etc.

# to think: dataclass or class?
# Note that some things are already computed (like area) by shapely
# like: polygon.area, polygon.centroid, etc.

# For now dataclass, if we need convert to regular class,
# init commented for now


class Geometry:
    '''Base class for a geometry object'''

    section_counter: t.ClassVar[int] = 0

    def __init__(
        self,
        name: t.Optional[str] = None,
        group_label: t.Optional[str] = None
    ) -> None:
        """Initializes a geometry object
        The name and grouplabel serve for filtering
        in a compound object. By default it
        creates a new name each time
        
        Arguments:
            name (Optional: default None): the name to be given to the object
            group_label (Optional: default None): a label for grouping several
            objects"""
        if name is not None:
            self._name = name
        else:
            counter = Geometry.return_global_counter_and_increase()
            self._name = f"Geometry_{counter}"
        self._group_label = group_label

    @property
    def name(self):
        """Returns the name of the Geometry"""
        return self._name
    
    @property
    def group_label(self):
        """Returns the group_label fo the Geometry"""
        return self._group_label

    @classmethod
    def _increase_global_counter(cls):
        '''Increases the global counter by one'''
        cls.section_counter += 1

    @classmethod
    def return_global_counter_and_increase(cls):
        '''Returns the current counter and increases it by one'''
        counter = cls.section_counter
        cls._increase_global_counter()
        return counter


class PointGeometry(Geometry):
    '''Class for a point geometry with material.

    Basically it is a wrapper for shapely Point
    including the material (and other parameters that
    may be needed)'''
    def __init__(
        self,
        point: t.Union[Point, ArrayLike],
        diameter: float,
        material: Material,
        name: t.Optional[str] = None
    ):
        super().__init__(name)
        # I check if point is a shapely Point or an ArrayLike object
        if not isinstance(point, Point):
            # It is an ArrayLike object -> create the Point given the
            # coordinates x and y (coordinates can be a List, Tuple, np.array, ...)
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
        if not isinstance(material, Material):
            raise TypeError('The material should be of Material class')
        self._point = point
        self._diameter = diameter
        self._material = material
        self._area = np.pi * diameter**2 / 4.0

    @property
    def diameter(self) -> float:
        '''Returns the point diameter'''
        return self._diameter

    @property
    def area(self) -> float:
        '''Returns the point area'''
        return self._area

    @property
    def material(self) -> Material:
        '''Returns the point material'''
        return self._material

    @property
    def x(self) -> float:
        '''Returns the x coordinate of the point'''
        return self._point.x

    @property
    def y(self) -> float:
        '''Returns the y coordinate of the point'''
        return self._point.y

    @property
    def point(self) -> Point:
        '''Returns the shapely Point object'''
        return self._point


def create_line_point_angle(
    point: t.Union[Point, t.Tuple[float, float]],
    theta: float,
    bbox: t.Tuple[float, float, float, float]
) -> LineString:
    '''Creates a line from point and angle within the bounding
    box'''
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
        # atg is almost zero
        atg = v[0] / v[1]
        y1 = bbox[1] - 1e-3
        y2 = bbox[2] + 1e-3
        x1 = point[0] + (y1 - point[1]) * atg
        x2 = point[0] + (y2 - point[1]) * atg
    # create the splitting line
    return LineString([(x1, y1), (x2, y2)])


class SurfaceGeometry:
    '''Class for a surface geometry with material.

    Basically it is a wrapper for shapely polygon
    including the material (and other parameters needed)
    As a shapely polygon it can contain one or more holes'''

    def __init__(
        self,
        poly: Polygon,
        mat: Material
    ) -> None:
        """Initializes a SurfaceGeometry object

        Args:
            polly (shapely.Polygon): Shapely polygon
            mat (Material): A Material class applied to the geometry
        """
        # Check if inputs are of the correct type, otherwise return error
        if not isinstance(poly, Polygon):
            raise ValueError(
                f'poly need to be a valid shapely.geometry.Polygon object. \
                {repr(poly)}'
            )
        if not isinstance(mat, Material):
            raise ValueError(
                f'mat should be a valid structuralcodes.base.Material object. \
                {repr(mat)}'
            )
        self.polygon = poly
        self.material = mat

    @property
    def area(self) -> float:
        """Returns the area of the geometry

        Returns:
            float: area of the geometry"""
        return self.polygon.area

    @property
    def centroid(self) -> t.Tuple[float, float]:
        """Returns the centroid of the geometry

        Returns:
            (float, float): x and y coordinates of centroid"""
        return self.polygon.centroid.coords[0]

    def split(
        self,
        line: t.Union[LineString, t.Tuple[t.Tuple[float, float], float]]
    ) -> t.Tuple[t.List[SurfaceGeometry], t.List[SurfaceGeometry]]:
        """Splits the geometry using a line

        Args:
            line: can be a LineString shapely object, or a tuple
                  (point, theta) where point is a tuple (x,y) and
                  theta is the angle respect the horizontal axis

        Returns:
            (above_polygons, below_polygons)
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
                if LinearRing((line.coords[0],
                               line.coords[1],
                               geom.centroid.coords[0])).is_ccw:
                    above_polygons.append(geom)
                else:
                    below_polygons.append(geom)
        else:
            # not intersecting, all the polygon is above or below the line
            geom = self.polygon
            if LinearRing((line.coords[0],
                           line.coords[1],
                           geom.centroid.coords[0])).is_ccw:
                above_polygons.append(geom)
            else:
                below_polygons.append(geom)

        return above_polygons, below_polygons

    def split_two_lines(
        self,
        lines: t.Union[t.Tuple[LineString, LineString], MultiLineString]
    ) -> Polygon:
        '''Docstrings'''
        if isinstance(lines, MultiLineString):
            multi_line = lines
        elif isinstance(lines, tuple):
            if len(lines) != 2:
                raise RuntimeError('Two lines must be input')
            multi_line = MultiLineString(lines)
        lines_polygon = multi_line.convex_hull

        # get the intersection
        return self.polygon.intersection(lines_polygon)

    # here we can also add static methods like:
    # from_points
    # from_points_and_facets
    # from_surface_geometry
    # from_dxf
    # from_ascii
    # ...
    # we could also add methods wrapping shapely function, like:
    # mirror, translation, rotation, etc.


class CompoundGeometry(Geometry):
    '''Class for a compound geometry

    it is basicaly a set of geometries, each one with its
    own materials and properties.
    '''
    def __init__(
        self,
        geometries: t.List[Geometry] | MultiPolygon,
        materials: t.Optional[t.List[Material] | Material] = None
    ) -> None:
        '''Creates a compound geometry

        Args:
            geometries: a list of SurfaceGeometry objects or a shapely
                        MultiPolygon object (in this case also a list of
                        materials should be given)
            materials (optional, default = None): a material (applied
                        to all polygons) or a list of materials. In this
                        case the number of polygons should match the number
                        of materials
        '''

        if isinstance(geometries, MultiPolygon):
            # a MultiPolygon is provided
            if isinstance(materials, Material):
                # one material is provided. Applied to all polygons
                pass
            elif isinstance(materials, list):
                # the list of materials is provided, one for each polygon
                pass
            raise NotImplementedError
            # x = np.asarray(x, dtype=np.float64)
            # y = np.asarray(y, dtype=np.float64)
            # if x.shape != y.shape:
            #     raise ValueError("x and y must be array-like with the same shape")
        elif isinstance(geometries, list):
            # a list of SurfaceGeometry is provided
            checked_geometries = []
            for geo in geometries:
                if isinstance(geo, SurfaceGeometry):
                    checked_geometries.append(geo)
                elif isinstance(geo, CompoundGeometry):
                    for g in geo.geometries:
                        checked_geometries.append(g)
            self.geometries = checked_geometries

    # we can add here static methods like
    # from_dxf
    # from_ascii
    # ...

    @property
    def area(self) -> float:
        """Return the area of the compund geometry"""
        area = 0
        for geo in self.geometries:
            area += geo.area
        return area

    # Add split method that call the split for each geometry

    # Add methods for translation and rotation of compound and single geoms

   
def add_reinforcement(
    geo: SurfaceGeometry | CompoundGeometry,
    coords: t.Tuple(float, float),
    diameter: float,
    material: Material
) -> CompoundGeometry:
    """Adds a reinforcement bar to the geometry

    Proposals:
        i. Geometry class (group_label -> to be use for filtering)
        ii. PointGeometry(Geometry) that is a point with a mat
        iii. SurfaceGeometry(Geometry) that is a polygon with a mat
        iv. CompoundGeometry(Geometry) that is a set of geometries

    For compound we could implement simple + method so that geo1 + geo2
    returns a compunt.
    So this method could become simply:

    bar = Polygon(xxxx)
    return geo + bar
    """
    raise NotImplementedError


class GenericSection(Section):
    """This is the implementation of the generic class section"""

    def __init__(
        self,
        geometry: SurfaceGeometry | CompoundGeometry,
        name: t.Optional[str] = None
    ) -> None:
        if name is None:
            name = 'GenericSection'
        super().__init__(name)
        self.section_analyzer: GenericSectionCalculator(self)
        self.gross_properties = (
            self.section_analyzer._calculate_gross_section_properties()
        )


class GenericSectionCalculator(SectionCalculator):
    '''Calculator class implementing analysis algortims for
    Code checks'''

    def __init__(self, sec: GenericSection) -> None:
        '''Initialize the GenericSectionCalculator
        Input:
        section (SectionCalculator): the section object'''
        super().__init__(section=sec)

    def _calculate_gross_section_properties(self) -> s_res.GrossProperties:
        '''Calculates the gross section properties of the GenericSection
        This function is private and called when the section is created
        It stores the result into the result object

        Returns:
        gross_section_properties (GrossSection)'''

        # It will use the algorithms for generic sections
        gp = s_res.GrossProperties()
        # gp.area = ...
        # ...
        return gp

    def calculate_bending_strength(
        self, theta=0, n=0
    ) -> s_res.UltimateBendingMomentResult:
        """Calculates the bending strength for given inclination of n.a.
        and axial load

        Input:
        theta (float, default = 0): inclination of n.a. respect to y axis
        n (float, default = 0): axial load applied to the section 
        (+: tension, -: compression)

        Return:
        ultimate_bending_moment_result (UltimateBendingMomentResult)"""
        # For now it returns an empty response. The proper algorithms 
        # for generic section will be here
        return s_res.UltimateBendingMomentResult()

    def calculate_moment_curvature(
        self, theta=0, n=0
    ) -> s_res.MomentCurvatureResults:
        """Calculates the moment-curvature relation for given inclination of
        n.a. and axial load

        Input:
        theta (float, default = 0): inclination of n.a. respect to y axis
        n (float, default = 0): axial load applied to the section
            (+: tension, -: compression)
        chi_incr (float, default = 1e-8): the curvature increment for the
            analysis

        Return:
        moment_curvature_result (MomentCurvatureResults)
        """
        # For now it returns an empty response. The proper algorithms for
        # generic section will be here
        return s_res.MomentCurvatureResults()


# Use examples:

# sec = Section(...xxx...)
# res = sec.section_analyzer.ultmateBendingStrength(n = -1000)
# stress = sec.section_analyzer.compute_stress_distribution(res)
