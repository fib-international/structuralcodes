""""Generic class section implemenetation"""
from __future__ import annotations

import typing as t
from math import sin, cos
from shapely.geometry import Polygon, LineString, LinearRing, MultiPolygon
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
        point: t.Tuple[float, float],
        theta: float
    ) -> t.Tuple[t.List[SurfaceGeometry], t.List[SurfaceGeometry]]:
        """Splits the geometry using a line

        Args:
            point (float,float): a point in the splitting line
            theta (float): angle (in radians) respect to horiztonal axis

        Returns:
            (above_polygons, below_polygons)
        """

        # get boundingbox of polygon
        bbox = self.polygon.bounds

        # create a unit vector defining the line
        v = (cos(theta), sin(theta))

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
        line = LineString([(x1, y1), (x2, y2)])

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

    # here we can also add static methods like:
    # from_point
    # from_points_and_facets
    # from_surface_geometry
    # from_dxf
    # from_ascii
    # ...


class CompoundGeometry(SurfaceGeometry):
    '''Class for a compound geometry

    it is basicaly a set of geometries, each one with its
    own materials and properties.
    '''
    def __init__(
        self,
        geometries: t.List[SurfaceGeometry] | MultiPolygon,
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
    1. we can create a SurfaceGeometry where the bar is
    of circular form. So we just need a method that
    creates a polygon discretizing a bar
    2. Or we can rethink the structure and add:
        i. Geometry class
        ii. PointGeometry(Geometry) that is a point with a mat
        iii. SurfaceGeometry(Geometry) that is a polygon with a mat
        iv. CompountGeometry(Geometry) that is a set of geometries

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
