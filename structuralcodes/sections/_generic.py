""""Generic class section implemenetation"""
import typing as t
from structuralcodes.core.base import Material
from structuralcodes.core.base import Section, SectionCalculator
import structuralcodes.sections._section_results as s_res

# Generic cross section implementation
# Characteristics:
# outer polygon (define a polygon class as a list of points)
# inner polygons (for holes)
# list of reinforcement positions, area and material (point reinforcement)
# list of reinforcement poistions, thickness and material (line reinforcement)?

# Useful classes and functions: where to put?????? (core? utility folder in sections? here in this file?)
from shapely import Polygon, Point
from dataclasses import dataclass

# Polygons, LineStrings, Points, MultiLyneStrings, MultiPolygons etc.

# to think: dataclass or class? 
# Note that some things are already computed (like area) by shapely
# ike: polygon.area, polygon.centroid, etc.

# For now dataclass, if we need convert to regular class, init commented for now
@dataclass
class SurfaceGeometry:
    '''Class for geometry. 
    Basically it is a wrapper for shapely polygon
    including the material (or the constitutive law directly?)'''

    polygon: Polygon
    material: Material
    # def __init__(
    #     self, poly: Polygon, mat: Material
    # ) -> None:
    #     if not isinstance(poly, Polygon):
    #         raise ValueError(
    #             f'poly need to be a valid shapely.geometry.Polygon object. {repr(poly)}'
    #         )
    #     if not isinstance(mat, Material):
    #         raise ValueError(
    #             f'mat should be a valid structuralcodes.base.Material object. {repr(mat)}'
    #         )
    #     self.polygon = poly
    #     self.material = mat

    # here we can add static methods like:
    # from_point
    # from_dxf
    # from_ascii
    # ...

# class CompoundGeometry:
# Here there will be a class for having sets of SurfaceGeometries, Points, Lines, etc.

# Add methods to add_bars or lines -> takes as input a geometry and returns a CompoundGeometry

# Geometry or CompoundGeometry are stored into the generic section

class GenericSection(Section):
    """This is the implementation of the generic class section"""

    def __init__(self, name: t.Optional[str] = None) -> None:
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
        n (float, default = 0): axial load applied to the section (+: tension, -: compression)

        Return:
        ultimate_bending_moment_result (UltimateBendingMomentResult)"""
        # For now it returns an empty response. The proper algorithms for generic section
        # will be here
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
