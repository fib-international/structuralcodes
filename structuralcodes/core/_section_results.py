"""Results."""

from __future__ import annotations  # To have clean hints of ArrayLike in docs

import math
import typing as t
from dataclasses import dataclass, field, fields

import numpy as np
from numpy.typing import ArrayLike, NDArray
from shapely import Point

from .base import Geometry, Section


@dataclass
class SectionProperties:
    """Simple dataclass for storing section properties."""

    # section areas
    area: float = field(default=0, metadata={'description': 'Total area'})
    area_reinforcement: float = field(
        default=0, metadata={'description': 'Reinforcement area'}
    )

    # axial rigidity
    ea: float = field(
        default=0, metadata={'description': 'Axial rigidity (E * A)'}
    )
    # ea = {'value': 0, 'description': 'Axial rigidity (E_A)'}

    # section mass
    mass: float = field(
        default=0, metadata={'description': 'Mass per unit length'}
    )

    # section perimeter
    perimeter: float = field(default=0, metadata={'description': 'Perimeter'})

    # first moments of area
    # global axes of section
    sy: float = field(
        default=0, metadata={'description': 'First moment of area (Sy)'}
    )
    sz: float = field(
        default=0, metadata={'description': 'Second moment of area (Sz)'}
    )

    # first moments of area * E
    # global axes of section
    e_sy: float = field(default=0, metadata={'description': 'E * Sy'})
    e_sz: float = field(default=0, metadata={'description': 'E * Sz'})

    # centroids
    cy: float = field(
        default=0, metadata={'description': 'Centroid y Coordinate'}
    )
    cz: float = field(
        default=0, metadata={'description': 'Centroid z Coordinate'}
    )

    # second moments of area
    # global axes of section
    iyy: float = field(
        default=0, metadata={'description': 'Second moment (Iyy)'}
    )
    izz: float = field(
        default=0, metadata={'description': 'Second moment (Izz)'}
    )
    iyz: float = field(
        default=0, metadata={'description': 'Product moment (Iyz)'}
    )
    # centroidal axes of section
    iyy_c: float = field(
        default=0, metadata={'description': 'Centroidal Second moment (Iyy_c)'}
    )
    izz_c: float = field(
        default=0, metadata={'description': 'Centroidal Second moment (Izz_c)'}
    )
    iyz_c: float = field(
        default=0,
        metadata={'description': 'Centroidal Product moment (Iyz_c)'},
    )
    # principal axes of section
    i11: float = field(
        default=0, metadata={'description': 'Principal Second Moment (I11)'}
    )
    i22: float = field(
        default=0, metadata={'description': 'Principal Second Moment (I22)'}
    )
    theta: float = field(
        default=0, metadata={'description': 'Principal axis angle (theta)'}
    )

    # section flexural rigidity
    # global axes of section
    e_iyy: float = field(default=0, metadata={'description': 'E * Iyy'})
    e_izz: float = field(default=0, metadata={'description': 'E * Izz'})
    e_iyz: float = field(default=0, metadata={'description': 'E * Iyz'})
    # centroidal axes of section
    e_iyy_c: float = field(default=0, metadata={'description': 'E * Iyy_c'})
    e_izz_c: float = field(default=0, metadata={'description': 'E * Izz_c'})
    e_iyz_c: float = field(default=0, metadata={'description': 'E * Iyz_c'})
    # principal axes of section
    e_i11: float = field(default=0, metadata={'description': 'E * I11'})
    e_i22: float = field(default=0, metadata={'description': 'E * I22'})
    e_theta: float = field(
        default=0, metadata={'description': 'Principal axis angle (theta)'}
    )

    def __format__(self, spec: str) -> str:
        """Defines the format for returning the string representation.

        Arguments:
            spec (str): The string specifying the format.
        """
        output_string = 'Section Properties:\n'
        for f in fields(self):
            value = getattr(self, f.name)
            description = f.metadata.get(
                'description', 'No description available'
            )
            output_string += f'{description}: {value:{spec}}\n'

        # etc. all other characteristics
        return output_string

    def __str__(self) -> str:
        """Returns the informal string representation of the section
        properties.
        """
        return f'{self}'

    def isclose(self, other, rtol=1e-5, atol=1e-8):
        """Check if two SectionProperties are close to each other.

        Arguments:
            other (SectionProperties): The other SectionProperties to compare.
            rtol (float): The relative tolerance.
            atol (float): The absolute tolerance.

        Returns:
            bool: True if the two SectionProperties are close, False otherwise.
        """
        if not isinstance(other, self.__class__):
            return NotImplemented

        a = np.array(list(vars(self).values()))
        b = np.array(list(vars(other).values()))
        return np.allclose(a, b, rtol=rtol, atol=atol)


# Utility functions for point response evaluation, used by both
# get_point_strain and get_point_stress.
# These are not meant to be public methods, but they are separated for clarity


def _matching_geometries(
    section: Section,
    name: t.Optional[str] = None,
    group_label: t.Optional[str] = None,
    case_sensitive: bool = True,
):
    """Return (surfaces, points) that match name/group filters."""
    geometries_name = section.geometry.name_filter(
        name, return_mode='split', case_sensitive=case_sensitive
    )
    geometries_group = section.geometry.group_filter(
        group_label, return_mode='split', case_sensitive=case_sensitive
    )

    surfaces = list(
        set(geometries_name['surfaces']) & set(geometries_group['surfaces'])
    )

    points = list(
        set(geometries_name['points']) & set(geometries_group['points'])
    )

    return surfaces, points


def _strain_from_kinematics(eps_a, chi_y, chi_z, y: float, z: float):
    """Get strain on points from section kinematics.
    Vectorized strain computation: works for scalars or arrays.
    """
    # strain = eps_a - chi_z * y + chi_y * z
    return eps_a - chi_z * y + chi_y * z


def _point_inside_geometry_surface(geo: Geometry, y: float, z: float) -> bool:
    """Check if a point is inside or on the boundary of a surface geometry."""
    p = Point(y, z)
    return geo.polygon.contains(p) or geo.polygon.touches(p)


def _point_matches_geometry_point(gp, y: float, z: float) -> bool:
    """Check if a point is close enough to the point geometry gp.

    Here we could implement in the future something that returns the
    strain at the point even if it is not exactly in the coordinates of
    the point geometries, but for now we will just check if it is close
    to any of the point geometries coordinates and return the strain if
    it is close enough.
    """
    return math.isclose(gp.x, y) and math.isclose(gp.y, z)


def _get_point_response(
    section: Section,
    eps_a: float,
    chi_y: float,
    chi_z: float,
    y: float,
    z: float,
    response_type: t.Literal['strain', 'stress'],
    name: t.Optional[str] = None,
    group_label: t.Optional[str] = None,
    case_sensitive: bool = True,
    all_results: bool = False,
):
    """Shared engine for get_point_strain / get_point_stress.
    Returns scalar/array (depending on eps_a/chi_y/chi_z) or dict if
    all_results=True.
    """
    if section is None:
        return None

    surfaces, points = _matching_geometries(
        section, name, group_label, case_sensitive
    )

    if all_results:
        response = {}

    # Surfaces: identify which surface(s) contain the point, then evaluate
    # response
    for geo in surfaces:
        if _point_inside_geometry_surface(geo, y, z):
            strain = _strain_from_kinematics(eps_a, chi_y, chi_z, y, z)
            value = (
                strain
                if response_type == 'strain'
                else geo.material.constitutive_law.get_stress(strain)
            )

            if all_results:
                response[(geo.name, geo.group_label)] = value
            else:
                return value

    # Points: check if we hit a point geometry exactly, then evaluate response
    for gp in points:
        if _point_matches_geometry_point(gp, y, z):
            strain = _strain_from_kinematics(eps_a, chi_y, chi_z, y, z)
            value = (
                strain
                if response_type == 'strain'
                else gp.material.constitutive_law.get_stress(strain)
            )

            if all_results:
                response[(gp.name, gp.group_label)] = value
            else:
                return value

    if all_results and len(response) > 0:
        return response
    return None


@dataclass
class MomentCurvatureResults:
    """Class for storing moment curvature results.

    The analysis will be done in general for a given inclination of n.a.
    """

    theta: float = 0  # the inclination of n.a.
    n: float = 0  # axial load - mantained constant during analysis
    chi_y: ArrayLike = None  # the curvatures
    chi_z: ArrayLike = None  # the curvatures
    eps_axial: ArrayLike = 0  # the axial strain (at section 0,0)
    m_y: ArrayLike = None  # the moment
    m_z: ArrayLike = None  # the moment

    section: Section = None

    detailed_result: SectionDetailedResultState = None
    seed: int = None
    current_step: int = None
    num_points: int = None

    def _create_detailed_result(self):
        self.detailed_result = SectionDetailedResultState(
            section=self.section,
            eps_a=self.eps_axial[self.current_step],
            chi_y=self.chi_y[self.current_step],
            chi_z=self.chi_z[self.current_step],
            n=self.n,
            m_y=self.m_y[self.current_step],
            m_z=self.m_z[self.current_step],
            num_points=self.num_points,
            seed=self.seed,
        )

    def create_detailed_result(self, num_points=1000):
        """Create the detailed result object.

        Arguments:
            num_points (int): Number of random points to sample for each
                surface geometry (default = 1000).
        """
        if self.seed is None:
            self.seed = np.random.randint(1, 100, 1)[0].item()

        self.current_step = 0
        self.num_points = num_points
        self._create_detailed_result()

    def next_step(self):
        """Advance to the next step in the detailed result."""
        if self.detailed_result is None:
            return
        if self.current_step < len(self.m_y) - 1:
            self.current_step += 1
            self._create_detailed_result()

    def previous_step(self):
        """Go back to the previous step in the detailed result."""
        if self.detailed_result is None:
            return
        if self.current_step > 0:
            self.current_step -= 1
            self._create_detailed_result()

    def set_step(self, step: int):
        """Set the detailed result to a specific step.

        Arguments:
            step (int): the step to set for the datailed_result object.
        """
        if self.detailed_result is None:
            return
        if 0 <= step < len(self.m_y):
            self.current_step = step
            self._create_detailed_result()

    def get_point_strain(
        self,
        y: float,
        z: float,
        name: t.Optional[str] = None,
        group_label: t.Optional[str] = None,
        case_sensitive: bool = True,
        all_results: bool = False,
    ) -> np.ndarray:
        """Return the strain at a given point (y,z).

        Arguments:
            y (float): The y-coordinate of the point.
            z (float): The z-coordinate of the point.
            name (str, optional): The name of the surface geometry to check.
            group_label (str, optional): The group label of the surface
                geometry to check.
            case_sensitive (bool, optional): If True (default) the matching is
                case sensitive.
            all_results (bool): If True, return the strain for all geometries
                that matches the filters, otherwise return the strain for the
                first geometry that matches the filters (default False).

        Returns:
            numpy.ndarray: The strain at the given point for all the steps, or
                None if the point is not within any of the geometries that
                match the filters.
        """
        return _get_point_response(
            section=self.section,
            eps_a=np.asarray(self.eps_axial),
            chi_y=np.asarray(self.chi_y),
            chi_z=np.asarray(self.chi_z),
            y=y,
            z=z,
            response_type='strain',
            name=name,
            group_label=group_label,
            case_sensitive=case_sensitive,
            all_results=all_results,
        )

    def get_point_stress(
        self,
        y: float,
        z: float,
        name: t.Optional[str] = None,
        group_label: t.Optional[str] = None,
        case_sensitive: bool = True,
        all_results: bool = False,
    ) -> np.ndarray:
        """Return the stress at a given point (y,z).

        Arguments:
            y (float): The y-coordinate of the point.
            z (float): The z-coordinate of the point.
            name (str, optional): The pattern for filtering the geometries by
                their name.
            group_label (str, optional): The pattern for filtering the
                geometries by their group_label.
            case_sensitive (bool, optional): If True (default) the matching is
                case sensitive.
            all_results (bool): If True, return the stress for all geometries
                that matches the filters, otherwise return the stress for the
                first geometry that matches the filters (default False).

        Returns:
            numpy.ndarray: The stress at the given point for all the steps, or
                None if the point is not within any of the geometries that
                match the filters.
        """
        return _get_point_response(
            section=self.section,
            eps_a=np.asarray(self.eps_axial),
            chi_y=np.asarray(self.chi_y),
            chi_z=np.asarray(self.chi_z),
            y=y,
            z=z,
            response_type='stress',
            name=name,
            group_label=group_label,
            case_sensitive=case_sensitive,
            all_results=all_results,
        )


class SectionDetailedResultState:
    """A class for storing section detailed results for specific state.

    This class stores in a specific data structure the results in terms of
    strain and stress fields for a specific state.

    The state is characterized by a strain plane (describing the section
    kinematics). This can correspond for instance to the plane strain
    computed when determining the ultimate strength of the section, or it can
    be a specific step during a moment-curvature analysis.
    """

    def __init__(
        self,
        section,
        eps_a,
        chi_y,
        chi_z,
        n,
        m_y,
        m_z,
        num_points=1000,
        seed=None,
    ):
        """Create the SectionDetailedResult.

        Arguments:
            section: The section object.
            eps_a: The axial strain at the section centroid.
            chi_y: The curvature about the y-axis.
            chi_z: The curvature about the z-axis.
            n: The axial force.
            m_y: The bending moment about the y-axis.
            m_z: The bending moment about the z-axis.
            num_points (int): Number of random points to sample in each surface
                geometry (default = 1000).
            seed (int): Seed for random number generator to ensure
                reproducibility of random points.
        """
        self._eps_a = eps_a
        self._chi_y = chi_y
        self._chi_z = chi_z
        self._n = n
        self._m_y = m_y
        self._m_z = m_z
        self._section = section
        self._seed = seed
        self._create_data_structure(num_points=num_points)

    def _create_data_structure(self, num_points=1000):
        """Create the data structure for storing the detailed results.

        This is designed to be easily used as input for DataFrames.
        """

        def _append_batch(cols: dict, batch: dict):
            """Helper function to append a batch of values to data."""
            # Determine batch length from first key
            first_key = next(iter(batch))
            n = len(batch[first_key])

            for k, v in batch.items():
                # Convert to list for extend
                if isinstance(v, np.ndarray):
                    cols[k].extend(v.tolist())
                elif isinstance(v, (list, tuple)):
                    cols[k].extend(v)
                else:
                    # this could be strings or materials
                    cols[k].extend([v] * n)

        # Data containers for surfaces and points
        surface_data = {
            'group_label': [],
            'name': [],
            'material': [],
            'y': [],
            'z': [],
            'strain': [],
            'stress': [],
        }
        point_data = {
            'group_label': [],
            'name': [],
            'material': [],
            'diameter': [],
            'area': [],
            'y': [],
            'z': [],
            'strain': [],
            'stress': [],
        }

        # SurfaceGeometry random points
        for surf in self.section.geometry.geometries:
            # Use surf.random_points_within to get random points
            y, z = surf.random_points_within(
                num_points=num_points, seed=self.seed
            )
            strain = self.eps_a - self.chi_z * y + self.chi_y * z
            stress = surf.material.constitutive_law.get_stress(strain)
            _append_batch(
                surface_data,
                {
                    'y': y,
                    'z': z,
                    'strain': strain,
                    'stress': stress,
                    'name': getattr(surf, 'name', None),
                    'group_label': getattr(surf, 'group_label', None),
                    'material': getattr(surf, 'material', None),
                },
            )

        # Get all point geometries
        point_geometries = self.section.geometry.group_filter(
            None, return_mode='split'
        )['points']

        if len(point_geometries) > 0:
            # We got some point geometries, so we will process them
            # For efficiency it is better to first group them by material,
            # so to use the vectorized get_stress of the constitutive law

            point_geometries_map = {}
            for pg in point_geometries:
                key = pg.material
                if key not in point_geometries_map:
                    point_geometries_map[key] = [pg]
                else:
                    point_geometries_map[key].append(pg)

            # Process each material group
            for material, pgs in point_geometries_map.items():
                y = np.array([pg.x for pg in pgs])
                z = np.array([pg.y for pg in pgs])
                strain = self.eps_a - self.chi_z * y + self.chi_y * z
                stress = material.constitutive_law.get_stress(strain)
                _append_batch(
                    point_data,
                    {
                        'y': y,
                        'z': z,
                        'strain': strain,
                        'stress': stress,
                        'name': [getattr(pg, 'name', None) for pg in pgs],
                        'group_label': [
                            getattr(pg, 'group_label', None) for pg in pgs
                        ],
                        'material': material,
                        'diameter': [
                            getattr(pg, 'diameter', None) for pg in pgs
                        ],
                        'area': [getattr(pg, 'area', None) for pg in pgs],
                    },
                )
        else:
            # If we don't have points return empty data structure
            point_data = None

        # Store results for later use (e.g., plotting)
        self._surface_data = surface_data
        self._point_data = point_data

    @property
    def section(self):
        """Return the section pointer."""
        return self._section

    @property
    def seed(self):
        """Return the seed for random sample of points."""
        return self._seed

    @property
    def n(self):
        """Return axial force."""
        return self._n

    @property
    def m_y(self):
        """Return bending moment m_y."""
        return self._m_y

    @property
    def m_z(self):
        """Return bending moment m_z."""
        return self._m_z

    @property
    def eps_a(self):
        """Return axial strain at (0, 0)."""
        return self._eps_a

    @property
    def chi_y(self):
        """Return curvature chi_y."""
        return self._chi_y

    @property
    def chi_z(self):
        """Return curvature chi_z."""
        return self._chi_z

    @property
    def strain(self):
        """Return the strain plane as a numpy array."""
        return np.array([self.eps_a, self.chi_y, self.chi_z])

    @property
    def surface_data(self):
        """Return the datastructure for surface geometries.

        The datastructure is a dictionary containing the following information
        for each point along the surfaces:

        - name: the name of the geometry
        - group_label: the group_label fo the geometry
        - material: the material of the geometry
        - y: the y coordinate of the point
        - z: the z coordinate of the point
        - strain: the strain of the point
        - stress: the stress of the point

        This dictionary can be easily given as input for creating a DataFrame.
        """
        return self._surface_data

    @property
    def point_data(self):
        """Return the datastructure for point geometries.

        The datastructure is a dictionary containing the following information
        for each point along the surfaces:

        - name: the name of the geometry
        - group_label: the group_label fo the geometry
        - material: the material of the geometry
        - diameter: the diameter of the geometry
        - area: the area of the geometry
        - y: the y coordinate of the point
        - z: the z coordinate of the point
        - strain: the strain of the point
        - stress: the stress of the point

        This dictionary can be easily given as input for creating a DataFrame.
        """
        return self._point_data


@dataclass
class UltimateBendingMomentResults:
    """Class for storing the ultimate bending moment computation for a given
    inclination of n.a. and axial load.
    """

    theta: float = 0  # the inclination of n.a.
    n: float = 0  # axial load - mantained constant during analysis
    m_y: float = 0  # the ultimate moment for given theta and n
    m_z: float = 0  # the ultimate moment for given theta and n
    chi_y: float = 0  # the curvature corresponding to the ultimate moment
    chi_z: float = 0  # the curvature corresponding to the ultimate moment
    eps_a: float = 0  # the axial strain at 0,0 corresponding to Mult

    section = None

    detailed_result: SectionDetailedResultState = None

    def create_detailed_result(self, num_points=1000):
        """Create the detailed result object.

        Arguments:
            num_points (int): Number of random points to sample for each
                surface geometry (default = 1000).
        """
        self.detailed_result = SectionDetailedResultState(
            section=self.section,
            eps_a=self.eps_a,
            chi_y=self.chi_y,
            chi_z=self.chi_z,
            n=self.n,
            m_y=self.m_y,
            m_z=self.m_z,
            num_points=num_points,
        )

    def get_point_strain(
        self,
        y: float,
        z: float,
        name: t.Optional[str] = None,
        group_label: t.Optional[str] = None,
        case_sensitive: bool = True,
        all_results: bool = False,
    ) -> float:
        """Return the strain at a given point (y,z).

        Arguments:
            y (float): The y-coordinate of the point.
            z (float): The z-coordinate of the point.
            name (str, optional): The name of the surface geometry to check.
            group_label (str, optional): The group label of the surface
                geometry to check.
            case_sensitive (bool, optional): If True (default) the matching is
                case sensitive.
            all_results (bool): If True, return the strain for all geometries
                that matches the filters, otherwise return the strain for the
                first geometry that matches the filters (default False).

        Returns:
            float: The strain at the given point, or None if the point is not
                within any of the geometries that match the filters.
        """
        return _get_point_response(
            section=self.section,
            eps_a=self.eps_a,
            chi_y=self.chi_y,
            chi_z=self.chi_z,
            y=y,
            z=z,
            response_type='strain',
            name=name,
            group_label=group_label,
            case_sensitive=case_sensitive,
            all_results=all_results,
        )

    def get_point_stress(
        self,
        y: float,
        z: float,
        name: t.Optional[str] = None,
        group_label: t.Optional[str] = None,
        case_sensitive: bool = True,
        all_results: bool = False,
    ) -> float:
        """Return the stress at a given point (y,z).

        Arguments:
            y (float): The y-coordinate of the point.
            z (float): The z-coordinate of the point.
            name (str, optional): The pattern for filtering the geometries by
                their name.
            group_label (str, optional): The pattern for filtering the
                geometries by their group_label.
            case_sensitive (bool, optional): If True (default) the matching is
                case sensitive.
            all_results (bool): If True, return the stress for all geometries
                that matches the filters, otherwise return the stress for the
                first geometry that matches the filters (default False).

        Returns:
            float: The strain at the given point, or None if the point is not
                within any of the geometries that match the filters.
        """
        return _get_point_response(
            section=self.section,
            eps_a=self.eps_a,
            chi_y=self.chi_y,
            chi_z=self.chi_z,
            y=y,
            z=z,
            response_type='stress',
            name=name,
            group_label=group_label,
            case_sensitive=case_sensitive,
            all_results=all_results,
        )


@dataclass(slots=True)
class StrainProfileResult:
    """Class for storing the results from calculate_strain_profile method."""

    # Solved generalized strains
    eps_a: float = 0.0  # the axial strain at 0, 0
    chi_y: float = 0.0  # the curvature respect y axes
    chi_z: float = 0.0  # the curvature respect z axes

    # target external loads assigned by user
    n_ext: float = 0.0  # Target axial load acting at 0, 0
    m_y_ext: float = 0.0  # Target moment My
    m_z_ext: float = 0.0  # Target moment Mz

    # Actual integrated loads at the converged strain state
    n: float = 0.0  # Axial load acting at 0, 0
    m_y: float = 0.0  # Bending moment My
    m_z: float = 0.0  # Bending moment Mz

    # Solver settings and diagnostics
    tolerance: float = 0.0
    max_iter: int = 0
    used_initial_tangent: bool = False
    iterations: int = 0
    converged: bool = False
    residual: NDArray[np.float64] = field(
        default_factory=lambda: np.zeros(3, dtype=float)
    )

    # Iteration history
    residual_history: list[NDArray[np.float64]] = field(default_factory=list)
    strain_history: list[NDArray[np.float64]] = field(default_factory=list)

    # For context store the section
    section: t.Any = None  # Note for future: if I want to type this I also have problem of circular import? #noqa E501
    # The detailed result data structure
    detailed_result: SectionDetailedResultState = None

    @property
    def residual_norm_history(self) -> t.List:
        """Returns the history of residual norm."""
        return [float(np.linalg.norm(x)) for x in self.residual_history]

    @property
    def delta_strain_history(self) -> t.List:
        """Returns as a list the history of delta_strain."""
        return [
            self.strain_history[i] - self.strain_history[i - 1]
            for i in range(1, len(self.strain_history))
        ]

    @property
    def delta_strain_norm_history(self) -> t.List:
        """Returns as a list the history of norm of delta_strain."""
        return [float(np.linalg.norm(x)) for x in self.delta_strain_history]

    @property
    def response_history(self) -> t.List:
        """Returns as a list the response (i.e. internal forces) history."""
        loads = np.array([self.n_ext, self.m_y_ext, self.m_z_ext])
        return [(loads - x) for x in self.residual_history]

    @property
    def strain_plane(self) -> NDArray[np.float64]:
        """Returns the strain profile as a numpy array."""
        return np.array([self.eps_a, self.chi_y, self.chi_z])

    @property
    def residual_norm(self) -> float:
        """Returns the norm of the residual at last iteration."""
        return float(np.linalg.norm(self.residual))

    def to_list(self) -> t.List:
        """Returns the strain profile coefficients in a list."""
        return [self.eps_a, self.chi_y, self.chi_z]

    def create_detailed_result(self, num_points=1000):
        """Create the detailed result object.

        Arguments:
            num_points (int): Number of random points to sample for each
                surface geometry (default = 1000).
        """
        self.detailed_result = SectionDetailedResultState(
            section=self.section,
            eps_a=self.eps_a,
            chi_y=self.chi_y,
            chi_z=self.chi_z,
            n=self.n,
            m_y=self.m_y,
            m_z=self.m_z,
            num_points=num_points,
        )

    def get_point_strain(
        self,
        y: float,
        z: float,
        name: t.Optional[str] = None,
        group_label: t.Optional[str] = None,
        case_sensitive: bool = True,
        all_results: bool = False,
    ) -> float:
        """Return the strain at a given point (y,z).

        Arguments:
            y (float): The y-coordinate of the point.
            z (float): The z-coordinate of the point.
            name (str, optional): The name of the surface geometry to check.
            group_label (str, optional): The group label of the surface
                geometry to check.
            case_sensitive (bool, optional): If True (default) the matching is
                case sensitive.
            all_results (bool): If True, return the strain for all geometries
                that matches the filters, otherwise return the strain for the
                first geometry that matches the filters (default False).

        Returns:
            float: The strain at the given point, or None if the point is not
                within any of the geometries that match the filters.
        """
        return _get_point_response(
            section=self.section,
            eps_a=self.eps_a,
            chi_y=self.chi_y,
            chi_z=self.chi_z,
            y=y,
            z=z,
            response_type='strain',
            name=name,
            group_label=group_label,
            case_sensitive=case_sensitive,
            all_results=all_results,
        )

    def get_point_stress(
        self,
        y: float,
        z: float,
        name: t.Optional[str] = None,
        group_label: t.Optional[str] = None,
        case_sensitive: bool = True,
        all_results: bool = False,
    ) -> float:
        """Return the stress at a given point (y,z).

        Arguments:
            y (float): The y-coordinate of the point.
            z (float): The z-coordinate of the point.
            name (str, optional): The pattern for filtering the geometries by
                their name.
            group_label (str, optional): The pattern for filtering the
                geometries by their group_label.
            case_sensitive (bool, optional): If True (default) the matching is
                case sensitive.
            all_results (bool): If True, return the stress for all geometries
                that matches the filters, otherwise return the stress for the
                first geometry that matches the filters (default False).

        Returns:
            float: The strain at the given point, or None if the point is not
                within any of the geometries that match the filters.
        """
        return _get_point_response(
            section=self.section,
            eps_a=self.eps_a,
            chi_y=self.chi_y,
            chi_z=self.chi_z,
            y=y,
            z=z,
            response_type='stress',
            name=name,
            group_label=group_label,
            case_sensitive=case_sensitive,
            all_results=all_results,
        )


@dataclass(slots=True)
class IntegrateStrainStiffnessResult:
    """Class for storing the results from integrating modulus."""

    # Input strain profile
    eps_a: float = 0.0  # the axial strain at 0, 0
    chi_y: float = 0.0  # the curvature respect y axes
    chi_z: float = 0.0  # the curvature respect z axes

    # Tangent
    tangent: NDArray[np.float64] = field(
        default_factory=lambda: np.zeros((3, 3), dtype=float)
    )

    def asarray(self, dtype=np.float64) -> NDArray:
        """Return an array representation of the tangent."""
        return np.asarray(self.tangent, dtype=dtype)


@dataclass(slots=True)
class IntegrateStrainForceResult:
    """Class for storing the results from integrating stresses."""

    # Input strain profile
    eps_a: float = 0.0  # the axial strain at 0, 0
    chi_y: float = 0.0  # the curvature respect y axes
    chi_z: float = 0.0  # the curvature respect z axes

    # Integrated loads
    n: float = 0.0  # Axial load acting at 0, 0
    m_y: float = 0.0  # Bending moment My
    m_z: float = 0.0  # Bending moment Mz

    # For context store the section
    section: t.Any = None  # Note for future: if I want to type this I also have problem of circular import? #noqa E501
    # The detailed result data structure
    detailed_result: SectionDetailedResultState = None

    def asarray(self, dtype=np.float64) -> NDArray:
        """Return an array representation of the forces."""
        return np.array([self.n, self.m_y, self.m_z], dtype=dtype)

    def astuple(self) -> tuple[float, float, float]:
        """Return a tuple representation of the forces."""
        return (self.n, self.m_y, self.m_z)

    def create_detailed_result(self, num_points=1000):
        """Create the detailed result object.

        Arguments:
            num_points (int): Number of random points to sample for each
                surface geometry (default = 1000).
        """
        self.detailed_result = SectionDetailedResultState(
            section=self.section,
            eps_a=self.eps_a,
            chi_y=self.chi_y,
            chi_z=self.chi_z,
            n=self.n,
            m_y=self.m_y,
            m_z=self.m_z,
            num_points=num_points,
        )

    def get_point_strain(
        self,
        y: float,
        z: float,
        name: t.Optional[str] = None,
        group_label: t.Optional[str] = None,
        case_sensitive: bool = True,
        all_results: bool = False,
    ) -> float:
        """Return the strain at a given point (y,z).

        Arguments:
            y (float): The y-coordinate of the point.
            z (float): The z-coordinate of the point.
            name (str, optional): The name of the surface geometry to check.
            group_label (str, optional): The group label of the surface
                geometry to check.
            case_sensitive (bool, optional): If True (default) the matching is
                case sensitive.
            all_results (bool): If True, return the strain for all geometries
                that matches the filters, otherwise return the strain for the
                first geometry that matches the filters (default False).

        Returns:
            float: The strain at the given point, or None if the point is not
                within any of the geometries that match the filters.
        """
        return _get_point_response(
            section=self.section,
            eps_a=self.eps_a,
            chi_y=self.chi_y,
            chi_z=self.chi_z,
            y=y,
            z=z,
            response_type='strain',
            name=name,
            group_label=group_label,
            case_sensitive=case_sensitive,
            all_results=all_results,
        )

    def get_point_stress(
        self,
        y: float,
        z: float,
        name: t.Optional[str] = None,
        group_label: t.Optional[str] = None,
        case_sensitive: bool = True,
        all_results: bool = False,
    ) -> float:
        """Return the stress at a given point (y,z).

        Arguments:
            y (float): The y-coordinate of the point.
            z (float): The z-coordinate of the point.
            name (str, optional): The pattern for filtering the geometries by
                their name.
            group_label (str, optional): The pattern for filtering the
                geometries by their group_label.
            case_sensitive (bool, optional): If True (default) the matching is
                case sensitive.
            all_results (bool): If True, return the stress for all geometries
                that matches the filters, otherwise return the stress for the
                first geometry that matches the filters (default False).

        Returns:
            float: The strain at the given point, or None if the point is not
                within any of the geometries that match the filters.
        """
        return _get_point_response(
            section=self.section,
            eps_a=self.eps_a,
            chi_y=self.chi_y,
            chi_z=self.chi_z,
            y=y,
            z=z,
            response_type='stress',
            name=name,
            group_label=group_label,
            case_sensitive=case_sensitive,
            all_results=all_results,
        )


@dataclass
class InteractionDomain:
    """Class for storing common data on all interaction domain results.

    Attributes:
        strains (numpy.Array): A numpy array with shape (n, 3) containing ea,
            ky and kz.
        forces (numpy.Array): A numpy array with shape (n, 3) containing n, my
            and mz.
        field_num (numpy.Array): a numpy array with shape (n,) containing a
            number between 1 and 6 indicating the failure field.
    """

    # array with shape (n,3) containing ea, ky, kz:
    strains: ArrayLike = None
    # array with shape(n,3) containing N, My, Mz
    forces: ArrayLike = None
    # array with shape(n,) containing the field number from 1 to 6
    field_num: ArrayLike = None

    @property
    def n(self):
        """Return axial force."""
        if self.forces is None:
            return None
        return self.forces[:, 0]

    @property
    def m_y(self):
        """Return my."""
        if self.forces is None:
            return None
        return self.forces[:, 1]

    @property
    def m_z(self):
        """Return mz."""
        if self.forces is None:
            return None
        return self.forces[:, 2]

    @property
    def e_a(self):
        """Return ea."""
        if self.strains is None:
            return None
        return self.strains[:, 0]

    @property
    def k_y(self):
        """Return ky."""
        if self.strains is None:
            return None
        return self.strains[:, 1]

    @property
    def k_z(self):
        """Return kz."""
        if self.strains is None:
            return None
        return self.strains[:, 2]


@dataclass
class NMMInteractionDomain(InteractionDomain):
    """Class for storing the NMM interaction domain results."""

    num_theta: int = 0  # number of discretizations along the angle
    num_axial: int = 0  # number of discretizations along axial load axis


@dataclass
class NMInteractionDomain(InteractionDomain):
    """Class for storing the NM interaction domain results."""

    theta: float = 0  # the inclination of n.a.
    num_axial: float = 0  # number of discretizations along axial load axis


@dataclass
class MMInteractionDomain(InteractionDomain):
    """Class for storing the MM interaction domain results."""

    num_theta: float = 0  # number of discretizations along the angle
    theta: ArrayLike = None  # Array with shape (n,) containing the angle of NA
