"""Base class for profiles."""

import numpy as np
from shapely import (
    LinearRing,
    LineString,
    Polygon,
)
from shapely.affinity import rotate, translate
from shapely.ops import split

from structuralcodes.core._marin_integration import marin_integration


class BaseProfile:
    """Base class representing a profile.

    Contains the common code for all profiles.
    """

    def __init__(self):
        """Creates an empty base profile."""
        self._polygon: Polygon = None
        self._A: float = None
        self._Iy: float = None
        self._Iz: float = None
        self._Icsi: float = None
        self._Ieta: float = None
        self._Iyz: float = None
        self._Wely: float = None
        self._Welz: float = None
        self._iy: float = None
        self._iz: float = None
        self._Wply: float = None
        self._Wplz: float = None

    def _check_polygon_defined(self):
        """Just checks if polygon attribute is defined.

        If the polygon is not defined (it should never happen), an exception
        is Raised.
        """
        # The polygon attribute should be already defined
        if self._polygon is None:
            raise RuntimeError(
                'The polygon for some reason was not correctly defined.'
            )

    def _find_plastic_neutral_axis_y(self) -> float:
        """Find posizion z of plastic neutral axes parallel to y.

        We use bisection algorithm within the section limits.
        """
        bounds = self.polygon.bounds
        zmin, zmax = bounds[1], bounds[3]

        zA = zmin
        zB = zmax

        daA = self._find_delta_area_above_minus_below(z=zA)

        ITMAX = 200
        it = 0

        while (it < ITMAX) and (abs(zB - zA) > (zmax - zmin) * 1e-10):
            zC = (zA + zB) / 2.0
            daC = self._find_delta_area_above_minus_below(z=zC)
            if abs(daC) < 1e-10:
                break
            if daA * daC < 0:
                # The solution is between A and C
                zB = zC
            else:
                # The solution is between C and B
                zA = zC
                daA = daC
            it += 1
        if it >= ITMAX:
            s = f'Last iteration reached a unbalance of {daC}'
            raise ValueError(f'Maximum number of iterations reached.\n{s}')

        return zC

    def _find_delta_area_above_minus_below(self, z: float) -> float:
        """Returns area difference between above and below parts.

        Above and below parts are computed splitting the polygon with a line
        parallel to Y axes at z coordinate.
        """
        bounds = self._polygon.bounds
        xmax = max(abs(bounds[0]), bounds[2])
        line = LineString([[-xmax * 1.05, z], [xmax * 1.05, z]])

        area_above = 0
        area_below = 0
        # divide polygons "above" and "below" line
        if line.intersects(self._polygon):
            result = split(self._polygon, line)
            # divide polygons "above" and "below" line
            for geom in result.geoms:
                if LinearRing(
                    (line.coords[0], line.coords[1], geom.centroid.coords[0])
                ).is_ccw:
                    area_above += geom.area
                else:
                    area_below += geom.area
        else:
            # not intersecting, all the polygon is above or below the line
            geom = self.polygon
            if LinearRing(
                (line.coords[0], line.coords[1], geom.centroid.coords[0])
            ).is_ccw:
                area_above += geom.area
            else:
                area_below += geom.area
        return area_above - area_below

    def _find_principals_direction_and_moments(self):
        """Computes principal direction and second area moments."""
        eigres = np.linalg.eig(
            np.array([[self.Iy, self.Iyz], [self.Iyz, self.Iz]])
        )
        max_idx = np.argmax(eigres[0])
        min_idx = 0 if max_idx == 1 else 1
        self._Icsi = eigres[0][max_idx]
        self._Ieta = eigres[0][min_idx]
        self._theta = np.arccos(
            np.dot(np.array([1, 0]), eigres[1][:, max_idx])
        )

    @property
    def A(self) -> float:
        """Returns area of profile."""
        if self._A is None:
            # Check if the polygon is defined
            self._check_polygon_defined()
            # Get the polygon coordinates:
            xy = self._polygon.exterior.coords.xy
            # Compute area
            self._A = marin_integration(xy[0], xy[1], 0, 0)
        return self._A

    @property
    def Iy(self) -> float:
        """Returns second moment of area around y axis."""
        if self._Iy is None:
            # Check if the polygon is defined
            self._check_polygon_defined()
            # Get the polygon coordinates:
            xy = self._polygon.exterior.coords.xy
            # Compute second moments of inertia
            self._Iy = marin_integration(xy[0], xy[1], 0, 2)
        return self._Iy

    @property
    def Iz(self) -> float:
        """Returns second moment of area around z axis."""
        if self._Iz is None:
            # Check if the polygon is defined
            self._check_polygon_defined()
            # Get the polygon coordinates:
            xy = self._polygon.exterior.coords.xy
            # Compute second moments of inertia
            self._Iz = marin_integration(xy[0], xy[1], 2, 0)
        return self._Iz

    @property
    def Iyz(self) -> float:
        """Returns product moment of inertia."""
        if self._Iyz is None:
            # Check if the polygon is defined
            self._check_polygon_defined()
            # Get the polygon coordinates:
            xy = self._polygon.exterior.coords.xy
            # Compute product moment of area
            self._Iyz = marin_integration(xy[0], xy[1], 1, 1)
        return self._Iyz

    @property
    def Icsi(self) -> float:
        """Returns second moment of area around principal csi axis.

        It is assumed that Icsi is maximum second moment, while Ieta is the
        minimum one.
        """
        if self._Icsi is None:
            self._find_principals_direction_and_moments()
        return self._Icsi

    @property
    def Ieta(self) -> float:
        """Returns second moment of area around principal eta axis.

        It is assumed that Icsi is maximum second moment, while Ieta is the
        minimum one.
        """
        if self._Ieta is None:
            self._find_principals_direction_and_moments()
        return self._Ieta

    @property
    def theta(self) -> float:
        """Returns angle between x and principal eta axis.

        It is assumed that Icsi is maximum second moment, while Ieta is the
        minimum one.

        Returns:
            float: The angle in radians.
        """
        if self._theta is None:
            self._find_principals_direction_and_moments()
        return self._theta

    def _compute_elastic_moduli(self):
        """Compute elastic moduli Wely and Welz."""
        # Check if the polygon is defined
        self._check_polygon_defined()
        # For computing section modulus get bounds
        bounds = self._polygon.bounds
        xmax = max(abs(bounds[0]), bounds[2])
        ymax = max(abs(bounds[1]), bounds[3])
        # Then compute section modulus
        self._Wely = self.Iy / ymax
        self._Welz = self.Iz / xmax

    @property
    def Wely(self) -> float:
        """Returns section modulus in y direction."""
        if self._Wely is None:
            # Compute elastic moduli
            self._compute_elastic_moduli()
        return self._Wely

    @property
    def Welz(self) -> float:
        """Returns section modulus in z direction."""
        if self._Welz is None:
            # Compute elastic moduli
            self._compute_elastic_moduli()
        return self._Welz

    @property
    def Wply(self) -> float:
        """Returns plastic section modulus in y direction."""
        if self._Wply is None:
            # Check if the polygon is defined
            self._check_polygon_defined()
            # For computing section modulus get bounds
            bounds = self._polygon.bounds
            xmax = max(abs(bounds[0]), bounds[2])
            # Compute plastic section modulus
            # find plastic neutral axis parallel to y
            self._Wply = 0
            z_pna = self._find_plastic_neutral_axis_y()
            poly = translate(self._polygon, xoff=0, yoff=-z_pna)
            result = split(
                poly,
                LineString([[-xmax * 1.05, 0], [xmax * 1.05, 0]]),
            )
            for poly in result.geoms:
                xy = poly.exterior.coords.xy
                self._Wply += abs(marin_integration(xy[0], xy[1], 0, 1))
        return self._Wply

    @property
    def Wplz(self) -> float:
        """Returns plastic section modulus in z direction."""
        if self._Wplz is None:
            # Check if the polygon is defined
            self._check_polygon_defined()
            # For computing section modulus get bounds
            bounds = self._polygon.bounds
            ymax = max(abs(bounds[1]), bounds[3])
            # Compute plastic section modulus
            # # find plastic neutral axis parallel to z
            self._polygon = rotate(geom=self._polygon, angle=90, origin=(0, 0))
            self._Wplz = 0
            y_pna = self._find_plastic_neutral_axis_y()
            poly = translate(self._polygon, xoff=0, yoff=-y_pna)
            result = split(
                poly,
                LineString([[-ymax * 1.05, 0], [ymax * 1.05, 0]]),
            )
            for poly in result.geoms:
                xy = poly.exterior.coords.xy
                self._Wplz += abs(marin_integration(xy[0], xy[1], 0, 1))
            self._polygon = rotate(
                geom=self._polygon, angle=-90, origin=(0, 0)
            )
        return self._Wplz

    @property
    def iy(self) -> float:
        """Returns radius of inertia of profile."""
        # Compute radius of inertia
        self._iy = self._iy or (self.Iy / self.A) ** 0.5
        return self._iy

    @property
    def iz(self) -> float:
        """Returns radius of inertia of profile."""
        # Compute radius of inertia
        self._iz = self._iz or (self.Iz / self.A) ** 0.5
        return self._iz
