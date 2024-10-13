""" Specific class section implementation"""

from __future__ import annotations

import typing as t
from shapely import Polygon
from structuralcodes.sections._generic import GenericSection
from structuralcodes.geometry import SurfaceGeometry
from structuralcodes.materials.concrete import Concrete

class RectangularSection(GenericSection): 
    """This is the specific implementation of a rectangular class section.

    The section is a 2D geometry where Y axis is horizontal while Z axis is
    vertical.

    The moments and curvatures around Y and Z axes are assumed positive
    according to RHR.

    Attributes: TODO
        geometry (Union(SurfaceGeometry, CompoundGeometry)): The geometry of
            the section.
        name (str): The name of the section.
        section_calculator (GenericSectionCalculator): The object responsible
            for performing different calculations on the section (e.g. bending
            strength, moment curvature, etc.).
    """

    def __init__(self, height: float, width: float, concrete: Concrete, name: t.Optional[str] = None, integrator: t.Literal['marin', 'fiber'] = 'marin') -> None:
        """Initialize a GenericSection.

        Arguments:
            geometry (Union(SurfaceGeometry, CompoundGeometry)): The geometry
                of the section.
            name (str): The name of the section.
            integrator (str): The name of the SectionIntegrator to use.
        """
        poly: Polygon = Polygon(((0, 0), (width, 0), (width, height), (0, height)))
        geometry: SurfaceGeometry = SurfaceGeometry(poly, concrete)
        geometry.translate(-width/2, -height/2)
        if name is None:
            name: str = 'RectangularSection'
        super().__init__(geometry, name, integrator)

        self.height: float = height
        self.width: float = width

        # Calculating gross section properties
        # Gross area
        area: float = height * width
        perimeter: float = 2 * (height + width)
        ea: float = area * geometry.material.get_tangent(eps=0)[0]
        if geometry.density is not None:
            mass = area * geometry.density * 1e-9
        else:
            mass = 0

        # TODO: more gross sectional properties, and store them in self._gross_properties

    @property
    def area(self) -> float:
        """Returns Area of the cross section.
        Only the gross area of concrete is computed.

        Returns:
            float: The area of concrete section.
        """
        return self.width * self.height


class RectangularRC(Section):
    """Class for RC rectangular cross section.
    The class assumes width and height respectively
    along local x and y axes.

    Input parameters:
    width: width of the section
    height: height of the section
    cover: net cover (distance from outer edget of the section to
    edge of the stirrup)
    reinforcement_data: ReinforcementData
    concrete: Concrete
    reinforcement: Reinforcement (TODO!!!!!! for now const law)
    confined_concrete: Concrete (default = None) When set to None the
    same concrete is adopted for concrete core. When using the string
    "auto" the confinement will be automatically computed
    stirrups_data: StirrupsData (default = None) to use when need
    automatic computation of concrete confinement
    name: optional string for identidying the section

    Axes definition:

                     z ▲
                       │
                       │
            ┌──────────┼───────────┐
            │          │           │
            │          │           │
            │          │           │
            │          │           │
            │          │           │
            │          │           │
            │          │           │
            │          └───────────┼─────►
            │                      │     y
            │                      │
            │                      │
            │                      │
            │                      │
            │                      │
            │                      │
            └──────────────────────┘

            RectangularRCSection.fromGeometry(Geometry(itPolygons,itMaterials))


    """

    def __init__(
        self,
        width: float,
        height: float,
        cover: float,
        reinforcement_data: ReinforcementData,
        concrete: Concrete,
        reinforcement: ElasticPlastic,
        confined_concrete: t.Union[Concrete, str, None] = None,
        stirrups_data: t.Optional[StirrupsData] = None,
        name: t.Optional[str] = None,
    ) -> None:
        if name is None:
            name = 'RectangularRCSection'
        super().__init__(name)

        # section geometry
        self.height = height
        self.width = width
        self.cover = cover

        # materials
        self.concrete = concrete
        self.reinforcement = reinforcement
        if confined_concrete is not None:
            if isinstance(confined_concrete, str):
                # provided a string, check if it auto
                if (
                    confined_concrete.lower() == 'automatic'
                    or confined_concrete.lower() == 'auto'
                ):
                    # TODO: compute automatic confinement
                    if stirrups_data is None:
                        raise ValueError(
                            'Asked for automatic confinement \
                            computation, but stirrups_data is None \
                            \nPlease provide a proper stirrups_data'
                        )
                    self.confined_concrete = concrete
                else:
                    raise ValueError(f'Value "{confined_concrete}" unknown')
            elif isinstance(confined_concrete, Concrete):
                # provided a Concrete material, used that one
                self.confined_concrete = confined_concrete
        else:
            self.confined_concrete = concrete

        # stirrups
        if stirrups_data is None:
            # If no stirrups are provided by default it is assumed the
            # diameter is 8 mm - The spacing is not significant
            stirrups_data = StirrupsData(8, 150)
        self.stirrups_data = stirrups_data

        # longitudinal reinforcement
        self.reinforcement_data = reinforcement_data
        # arrange in the section the reinforcement
        self._arrange_reinforcement_in_section()
        # For now discretize in fibers here
        self._discretize_fibers(10, 25)

 

    def _discretize_fibers(self, nx: int, ny: int):
        """Function for creating a vector of all fibers of concrete."""

        def discretize_rectangle(x1, y1, x2, y2, nx, ny):
            if x2 < x1 or y2 < y1:
                raise ValueError('x2 and y2 cannot be lower than x1 and y1')
            w = x2 - x1
            h = y2 - y1
            dx = w / nx
            dy = h / ny
            x = []
            y = []
            Ac = []
            for i in range(ny):
                for j in range(nx):
                    x.append(x1 + dx * j + dx / 2)
                    y.append(y1 + dy * i + dy / 2)
                    Ac.append(dx * dy)
            return x, y, Ac

        # TODO confinement
        x, y, A = discretize_rectangle(
            -self.width / 2,
            -self.height / 2,
            self.width / 2,
            self.height / 2,
            nx,
            ny,
        )

        self.xc = np.array(x)
        self.yc = np.array(y)
        self.Ac = np.array(A)

    def max_compression_force(self):
        """Compute maximum compression force that can be
        taken by the section.
        """
        N = 0.0
        # Contribution from concrete
        # We need fcd actually: TODO
        N -= np.sum(self.Ac * 0.8 * self.concrete.fck)
        # Contribution from steel
        # Trick: I have not the material class for steel yet.
        N -= np.sum(self.As * self.reinforcement._fy)

        return N

    def max_tension_force(self):
        """Compute maximum tensile force that can be
        taken by the section.
        """
        return np.sum(self.As * self.reinforcement._fy)

    def bending_strength_xp(self, N: float = 0) -> float:
        """Returns the beding strength in x+ direction for a given
        value of axial force (+: tension, -: compression).
        """
        return 0.0

    def bending_strength_xn(self, N: float = 0) -> float:
        """Returns the beding strength in x- direction for a given
        value of axial force (+: tension, -: compression).
        """
        if self.max_tension_force() < N:
            raise ValueError('Too much tension applied to the section')
        if self.max_compression_force() > N:
            raise ValueError('Too much compression is applied to the section')

        # Find limit condition starting from contemporary failure
        # of steel and concrete and then perform a bisection method
        # for finding the NA position that equilibrates external N
        # Add this property to concrete material: TODO
        eps_cu = -0.0035
        p1 = (self.height / 2, -abs(eps_cu))
        p2 = (self.ymin, self.reinforcement._eps_su)
        strain_c = p2[1] - (p2[1] - p1[1]) / (p2[0] - p1[0]) * (
            p2[0] - self.yc
        )
        strain_s = p2[1] - (p2[1] - p1[1]) / (p2[0] - p1[0]) * (
            p2[0] - self.ys
        )
        chi = (p2[1] - p1[1]) / (p2[0] - p1[0])
        # Integrate internal axial force
        Nint = np.sum(
            self.concrete._stress_strain.get_stress(strain_c) * self.Ac
        ) + np.sum(self.reinforcement.get_stress(strain_s) * self.As)

        if Nint > N:
            # Too much tension, lowering NA
            # print('Too much compression, lowering NA')
            chi_a = chi
            dN_a = Nint - N

            chi_b = 0
            strain_c = p1[1] + chi_b * (self.yc - p1[0])
            strain_s = p1[1] + chi_b * (self.ys - p1[0])
            Nint = np.sum(
                self.concrete._stress_strain.get_stress(strain_c) * self.Ac
            ) + np.sum(self.reinforcement.get_stress(strain_s) * self.As)
            dN_b = Nint - N
            if dN_a * dN_b > 0:
                raise ValueError(
                    'Same sign on the interval, bisection will not work'
                )
            IT = 1
            ITMAX = 100
            while (abs(dN_b - dN_a) > 100) and (IT < ITMAX):
                chi_c = (chi_a + chi_b) / 2.0
                strain_c = p1[1] + chi_c * (self.yc - p1[0])
                strain_s = p1[1] + chi_c * (self.ys - p1[0])
                Nint = np.sum(
                    self.concrete._stress_strain.get_stress(strain_c) * self.Ac
                ) + np.sum(self.reinforcement.get_stress(strain_s) * self.As)
                dN_c = Nint - N
                if dN_c * dN_a < 0:
                    chi_b = chi_c
                    dN_b = dN_c
                else:
                    chi_a = chi_c
                    dN_a = dN_c
                IT += 1
        else:
            # Too much compression, rising NA
            # print('Too much compression, rising NA')
            chi_a = chi
            dN_a = Nint - N

            chi_b = 0
            strain_c = p2[1] + chi_b * (self.yc - p2[0])
            strain_s = p2[1] + chi_b * (self.ys - p2[0])
            Nint = np.sum(
                self.concrete._stress_strain.get_stress(strain_c) * self.Ac
            ) + np.sum(self.reinforcement.get_stress(strain_s) * self.As)
            dN_b = Nint - N
            print(f'chi_b = {chi_b}, dN_b = {dN_b}')
            if dN_a * dN_b > 0:
                raise ValueError(
                    'Same sign on the interval, bisection will not work'
                )
            IT = 1
            ITMAX = 100
            while (abs(dN_b - dN_a) > 100) and (IT < ITMAX):
                chi_c = (chi_a + chi_b) / 2.0
                strain_c = p2[1] + chi_c * (self.yc - p2[0])
                strain_s = p2[1] + chi_c * (self.ys - p2[0])
                Nint = np.sum(
                    self.concrete._stress_strain.get_stress(strain_c) * self.Ac
                ) + np.sum(self.reinforcement.get_stress(strain_s) * self.As)
                dN_c = Nint - N
                if dN_c * dN_a < 0:
                    chi_b = chi_c
                    dN_b = dN_c
                else:
                    chi_a = chi_c
                    dN_a = dN_c
                IT += 1
        # print(f'Found equilibrium after {IT} iterations')

        # For now I am saving these for plotting, will think if they are needed
        self.strain_c = strain_c
        self.strain_s = strain_s
        self.stress_c = self.concrete._stress_strain.get_stress(self.strain_c)
        self.stress_s = self.reinforcement.get_stress(self.strain_s)
        return np.sum(
            self.concrete._stress_strain.get_stress(strain_c)
            * self.yc
            * self.Ac
        ) + np.sum(self.reinforcement.get_stress(strain_s) * self.ys * self.As)

    def bending_strength_yp(self, N: float = 0) -> float:
        """Returns the beding strength in y+ direction for a given
        value of axial force (+: tension, -: compression).
        """
        return 0.0

    def bending_strength_yn(self, N: float = 0) -> float:
        """Returns the beding strength in y- direction for a given
        value of axial force (+: tension, -: compression).
        """
        return 0.0

    def _arrange_reinforcement_in_section(self):
        """The method creates the arrays for reinforcements
        starting from reinforcement_data.
        """
        delta = 2 * (self.cover + self.stirrups_data.diameter / 2.0)
        hc = self.height - delta
        wc = self.width - delta
        delta = (
            self.reinforcement_data.diameter_corner
            + self.stirrups_data.diameter
        )
        hcc = hc - delta
        wcc = wc - delta
        xmin = -wcc / 2.0
        xmax = wcc / 2.0
        ymin = -hcc / 2.0
        ymax = hcc / 2.0
        self.xs = [xmin, xmax, xmax, xmin]
        self.ys = [ymin, ymin, ymax, ymax]
        self.d = [self.reinforcement_data.diameter_corner] * 4
        # bottom bars
        for i in range(self.reinforcement_data.number_bottom):
            self.xs.append(
                xmin
                + (xmax - xmin)
                / (self.reinforcement_data.number_bottom + 1)
                * (i + 1)
            )
            self.ys.append(ymin)
            self.d.append(self.reinforcement_data.diameter_bottom)
        # top bars
        for i in range(self.reinforcement_data.number_top):
            self.xs.append(
                xmin
                + (xmax - xmin)
                / (self.reinforcement_data.number_top + 1)
                * (i + 1)
            )
            self.ys.append(ymax)
            self.d.append(self.reinforcement_data.diameter_top)
        # left bars
        for i in range(self.reinforcement_data.number_left):
            self.xs.append(xmin)
            self.ys.append(
                ymin
                + (ymax - ymin)
                / (self.reinforcement_data.number_left + 1)
                * (i + 1)
            )
            self.d.append(self.reinforcement_data.diameter_left)
        # right bars
        for i in range(self.reinforcement_data.number_right):
            self.xs.append(xmax)
            self.ys.append(
                ymin
                + (ymax - ymin)
                / (self.reinforcement_data.number_right + 1)
                * (i + 1)
            )
            self.d.append(self.reinforcement_data.diameter_right)

        self.xs = np.array(self.xs)
        self.ys = np.array(self.ys)
        self.d = np.array(self.d)
        # compute Areas for each reinforcement
        self.As = np.pi * (self.d**2) / 4

        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax