"""Results."""

from __future__ import annotations  # To have clean hints of ArrayLike in docs

import typing as t
from dataclasses import dataclass, field, fields

import matplotlib.pyplot as plt
from matplotlib.axis import Axis
from matplotlib.figure import Figure
from numpy.typing import ArrayLike


@dataclass
class GrossProperties:
    """Simple dataclass for storing gross section properties."""

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
        output_string = 'Gross Concrete Section Properties:\n'
        for f in fields(self):
            value = getattr(self, f.name)
            description = f.metadata.get(
                'description', 'No description available'
            )
            output_string += f'{description}: {value:{spec}}\n'

        # etc. all other characteristics
        return output_string

    def __str__(self) -> str:
        """Returns the informal string representation of the gross concrete
        section properties.
        """
        return f'{self}'


@dataclass
class CrackedProperties:
    """Simple dataclass for storing cracked section properties."""

    # second moments of area
    i_yy: float = 0
    i_zz: float = 0
    i_yz: float = 0

    # section cracked flexural rigidity
    ei_yy: float = 0
    ei_zz: float = 0


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

    # The strains can be reconstructed at each step from chi and eps_axial
    # The stresses can be recomputed if needed on the fly? Or storing them?


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


@dataclass
class NMMInteractionDomain:
    """Class for storing the NMM interaction domain results."""

    num_theta: int = 0  # number of discretizations along the angle
    num_axial: int = 0  # number of discretizations along axial load axis

    strains: ArrayLike = None  # array with shape (n,3) containing strains
    forces: ArrayLike = None  # array with shape(n,3) containing N, My, Mz


@dataclass
class NMInteractionDomain:
    """Class for storing the NM interaction domain results."""

    theta: float = 0  # the inclination of n.a.
    num_axial: float = 0  # number of discretizations along axial load axis

    n: ArrayLike = None  # Axial loads
    m_y: ArrayLike = None  # Moments My
    m_z: ArrayLike = None  # Moments Mz

    strains: ArrayLike = None

    def get_2d_plot(
        self,
        horizontal_axis: t.Literal['n', 'my', 'mz'] = 'n',
        vertical_axis: t.Literal['n', 'my', 'mz'] = 'my',
        ax: t.Optional[Axis] = None,
    ) -> t.Tuple[Figure, Axis]:
        """Retuns figure and axes handlers for a matplotlib plot."""
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.get_figure()
        if horizontal_axis == 'n':
            x = self.n
        elif horizontal_axis == 'my':
            x = self.m_y
        elif horizontal_axis == 'mz':
            x = self.m_z
        if vertical_axis == 'n':
            y = self.n
        elif vertical_axis == 'my':
            y = self.m_y
        elif vertical_axis == 'mz':
            y = self.m_z
        ax.plot(x, y)
        ax.set(xlabel=horizontal_axis, ylabel=vertical_axis)
        return fig, ax

    def get_3d_plot(
        self,
        x_axis: t.Literal['n', 'my', 'mz'] = 'n',
        y_axis: t.Literal['n', 'my', 'mz'] = 'my',
        z_axis: t.Literal['n', 'my', 'mz'] = 'mz',
        ax: t.Optional[Axis] = None,
    ) -> t.Tuple[Figure, Axis]:
        """Retuns figure and axes handlers for a matplotlib 3d plot."""
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        else:
            fig = ax.get_figure()
        if x_axis == 'n':
            x = self.n
        elif x_axis == 'my':
            x = self.m_y
        elif x_axis == 'mz':
            x = self.m_z
        if y_axis == 'n':
            y = self.n
        elif y_axis == 'my':
            y = self.m_y
        elif y_axis == 'mz':
            y = self.m_z
        if z_axis == 'n':
            z = self.n
        elif z_axis == 'my':
            z = self.m_y
        elif z_axis == 'mz':
            z = self.m_z
        ax.plot(x, y, z)
        ax.set(xlabel=x_axis, ylabel=y_axis, zlabel=z_axis)
        return fig, ax


@dataclass
class MMInteractionDomain:
    """Class for storing the MM interaction domain results."""

    num_theta: float = 0  # number of discretizations along the angle
    n: float = 0  # axial load

    theta: ArrayLike = None  # Angle theta respect axis Y
    m_y: ArrayLike = None  # Moments My
    m_z: ArrayLike = None  # Moments Mz
