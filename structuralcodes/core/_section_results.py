"""Results."""

from __future__ import annotations  # To have clean hints of ArrayLike in docs

from dataclasses import dataclass, field, fields

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


@dataclass
class MMInteractionDomain:
    """Class for storing the MM interaction domain results."""

    num_theta: float = 0  # number of discretizations along the angle
    n: float = 0  # axial load

    theta: ArrayLike = None  # Angle theta respect axis Y
    m_y: ArrayLike = None  # Moments My
    m_z: ArrayLike = None  # Moments Mz
