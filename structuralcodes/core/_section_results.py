"""Results."""

# import typing as t
from dataclasses import dataclass

from numpy.typing import ArrayLike


@dataclass
class GrossProperties:
    """Simple dataclass for storing gross section properties."""

    # section areas
    area: float = 0
    area_concrete: float = 0
    area_reinforcement: float = 0

    # axial rigidity
    ea: float = 0

    # section mass
    mass: float = 0

    # section perimeter
    perimeter: float = 0

    # first moments of area
    sy: float = 0
    sz: float = 0

    # centroids
    cy: float = 0
    cz: float = 0

    # second moments of area
    i_yy: float = 0  # strong axis
    i_zz: float = 0  # weak axis
    i_yz: float = 0
    # TODO: principal axes i_11 i_22 and angle

    # section flexural rigidity
    ei_yy: float = 0
    ei_zz: float = 0

    def __format__(self, spec: str) -> str:
        """Defines the format for returning the string representation.

        Arguments:
        spec: the string specifying the format
        """
        output_string = 'Gross Concrete Section Properties:\n'
        output_string += f'Total area: {self.area:{spec}}'
        # etc. all other characteristics
        return output_string

    def __str__(self) -> str:
        """Returns the informal string representation.
        Returns the informal string representation of the gross concrete
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
    """class for storing moment curvature results
    the analysis will be done in general for a given inclination
    of n.a.
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
    """class for storing the ultimate bending moment computation
    for a given inclination of n.a. and axial load.
    """

    theta: float = 0  # the inclination of n.a.
    n: float = 0  # axial load - mantained constant during analysis
    m_y: float = 0  # the ultimate moment for given theta and n
    m_z: float = 0  # the ultimate moment for given theta and n
    chi_y: float = 0  # the curvature corresponding to the ultimate moment
    chi_z: float = 0  # the curvature corresponding to the ultimate moment
    eps_a: float = 0  # the axial strain at 0,0 corresponding to Mult


@dataclass
class NmmInteractionDomain:
    """Class for storing the NMM interaction domain results."""

    n_theta: int = 0  # number of discretizations along the angle
    n_axial: int = 0  # number of discretizations along axial load axis


@dataclass
class NmInteractionDomain:
    """Class for storing the NM interaction domain results."""

    theta: float = 0  # the inclination of n.a.
    n_axial: float = 0  # number of discretizations along axial load axis

    n: ArrayLike = None  # Axial loads
    m_y: ArrayLike = None  # Moments My
    m_z: ArrayLike = None  # Moments Mz


@dataclass
class MmInteractionDomain:
    """Class for storing the MM interaction domain results."""

    n_theta: float = 0  # number of discretizations along the angle
    n: float = 0  # axial load
