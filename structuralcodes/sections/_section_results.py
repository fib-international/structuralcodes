''' Results '''
# import typing as t
from dataclasses import dataclass


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
    sx: float = 0
    sy: float = 0

    # centroids
    cx: float = 0
    cy: float = 0

    # second moments of area
    i_yy: float = 0  # strong axis
    i_zz: float = 0  # weak axis
    i_yz: float = 0
    # TODO: principal axes i_11 i_22 and angle

    # section flexural rigidity
    ei_xx: float = 0
    ei_yy: float = 0

    def __repr__(self, fmt: str = ".3f") -> str:
        """ Prints the gross concrete section properties.
        TODO: optionally give a file for printing into?

        Arguments:
         fmt: Number format"""

        output_string = 'Gross Concrete Section Properties:\n'
        output_string += f'Total area: {self.area:{fmt}}'
        # etc...
        return output_string


@dataclass
class CrackedProperties:
    """Simple dataclass for storing cracked section properties."""

    # second moments of area
    i_xx: float = 0
    i_yy: float = 0
    i_xy: float = 0

    # section cracked flexural rigidity
    ei_xx: float = 0
    ei_yy: float = 0


@dataclass
class MomentCurvatureResults:
    '''class for storing moment curvature results
    the analysis will be done in general for a given inclination
    of n.a.'''

    theta: float = 0  # the inclination of n.a.
    n: float = 0  # axial load - mantained constant during analysis
    chi: float = 0  # the curvatures
    eps_axial: float = 0  # the axial strain (at csection entroid)
    m: float = 0  # the moment

    # The strains can be reconstructed at each step from chi and eps_axial
    # The stresses can be recomputed if needed on the fly? Or storing them?


@dataclass
class UltimateBendingMomentResult:
    '''class for storing the ultimate bending moment computation
    for a given inclination of n.a. and axial load'''

    theta: float = 0  # the inclination of n.a.
    n: float = 0  # axial load - mantained constant during analysis
    m: float = 0  # the ultimate moment for given theta and n
    chi: float = 0  # the curvature corresponding to the ultimate moment
    eps_a: float = 0  # the axial strain at the centroid corresponding to Mult