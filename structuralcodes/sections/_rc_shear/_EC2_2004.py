import typing as t

import numpy as np

from ...codes.ec2_2004 import (
    Asw_max,
    Asw_s_required,
    VEdmax_unreinf,
    VRdc,
    VRdc_prin_stress,
    VRdmax,
    VRds,
)
from ...geometry import SurfaceGeometry
from ...materials.reinforcement import ReinforcementEC2_2004
from ...sections import BeamSection


class ShearReinforcement:
    """Shear Reinforcement implementation."""

    def __init__(
        self,
        diameter: float,
        s: float,
        material: ReinforcementEC2_2004,
        n: t.Optional[int] = 2,
        alpha: t.Optional[float] = 90,
    ):
        """Initializes a new instance of Shear Reinforcement.

        Arguments:
            diameter (float): The diameter of the shear reinforcement
            s (float): The centre-to-centre distance of the shear reinforcement
                in mm.
            material (ReinforcementEC2_2004): A material for the shear
                reinforcement.
            n (Optional(int)): The number of legs in the shear reinforcement
                stirrups. Default value is 2.
            alpha (Optional(float)): The angle of the shear reinforcement with
                respect to the neutral axis in degrees. Default value is 90
                degrees.
        """
        Asw = n * diameter**2 / 4 * np.pi

        self.diameter = diameter
        self.Asw = Asw
        self.Asw_s = Asw / s
        self.s = s
        self.material = material
        self.n = n
        self.alpha = alpha


def shearcap_rectangular_section(
    section: BeamSection,
    NEd: float = 0,
    k1: float = 0.15,
):
    """Compute the design strength of the shear resistance and the maximum
        allowable shear force for rectangular cross-sections without shear
        reinforcement.

    EN 1992-1-1 (2005), Eq. (6.2) and (6.5).

    Arguments:
        section (BeamSection): The section to use as basis for the
            calculation.

    Keyword Args:
        NEd (float): The normal force in the cross-section due to loading or
            prestress (NEd > 0 for compression) in N. Default value is 0 N.
        k1 (float): Factor used to include the effect of the normal stress
            into the shear resistance of the concrete. Default value = 0.15,
            value might differ between National Annexes.

    Returns:
        float: The concrete shear resistance in N.
        float: The maximum allowable shear force in the cross-section in N.
            When a reduced shear force may be considered for the calculations,
            the unreduced shear force has to comply yo this value.
    """
    # Parameters from section geometry
    srf_geoms = section.geometry.geometries[0]
    concretebounds = srf_geoms.polygon.bounds
    h = concretebounds[3] - concretebounds[1]
    bw = concretebounds[2] - concretebounds[0]
    # Assume tensile reinforcement is in the lower half of the section
    reinf_bars = []
    for bar in section.geometry.point_geometries:
        if abs(concretebounds[1] - bar.point.bounds[1]) < abs(
            concretebounds[3] - bar.point.bounds[1]
        ):
            reinf_bars.append(bar)
    # Calculate d as the average height of the tensile reinforcement
    avg_reinf_height = 0
    for bar in reinf_bars:
        avg_reinf_height += bar.point.bounds[1]
    avg_reinf_height /= len(reinf_bars)
    d = h + concretebounds[1] - avg_reinf_height

    Ac = section.gross_properties.area
    Asl = 0
    for bar in reinf_bars:
        Asl += bar.area

    # Concrete parameters
    fck = srf_geoms.material.fck
    fcd = srf_geoms.material.fcd()
    gamma_c = srf_geoms.material.gamma_c

    design_value = VRdc(
        fck,
        d,
        Asl,
        bw,
        NEd,
        Ac,
        fcd,
        k1=k1,
        gamma_c=gamma_c,
    )

    max_allowable_shearforce = VEdmax_unreinf(
        bw,
        d,
        fck,
        fcd,
    )

    return design_value, max_allowable_shearforce


def shearcap_reinf_rectangular_section(
    section: BeamSection,
    shear_reinf: ShearReinforcement,
    NEd: float = 0,
    theta: float = 21.8,  # cot(theta) = 2.5
    limit_fyd: bool = False,
):
    """Calculate the shear resistance of vertical shear reinforcement and the
    maximum shear strength of the compression strut for rectangular
    cross-sections.

    EN 1992-1-1 (2005), Eq. (6.8) and (6.9).

    Argsuments:
        section (BeamSection): The section to use as basis for the
            calculation.
        shear_reinf (ShearReinforcement): The shear reinforcement to use as
            basis for the calculation.

    Keyword Args:
        NEd (float): The normal force in the cross-section due to loading or
            prestress (NEd > 0 for compression) in N. Default value is 0 N.
        theta (float): The angle of the compression strut in degrees. Default
            value is 21.8 degrees.
        limit_fyd (bool): Flag to indicate if the design yield stress is
            limited to 0.8 * fyk or not. This controls whether the stress
            reduction factor of concrete is given by Eq. (6.6) (False) or
            (6.10) (True). Default value is False.

    Returns:
        float: The shear resistance of the shear reinforcement in N.
        float: The shear strength of the compression strut in N.

    Raises:
        ValueError: When theta < 21.8 degrees or theta > 45 degrees.
        ValueError: The applied prestress exceeds the concrete design strength.
    """
    # Parameters from section geometry
    srf_geoms = section.geometry.geometries[0]
    concretebounds = srf_geoms.polygon.bounds
    h = concretebounds[3] - concretebounds[1]
    bw = concretebounds[2] - concretebounds[0]
    # Assume tensile reinforcement is in the lower half of the section
    reinf_bars = []
    for bar in section.geometry.point_geometries:
        if abs(concretebounds[1] - bar.point.bounds[1]) < abs(
            concretebounds[3] - bar.point.bounds[1]
        ):
            reinf_bars.append(bar)
    # Calculate d as the average height of the tensile reinforcement
    avg_reinf_height = 0
    for bar in reinf_bars:
        avg_reinf_height += bar.point.bounds[1]
    avg_reinf_height /= len(reinf_bars)
    d = h + concretebounds[1] - avg_reinf_height
    z = 0.9 * d

    Ac = section.gross_properties.area
    Asl = 0
    for bar in reinf_bars:
        Asl += bar.area

    # Concrete parameters
    fck = srf_geoms.material.fck
    fcd = srf_geoms.material.fcd()

    # Reinforcement parameters
    gamma_s = shear_reinf.material.gamma_s
    fyk = shear_reinf.material.fyk

    # Shear reinforcment parameters
    Asw = shear_reinf.Asw
    s = shear_reinf.s
    alpha = shear_reinf.alpha

    shear_resistance_reinforcement = VRds(
        Asw,
        s,
        z,
        theta,
        fyk,
        alpha,
        gamma_s,
    )

    shear_strength_compression_strut = VRdmax(
        bw,
        z,
        fck,
        theta,
        NEd,
        Ac,
        fcd,
        alpha,
        limit_fyd,
    )

    return shear_resistance_reinforcement, shear_strength_compression_strut


def required_shear_reinf(
    section: BeamSection,
    material: ReinforcementEC2_2004,
    VEd: float,
    theta: float = 21.8,  # Gir cot theta = 2.5
    alpha: float = 90,
    diameter: float = 10,
    n: int = 2,
):
    """Calculates the required shear reinforcement.

    EN 1992-1-1 (2005), Eq. (6.13).

    Arguments:
        section (BeamSection): The section to use as basis for the
            calculation.
        material (ReinforcementEC2_2004): The material of the shear
            reinforcement.
        VEd (float): The shear force in N.

    Keyword Args:
        theta (float): The angle of the compression strut in degrees. Default
            value is 21.8 degrees.
        alpha (float): The angle of the shear reinforcement with respect to the
            neutral axis in degrees. Default value is 90 degrees.
        diameter (float): The diameter of the shear reinforcement stirrups in
            mm. Default value is 10 mm.
        n (int): The number of legs in the shear reinforcement
                stirrups. Default value is 2.

    Returns:
        (ShearReinforcement): The required shear reinforcement expressed in the
            ShearReinforcement class.

    Raises:
        ValueError: When theta < 21.8 degrees or theta > 45 degrees.
    """
    # Parameters from section geometry
    srf_geoms = section.geometry.geometries[0]
    concretebounds = srf_geoms.polygon.bounds
    h = concretebounds[3] - concretebounds[1]
    # Assume tensile reinforcement is in the lower half of the section
    reinf_bars = []
    for bar in section.geometry.point_geometries:
        if abs(concretebounds[1] - bar.point.bounds[1]) < abs(
            concretebounds[3] - bar.point.bounds[1]
        ):
            reinf_bars.append(bar)
    # Calculate d as the average height of the tensile reinforcement
    avg_reinf_height = 0
    for bar in reinf_bars:
        avg_reinf_height += bar.point.bounds[1]
    avg_reinf_height /= len(reinf_bars)
    d = h + concretebounds[1] - avg_reinf_height
    z = 0.9 * d

    # Shear reinforcement parameters
    gamma_s = material.gamma_s
    fywk = material.fyk
    fywd = fywk / gamma_s

    req_reinf = Asw_s_required(
        VEd,
        z,
        theta,
        fywd,
        alpha,
    )

    s = (n * diameter**2 / 4 * np.pi) / req_reinf

    return ShearReinforcement(
        diameter=diameter, s=s, material=material, n=n, alpha=alpha
    )


def shearcap_rectangular_uncracked_prestressed(
    section: BeamSection,
    NEd: float,
    alpha_ct: float = 1.0,
    L_x: float = None,
    L_pt2: float = None,
):
    """Calculate the shear resistance in rectangular, uncracked, prestressed
    elements without shear reinforcement, value is determined via Mohr's
    circle.

    The maximal value of the principle tensile stress does not necessarily lay
    at the centre of gravity. If this is the case the minimum value of the
    shear resistance and corresponding stress needs to be found at the relevant
    location.

    EN 1992-1-1 (2005), Eq. (6.4).

    Arguments:
        section (BeamSection): The section to use as basis for the
            calculation.
        NEd (float): The normal force in the cross-section due to loading or
            prestress (NEd > 0 for compression) in N.

    Keyword Args:
        alpcha_ct (float): Coefficient for taking account of long term effects
            on the tensile strength and of unfavourable effects, resulting
            from the way the load is applied. Default value is 1.
        L_x (float): Distance from the considered cross-section until the
            starting point of the transference length of the prestress steel.
            This value should be provided when the prestressing steel is
            prestreched. Default value is None.
        L_pt2 (float): Maximum value of the transference length of the
            prestress steel, according to Eq. (8.18). This value should be
            provided when the prestressing steel is prestreched. Default value
            is None.

    Returns:
        float: The maximum allowable shear force in N for an uncracked,
            prestressed element without shear reinfordement, determined from
            maximum allowable principle stress.
    """
    # Parameters from section geometry
    srf_geoms = section.geometry.geometries[0]
    concretebounds = srf_geoms.polygon.bounds
    bw = concretebounds[2] - concretebounds[0]
    Ac = section.gross_properties.area
    # Translate section so that coord (0,0) is aligned with centroid
    c = srf_geoms.centroid
    srf_geoms = srf_geoms.translate(-c[0], -c[1])
    section = BeamSection(srf_geoms)
    # Split section to obtain correct value for S
    split_poly = srf_geoms.split(((0, 0), 0))[0][0]
    split_geo = SurfaceGeometry(split_poly, srf_geoms.material)
    split_sec = BeamSection(split_geo)

    Iy = section.gross_properties.iyy
    S = split_sec.gross_properties.sy

    # Concrete properties
    gamma_c = srf_geoms.material.gamma_c
    fctk_5 = srf_geoms.material.fctk_5
    fctd = alpha_ct * fctk_5 / gamma_c  # Should be taken from concrete class

    return VRdc_prin_stress(
        Iy,
        bw,
        S,
        fctd,
        NEd,
        Ac,
        L_x,
        L_pt2,
    )


def max_area_shear_reinf(
    section: BeamSection,
    shear_reinf: ShearReinforcement,
    NEd: float = 0,
):
    """Calculate the maximum cross-sectional area of the shear reinforcement
    based on hte assumption 1/tan(theta) == 1.

    EN 1992-1-1 (2005), Eq. (6.13)

    Arguments:
        section (BeamSection): The section to use as basis for the
            calculation.
        shear_reinf (ShearReinforcement): The shear reinforcement to use as
            basis for the calculation.

    Keyword Args:
        NEd (float): The normal force in the cross-section due to loading or
            prestress (NEd > 0 for compression) in N. Default value is 0.

    Returuns:
        float: The maximum allowable cross-sectional area of the shear
            reinforcement in mm2.

    Raises:
        ValueError: The applied prestress exceeds the concrete design strength.
    """
    # Parameters from section geometry
    srf_geoms = section.geometry.geometries[0]
    concretebounds = srf_geoms.polygon.bounds
    bw = concretebounds[2] - concretebounds[0]
    Ac = section.gross_properties.area

    # Concrete parameters
    fck = srf_geoms.material.fck
    fcd = srf_geoms.material.fcd()

    # Shear reinforcement parameters
    s = shear_reinf.s
    fywd = shear_reinf.material.fyd()
    alpha = shear_reinf.alpha

    return Asw_max(
        fcd,
        fck,
        bw,
        s,
        fywd,
        NEd,
        Ac,
        alpha,
    )
