"""Functions from Section 8.4 of EN 1992-1-1:2023."""

import math
from typing import Literal, Optional


def dv(dvx: float, dvy: float) -> float:
    """Calculate the shear-resisting effective depth of a slab.

    EN1992-1-1:2023 Eq. (8.91)

    The shear-resisting effective depth (dv) is calculated
        as the average of dvx and dvy.

    Args:
        dvx (float): Nominal value in the x-direction in mm.
        dvy (float): Nominal value in the y-direction in mm.

    Returns:
        float: Shear-resisting effective depth [mm].

    Raises:
        ValueError: If dvx or dvy is negative.
    """
    if dvx < 0:
        raise ValueError(f'dvx must not be negative. Got {dvx}')
    if dvy < 0:
        raise ValueError(f'dvy must not be negative. Got {dvy}')

    return (dvx + dvy) / 2


def beta_e(
    support_type: Literal[
        'internal columns',
        'edge columns',
        'corner columns',
        'ends of walls',
        'corners of walls',
    ],
    ebx: Optional[float] = None,
    eby: Optional[float] = None,
    bb_min: Optional[float] = None,
    bb_max: Optional[float] = None,
    refined: bool = False,
) -> float:
    """Calculate the coefficient βe accounting for
        concentrations of shear forces.

    EN1992-1-1:2023 Table 8.3

    The coefficient βe is determined based on the support type
    and whether a refined calculation is required.
    The refined calculation considers the eccentricities ebx, eby,
    and the geometric mean of the control perimeter widths bb.

    Args:
        support_type (Literal): Type of support
            ('internal columns', 'edge columns',
            'corner columns', 'ends of walls', 'corners of walls').
        ebx (float, optional): Eccentricity in the x-direction in mm.
        eby (float, optional): Eccentricity in the y-direction in mm.
        bb_min (float, optional): Minimum width of the control perimeter in mm.
        bb_max (float, optional): Maximum width of the control perimeter in mm.
        refined (bool, optional): Whether to use the refined calculation.
            Default is False.

    Returns:
        float: Coefficient βe.

    Raises:
        ValueError: If bb_min or bb_max is negative, or
            if bb_min is greater than bb_max.

    Notes:
        - ebx, eby are the eccentricities of the line
            of action of the support forces.
        - bb is the geometric mean of the minimum and
            maximum widths of the control perimeter.
    """
    if refined and support_type in (
        'internal columns',
        'edge columns',
        'corners of walls',
    ):
        if bb_min is None or bb_max is None or bb_min < 0 or bb_max < 0:
            raise ValueError(
                'bb_min and bb_max must not be negative. '
                + f' Got bb_min={bb_min}, bb_max={bb_max}'
            )

        if ebx is None or eby is None:
            raise ValueError('ebx and eby must not be None. ')

        if bb_min > bb_max:
            raise ValueError(
                'bb_min must not be greater than bb_max. '
                + f' Got bb_min={bb_min}, bb_max={bb_max}'
            )

        bb = math.sqrt(bb_min * bb_max)
        if support_type == 'internal columns':
            eb = math.sqrt(ebx**2 + eby**2)
        elif support_type == 'edge columns':
            eb = 0.5 * abs(ebx) + abs(eby)
        elif support_type == 'corners of walls':
            eb = 0.27 * (abs(ebx) + abs(eby))

        return max(1.05, 1 + 1.1 * eb / bb)

    if support_type == 'internal columns':
        return 1.15
    if support_type == 'edge columns':
        return 1.4
    if support_type == 'corner columns':
        return 1.5
    if support_type == 'ends of walls':
        return 1.4

    return 1.2


def tau_Ed_punch(
    V_ed: float,
    beta_e: float,
    b_0_5: float,
    d_v: float,
) -> float:
    """Calculate the punchinh design shear stress (τEd).

    EN1992-1-1:2023 Eq. (8.92)

    Args:
        V_ed (float): Design shear force at the control perimeter in kN.
        beta_e (float): Coefficient accounting for concentrations of
            the shear forces.
        b_0_5 (float): Length of the control perimeter in mm.
        d_v (float): Effective depth in mm.

    Returns:
        float: Design shear stress τEd in MPa.

    Raises:
        ValueError: If any input is negative or beta_e is not greater than 0.
    """
    if beta_e <= 0:
        raise ValueError(f'beta_e must be greater than 0. Got {beta_e}')
    if b_0_5 < 0:
        raise ValueError(f'b_0_5 must not be negative. Got {b_0_5}')
    if d_v < 0:
        raise ValueError(f'd_v must not be negative. Got {d_v}')

    return beta_e * V_ed * 1000 / (b_0_5 * d_v)


def tau_Ed_punch_d(v_ed: float, d_v: float) -> float:
    """Calculate the design shear stress (τEd) directly from a
        detailed analysis of the shear stress distribution.

    EN1992-1-1:2023 Eq. (8.93).

    Args:
        v_ed (float): Shear force per unit width in kN/m
        d_v (float): Effective depth mm

    Returns:
        float: Design shear stress τEd in MPa.

    Raises:
        ValueError: If any input is negative.
    """
    if d_v < 0:
        raise ValueError(f'd_v must not be negative. Got {d_v}')

    return v_ed / d_v


def tau_Rd_punch(
    gamma_v: float,
    rho_l_x: float,
    rho_l_y: float,
    f_ck: float,
    d_v: float,
    d_dg: float,
    b_0: float,
    b_0_5: float,
) -> float:
    """Calculate the design punching shear stress resistance.

    EN1992-1-1:2023 Eq. (8.94), (8.95), (8.96)

    Args:
        gamma_v (float): Partial factor for shear, unitless.
        rho_l_x (float): Reinforcement ratio in the x-direction, unitless.
        rho_l_y (float): Reinforcement ratio in the y-direction, unitless.
        f_ck (float): Characteristic compressive cylinder strength
            of concrete, in MPa.
        d_v (float): Effective depth of the slab, in mm.
        d_dg (float): Maximum aggregate size, in mm.
        b_0 (float): Length of the perimeter at the face of the
            supporting area, in mm.
        b_0_5 (float): Length of the perimeter at 0.5d from the
            face of the support, in mm.

    Returns:
        float: Design punching shear stress resistance, τ_Rd,c, in MPa.

    Raises:
        ValueError: If any input value is negative.
    """
    # Input validation
    if gamma_v <= 0:
        raise ValueError(f'gamma_v must be positive. Got {gamma_v}')
    if rho_l_x < 0 or rho_l_y < 0:
        raise ValueError(
            'Reinforcement ratios must not be negative. '
            + f' Got rho_l_x={rho_l_x}, rho_l_y={rho_l_y}'
        )
    if f_ck <= 0:
        raise ValueError(f'f_ck must be positive. Got {f_ck}')
    if d_v <= 0:
        raise ValueError(f'd_v must be positive. Got {d_v}')
    if d_dg <= 0:
        raise ValueError(f'd_dg must be positive. Got {d_dg}')
    if b_0 <= 0:
        raise ValueError(f'b_0 must be positive. Got {b_0}')
    if b_0_5 <= 0:
        raise ValueError(f'b_0_5 must be positive. Got {b_0_5}')

    # Calculate reinforcement ratio
    rho_l = math.sqrt(rho_l_x * rho_l_y)

    # Calculate k_pb
    k_pb = 3.6 * (1 - b_0 / b_0_5)
    k_pb = max(1, min(k_pb, 2.5))

    # Calculate τ_Rd,c
    tau_Rd_c = (
        (0.6 / gamma_v) * k_pb * ((100 * rho_l * f_ck * d_dg / d_v) ** (1 / 3))
    )

    # Limit τ_Rd,c
    return min(tau_Rd_c, (0.5 / gamma_v) * math.sqrt(f_ck))


def a_pd(a_p_x: float, a_p_y: float, d_v: float) -> float:
    """Calculate the modified effective depth for distances
    between the center of the support area and the point of contraflexure.

    EN1992-1-1:2023 Eq. (8.97), (8.98)

    Args:
        a_p_x (float): Maximum distance from the centroid of the
            control perimeter to the point where the bending moment
            is zero on the x-axis, in mm.
        a_p_y (float): Maximum distance from the centroid of
            the control perimeter to the point where
            the bending moment is zero on the y-axis, in mm.
        d_v (float): Effective depth of the slab, in mm.

    Returns:
        float: Modified effective depth, a_pd, in mm.

    Raises:
        ValueError: If any input value is negative.
    """
    # Input validation
    if a_p_x < 0:
        raise ValueError(f'a_p_x must not be negative. Got {a_p_x}')
    if a_p_y < 0:
        raise ValueError(f'a_p_y must not be negative. Got {a_p_y}')
    if d_v <= 0:
        raise ValueError(f'd_v must be positive. Got {d_v}')

    a_p = max(math.sqrt(a_p_x * a_p_y), d_v)

    return math.sqrt(a_p * d_v / 8)


def kpp(
    sigma_d: float,
    f_ck: float,
    b_0: float,
    b_0_5: float,
    is_compressive: bool,
    e_p: float = 0,
    d: float = 0,
) -> float:
    """Calculate the coefficient kpp for slabs with
        axial forces or prestressed slabs.

    EN1992-1-1:2023 Eq. (8.99), (8.100), (8.101), (8.103)

    Args:
        sigma_d (float): Average normal stress over the width in MPa.
        f_ck (float): Characteristic compressive cylinder strength of concrete
            in MPa.
        b_0 (float): Length of the perimeter at the face of the supporting
            area in mm.
        b_0_5 (float): Length of the perimeter at 0.5d from the face
            of the support in mm.
        is_compressive (bool): True for compressive axial forces,
            False for tensile forces.
        e_p (float, optional): Eccentricity of tendons, in mm. Default is 0.
        d (float, optional): Effective depth of the slab, in mm. Default is 0.

    Returns:
        float: Coefficient kpp, unitless.

    Raises:
        ValueError: If any input value is invalid.
    """
    # Input validation
    if sigma_d < 0:
        raise ValueError(f'sigma_d must not be negative. Got {sigma_d}')
    if f_ck <= 0:
        raise ValueError(f'f_ck must be positive. Got {f_ck}')
    if b_0 <= 0:
        raise ValueError(f'b_0 must be positive. Got {b_0}')
    if b_0_5 <= 0:
        raise ValueError(f'b_0_5 must be positive. Got {b_0_5}')
    if e_p < 0:
        raise ValueError(f'e_p must not be negative. Got {e_p}')
    if d < 0:
        raise ValueError(f'd must not be negative. Got {d}')

    # Calculate k_N
    k_N = math.sqrt(
        1
        + 0.47
        / (1 - b_0 / b_0_5)
        * abs(sigma_d)
        / math.sqrt(f_ck)
        * (1 + 6 * e_p / d)
    )

    # Calculate k_pp based on the type of axial force
    return k_N if is_compressive else 1 / k_N


def kpp_xy(kpp_x: float, kpp_y: float) -> float:
    """Comptues the geometric value for kpp when stresses in both directions.

    EN1992-1-1:2023 Eq. (8.102)

    Args:
        kpp_x (float): kpp coefficient in the x-direction, non dimensional.
        kpp_y (float): kpp coefficient in the y-direction, non dimensional.

    Returns:
        float: the computed kpp_xy coefficient, non-dimensional
    """
    if kpp_x < 0:
        raise ValueError(f'kpp_x must not be negative. Got {kpp_x}')
    if kpp_y < 0:
        raise ValueError(f'kpp_y must not be negative. Got {kpp_y}')

    return math.sqrt(kpp_x * kpp_y)


def tau_Rd_cs(
    tau_Rd_c: float,
    tau_Ed: float,
    rho_w: float,
    f_ywd: float,
    d_v: float,
    phi_w: float,
    d_dg: float,
    k_pb: float,
) -> float:
    """Calculate the punching shear resistance of
        slabs with shear reinforcement.

    EN1992-1-1:2023 Eq. (8.104), (8.105), (8.106)

    Args:
        tau_Rd_c (float): Punching shear stress resistance of slabs
            without shear reinforcement in MPa.
        tau_Ed (float): Design value of the shear stress in MPa.
        rho_w (float): Vertical shear reinforcement ratio (dimensionless).
        f_ywd (float): Yield strength of shear reinforcement in MPa.
        d_v (float): Effective depth of the slab in mm.
        phi_w (float): Diameter of the shear reinforcement in mm.
        d_dg (float): Maximum aggregate size in mm.
        k_pb (float): Coefficient related to the punching shear
            resistance (dimensionless).

    Returns:
        float: Punching shear resistance of slabs with shear
            reinforcement in MPa.

    Raises:
        ValueError: If any input is negative.
    """
    # Validate inputs
    tau_Rd_c = abs(tau_Rd_c)
    tau_Ed = abs(tau_Ed)
    if any(x < 0 for x in [rho_w, f_ywd, d_v, phi_w, d_dg, k_pb]):
        raise ValueError('Input values must not be negative.')

    # Calculate ηc
    eta_c = tau_Rd_c / tau_Ed

    # Calculate ηs
    eta_s = d_v / (150 * phi_w) + math.sqrt(15 * d_dg / d_v) * (
        1 / (eta_c * k_pb)
    ) ** (3 / 2)
    eta_s = min(eta_s, 0.8)

    # Calculate τRd,cs
    tau_Rd_cs = eta_c * tau_Rd_c + eta_s * rho_w * f_ywd
    return max(tau_Rd_cs, rho_w * f_ywd)


def rho_w(Asw: float, sr: float, st: float) -> float:
    """Calculate the vertical shear reinforcement ratio.

    EN1992-1-1:2023 Eq. (8.107)

    Args:
        Asw (float): Area of one leg of shear reinforcement in mm2.
        sr (float): Radial spacing of shear reinforcement in mm.
        st (float): Average tangential spacing of perimeters of
            shear reinforcement in mm.

    Returns:
        float: Vertical shear reinforcement ratio (dimensionless).

    Raises:
        ValueError: If any input is negative.
    """
    # Validate inputs
    if any(x < 0 for x in [Asw, sr, st]):
        raise ValueError('Input values must not be negative.')

    # Calculate ρw
    return Asw / (sr * st)


def dv_out(d_x: float, d_y: float, c_v: float) -> float:
    """Calculate the shear resisting effective depth of the
        outer shear reinforcement.

    EN1992-1-1:2023 Eq. (8.108)

    Args:
        d_x (float): Effective depth in the x-direction in mm.
        d_y (float): Effective depth in the y-direction in mm.
        c_v (float): Concrete cover to the shear reinforcement in mm.

    Returns:
        float: Effective depth of the outer shear reinforcement in mm.

    Raises:
        ValueError: If any input is negative.
    """
    # Validate inputs
    if any(x < 0 for x in [d_x, d_y, c_v]):
        raise ValueError('Input values must not be negative.')

    # Calculate dv,out
    return (d_x + d_y) / 2 - c_v


def tau_Rd_max_punch(
    tau_Rd_c: float,
    b_0: float,
    d_v: float,
    reinforcement_type: Literal['studs', 'links and stirrups'],
) -> float:
    """Calculate the maximum punching shear resistance.

    EN1992-1-1:2023 Eq. (8.109), (8.110), (8.111)

    Args:
        tau_Rd_c (float): Punching shear stress resistance of
            slabs without shear reinforcement in MPa.
        b_0 (float): Control perimeter at a distance dv,out/2 in mm.
        d_v (float): Effective depth of the slab in mm.
        reinforcement_type (str): Type of shear reinforcement,
            either 'studs' or 'links_and_stirrups'.

    Returns:
        float: Maximum punching shear resistance in MPa.

    Raises:
        ValueError: If any input is negative or reinforcement_type is invalid.
    """
    # Validate inputs
    if any(x < 0 for x in [tau_Rd_c, b_0, d_v]):
        raise ValueError('Input values must not be negative.')
    if reinforcement_type not in ['studs', 'links_and_stirrups']:
        raise ValueError(
            "reinforcement_type must be either 'studs' or 'links_and_stirrups'"
        )

    # Calculate ηsys
    if reinforcement_type == 'studs':
        eta_sys = max(0.70 + 0.63 * (b_0 / d_v) ** (1 / 4), 1)
    else:  # links_and_stirrups
        eta_sys = max(0.50 + 0.63 * (b_0 / d_v) ** (1 / 4), 1)

    # Calculate τRd,max
    return eta_sys * tau_Rd_c


def b0_5_out(b0_5: float, d_v: float, dv_out: float, eta_c: float) -> float:
    """Calculate the outer control perimeter at which
        shear reinforcement is not required.

    EN1992-1-1:2023 Eq. (8.112)

    Args:
        b0_5 (float): Control perimeter located at a distance dv/2
            from the face of the supporting area in mm.
        d_v (float): Effective depth of the slab in mm.
        dv_out (float): Shear-resisting effective depth
            of the outer shear reinforcement in mm.
        eta_c (float): Coefficient defined in Formula (8.105) (dimensionless).

    Returns:
        float: Outer control perimeter where shear reinforcement
            is not required in mm.

    Raises:
        ValueError: If any input is negative.
    """
    # Validate inputs
    if any(x < 0 for x in [b0_5, d_v, dv_out, eta_c]):
        raise ValueError('Input values must not be negative.')

    # Calculate b0,5,out
    return b0_5 * ((d_v / dv_out) * (1 / eta_c)) ** 2
