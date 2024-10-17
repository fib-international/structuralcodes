"""Functions from Section 11 of EN 1992-1-1:2023."""

import math
from typing import Literal, Union

from scipy.interpolate import RegularGridInterpolator


def min_clear_distance_between_bars(phi: float, D_upper: float) -> float:
    """Calculate the clear distance between individual parallel bars.

    EN1992-1-1:2023 Eq. (11.2-2)

    Args:
        phi (float): Diameter of the bar in mm.
        D_upper (float): Distance to the upper bar layer in mm.

    Returns:
        float: The minimum clear distance in mm.

    Raises:
        ValueError: If any input is negative.
    """
    if phi < 0:
        raise ValueError(f'phi must not be negative. Got {phi}')
    if D_upper < 0:
        raise ValueError(f'D_upper must not be negative. Got {D_upper}')

    return max(phi, D_upper + 5, 20)


def min_clear_distance_to_poured_concrete(
    surface_roughness: bool, min_bond_distance: float
) -> float:
    """Calculate the clear distance between the face
        of already poured concrete and a parallel bar.

    EN1992-1-1:2023 Eq. (11.2-4)

    Args:
        surface_roughness (bool): Whether the surface is at least rough.
        min_bond_distance (float): Minimum distance required
            for bond according to 6.5.2.3(1) in mm.

    Returns:
        float: The minimum clear distance in mm.

    Raises:
        ValueError: If min_bond_distance is negative.
    """
    if min_bond_distance < 0:
        raise ValueError(
            f'min_bond_distance must not be negative. Got {min_bond_distance}'
        )

    if surface_roughness:
        return 5.0
    return min_bond_distance


def min_mandrel_phi(phi: float) -> float:
    """Calculate the minimum mandrel diameter for bending bars.

    EN1992-1-1:2023 Eq. (11.3-2)

    Args:
        phi (float): Diameter of the bar in mm.

    Returns:
        float: The minimum mandrel diameter in mm.

    Raises:
        ValueError: If phi is negative.
    """
    if phi < 0:
        raise ValueError(f'phi must not be negative. Got {phi}')

    if phi <= 16:
        return 4 * phi
    return 7 * phi


def k_bend(alpha_bend: float) -> float:
    """Calculate the bend parameter k_bend based on the bend angle alpha_bend.

    EN1992-1-1:2023 Eq. (11.3-4)

    Args:
        alpha_bend (float): Bend angle in degrees.

    Returns:
        float: The bend parameter k_bend.

    Raises:
        ValueError: If alpha_bend is not greater than zero.
    """
    if alpha_bend <= 0:
        raise ValueError(
            f'alpha_bend must be greater than zero. Got {alpha_bend}'
        )

    return 32 * (45 / alpha_bend)


def sigma_sd_bend_concrete(
    fck: float,
    gamma_c: float,
    d_g: float,
    phi: float,
    phi_mand: float,
    c_d: float,
    k_bend: float,
) -> float:
    """Calculate the steel stress limit to avoid concrete
        failure inside the bend.

    EN1992-1-1:2023 Eq. (11.1)

    Args:
        fck (float): Characteristic compressive strength of concrete in MPa.
        gamma_c (float): Partial safety factor for concrete.
        d_g (float): Maximum aggregate size in mm.
        phi (float): Diameter of the bar in mm.
        phi_mand (float): Mandrel diameter in mm.
        c_d (float): Minimum clear distance from the
            edge or between bars in mm.
        k_bend (float): Bend parameter considering bend angle alpha_bend.

    Returns:
        float: The maximum permissible steel stress σ_sd in MPa.

    Raises:
        ValueError: If any input is invalid.
    """
    if fck < 0:
        raise ValueError(f'fck must not be negative. Got {fck}')
    if gamma_c <= 0:
        raise ValueError(f'gamma_c must be greater than zero. Got {gamma_c}')
    if d_g < 0:
        raise ValueError(f'd_g must not be negative. Got {d_g}')
    if phi < 0:
        raise ValueError(f'phi must not be negative. Got {phi}')
    if phi_mand < 0:
        raise ValueError(f'phi_mand must not be negative. Got {phi_mand}')
    if c_d < 0:
        raise ValueError(f'c_d must not be negative. Got {c_d}')
    if k_bend < 0:
        raise ValueError(f'k_bend must not be negative. Got {k_bend}')

    a = 0.65 * fck / gamma_c * phi_mand / phi
    b = math.sqrt(fck) / gamma_c * ((d_g / phi) ** (1 / 3))
    c = c_d / phi + 1 / 2
    d = k_bend + 0.7 * phi_mand / phi

    return a + b * c * d


def k_trans(
    n_trans: int,
    phi: float,
    phi_mand: float,
    phi_trans: float,
    alpha_bend: float,
) -> float:
    """Calculate the increase factor for steel stress
        limit when transverse bars are within the bend.

    EN1992-1-1:2023 Eq. (11.2)

    Args:
        n_trans (int): Number of transverse bars within the bend.
        phi (float): Diameter of the main bar in mm.
        phi_mand (float): Mandrel diameter in mm.
        phi_trans (float): Diameter of the transverse bars in mm.
        alpha_bend (float): Bend angle in degrees. Must be greater than zero.

    Returns:
        float: The increase factor for steel stress limit.

    Raises:
        ValueError: If any input is invalid.

    Reference:
        EN1992-1-1:2023 Eq. (11.2)
    """
    if n_trans < 0:
        raise ValueError(f'n_trans must not be negative. Got {n_trans}')
    if phi < 0:
        raise ValueError(f'phi must not be negative. Got {phi}')
    if phi_mand < 0:
        raise ValueError(f'phi_mand must not be negative. Got {phi_mand}')
    if phi_trans < 0:
        raise ValueError(f'phi_trans must be non-negative. Got {phi_trans}')
    if alpha_bend <= 0:
        raise ValueError(
            f'alpha_bend must be greater than zero. Got {alpha_bend}'
        )

    phi_trans = min(phi_trans, 1.35 * phi)
    return 1 + 4 * n_trans * (phi / phi_mand) * pow(phi_trans / phi, 2) * (
        45 / alpha_bend
    )


def lbd_simple(
    fck: float,
    phi: float,
    bond_condition: Literal['good', 'poor'],
    poor_bond_multiplier: float = 1.2,
) -> float:
    """Calculate the design anchorage length for
        straight bars using linear interpolation for fck.

    EN1992-1-1:2023 Table. (11.1)

    Args:
        fck (float): Characteristic compressive strength of concrete in MPa.
        phi (float): Diameter of the reinforcing bar in mm.
        bond_condition (str): Bond condition ('good' or 'poor').
        poor_bond_multiplier (float, optional): Multiplier for
            poor bond conditions. Defaults to 1.2.

    Returns:
        float: Design anchorage length in mm.

    Raises:
        ValueError: If `phi` is not supported or if `fck` is outside
            the supported range for the given `phi`.
    """
    if phi <= 0:
        raise ValueError(f'Bar diameter must be positive. Got {phi}')

    # Good bond conditions anchorage lengths table from Table 11.1 (NDP)
    good_bond_length_table = {
        12: {20: 47, 25: 42, 30: 38, 35: 36, 40: 33, 45: 31, 50: 30, 60: 27},
        14: {20: 50, 25: 44, 30: 41, 35: 38, 40: 35, 45: 33, 50: 31, 60: 29},
        16: {20: 52, 25: 46, 30: 42, 35: 39, 40: 37, 45: 35, 50: 33, 60: 30},
        20: {20: 56, 25: 50, 30: 46, 35: 42, 40: 40, 45: 37, 50: 35, 60: 32},
        25: {20: 60, 25: 54, 30: 49, 35: 46, 40: 43, 45: 40, 50: 38, 60: 35},
        28: {20: 63, 25: 56, 30: 51, 35: 47, 40: 44, 45: 42, 50: 40, 60: 36},
        32: {20: 65, 25: 58, 30: 53, 35: 49, 40: 46, 45: 44, 50: 41, 60: 38},
    }

    # Get the unique sorted lists of phis and fck values from the table
    phi_values = sorted(good_bond_length_table.keys())
    fck_values = sorted(good_bond_length_table[phi_values[0]].keys())

    # Create a matrix of lbd/phi values for interpolation
    anchorage_matrix = [
        [good_bond_length_table[p][fck] for fck in fck_values]
        for p in phi_values
    ]

    # Create an interpolation function
    interpolation_function = RegularGridInterpolator(
        (phi_values, fck_values), anchorage_matrix, method='linear'
    )

    # Validate fck and phi are within the interpolation range
    if fck < fck_values[0] or fck > fck_values[-1]:
        raise ValueError(
            f'fck value {fck} is outside the supported'
            + f' range ({fck_values[0]}, {fck_values[-1]})'
        )
    if phi < phi_values[0] or phi > phi_values[-1]:
        raise ValueError(
            f'phi value {phi} is outside the supported'
            + f' range ({phi_values[0]}, {phi_values[-1]})'
        )
    lbd_phi = interpolation_function((phi, fck))

    if bond_condition.lower() == 'poor':
        lbd_phi *= poor_bond_multiplier

    return lbd_phi * phi


def lbd(
    phi: float,
    fck: float,
    sigma_sd: float,
    cd: float,
    c: float,
    bond_condition: Literal['good', 'poor', 'bentonite'],
    design_situation: Literal['persistent', 'accidental'] = 'persistent',
) -> float:
    """Calculate the design anchorage length for straight bars using.

    EN1992-1-1:2023 Eq. (11.3)

    Args:
        phi (float): Diameter of the reinforcing bar in mm.
        fck (float): Characteristic compressive strength of concrete in MPa.
        sigma_sd (float): Design value of the tensile stress in the bar in MPa.
        cd (float): Minimum nominal cover to the reinforcement in mm in X or Y,
            whatever is smaller.
        cx (float): Minimum separation between bars in mm in X or Y,
            whatever is smaller.
        bond_condition (str): Bond condition ('good', 'poor', 'bentonite').
        design_situation (str): Design situation ('persistent', 'accidental').

    Returns:
        float: Design anchorage length in mm.

    Raises:
        ValueError: If any of the input values are
            outside the specified limits.

    """
    # Validate inputs
    if phi <= 0:
        raise ValueError(f'Bar diameter must be positive. Got {phi}')
    if fck <= 0:
        raise ValueError(
            f'Concrete compressive strength must be positive. Got {fck}'
        )
    if sigma_sd <= 0:
        raise ValueError(f'Tensile stress must be positive. Got {sigma_sd}')
    if cd <= 0:
        raise ValueError(f'Nominal cover must be positive. Got {cd}')

    # Define coefficients
    if design_situation.lower() == 'persistent':
        klb = 50
        n_sigma = 3 / 2
    elif design_situation.lower() == 'accidental':
        klb = 35
        n_sigma = 3 / 2
    else:
        raise ValueError(f'Unknown design situation: {design_situation}')

    if bond_condition.lower() == 'good':
        kcp = 1.0
    elif bond_condition.lower() == 'poor':
        kcp = 1.2
    elif bond_condition.lower() == 'bentonite':
        kcp = 1.4
    else:
        raise ValueError(f'Unknown bond condition: {bond_condition}')

    # Calculate minimum cover cd
    min_cd = min(0.5 * cd, c, 3.75 * phi)

    # Calculate the ratios with limits
    phi_ratio = max(phi / 20, 0.6)
    fck_ratio = max(25 / fck, 0.3)

    # Calculate the design anchorage length
    lbd = (
        klb
        * kcp
        * phi
        * ((sigma_sd / 435) ** n_sigma)
        * (fck_ratio**0.5)
        * (phi_ratio ** (1 / 3))
        * ((1.5 * phi) / min_cd) ** 0.5
    )

    # Ensure the minimum anchorage length
    return max(lbd, 10 * phi)


def calculate_bond_condition(
    inclination_angle: float,
    distance_from_bottom: float,
    distance_from_surface: float,
) -> str:
    """Determine bond condition based on bar inclination
        and positioning during concreting.

    EN-1992-1-1:2023 11.4.2(4)

    Args:
        inclination_angle (float): Inclination angle of the
            bar to the horizontal (degrees).
        distance_from_bottom (float): Distance of the bar from
            the bottom of the formwork in mm.
        distance_from_surface (float): Distance of the bar from
            the free surface during concreting in mm.

    Returns:
        str: 'good' for good bond conditions, 'poor' for poor bond conditions.
    """
    if inclination_angle < 0 or inclination_angle > 90:
        raise ValueError(
            'Inclination angle must be between '
            + f' 0 and 90 degrees. Got {inclination_angle}'
        )

    # Determine bond condition based on provided rules
    if inclination_angle >= 45:
        return 'good'
    if distance_from_bottom <= 300 or distance_from_surface >= 300:
        return 'good'
    return 'poor'


def cd_conf(
    cx: float,
    cy: float,
    phi_t: float,
    st: float,
    cs: float,
    phi: float,
    rho_conf: float,
    sigma_ccd: float,
    fck: float,
) -> float:
    """Calculate the design anchorage length reduction parameter `cd,conf`.

    EN-1992-1-1:2023 Eq. 11.4

    Args:
        cx (float): minimum cover to longitudinal
            reinforcement along x-axis in mm.
        cy (float): minimum cover to longitudinal
            reinforcement along y-axis in mm.
        phi_t (float): Diameter of transverse reinforcement in mm.
        st (float): Spacing of transverse reinforcement along the bar in mm.
        cs (float): Clear spacing between bars or bundles in mm.
        phi (float): Diameter of the bar to be anchored or spliced in mm.
        rho_conf (float): Ratio of the confinement reinforcement
        sigma_ccd (float): Design value of mean compression stress
            perpendicular to the potential splitting failure in MPa.
        fck (float): Characteristic compressive cylinder
            strength of concrete in MPa.

    Returns:
        float: Reduced design anchorage length in mm.

    Raises:
        ValueError: If any input dimension is negative.
    """
    fck = abs(fck)
    sigma_ccd = abs(sigma_ccd)
    if any(x < 0 for x in [cx, cy, phi_t, st, cs, phi, rho_conf]):
        raise ValueError('Input values must not be negative.')

    delta_cd = (70 * rho_conf + 12 * sigma_ccd / (fck**0.5)) * phi
    cd_conf = min(cx, cy + (25 * phi_t**2 / st), cs / 2, 3.75 * phi) + delta_cd

    return min(cd_conf, 6 * phi)


def rho_conf(nc: int, phi_c: float, nb: int, phi: float, sc: float) -> float:
    """Calculate the ratio of confinement reinforcement.

    EN1992-1-1:2023 Eq. (11.5)

    Args:
        nc (int): Number of legs of confinement reinforcement
        crossing the potential splitting surface.
        phi_c (float): Diameter of the confinement
            reinforcement in mm.
        nb (int): Number of anchored bars or pairs of lapped bars
            in the potential splitting surface.
        phi (float): Diameter of the bar to be
            anchored or spliced in mm.
        sc (float): Spacing of the confinement
            reinforcement along the bar in mm.

    Returns:
        float: Ratio of the confinement reinforcement.

    Raises:
        ValueError: If any input dimension is negative.
    """
    if any(x < 0 for x in [nc, phi_c, nb, phi, sc]):
        raise ValueError('Input values must not be negative.')

    return (nc * math.pi * (phi_c**2) / 4) / (nb * phi * sc)


def additional_transverse_confinement_reinf(
    phi: float, n1: int, n2: int
) -> Union[float, float]:
    """Calculate additional transverse and confinement reinforcement.

    EN1992-1-1:2023 (11.5.2(6))

    Args:
        phi (float): Diameter of bars in mm.
        n1 (int): Number of layers with bars anchored at the
        same point in the member.
        n2 (int): Number of bars anchored in each layer.

    Returns:
        tuple: Additional transverse and confinement reinforcement in mm2.

    Raises:
        ValueError: If any input dimension is negative.
    """
    if any(x < 0 for x in [phi, n1, n2]):
        raise ValueError('Input values must not be negative.')

    Ast = 0.20 * (phi**2) * n1
    Asc = 0.20 * (phi**2) * n2

    return Ast, Asc


def phi_b_anchor(area_bars: float) -> float:
    """Calculate the equivalent diameter of a bundle of bars.

    EN1992-1-1:2023 Eq. (11.6)


    Args:
        area_bars (float): Total area of all bars contained
            in the bundle in mm2.

    Returns:
        float: Equivalent diameter of the bundle in mm.

    Raises:
        ValueError: If the area is negative.
    """
    if area_bars < 0:
        raise ValueError('Area of bars must not be negative.')

    return (4 * area_bars / math.pi) ** 0.5


def cd(
    cs: float,
    cx: float,
    cy: float,
    cyb: float,
) -> float:
    """Calculate the nominal cover cd as defined in Figure 11.6c).

    EN-1992-1-1:2023 (11.4.4(Fig 11.6c))

    Args:
        cs (float): Minimum cover to the bar in mm.
        cx (float): Cover in x-direction in mm.
        cy (float): Cover in y-direction in mm.
        cyb (float): Cover to bottom surface in mm.

    Returns:
        float: Nominal cover cd in mm.

    Raises:
        ValueError: If any input value is negative.
    """
    if any(x < 0 for x in (cs, cx, cy, cyb)):
        raise ValueError('Cover dimensions must not be negative.')

    return min(cs / 2, cx, cy, cyb)


def lb_anchor_tension(lbd: float, phi: float) -> float:
    """Calculate the design anchorage length for bars
        with bends and hooks in tension.

    EN-1992-1-1:2023 (11.4.4(1))

    Args:
        lbd (float): Basic design anchorage length in mm.
        phi (float): Diameter of the bar in mm.

    Returns:
        float: Total design anchorage length lbd,tot [mm].

    Raises:
        ValueError: If lbd or phi is negative.

    Reference:
        EN1992-1-1:2023 Eq. (XX), Figure 11.7a)
    """
    if lbd < 0 or phi < 0:
        raise ValueError(
            f'lbd and phi must not be negative. Got lbd={lbd}, phi={phi}.'
        )

    lbd_tot = lbd - 15 * phi
    return max(lbd_tot, 10 * phi)


def lb_anchor_compression(
    lbd: float,
    phi: float,
    d_p: float,
) -> float:
    """Calculate the design anchorage length for
        bars with bends and hooks in compression.

    Args:
        lbd (float): Basic design anchorage length in mm.
        phi (float): Diameter of the bar in mm.
        d_p (float): Distance of free surfaces
            perpendicular to the bar in mm.

    Returns:
        float: Total design anchorage length in mm.

    Raises:
        ValueError: If lbd, phi, or d_p is negative.

    Reference:
        EN1992-1-1:2023 Figure 11.7b)
    """
    if lbd < 0 or phi < 0 or d_p < 0:
        raise ValueError(
            'lbd, phi, and perpendicular_distance must not be negative. '
            + f'Got lbd={lbd}, phi={phi}, d_p={d_p}.'
        )

    if d_p >= 3.5 * phi:
        return max(lbd - 15 * phi, 10 * phi)

    return lbd


def lb_anchor_weld_trans(
    lbd: float,
    phi: float,
    phi_t: float,
    s: float,
    number_of_transverse_bars: int,
) -> float:
    """Calculate the design anchorage length for bars
        with welded transverse reinforcement in tension and compression.

    EN1992-1-1:2023 (11.4.5)

    Args:
        lbd (float): Basic design anchorage length in mm.
        phi (float): Diameter of the main bar in mm.
        phi_t (float): Diameter of the transverse bar in mm.
        s (float): Spacing of transverse bars in mm.
        number_of_transverse_bars (int): Number of transverse
            bars within the anchorage length.

    Returns:
        float: Total design anchorage length in mm.

    Raises:
        ValueError: If lbd, phi, phi_t, s are negative or if
            the number_of_transverse_bars is not a positive integer.
    """
    if lbd < 0 or phi < 0 or phi_t < 0 or s < 0:
        raise ValueError(
            'lbd, phi, phi_t, and s must not be negative. Got '
            + f'lbd={lbd}, phi={phi}, phi_t={phi_t}, s={s}.'
        )
    if number_of_transverse_bars <= 0:
        raise ValueError(
            'number_of_transverse_bars must be a positive integer. '
            + f'Got {number_of_transverse_bars}.'
        )

    if (phi_t >= 0.6 * phi and number_of_transverse_bars < 1) or (
        phi_t < 0.6
        and (not (50 <= s <= 100) or number_of_transverse_bars < 2 or phi > 16)
    ):
        return lbd

    # Reduce design anchorage length by 15*phi, but not less than 5*phi
    lbd_tot = lbd - 15 * phi
    return max(lbd_tot, 5 * phi)


def lb_anchor_ubar_loops_tension(lbd: float, phi: float) -> float:
    """Calculate the design anchorage length for U-bar loops in pure tension.

    EN1992-1-1:2023 (11.4.6)

    Args:
        lbd (float): Basic design anchorage length in mm.
        phi (float): Diameter of the bar in mm.

    Returns:
        float: Total design anchorage length in mm.

    Raises:
        ValueError: If lbd or phi is negative.

    Reference:

    """
    if lbd < 0 or phi < 0:
        raise ValueError(
            f'lbd and phi must not be negative. Got lbd={lbd}, phi={phi}.'
        )

    # Reduce design anchorage length by 20*phi, but not less than 10*phi
    lbd_tot = lbd - 20 * phi
    return max(lbd_tot, 10 * phi)


def phi_h(Ah: Union[float, int]) -> float:
    """Calculate the equivalent diameter of a circular head
         based on the head area Ah.

    EN1992-1-1:2023 Eq. (11.7)

    Args:
        Ah (float): Total area of the head in mm2.

    Returns:
        float: Equivalent diameter of the circular head in mm.

    Raises:
        ValueError: If Ah is negative.
    """
    if Ah < 0:
        raise ValueError(f'Ah must not be negative. Got Ah={Ah}.')

    return 2 * math.sqrt(Ah / math.pi)


def check_anchorage_headed_bars_in_tension(
    fck: float,
    phi: float,
    phi_h: float,
    th: float,
    ay: float,
    ax: float,
    sx: float,
    Ah: float,
    concrete_cracked: bool,
) -> bool:
    """Check if a headed bar in tension satisfies
         requirements for anchorage without additional length.

    EN-1992-1-1:2023 11.4.7

    Args:
        fck (float): Characteristic compressive cylinder
            strength of concrete in MPa.
        phi (float): Diameter of the reinforcing bar in mm.
        phi_h (float): Diameter of the circular head or equivalent
            circular diameter in mm.
        th (float): Thickness of the head in mm.
        ay (float): Distance between the bar axis
            and the nearest edge in mm.
        ax (float): Minimum distance between
            the bar axis and a corner in mm.
        sx (float): Bar spacing of a group of
            bars along the considered edge in mm.
        Ah (float): Total area of the head in mm2.
        concrete_cracked (bool): True if the concrete
            is cracked, False if uncracked.

    Returns:
        bool: True if the configuration satisfies
            the conditions, False otherwise.

    Raises:
        ValueError: If any input value is negative.
    """
    if any(x < 0 for x in (fck, phi, phi_h, th, ay, ax, sx, Ah)):
        raise ValueError('Input values must not be negative.')

    # Check the given conditions for headed bars in tension
    if fck < 25 or phi > 25 or (phi_h < 3 * phi):
        return False

    if concrete_cracked and ay < 4 * phi:
        return False

    if ay < 3 * phi:
        return False

    if not (ax >= 2 * ay + 1.2 * phi_h):
        return False

    if not (sx >= 4 * ay):
        return False

    return True


def kh_A(phi_h: float, phi: float) -> float:
    """Calculate the ratio between the net area of the head and the
        cross-sectional area of the reinforcement.

    EN-1992-1-1:2023 Eq. (11.9)

    Args:
        phi_h (float): Diameter of the circular head in mm.
        phi (float): Diameter of the reinforcing bar in mm.

    Returns:
        float: Ratio kh,A.

    Raises:
        ValueError: If phi_h or phi is negative.
    """
    if phi_h < 0 or phi < 0:
        raise ValueError(
            'phi_h and phi must not be negative.'
            + f' Got phi_h={phi_h}, phi={phi}.'
        )

    return (phi_h / phi) ** 2 - 1


def a_d(
    ay: float,
    ax: float,
    sx: float,
    phi: float,
    phi_h: float,
) -> float:
    """Calculate the nominal value of the distance
        a_d between the bar and a free surface.

    EN-1992-1-1:2023 Eq. (11.10)

    Args:
        ay (float): Distance between the bar axis and the nearest edge in mm.
        ax (float): Minimum distance between the bar axis and a corner in mm.
        sx (float): Bar spacing of a group of bars along the
            considered edge in mm.
        phi (float): Diameter of the reinforcing bar in mm.
        phi_h (float): Diameter of the circular head in mm.

    Returns:
        float: Nominal distance a_d in mm.

    Raises:
        ValueError: If any input value is negative.
    """
    if any(x < 0 for x in (ay, ax, sx, phi, phi_h)):
        raise ValueError('Input values must not be negative.')

    if ax >= 2 * ay + 1.2 * phi_h and sx >= 4 * ay:
        # Single bar near an edge or group of bars with spacing sx >= 4ay
        return ay
    if ax < 2 * ay + 1.2 * phi_h:
        # Single bar near a corner
        return 0.5 * ay + 0.25 * ax - 0.3 * phi_h

    # Group of bars with spacing sx < 4ay
    a = ay * (sx - phi_h) / (4 * ay - phi_h)
    b = 2.3 * (ay - phi_h / 2)
    c = (4 * ay - sx) / (4 * ay - phi_h)
    d = 1 - 1 / (phi_h / phi) ** 2
    return a + b * c * d


def sigma_sd_prime(
    fck: float,
    gamma_c: float,
    phi: float,
    phi_h: float,
    a_d: float,
    kh_A: float,
    nu_part: float,
    ddg: float,
) -> float:
    """Calculate the maximum tensile stress in the
        reinforcing steel developed by the head.

    EN1992-1-1:2023 Eq. (11.8)

    Args:
        fck (float): Characteristic compressive cylinder
            strength of concrete in MPa.
        gamma_c (float): safety coefficient for concrete.
        phi (float): Diameter of the reinforcing bar in mm.
        phi_h (float): Diameter of the circular head in mm.
        a_d (float): Nominal value of the distance
            between the bar and a free surface in mm.
        kh_A (float): Ratio kh_A.
        nu_part (float): Coefficient for cracked or uncracked concrete.
            8.0 for concrete cracked region of the head and 11.0 for
            uncracked concrete in the region of the head.
        ddg (float): Maximum aggregate size in mm.

    Returns:
        float: Maximum tensile stress in the reinforcing steel in MPa.

    Raises:
        ValueError: If any input value is negative.

    Reference:

    """
    if any(x < 0 for x in (fck, phi, phi_h, a_d, fck, ddg, gamma_c)):
        raise ValueError('Input values must not be negative.')

    sigma_sd = kh_A * fck / gamma_c + (
        nu_part
        * (math.sqrt(fck) / gamma_c)
        * (a_d / phi)
        * (phi_h / phi) ** (5 / 6)
        * (ddg / phi) ** (1 / 3)
    )

    return min(sigma_sd, kh_A * nu_part * fck / gamma_c)


def lbd_bar_tension(
    lbd_sigma_sd: Union[float, int], lbd_sigma_sd_prime: Union[float, int]
) -> float:
    """Calculate the design length in the reinforcing
        bar lbd to develop the remaining stress.

    EN1992-1-1:2023 Eq. (11.11)

    Args:
        lbd_sigma_sd (float): Design length corresponding
            to the stress sigma_sd, in mm.
        lbd_sigma_sd_prime (float): Design length corresponding
            to the stress sigma_sd_prime, in mm.

    Returns:
        float: Design length lbd in mm.

    Raises:
        ValueError: If lbd_sigma_sd or lbd_sigma_sd_prime is negative.
    """
    if lbd_sigma_sd < 0 or lbd_sigma_sd_prime < 0:
        raise ValueError(
            f'lbd_sigma_sd and lbd_sigma_sd_prime must not be negative. '
            f'Got lbd_sigma_sd={lbd_sigma_sd}, '
            + f'lbd_sigma_sd_prime={lbd_sigma_sd_prime}.'
        )

    return 1.1 * (lbd_sigma_sd - lbd_sigma_sd_prime)


def cmin_b(
    phi: float,
    lbd_pi: float,
    drilling_method: Literal['rotary_percussion', 'diamond', 'compressed_air'],
    use_drilling_aid: bool,
) -> float:
    """Calculate the minimum concrete cover cmin,b
        for post-installed reinforcing steel bars.

    EN1992-1-1:2023 Table (11.2)

    Args:
        phi (float): Diameter of the reinforcing bar in mm.
        lbd_pi (float): Design anchorage length for post-installed bars in mm.
        drilling_method (str): Drilling method used
            ('rotary_percussion', 'diamond', 'compressed_air').
        use_drilling_aid (bool): True if a drilling
            aid is used, False otherwise.

    Returns:
        float: Minimum concrete cover cmin_b in mm.

    Raises:
        ValueError: If any input value is invalid.
    """
    if phi < 0 or lbd_pi < 0:
        raise ValueError(
            'phi and lbd_pi must not be negative. '
            + f'Got phi={phi}, lbd_pi={lbd_pi}.'
        )

    if drilling_method in ('rotary_percussion', 'diamond'):
        if phi < 25:
            base_cover = 30
            factor = 0.06 if not use_drilling_aid else 0.02
        else:
            base_cover = 40
            factor = 0.06 if not use_drilling_aid else 0.02
    elif drilling_method == 'compressed_air':
        if phi < 25:
            base_cover = 50
            factor = 0.08 if not use_drilling_aid else 0.02
        else:
            base_cover = 60
            factor = 0.08 if not use_drilling_aid else 0.02

    # Calculate cmin,b
    return base_cover + factor * lbd_pi


def lbd_pi(
    lbd: float,
    kb_pi: float,
    phi: float,
    alpha_lb: float = 1.5,
) -> float:
    """Calculate the design anchorage length lbd,pi
        of post-installed reinforcing steel bars in tension.

    EN1992-1-1:2023 Eq. (11.12)

    Args:
        lbd (float): Design anchorage length calculated
            according to 11.4.2 in mm.
        kb_pi (float): Bond efficiency factor.
        phi (float): Diameter of the reinforcing bar in mm.
        alpha_lb (float): Factor accounting for cracks along the bar.

    Returns:
        float: Design anchorage length lbd,pi in mm.

    Raises:
        ValueError: If any input value is invalid.

    """
    if lbd < 0 or phi < 0 or kb_pi <= 0:
        raise ValueError(
            f'lbd, phi must not be negative, and kb_pi must be positive. '
            f'Got lbd={lbd}, phi={phi}, kb_pi={kb_pi}.'
        )

    lbd_pi = lbd / kb_pi
    return max(lbd_pi, 10 * phi * alpha_lb)


def lsd(
    type_of_lap: Literal[
        'straight_bars',
        'bends_and_hooks',
        'loops',
        'headed_bars',
        'intermeshed_fabric',
        'layered_fabric',
        'bonded_postinstalled',
    ],
    state: Literal['tension', 'compression'],
    lbd: float,
    phi: float,
    kls: float = 1.2,
    lbd_pi: float = 0,
    alpha_lb: float = 1.5,
    phi_mand: float = 0,
) -> float:
    """Calculate the design lap length lsd based on the type of lap splice.

    EN1992-1-1:2023 Table (11.3)

    Args:
        type_of_lap (str): Type of lap splice
            ('straight_bars', 'bends_and_hooks', 'loops', 'headed_bars',
            'intermeshed_fabric', 'layered_fabric', 'bonded_postinstalled').
        state (str): where the lap is in tension or compression.
        lbd (float): Basic design anchorage length in mm.
        phi (float): Diameter of the reinforcing bar in mm.
        kls (float): Coefficient for lap length, default is
            1.2 unless specified otherwise.
        lbd_pi (float): Design anchorage length
            for post-installed bars in mm.
        alpha_lb (float): Factor accounting for cracks along the bar.
            Normally the value is 1.5.
        phi_mand (float or int): Mandrel diameter for loops in mm.

    Returns:
        float: Design lap length lsd in mm.

    Raises:
        ValueError: If any input value is invalid or the
            result is incompatible according to the Table 11.3.
    """
    if (
        phi < 0
        or lbd < 0
        or lbd_pi < 0
        or phi_mand < 0
        or kls < 0
        or alpha_lb < 0
    ):
        raise ValueError(
            'Invalid input values. phi, lbd, lbd_pi, alpha_lb, kls'
            + 'and phi_mand must not be negative.'
        )
    if (
        type_of_lap in ('straight_bars', 'bends_and_hooks')
        and state == 'compression'
    ):
        raise ValueError(f'Not possible to have {type_of_lap} laps in {state}')

    if type_of_lap in ('straight_bars', 'bends_and_hooks'):
        return max(kls * lbd, 15 * phi)
    if type_of_lap == 'loops':
        return phi_mand + 4 * phi  # Calculated according to 11.5.4
    if type_of_lap == 'headed_bars':
        return lbd  # Calculated according to 11.5.5
    if type_of_lap in ('intermeshed_fabric', 'layered_fabric'):
        return max(kls * lbd, 15 * phi, 250)

    # If bonded_postinstalled
    return max(kls * lbd_pi, 15 * phi * alpha_lb)


def cs_laps(
    phi: float,
) -> float:
    """Calculate the minimum clear distance between adjacent laps.

    EN1992-1-1:2023 Section 11.5.2 (Key 1)

    Args:
        phi (float): Diameter of the reinforcing bar in mm.

    Returns:
        float: Minimum clear distance between adjacent laps in mm.

    Raises:
        ValueError: If phi is negative.
    """
    if phi < 0:
        raise ValueError('phi must not be negative.')

    return max(2 * phi, 20)


def cs_bars(phi: float) -> float:
    """Calculate the clear distance between lapping bars.

    EN1992-1-1:2023 Section 11.5.2 (Key 2)

    Args:
        phi (float): Diameter of the reinforcing bar in mm.

    Returns:
        float: Clear distance between lapping bars in mm.

    Raises:
        ValueError: If phi is negative.
    """
    if phi < 0:
        raise ValueError('phi must not be negative.')

    return min(4 * phi, 50)


def min_clear_distante_laps(phi: float) -> float:
    """Calculate the minimum clear distance between lapped bars.

    EN1992-1-1:2023 Section 11.5.2(7)

    Args:
        phi (float): Diameter of the reinforcing bar in mm.

    Returns:
        float: Minimum clear distance between lapped bars in mm.

    Raises:
        ValueError: If phi is negative.
    """
    if phi < 0:
        raise ValueError('phi must not be negative.')

    return min(4 * phi, 50)


def Ac_u_bars(phi_mand: float, phi: float, lsd: float) -> float:
    """Calculate the total effective concrete area
        within the curved parts of the overlapping U-bars.

    EN1992-1-1:2023 Eq. (11.14)

    Args:
        phi_mand (float): Mandrel diameter for loops in mm.
        phi (float): Diameter of the reinforcing bar in mm.
        lsd (float): Lap length in mm.

    Returns:
        float: Total effective concrete area Ac in mm2.

    Raises:
        ValueError: If any input value is negative.
    """
    if phi_mand < 0 or phi < 0 or lsd < 0:
        raise ValueError(
            'Invalid input values. phi_mand, phi, and '
            + 'lsd must not be negative.'
        )

    return (phi_mand + phi) * (lsd - 0.21 * (phi_mand + phi))


def kst_u_bar_loops(
    Ast: float,
    fyd: float,
    Ac: float,
    fcd: float,
    lsd: float,
    ddg: float,
) -> float:
    """Calculate the resistance factor of the confinement reinforcement kst.

    EN1992-1-1:2023 Eq. (11.15)

    Args:
        Ast (float): Total area of the fully
        anchored confinement reinforcement positioned within Ac in mm2.
        fyd (float): Design yield strength
            of the confinement reinforcement in MPa.
        Ac (float): total effective concrete area within curved parts
            of overlapping U-bars in mm2.
        fcd (float): Design compressive strength of concrete in MPa.
        lsd (float): lap length in mm.
        ddg (float): Coefficient that takes into account
            the concrete type and its aggregate properties.

    Returns:
        float: Resistance factor kst.

    Raises:
        ValueError: If any input value is negative.
    """
    if Ast < 0 or fyd < 0 or fcd < 0 or ddg < 0 or Ac < 0:
        raise ValueError(
            'Invalid input values. All parameters must not be negative.'
        )

    w = (Ast * fyd) / (0.85 * pow((ddg / lsd), 1 / 3) * Ac * fcd)

    if w >= 0.5:
        return 1.0

    return 4 * w * (1 - w)


def TRd_c_u_bar_loops(
    fcd: float,
    lsd: float,
    ddg: float,
    kst: float,
    cs: float,
    Ac: float,
) -> float:
    """Calculate the resistance of a single lap splice TRd_c.

    EN1992-1-1:2023 Eq. (11.13)

    Args:
        fcd (float): Design compressive strength of concrete in MPa.
        lsd (float): Lap length in mm.
        ddg (float): Coefficient that takes into account the
            concrete type and its aggregate properties.
        kst (float): Resistance factor of the confinement reinforcement.
        cs (float): Clear spacing of U-bars in mm.
        Ac (float): total effective concrete area within curved parts
            of overlapping U-bars in mm2.

    Returns:
        float: Resistance of a single lap splice TRd,c in kN.

    Raises:
        ValueError: If any input value is negative.
    """
    if fcd < 0 or lsd < 0 or ddg < 0 or cs < 0 or Ac < 0:
        raise ValueError(
            'Invalid input values. All parameters must not be negative.'
        )

    trd_c = (
        (0.2 * fcd * Ac)
        * pow((ddg / lsd), 1 / 3)
        * (math.sqrt(kst + (cs / lsd) ** 2) - cs / lsd)
    )
    return trd_c / 1000  # To convert to kN


def Ast_min_u_bar_loops(fck: float, Ac: float, fyk: float) -> float:
    """Calculate the minimum amount of confinement reinforcement
        in a double symmetric configuration within Ac to avoid
        brittle behaviour.

    EN1992-1-1:2023 Eq. (11.16)

    Args:
        fck (float): Characteristic compressive
            strength of concrete in MPa.
        Ac (float): Total effective concrete area in mm2.
        fyk (float): Characteristic yield strength of
            the confinement reinforcement in MPa.

    Returns:
        float: Minimum amount of confinement reinforcement Ast in mm2.

    Raises:
        ValueError: If any input value is negative.
    """
    if fck < 0 or Ac < 0 or fyk < 0:
        raise ValueError(
            'Invalid input values. All parameters must not be negative.'
        )

    return 0.5 * math.sqrt(fck) * Ac / fyk


def Ac_headed_laps(lsd: float, phi: float, bh1: float) -> float:
    """Calculate the effective concrete area within
        the heads of the overlapping bars.

    EN1992-1-1:2023 Eq. (11.18)

    Args:
        lsd (float): Design anchorage length in mm.
        phi (float): Diameter of the reinforcing bar in mm.
        bh1 (float): Effective width of the head
            perpendicular to the plane of the lap in mm.

    Returns:
        float: Effective concrete area Ac in mm2.

    Raises:
        ValueError: If any input value is negative.
    """
    if lsd < 0 or phi < 0 or bh1 < 0:
        raise ValueError(
            'Invalid input values. lbd, phi, and bh1 must not be negative.'
        )

    return (lsd - 2 * phi) * bh1


def bh1_headed_laps(phi_h: float) -> float:
    """Calculate the effective width of the head
        perpendicular to the plane of the lap for circular heads.

    EN1992-1-1:2023 Eq. (11.19)

    Args:
        phi_h (float): Diameter of the circular head in mm.

    Returns:
        float: Effective width bh1 in mm.

    Raises:
        ValueError: If phi_h is negative.
    """
    if phi_h < 0:
        raise ValueError('Invalid input value. phi_h must not be negative.')

    return 0.5 * phi_h * math.sqrt(math.pi)


def kst_headed_bars(
    Ast: float,
    fyd: float,
    fcd: float,
    Ac: float,
    lsd: float,
    phi: float,
    ddg: float,
) -> float:
    """Calculate the resistance factor of the transverse reinforcement kst.

    EN-1992-1-1:2023 Eq. (11.20)

    Args:
        Ast (float): Total area of the fully anchored
            transverse reinforcement positioned within Ac in mm2.
        fyd (float): Design yield strength of the
            transverse reinforcement in MPa.
        fcd (float): Design compressive
            strength of concrete in MPa.
        phi (float): Diameter of the reinforcing bar in mm.
        Ac (float):  Total effective concrete area in mm2.
        lsd (float): Design anchorage length in mm.
        ddg (float): Coefficient that takes into
            account the concrete type and its aggregate properties.

    Returns:
        float: Resistance factor kst.

    Raises:
        ValueError: If any input value is negative.
    """
    if Ast < 0 or fyd < 0 or fcd < 0 or phi < 0 or ddg < 0 or lsd < 0:
        raise ValueError(
            'Invalid input values. All parameters must not be negative.'
        )

    w = (Ast * fyd) / (1.3 * pow((ddg / (lsd - 2 * phi)), 1 / 3) * Ac * fcd)

    if w >= 0.5:
        return 1.0
    return 4 * w * (1 - w)


def TRd_c_headed(
    fcd: float,
    lsd: float,
    ddg: float,
    kst: float,
    cs: float,
    phi: float,
    bh1: float,
    Ac: float,
) -> float:
    """Calculate the resistance of a single headed bar lap TRd_c.

    EN-1992-1-1:2023 Eq. (11.17)

    Args:
        fcd (float): Design compressive strength of concrete in MPa.
        lsd (float): Design anchorage length in mm.
        ddg (float): Coefficient that takes into account
            the concrete type and its aggregate properties.
        kst (float): Resistance factor of the transverse reinforcement.
        cs (float): Clear spacing of headed bars in mm.
        phi (float): Diameter of the reinforcing bar in mm.
        bh1 (float): Effective width of the head perpendicular
            to the plane of the map in mm2.
        Ac (float): Effective concrete area in mm2.

    Returns:
        float: Resistance of a single headed bar lap TRd_c in kN.

    Raises:
        ValueError: If any input value is negative.
    """
    if fcd < 0 or lsd < 0 or ddg < 0 or cs < 0 or phi < 0 or bh1 < 0:
        raise ValueError(
            'Invalid input values. All parameters must not be negative.'
        )

    tRd_c = (
        (0.6 * fcd * Ac)
        * pow((ddg / (lsd - 2 * phi)), 1 / 3)
        * (math.sqrt(kst + (cs / (lsd - 2 * phi)) ** 2) - cs / (lsd - 2 * phi))
    )
    return tRd_c / 1000  # Convert to kN


def Ast_min_headed_bars(
    fck: float,
    Ac: float,
    fyk: float,
    phi: float,
) -> float:
    """Calculate the minimum amount of transverse reinforcement.

    EN-1992-1-1:2023 Eq. (11.21)

    Args:
        fck (float): Characteristic compressive strength of concrete in MPa.
        Ac (float): Effective concrete area in mm2.
        fyk (float): Characteristic yield strength
            of the transverse reinforcement in MPa.
        phi (float): Diameter of the reinforcing bar in mm.

    Returns:
        float: Minimum amount of transverse reinforcement Ast in mm2.

    Raises:
        ValueError: If any input value is negative.
    """
    if fck < 0 or Ac < 0 or fyk < 0 or phi < 0:
        raise ValueError(
            'Invalid input values. All parameters must not be negative.'
        )

    return max(0.75 * math.sqrt(fck) / fyk * Ac, (math.pi / 8) * phi**2)


def Ast_min_headed_single_bars(phi: float) -> float:
    """Calculate the minimum area of tie down reinforcement.

    EN-1992-1-1:2023 Eq. (11.22)

    Args:
        phi (float): Diameter of the reinforcing bar in mm.

    Returns:
        float: Minimum area of tie down reinforcement Astd in mm2.

    Raises:
        ValueError: If phi is negative.

    """
    if phi < 0:
        raise ValueError('Invalid input value. phi must not be negative.')

    return 0.12 * phi**2


def min_clear_distante_laps_couplers(
    Dupper: float, max_bar_diameter: float
) -> float:
    """Calculate the minimum clear distance (horizontal and vertical)
        between couplers and between couplers and adjacent bars.

    EN1992-1-1:2023 Section 11.5.6(1)

    Args:
        Dupper (float): Upper diameter of the coupler in mm.
        max_bar_diameter (float): Maximum diameter of the bars in mm.

    Returns:
        float: Minimum clear distance in mm.

    Raises:
        ValueError: If any input value is negative.
    """
    if Dupper < 0 or max_bar_diameter < 0:
        raise ValueError(
            'Invalid input values. Dupper and max_bar_diameter '
            + 'must not be negative.'
        )

    return max(Dupper + 5, max_bar_diameter)


def csx_duct(Dupper: float, phi_duct: float) -> float:
    """Calculate the minimum horizontal spacing
        csx for post-tensioning tendons.

    EN1992-1-1:2023 11.6.2(1)

    Args:
        Dupper (float): Upper diameter of the duct in mm.
        phi_duct (float): Diameter of the duct in mm.

    Returns:
        float: Minimum horizontal spacing csx in mm.

    Raises:
        ValueError: If any input value is negative.
    """
    if Dupper < 0 or phi_duct < 0:
        raise ValueError(
            'Invalid input values. Dupper and phi_duct must not be negative.'
        )

    return max(Dupper + 5, phi_duct, 50)


def csy_duct(Dupper: float, phi_duct: float) -> float:
    """Calculate the minimum vertical spacing csy
        for post-tensioning tendons.

    EN1992-1-1:2023 Section 11.6.2(1)

    Args:
        Dupper (float): Upper diameter of the duct in mm.
        phi_duct (float): Diameter of the duct in mm.

    Returns:
        float: Minimum vertical spacing csy in mm.

    Raises:
        ValueError: If any input value is negative.
    """
    if Dupper < 0 or phi_duct < 0:
        raise ValueError(
            'Invalid input values. Dupper and phi_duct must not be negative.'
        )

    return max(Dupper, phi_duct, 40)


def min_spacing_bundled_tendons() -> float:
    """Calculate the minimum spacing between bundled tendons.

    EN1992-1-1:2023 Section 11.6.2(2)

    Returns:
        float: Minimum spacing between bundled tendons in mm.
    """
    return 100.0


def Rmin_unconf_concrete(sigma_pd: float, pRd: float, Ap: float) -> float:
    """Calculate the minimum radius of curvature Rmin
        for tendons to prevent damage of the unconfined concrete.

    EN-1992-1-1:2023 Eq. (11.23)

    Args:
        sigma_pd (float): Tendon design stress in MPa.
        pRd (float): Maximum transverse bearing
            stress on the prestressing tendon in MPa.
        Ap (float): Cross-sectional area
            of the prestressing tendon in mm2.

    Returns:
        float: Minimum radius of curvature Rmin in mm.

    Raises:
        ValueError: If any input value is negative.
    """
    if sigma_pd < 0 or pRd < 0 or Ap < 0:
        raise ValueError(
            'Invalid input values. sigma_pd, pRd, and Ap must not be negative.'
        )

    return sigma_pd / pRd * math.sqrt(Ap)


def max_pRd(
    case: Literal[
        'internal_corrugated', 'internal_U_shape', 'external_smooth'
    ],
    fcd: float,
) -> float:
    """Get the maximum transverse bearing stress pRd
        based on the case and concrete design compressive strength.

    EN1992-1-1:2023 Table (11.4)

    Args:
        case (str): Case of the tendon ('internal_corrugated',
            'internal_U_shape', 'external_smooth').
        fcd (float): Design compressive strength of concrete in MPa.

    Returns:
        float: Maximum transverse bearing stress pRd in MPa.

    Raises:
        ValueError: If case is invalid or fcd is negative.
    """
    if fcd < 0:
        raise ValueError('Invalid input value. fcd must not be negative.')

    if case == 'internal_corrugated':
        return min(0.75 * fcd, 15)
    if case == 'internal_U_shape':
        return 70
    # if external_smooth
    return 30


def cu(
    cs: float,
    cy: float,
    phi: float,
) -> float:
    """Calculate the parameter cu.

    EN1992-1-1:2023 (11.7)

    Args:
        cs (float): Clear spacing of reinforcing bars or ducts in mm.
        cy (float): Distance from the center of the bar to the
            concrete surface in mm.
        phi (float): Diameter of the reinforcing bar or duct in mm.
            Use phi_duct for post-tensioning tendons.
        is_post_tensioning (bool): True if the calculation is for
            post-tensioning tendons, False otherwise.

    Returns:
        float: Parameter cu in mm.

    Raises:
        ValueError: If any input value is negative.
    """
    if cs < 0 or cy < 0 or phi < 0:
        raise ValueError(
            'Invalid input values. cs, cy, and phi must not be negative.'
        )

    return min(cs, 2 * math.sqrt(3) * (cy + 0.5 * phi))


def verify_deviation_forces(
    Ftd: float,
    r: float,
    cu: float,
    fck: float,
    gamma_c: float,
) -> bool:
    """Verify if the deviation forces can be resisted by
         the concrete in tension.

        EN1992-1-1:2023 Eq. (11.24)

    Args:
        Ftd (float): Deviation force in kN.
        r (float): Radius of curvature in mm.
        cu (float): Parameter cu in mm.
        fck (float): Characteristic compressive strength of concrete in MPa.
        gamma_c (float): Partial safety factor for concrete.

    Returns:
        bool: True if the deviation forces can be
            resisted by the concrete in tension, False otherwise.

    Raises:
        ValueError: If any input value is negative.
    """
    if Ftd < 0 or r < 0 or cu < 0 or fck < 0 or gamma_c <= 0:
        raise ValueError(
            'Invalid input values. Ftd, r, cu, and fck must not be '
            + 'negative, and gamma_c must be positive.'
        )

    return (Ftd * 1000 / (r * cu)) <= (0.125 / gamma_c) * math.sqrt(fck)


def verify_interaction_bond(
    Ftd: float,
    r: float,
    cu: float,
    fck: float,
    gamma_c: float,
    lsd: float,
    ls: float,
) -> bool:
    """Verify the interaction between bond and transverse
        tensile stresses due to deviation forces.

    EN1992-1-1:2023 Eq. (11.25)

    Args:
        Ftd (float): Deviation force in kN.
        r (float): Radius of curvature in mm.
        cu (float): Parameter cu in mm.
        fck (float): Characteristic compressive strength of concrete in MPa.
        gamma_c (float): Partial safety factor for concrete.
        lsd (float): Design lap length in mm.
        ls (float): Actual lap length in mm.

    Returns:
        bool: True if the interaction criterion is satisfied, False otherwise.

    Raises:
        ValueError: If any input value is negative.
    """
    if (
        Ftd < 0
        or r < 0
        or cu < 0
        or fck < 0
        or gamma_c <= 0
        or lsd < 0
        or ls < 0
    ):
        raise ValueError(
            'Invalid input values. Ftd, r, cu, fck, lsd, and ls must not be '
            + 'negative, and gamma_c must be positive.'
        )

    return (Ftd * 1000 * 8 * gamma_c) / (r * cu * math.sqrt(fck)) + (
        lsd / ls
    ) <= 1
