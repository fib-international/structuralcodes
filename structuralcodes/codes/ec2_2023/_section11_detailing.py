"""Functions from Section 11 of EN 1992-1-1:2023."""

import math
from typing import Literal

from scipy.interpolate import interp2d


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
    interpolation_function = interp2d(
        fck_values, phi_values, anchorage_matrix, kind='linear'
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
    lbd_phi = interpolation_function(fck, phi)[0]

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
        ValueError: If any of the input values are outside the specified limits.

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


def calculate_cd_conf(
    cx: float,
    cy: float,
    phi_t: float,
    st: float,
    cs: float,
    phi: float,
    rho_conf: float,
    sigma_ccd: float,
    f_ck: float,
) -> float:
    """Calculate the reduced design anchorage length parameter cd,conf.

    EN1992-1-1:2023 Eq. (11.4)

    Args:
        cx (float): Distance to the free edge in the x-direction in mm.
        cy (float): Distance to the free edge in the y-direction in mm.
        phi_t (float): Diameter of the transverse reinforcement in mm.
        st (float): Spacing of the transverse reinforcement in mm.
        cs (float): Spacing of the confinement reinforcement in mm.
        phi (float): Diameter of the bar to be anchored in mm.
        rho_conf (float): Ratio of the confinement reinforcement (dimensionless).
        sigma_ccd (float): Design value of the mean compression stress in MPa.
        f_ck (float): Characteristic compressive strength of concrete in MPa.

    Returns:
        float: The reduced design anchorage length parameter cd,conf in mm.

    Raises:
        ValueError: If any input is negative.

    References:
        EN1992-1-1:2023 Eq. (11.4)
    """
    # Validate inputs
    if any(
        param < 0
        for param in [cx, cy, phi_t, st, cs, phi, rho_conf, sigma_ccd, f_ck]
    ):
        raise ValueError('Inputs must not be negative.')

    # Calculate delta_cd
    delta_cd = (70 * rho_conf + 12 * sigma_ccd / f_ck) * phi

    # Calculate cd,conf
    cd_conf = min(cx, cy + 25 * (phi_t**2) / st, cs / 2, 3.75 * phi) + delta_cd

    # Limit cd_conf
    return min(cd_conf, 6 * phi)
