"""fib MC2020 Chapter 14.6.1."""

import math
from typing import List, Literal, Optional, Tuple


def fcm(fck: float, delta_f: float = 8) -> float:
    """Calculate the mean compressive strength fcm based on the
        characteristic compressive strength fck.

    fib Model Code 2020, eq. (14.6-1)

    Args:
        fck (float): Characteristic compressive strength fck in MPa.
        delta_f (float): increment in MPa. Defaults to 8MPa.

    Returns:
        float: Mean compressive strength fcm in MPa.

    Raises:
        ValueError: If fck is negative or zero.
    """
    if fck <= 0:
        raise ValueError(
            'Characteristic compressive strength fck must be positive.'
        )

    return fck + delta_f


def flcm(flck: float, delta_f: float = 8) -> float:
    """Calculate the mean compressive strength flcm for lightweight
        aggregate concrete based on the characteristic
        compressive strength flck.

    fib Model Code 2020, eq. (14.6-2)

    Args:
        flck (float): Characteristic compressive strength flck in MPa.

    Returns:
        float: Mean compressive strength flcm in MPa.
        delta_f (float): increment in MPa. Defaults to 8MPa.

    Raises:
        ValueError: If flck is negative or zero.
    """
    if flck <= 0:
        raise ValueError(
            'Characteristic compressive strength flck must be positive.'
        )

    return flck + delta_f


# Dictionary for normal concrete grades and
# their respective fck and fck,cube values
_concrete_strengths = {
    'C12': {'fck': 12, 'fck_cube': 15},
    'C16': {'fck': 16, 'fck_cube': 20},
    'C20': {'fck': 20, 'fck_cube': 25},
    'C25': {'fck': 25, 'fck_cube': 30},
    'C30': {'fck': 30, 'fck_cube': 37},
    'C35': {'fck': 35, 'fck_cube': 45},
    'C40': {'fck': 40, 'fck_cube': 50},
    'C45': {'fck': 45, 'fck_cube': 55},
    'C50': {'fck': 50, 'fck_cube': 60},
    'C55': {'fck': 55, 'fck_cube': 67},
    'C60': {'fck': 60, 'fck_cube': 75},
    'C70': {'fck': 70, 'fck_cube': 85},
    'C80': {'fck': 80, 'fck_cube': 95},
    'C90': {'fck': 90, 'fck_cube': 105},
    'C100': {'fck': 100, 'fck_cube': 115},
    'C110': {'fck': 110, 'fck_cube': 130},
    'C120': {'fck': 120, 'fck_cube': 140},
    'LC8': {'fck': 8, 'fck_cube': 9},
    'LC12': {'fck': 12, 'fck_cube': 13},
    'LC16': {'fck': 16, 'fck_cube': 18},
    'LC20': {'fck': 20, 'fck_cube': 22},
    'LC25': {'fck': 25, 'fck_cube': 28},
    'LC30': {'fck': 30, 'fck_cube': 33},
    'LC35': {'fck': 35, 'fck_cube': 38},
    'LC40': {'fck': 40, 'fck_cube': 44},
    'LC45': {'fck': 45, 'fck_cube': 50},
    'LC50': {'fck': 50, 'fck_cube': 55},
    'LC55': {'fck': 55, 'fck_cube': 60},
    'LC60': {'fck': 60, 'fck_cube': 66},
    'LC70': {'fck': 70, 'fck_cube': 77},
    'LC80': {'fck': 80, 'fck_cube': 88},
}


def fck(
    grade: Literal[
        'C12',
        'C16',
        'C20',
        'C25',
        'C30',
        'C35',
        'C40',
        'C45',
        'C50',
        'C55',
        'C60',
        'C70',
        'C80',
        'C90',
        'C100',
        'C110',
        'C120',
        'LC8',
        'LC12',
        'LC16',
        'LC20',
        'LC25',
        'LC30',
        'LC35',
        'LC40',
        'LC45',
        'LC50',
        'LC55',
        'LC60',
        'LC70',
        'LC80',
    ],
) -> float:
    """Retrieve the characteristic strength values (fck)
        for a given con.

    fib Model Code 2020, Tables 14.6-1, 14.6-2

    Args:
        grade (str): Concrete grade (e.g., 'C20', 'C30', 'LC12').

    Returns:
        dict: fck in MPa.

    Raises:
        ValueError: If the grade is not found.

    Reference:

    """
    if grade not in _concrete_strengths:
        raise ValueError(
            f'Invalid concrete grade: {grade}. '
            + f'Available grades are {list(_concrete_strengths.keys())}.'
        )

    return _concrete_strengths[grade]['fck']


def fck_cube(
    grade: Literal[
        'C12',
        'C16',
        'C20',
        'C25',
        'C30',
        'C35',
        'C40',
        'C45',
        'C50',
        'C55',
        'C60',
        'C70',
        'C80',
        'C90',
        'C100',
        'C110',
        'C120',
        'LC8',
        'LC12',
        'LC16',
        'LC20',
        'LC25',
        'LC30',
        'LC35',
        'LC40',
        'LC45',
        'LC50',
        'LC55',
        'LC60',
        'LC70',
        'LC80',
    ],
) -> float:
    """Retrieve the characteristic cube strength values (fck_cube)
        for a given con.

    fib Model Code 2020, Tables 14.6-1, 14.6-2

    Args:
        grade (str): Concrete grade (e.g., 'C20', 'C30', 'LC12').

    Returns:
        dict: fck_cube in MPa.

    Raises:
        ValueError: If the grade is not found.

    Reference:

    """
    if grade not in _concrete_strengths:
        raise ValueError(
            f'Invalid concrete grade: {grade}. '
            + f'Available grades are {list(_concrete_strengths.keys())}.'
        )

    return _concrete_strengths[grade]['fck_cube']


def fctm(fck: float) -> float:
    """Calculate the mean tensile strength fctm for normal weight concrete.

    fib Model Code 2020, eq. (14.6-3)

    Args:
        fck (float): Characteristic compressive strength fck in MPa.

    Returns:
        float: Mean tensile strength fctm in MPa.

    Raises:
        ValueError: If fck is negative or zero.
    """
    if fck <= 0:
        raise ValueError(
            'Characteristic compressive strength fck must be positive.'
        )

    return 1.8 * math.log(fck) - 3.1


def fctk_min(fctm: float) -> float:
    """Calculate the lower bound characteristic tensile strength fctk,min.

    fib Model Code 2020, eq. (14.6-4)

    Args:
        fctm (float): Mean tensile strength fctm in MPa.

    Returns:
        float: Lower bound characteristic tensile strength fctk,min in MPa.
    """
    return 0.7 * fctm


def fctk_max(fctm: float) -> float:
    """Calculate the upper bound characteristic tensile strength fctk,max.

    fib Model Code 2020, eq. (14.6-5)

    Args:
        fctm (float): Mean tensile strength fctm in MPa.

    Returns:
        float: Upper bound characteristic tensile strength fctk,max in MPa.
    """
    return 1.3 * fctm


def flctm(fctm: float, density: float) -> float:
    """Calculate the mean tensile strength flctm for
        lightweight aggregate concrete.

    fib Model Code 2020, eqs. (14.6-6a) and (14.6-6b)

    Args:
        fctm (float): Mean tensile strength for normal weight concrete in MPa.
        density (float): Oven-dry density of the lightweight
            aggregate concrete in kg/m3.

    Returns:
        float: Mean tensile strength flctm in MPa.

    Raises:
        ValueError: If density is outside the
            reasonable range for lightweight aggregate concrete.
    """
    if density < 1000 or density > 2200:
        raise ValueError(
            'Density must be between 1000 and 2200 kg/m³ '
            + 'for lightweight aggregate concrete.'
        )

    eta_l = 0.4 + 0.6 * (density / 2200)
    return eta_l * fctm


def flctk_min(flctm: float) -> float:
    """Calculate the lower bound characteristic
        tensile strength flctk,min for lightweight aggregate concrete.

    fib Model Code 2020, similar to eq. (14.6-4)

    Args:
        flctm (float): Mean tensile strength flctm in MPa.

    Returns:
        float: Lower bound characteristic tensile strength flctk,min in MPa.
    """
    return 0.7 * flctm


def flctk_max(flctm: float) -> float:
    """Calculate the upper bound characteristic tensile strength
        flctk,max for lightweight aggregate concrete.

    fib Model Code 2020, similar to eq. (14.6-5)

    Args:
        flctm (float): Mean tensile strength flctm in MPa.

    Returns:
        float: Upper bound characteristic tensile strength flctk,max in MPa.
    """
    return 1.3 * flctm


def fctm_from_splitting(fctm_sp: float, alpha_sp: float = 1.0) -> float:
    """Calculate the mean uniaxial tensile strength
        fctm from the mean splitting tensile strength fctm,sp.

    fib Model Code 2020, eq. (14.6-7)

    Args:
        fctm_sp (float): Mean splitting tensile strength fctm,sp in MPa.
        alpha_sp (float): Conversion factor (default is 1.0).

    Returns:
        float: Mean uniaxial tensile strength fctm iun MPa.
    """
    return alpha_sp * fctm_sp


def fctm_from_flexural(fctm_fl: float, hb: float) -> float:
    """Calculate the mean uniaxial tensile strength
        fctm from the mean flexural tensile strength fctm,fl.

    fib Model Code 2020, eqs. (14.6-8a) and (14.6-8b)

    Args:
        fctm_fl (float): Mean flexural tensile strength fctm,fl in MPa.
        hb (float): Beam depth in mm.

    Returns:
        float: Mean uniaxial tensile strength fctm in MPa.

    Raises:
        ValueError: If hb is non-positive.
    """
    if hb <= 0:
        raise ValueError(f'Beam depth hb must be positive. Got {hb}')

    alpha_fl = 0.06 * hb**0.7 / (1 + 0.06 * hb**0.7)
    return alpha_fl * fctm_fl


def GF(fck: float) -> float:
    """Calculate the fracture energy GF for normal weight concrete.

    fib Model Code 2020, eq. (14.6-9)

    Args:
        fck (float): Characteristic compressive strength fck in MPa.

    Returns:
        float: Fracture energy GF in N/m.

    Raises:
        ValueError: If fck is negative or zero.
    """
    if fck <= 0:
        raise ValueError(
            'Characteristic compressive strength fck must be positive.'
        )

    return 85 * fck**0.15


def GF_l(flctm: float, use_normal_weight_sand: bool = True) -> float:
    """Calculate the fracture energy GF,l for lightweight aggregate concrete.

    fib Model Code 2020, eq. (14.6-10)

    Args:
        flctm (float): Mean tensile strength flctm in MPa.
        use_normal_weight_sand (bool): True if using normal
            weight sand, False if using lightweight sand.

    Returns:
        float: Fracture energy GF,l in N/m.

    Raises:
        ValueError: If flctm is negative or zero.
    """
    if flctm <= 0:
        raise ValueError('Mean tensile strength flctm must be positive.')

    G_FoA = 24 if use_normal_weight_sand else 0
    return G_FoA + 16 * flctm


def multiaxial_stress_equation(
    sigma_1: float,
    sigma_2: float,
    sigma_3: float,
    fcm: float,
    fctm: float,
    fc2cm: float,
    sigma_com: Optional[float] = None,
    tau_com: Optional[float] = None,
    is_lightweight: bool = False,
) -> float:
    """Returns the symbolic equation for the mean
        strength under multiaxial states of stress. The result
        of this evaluation must be equal to zero to
        fulfill the failure criterion. So a custom
        solver defined by the user can be used.

    fib Model Code 2020, eq. (14.6-11 to 14.6-19.b)

    Args:
        sigma_1 (float): first principal stress in MPa.
        sigma_2 (float): second principal stress in MPa.
        sigma_3 (float): third principal stress in MPa.
        fcm (float): mean compressive strength of concrete in MPa. Use
            flcm for lightweight concrete.
        fctm (float): mean tensile strength of concrete in MPa. Use flctm
            for lightweight concrete.
        fc2cm (float): biaxial compressive strength in MPa.
        sigma_com (float): triaxial compressive strength at on point
            on the compressive meridian in MPa. (sigma_1 = sigma_2 > sigma_3).
            If None, it computes the recommended value.
        tau_com (float): triaxial comrpessive shear strength at one point
            on the compressive meridian in MPa. (sigma_1 = sigma_2 > sigma_3).
            If None, it computes the recommended value.
        is_lightweight (bool): True if lightweight concrete. False otherwise.

    Returns:
        float: the result of the evaluation of the multiaxial stress.
            When this result is equal to zero, the equilibrium is met.
            A custom solver/standard solver can be used passing the 3
            principal stresses (sigma_1, sigma_2, sigma_3).
    """
    # Compute the invariants
    I1 = sigma_1 + sigma_2 + sigma_3
    sigma_m = I1 / 3
    J2 = (1 / 6) * (
        (sigma_1 - sigma_2) ** 2
        + (sigma_2 - sigma_3) ** 2
        + (sigma_3 - sigma_1) ** 2
    )
    J3 = (sigma_1 - sigma_m) * (sigma_2 - sigma_m) * (sigma_3 - sigma_m)

    # Compute sigma_com and tau_com if not specified
    if sigma_com is None:
        sigma_com = -240 if not is_lightweight else -60
    if tau_com is None:
        if is_lightweight:
            tau_com = (
                185
                - 180 * (fcm / 100)
                + 260 * (fcm / 100) ** 2
                - 84 * (fcm / 100) ** 3
            )
        else:
            tau_com = (
                250 * (fcm / 100)
                - 460 * (fcm / 100) ** 2
                + 310 * (fcm / 100) ** 3
            )

    # Compute k, f2c, x, y, h, alpha and beta
    k = fctm / fcm
    f2c = fc2cm / fcm
    x = sigma_com / fcm
    y = tau_com / fcm
    h = -(math.sqrt(2) * x + y) / ((y / math.sqrt(2)) - (1 / 3))
    beta = (math.sqrt(2) - (3 * y) / (k * f2c)) / (h - (9 * y) / (f2c - k))
    alpha = (h * beta - math.sqrt(2)) / y

    # Compute lambda_c and lambda_t
    lambda_c = (
        (1 - h / (3 * y)) * math.sqrt(3) * beta
        + math.sqrt(3)
        + math.sqrt(2) / (math.sqrt(3) * y)
    )
    lambda_t = (
        2 * math.sqrt(3)
        - (f2c * h) / (math.sqrt(3) * y) * beta
        + math.sqrt(3) / f2c
        + (math.sqrt(2) * f2c) / (math.sqrt(3) * y)
    )
    lambda_r = lambda_c / lambda_t

    # Compute theta
    cos_3_theta = (3 * math.sqrt(3) / 2) * J3 / (J2 ** (3 / 2))
    theta = math.acos(cos_3_theta) / 3
    cos_theta = math.cos(theta)

    # Compute c1 and c2
    if lambda_r <= 1 / 2:
        c2 = 1
    else:
        c2 = math.cos(3 * math.atan((2 * lambda_r - 1) / (math.sqrt(3))))

    if lambda_r <= 1 / 2:
        c1 = (2 * cos_theta - 1) * lambda_t + 4 * (1 - cos_theta) * lambda_c
    else:
        c1 = lambda_c / (math.cos(math.pi / 3 - 1 / 3 * math.acos(c2)))

    # Compute lambda
    _lambda = c1 * math.cos(1 / 3 * math.acos(c2 * cos_3_theta))

    # Return the total result
    return (
        (alpha * J2 / fcm**2)
        + (_lambda * math.sqrt(J2) / fcm)
        + (beta * I1 / fcm)
        - 1
    )


def E_ci(fcm: float, alpha_E: float = 1.0, Ec0: float = 21500) -> float:
    """Calculate the modulus of elasticity Eci for normal
        weight concrete at 28 days.

    fib Model Code 2020, eq. (14.6-20)

    Args:
        fcm (float): Mean compressive strength fck in MPa.
        alpha_E (float): Coefficient depending on
            the aggregate type (default is 1.0 for quartzite).
        Ec0 (float): base elastic modulus value in MPa.

    Returns:
        float: Modulus of elasticity E_ci in MPa.

    Raises:
        ValueError: If fck or alpha_E is invalid.
    """
    if fcm <= 0:
        raise ValueError(
            'Characteristic compressive strength fck must be positive.'
        )
    if not (0.7 <= alpha_E <= 1.2):
        raise ValueError('Coefficient alpha_E must be between 0.7 and 1.2.')

    return alpha_E * Ec0 * (fcm / 10) ** (1 / 3)


def El_ci(
    fcm: float, density: float, alpha_E: float = 1.0, Ec0: float = 21500
) -> float:
    """Calculate the modulus of elasticity Elci for
        lightweight aggregate concrete.

    fib Model Code 2020, eq. (14.6-22)

    Args:
        fcm (float): Mean compressive strength fcm in MPa.
        density (float): Oven-dry density of the lightweight
            aggregate concrete in kg/m3.
        alpha_E (float): Coefficient depending on the aggregate
            type (default is 1.0 for quartzite).
        Ec0 (float): base elastic modulus value in MPa.

    Returns:
        float: Modulus of elasticity El_ci in MPa.

    Raises:
        ValueError: If fcm, density, or alpha_E is invalid.
    """
    if fcm <= 0:
        raise ValueError('Mean compressive strength fcm must be positive.')
    if density <= 0:
        raise ValueError('Density must be positive.')
    if not (0.7 <= alpha_E <= 1.2):
        raise ValueError('Coefficient alpha_E must be between 0.7 and 1.2.')

    eta_E = (density / 2200) ** 2

    Eci = alpha_E * Ec0 * (fcm / 10) ** (1 / 3)
    return eta_E * Eci


def E_ci_red_el(Eci: float, fcm: float) -> float:
    """Calculate the reduced modulus of elasticity Ec for concrete
        for only static analysis.

    fib Model Code 2020, eq. (14.6-23) (14.6-24)

    Args:
        Eci (float): Modulus of elasticity Eci in MPa.
        fcm (float): Mean compressive strength fcm in MPa.

    Returns:
        float: Reduced modulus of elasticity Ec in MPa.

    Raises:
        ValueError: If Eci or fcm is negative or zero.
    """
    if Eci < 0:
        raise ValueError('Modulus of elasticity Eci must be positive.')

    alpha_i = min(0.8 + 0.2 * fcm / 88, 1)
    return alpha_i * Eci


def nu_c() -> float:
    """Get the Poisson's ratio νc for concrete for design purposes.

    fib Model Code 2020 (14.6.1.4.3)

    Returns:
        float: Poisson's ratio νc for concrete,
            typically 0.20 for design purposes.
    """
    return 0.20


def sigma_c(
    eps_c: float,
    eps_c1: float,
    eps_lim: float,
    Eci: float,
    Ec1: float,
    fcm: float,
    k: float,
) -> float:
    """Calculate the compressive stress σc for short-term
        uniaxial compression of concrete.

    fib Model Code 2020, eq. (14.6-26)

    Args:
        eps_c (float): Compressive strain εc.
        eps_c1 (float): Strain at maximum compressive stress εc1.
        eps_lim (float): Limit strain.
        Eci (float): Modulus of elasticity of concrete in MPa.
        Ec1 (float): Secant modulus of elasticity from origin
            to peak stress Ec1 in MPa.
        fcm (float): Mean compressive strength fcm in MPa.
        k (float): Plasticity number k.

    Returns:
        float: Compressive stress σc in MPa.

    Raises:
        ValueError: If input values are not positive or within valid ranges.
    """
    if eps_c < 0 or eps_c1 < 0 or Ec1 < 0 or fcm < 0 or k < 0:
        raise ValueError('All input values must be positive.')
    if eps_c > eps_lim:
        raise ValueError('eps_c must be equal or less than eps_lim')

    # Calculate η (eta)
    eta = eps_c / eps_c1
    k = Eci / Ec1

    # Calculate compressive stress σc using the given formula
    return -(k * eta - eta**2) / (1 + (k - 2) * eta)


# Dictionary for concrete properties based on grade
concrete_properties = {
    'C12': {
        'Eci': 27.1,
        'Ec1': 11.1,
        'epsilon_c1': 1.9,
        'epsilon_c_lim': 3.5,
        'k': 2.44,
    },
    'C16': {
        'Eci': 28.8,
        'Ec1': 12.2,
        'epsilon_c1': 2.0,
        'epsilon_c_lim': 3.5,
        'k': 2.36,
    },
    'C20': {
        'Eci': 30.3,
        'Ec1': 13.3,
        'epsilon_c1': 2.1,
        'epsilon_c_lim': 3.5,
        'k': 2.28,
    },
    'C25': {
        'Eci': 32.0,
        'Ec1': 14.9,
        'epsilon_c1': 2.2,
        'epsilon_c_lim': 3.5,
        'k': 2.15,
    },
    'C30': {
        'Eci': 33.6,
        'Ec1': 16.5,
        'epsilon_c1': 2.3,
        'epsilon_c_lim': 3.5,
        'k': 2.04,
    },
    'C35': {
        'Eci': 35.0,
        'Ec1': 18.2,
        'epsilon_c1': 2.3,
        'epsilon_c_lim': 3.5,
        'k': 1.92,
    },
    'C40': {
        'Eci': 36.3,
        'Ec1': 20.0,
        'epsilon_c1': 2.4,
        'epsilon_c_lim': 3.5,
        'k': 1.82,
    },
    'C45': {
        'Eci': 37.5,
        'Ec1': 21.6,
        'epsilon_c1': 2.5,
        'epsilon_c_lim': 3.5,
        'k': 1.74,
    },
    'C50': {
        'Eci': 38.6,
        'Ec1': 23.2,
        'epsilon_c1': 2.6,
        'epsilon_c_lim': 3.4,
        'k': 1.66,
    },
    'C55': {
        'Eci': 39.7,
        'Ec1': 24.7,
        'epsilon_c1': 2.6,
        'epsilon_c_lim': 3.4,
        'k': 1.61,
    },
    'C60': {
        'Eci': 40.7,
        'Ec1': 26.2,
        'epsilon_c1': 2.7,
        'epsilon_c_lim': 3.3,
        'k': 1.55,
    },
    'C70': {
        'Eci': 42.6,
        'Ec1': 28.9,
        'epsilon_c1': 2.7,
        'epsilon_c_lim': 3.2,
        'k': 1.47,
    },
    'C80': {
        'Eci': 44.4,
        'Ec1': 31.4,
        'epsilon_c1': 2.8,
        'epsilon_c_lim': 3.1,
        'k': 1.41,
    },
    'C90': {
        'Eci': 46.0,
        'Ec1': 33.8,
        'epsilon_c1': 2.9,
        'epsilon_c_lim': 3.0,
        'k': 1.36,
    },
    'C100': {
        'Eci': 47.5,
        'Ec1': 36.0,
        'epsilon_c1': 3.0,
        'epsilon_c_lim': 3.0,
        'k': 1.32,
    },
    'C110': {
        'Eci': 48.9,
        'Ec1': 39.3,
        'epsilon_c1': 3.0,
        'epsilon_c_lim': 3.0,
        'k': 1.24,
    },
    'C120': {
        'Eci': 50.3,
        'Ec1': 42.7,
        'epsilon_c1': 3.0,
        'epsilon_c_lim': 3.0,
        'k': 1.18,
    },
}


def Ec1(concrete_grade: str) -> float:
    """Get the secant modulus of elasticity Ec1 for a given concrete grade.

    Args:
        concrete_grade (str): Concrete grade (e.g., 'C20', 'C30').

    Returns:
        float: Secant modulus of elasticity Ec1 in MPa.

    Raises:
        ValueError: If the concrete grade is not recognized.
    """
    if concrete_grade not in concrete_properties:
        raise ValueError(
            f'Invalid concrete grade: {concrete_grade}. '
            + f'Available grades are {list(concrete_properties.keys())}.'
        )

    return concrete_properties[concrete_grade]['Ec1'] * 1000


def eps_c1(concrete_grade: str) -> float:
    """Get the strain at maximum compressive stress
        εc1 for a given concrete grade.

    Args:
        concrete_grade (str): Concrete grade (e.g., 'C20', 'C30').

    Returns:
        float: Strain at maximum compressive stress εc1.

    Raises:
        ValueError: If the concrete grade is not recognized.
    """
    if concrete_grade not in concrete_properties:
        raise ValueError(
            f'Invalid concrete grade: {concrete_grade}. '
            + f'Available grades are {list(concrete_properties.keys())}.'
        )

    return concrete_properties[concrete_grade]['epsilon_c1'] / 1000


def eps_c_lim(concrete_grade: str) -> float:
    """Get the strain limit εc,lim for a given concrete grade.

    Args:
        concrete_grade (str): Concrete grade (e.g., 'C20', 'C30').

    Returns:
        float: Strain limit εc,lim.

    Raises:
        ValueError: If the concrete grade is not recognized.
    """
    if concrete_grade not in concrete_properties:
        raise ValueError(
            f'Invalid concrete grade: {concrete_grade}. '
            + f'Available grades are {list(concrete_properties.keys())}.'
        )

    return concrete_properties[concrete_grade]['epsilon_c_lim'] / 1000


def eps_lc1(
    flck: float, Elc: float, sand: Literal['light', 'natural']
) -> float:
    """Calculate the strain εlc1 for lightweight aggregate concrete.

    fib Model Code 2020, eq. (14.6-27)

    Args:
        flck (float): Characteristic compressive strength flck in MPa.
        E_lc (float): Modulus of elasticity E_lc in MPa.
        sand (float): Specifies the type of sand used.

    Returns:
        float: Strain εlc1.

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if flck < 0 or Elc <= 0:
        raise ValueError(
            'Characteristic compressive strength f_lck '
            + 'and modulus of elasticity E_lc must be positive.'
        )

    kappa_lc = 1.1 if sand == 'light' else 1.3
    return kappa_lc * (flck + 8) / Elc


def sigma_ct_uncracked(eps_ct: float, E_ci: float, fctm: float) -> float:
    """Calculate the tensile stress σct for uncracked
        normal weight concrete subjected to tension.

    fib Model Code 2020, eqs. (14.6-29) and (14.6-30)

    Args:
        eps_ct (float): Tensile strain εct.
        E_ci (float): Tangent modulus of elasticity Eci in MPa.
        fctm (float): Mean tensile strength fctm in MPa.

    Returns:
        float: Tensile stress σct in MPa.

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if E_ci <= 0 or fctm < 0:
        raise ValueError(
            'Modulus of elasticity Eci and tensile '
            + 'strength fctm must be positive.'
        )

    sigma_ct_linear = E_ci * eps_ct

    if sigma_ct_linear <= 0.9 * fctm:
        return sigma_ct_linear
    if sigma_ct_linear <= fctm:
        return fctm * (
            1 - 0.1 * (0.00015 - eps_ct) / (0.00015 - 0.9 * fctm / E_ci)
        )
    raise ValueError('Stress is larger than mean concrete tensile stress.')


def sigma_ct_cracked(w: float, GF: float, fctm: float) -> float:
    """Calculate the tensile stress σct for cracked normal
        weight concrete subjected to tension.

    fib Model Code 2020, eqs. (14.6-31) and (14.6-32)

    Args:
        w (float): Crack opening w in mm.
        GF (float): Fracture energy GF in N/mm.
        fctm (float): Mean tensile strength fctm in MPa.

    Returns:
        float: Tensile stress σct in MPa.

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if GF <= 0 or fctm <= 0:
        raise ValueError(
            'Fracture energy GF and tensile strength fctm must be positive.'
        )

    w1 = GF / fctm
    wc = 5 * GF / fctm

    if w <= w1:
        return fctm * (1 - 0.8 * w / w1)
    if w <= wc:
        return fctm * (0.25 - 0.05 * w / w1)

    raise ValueError('w cannot be larger than wc')


# Chapter 14.6.1.5.3 Multiaxial states of stress skipped


def tau_crack_friction(
    delta: float, w: float, fcm: float, Cf: float = 1.0
) -> float:
    """Calculate the mean shear stress τ in an open
        crack subjected to shear displacements.

    fib Model Code 2020, eq. (14.6-43)

    Args:
        delta (float): Shear displacement δ in mm.
        w (float): Crack width w in mm.
        fcm (float): Mean compressive strength fcm in MPa.
        Cf (float): Constant C (default is 1.0).

    Returns:
        float: Mean shear stress τ in MPa.

    Raises:
        ValueError: If input values are not within valid ranges.

    """
    if delta <= 0 or w <= 0 or fcm <= 0:
        raise ValueError(
            'Shear displacement δ, crack width w, '
            + 'and mean compressive strength fcm must be positive.'
        )

    return Cf * (
        -0.04 * fcm + (1.8 * w**-0.8 + (0.292 * w**-0.7 - 0.25) * fcm) * delta
    )


def sigma_crack_friction(
    delta: float, w: float, fcm: float, Cf: float = 1.0
) -> float:
    """Calculate the mean normal stress σ in an open
        crack subjected to shear displacements.

    fib Model Code 2020, eq. (14.6-44)

    Args:
        delta (float): Shear displacement δ in mm.
        w (float): Crack width w in mm.
        fcm (float): Mean compressive strength fcm in MPa.
        Cf (float): Constant C (default is 1.0).

    Returns:
        float: Mean normal stress σ in MPa.

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if delta <= 0 or w <= 0 or fcm <= 0:
        raise ValueError(
            'Shear displacement δ, crack width w, and '
            + 'mean compressive strength fcm must be positive.'
        )

    return Cf * (
        -0.06 * fcm
        + (1.35 * w**-0.63 + (0.242 * w**-0.55 - 0.19) * fcm) * delta
    )


SC_dict = {'CR': [0.3, 0.2, 0.1], 'CN': [0.5, 0.4, 0.3], 'CS': [0.6, 0.5, 0.4]}


def sC(fck: float, strength_class: Literal['CR', 'CN', 'CS']) -> float:
    """Get the coefficient sC based on concrete strength class and fck.

    fib Model Code 2020, Table (14.6-7)

    Args:
        f_ck (float): Characteristic compressive strength fck in MPa.
        strength_class (str): Strength development class ('CR', 'CN', or 'CS').

    Returns:
        float: Coefficient sC.

    Raises:
        ValueError: If the input values are not valid.
    """
    if fck <= 0:
        raise ValueError(
            'Characteristic compressive strength fck must be positive.'
        )

    if strength_class not in ['CR', 'CN', 'CS']:
        raise ValueError("Strength class must be 'CR', 'CN', or 'CS'.")

    t = SC_dict[strength_class]
    if fck <= 35:
        return t[0]
    if fck < 60:
        return t[1]
    return t[2]


def beta_cc(t: float, t_ref: float, sC: float) -> float:
    """Calculate the strength development function βcc(t) for concrete.

    fib Model Code 2020, eq. (14.6-46)


    Args:
        t (float): Concrete age t in days.
        t_ref (float): Reference age tref in days.
        sC (float): Coefficient sC based on strength development class.

    Returns:
        float: Strength development function βcc(t).

    Raises:
        ValueError: If the input values are not valid.
    """
    if t <= 0 or t_ref <= 0:
        raise ValueError(
            'Concrete age t and reference age tref must be positive.'
        )

    return math.exp(sC * (1 - (t_ref / t) ** 0.5) * ((28 / t_ref) ** 0.5))


def fcm_t(fcm: float, beta_cc: float) -> float:
    """Calculate the mean compressive strength fcm(t) at a given age t.

    fib Model Code 2020, eq. (14.6-45)

    Args:
        fcm (float): Mean compressive strength at reference age fcm in MPa.
        beta_cc (float): Strength development function βcc(t).

    Returns:
        float: Mean compressive strength fcm(t) in MPa.

    Raises:
        ValueError: If the input values are not valid.
    """
    if fcm <= 0 or beta_cc <= 0:
        raise ValueError('Mean compressive strength fcm must be positive.')

    return beta_cc * fcm


def slC(aggregate_strength: Literal['high', 'low']) -> float:
    """Get the coefficient slC based on the strength of lightweight aggregate.

    fib Model Code 2020, eq. (14.6-47)

    Args:
        aggregate_strength (str): Strength of
            lightweight aggregate ('high' or 'low').

    Returns:
        float: Coefficient slC.

    Raises:
        ValueError: If the input values are not valid.
    """
    if aggregate_strength == 'high':
        return 0.05
    # if aggregate_strength == 'low':
    return 0.25


def beta_lcc(t: float, t_ref: float, slC: float) -> float:
    """Calculate the strength development function βlcc(t)
        for lightweight aggregate concrete.

    fib Model Code 2020, eq. (14.6-47)

    Args:
        t (float): Concrete age t in days.
        t_ref (float): Reference age tref ins days (typically 28 days).
        slC (float): Coefficient slc based on aggregate strength.

    Returns:
        float: Strength development function βlcc(t).

    Raises:
        ValueError: If the input values are not valid.
    """
    if t <= 0 or t_ref <= 0:
        raise ValueError(
            'Concrete age t and reference age tref must be positive.'
        )

    return math.exp(slC * (1 - (t_ref / t) ** 0.5) * ((28 / t_ref) ** 0.5))


def flcm_t(f_lcm: float, beta_lcc: float) -> float:
    """Calculate the mean compressive strength flcm(t) at
        a given age t for lightweight aggregate concrete.

    fib Model Code 2020, eq. (14.6-47)

    Args:
        f_lcm (float): Mean compressive strength at reference age flcm in MPa.
        beta_lcc (float): Strength development function βlcc(t).

    Returns:
        float: Mean compressive strength flcm(t) in MPa.

    Raises:
        ValueError: If the input values are not valid.
    """
    if f_lcm <= 0 or beta_lcc <= 0:
        raise ValueError('Mean compressive strength flcm must be positive.')

    return beta_lcc * f_lcm


def beta_E(beta_cc: float) -> float:
    """Calculate the modulus of elasticity development
        coefficient βE(t) for concrete.

    fib Model Code 2020, eq. (14.6-49)

    Args:
        beta_cc (float)

    Returns:
        float: Modulus of elasticity development coefficient βE(t).

    Raises:
        ValueError: If the input values are not valid.
    """
    return beta_cc**0.33


def E_ci_t(E_ci: float, beta_E: float) -> float:
    """Calculate the modulus of elasticity Eci(t) at a given age t.

    fib Model Code 2020, eq. (14.6-48)

    Args:
        E_ci (float): Modulus of elasticity at reference age Eci in MPa.
        beta_E (float): Modulus of elasticity development coefficient βE(t).

    Returns:
        float: Modulus of elasticity Eci(t) in MPa.

    Raises:
        ValueError: If the input values are not valid.
    """
    if E_ci <= 0:
        raise ValueError('Modulus of elasticity Eci must be positive.')

    return beta_E * E_ci


def beta_t0(t0: float) -> float:
    """Calculate the parameter βt0(t0) which
        considers the maturity of the concrete at loading.

    fib Model Code 2020, eq. (14.6-51)

    Args:
        t0 (float): Age of concrete at loading t0 in days.

    Returns:
        float: Parameter βt0(t0).

    Raises:
        ValueError: If the input values are not valid.
    """
    if t0 < 7:
        raise ValueError(
            'Age of concrete at loading t0 must be at least 7 days.'
        )

    return 0.64 + 0.01 * math.log(t0)


def beta_c_sus(t: float, t0: float, beta_t0: float) -> float:
    """Calculate the function βc,sus(t, t0) describing the
        decrease of strength with time under high sustained load.

    fib Model Code 2020, eq. (14.6-51)

    Args:
        t (float): Concrete age t in days.
        t0 (float): Age of concrete at loading t0 in days.
        beta_t0 (float): Value of evaluated function beta_T0.

    Returns:
        float: Function βc,sus(t, t0).

    Raises:
        ValueError: If the input values are not valid.
    """
    if t <= t0:
        raise ValueError(
            'Concrete age t must be greater than age at loading t0.'
        )

    return beta_t0 + (1 - beta_t0) * (1 + 10000 * (t - t0) / t0) ** -0.1


def fcm_sus(fcm: float, beta_cc_t: float, beta_c_sus: float) -> float:
    """Calculate the mean compressive strength
        fcm,sus(t, t0) under sustained load.

    fib Model Code 2020, eq. (14.6-50)

    Args:
        f_cm (float): Mean compressive strength at reference age fcm in MPa.
        beta_cc_t (float): Strength development function βcc(t).
        beta_c_sus_t_t0 (float): Function βc,sus(t, t0) describing
            the decrease of strength under high sustained load.

    Returns:
        float: Mean compressive strength fcm,sus(t, t0) in MPa.

    Raises:
        ValueError: If the input values are not valid.
    """
    if fcm <= 0:
        raise ValueError('Mean compressive strength fcm must be positive.')

    return beta_cc_t * beta_c_sus * fcm


def fctk_sus(
    fctk: float, concrete: Literal['normal', 'high-strength']
) -> float:
    """Calculate the sustained tensile strength fctk,sus
        for concrete under sustained loading.

    fib Model Code 2020, eq. (14.6-52)

    Args:
        fctk (float): Short-term tensile strength fctk in MPa.
        alpha (float): Reduction factor for sustained loading
            (default is 0.60).

    Returns:
        float: Sustained tensile strength fctk,sus in MPa.

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if fctk <= 0:
        raise ValueError('Short-term tensile strength f_ctk must be positive.')

    alpha = 0.60 if concrete == 'normal' else 0.75
    return alpha * fctk


def eps_cc(sigma_c_t0: float, phi_t_t0: float, E_ci: float) -> float:
    """Calculate the creep strain εcc(t, t0) at time
        t for a concrete member under sustained loading.

    fib Model Code 2020, eq. (14.6-55)

    Args:
        sigma_c_t0 (float): Constant stress applied at time t0, σc(t0) in MPa.
        phi_t_t0 (float): Creep coefficient φ(t, t0).
        E_ci (float): Modulus of elasticity at the age of 28 days, Eci in MPa.

    Returns:
        float: Creep strain εcc(t, t0).

    Raises:
        ValueError: If input values are not within valid ranges
    """
    if E_ci <= 0:
        raise ValueError('Modulus of elasticity Eci must be positive.')

    return (sigma_c_t0 / E_ci) * phi_t_t0


def eps_c_sigma(
    sigma_c_t0: float, E_ci_t0: float, E_ci: float, phi_t_t0: float
) -> float:
    """Calculate the stress-dependent strain εcσ(t, t0) at time t
        for a concrete member under sustained loading.

    fib Model Code 2020, eq. (14.6-56)

    Args:
        sigma_c_t0 (float): value of the stress at time t0 in MPa.
        E_ci_t0 (float): modulus of elasticity at time t0 in MPa.
        E_ci (float): modulus of elasticity of concrete in MPa.
        phi_t_t0 (float): Creep coefficient φ(t, t0).

    Returns:
        float: Stress-dependent strain εcσ(t, t0).

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if E_ci_t0 <= 0 or sigma_c_t0 < 0 or E_ci <= 0 or phi_t_t0 < 0:
        raise ValueError('Modulus of elasticity Eci(t0) must be positive.')

    return sigma_c_t0 * (1 / E_ci_t0 + phi_t_t0 / E_ci)


def phi_t_t0(phi_bc_t_t0: float, phi_dc_t_t0: float) -> float:
    """Calculate the total creep coefficient φ(t, t0)
        for concrete under sustained loading.

    fib Model Code 2020, eq. (14.6-58)

    Args:
        phi_bc_t_t0 (float): basic creep coefficient.
        phi_dc_t_t0 (float): drying creep coefficient.

    Returns:
        float: Total creep coefficient φ(t, t0).

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    return phi_bc_t_t0 + phi_dc_t_t0


def beta_bc_fcm(fcm: float) -> float:
    """Calculate the mean compressive strength factor
        (fcm) for drying creep coefficient.

    fib Model Code 2020, eq. (14.6-60)

    Args:
        fcm (float): Mean compressive strength at an age of 28 days fcm in MPa.

    Returns:
        float: Mean compressive strength factor β_bc_fcm(f_cm).

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if fcm <= 0:
        raise ValueError('Mean compressive strength fcm must be positive.')

    return 1.8 / fcm**0.7


def beta_bc_t_t0(t: float, t0: float, t0_adj: float) -> float:
    """Calculte the beta_bc coefficient evaluated for
        t,t0.

    fib Model Code 2020, eq. (14.6-61)

    Args:
        t (float): Age of concrete at the moment considered in days.
        t0 (float): Age of concrete at loading in days.
        t0_adj (float): Adjusted age at loading t0,adj in dats.

    Returns:
        float: the value of beta_bc.

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if t <= t0_adj:
        raise ValueError(
            'Age of concrete t must be greater than'
            + ' adjusted age at loading t0_adj.'
        )
    if t0_adj <= 0:
        raise ValueError(
            'Adjusted age at loading t0_adj and mean'
            + ' compressive strength f_cm must be positive.'
        )
    return math.log((30 / t0_adj + 0.0035) ** 2 * (t - t0) + 1)


def phi_bc_t_t0(beta_bc_t_t0: float, beta_bc_fcm: float) -> float:
    """Calculate the basic creep coefficient φbc(t, t0)
        for concrete under sustained loading.

    fib Model Code 2020, eqs. (14.6-59)

    Args:
        beta_bc_t_t0 (float): Value of beta_bc_t_t0.
        beta_bc_fcm (float): value of beta_bc_fcm.

    Returns:
        float: Basic creep coefficient φbc(t, t0).

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if beta_bc_t_t0 <= 0 or beta_bc_fcm <= 0:
        raise ValueError('beta_bc_t_t0 and beta_bc_fcm should be positive.')

    return beta_bc_fcm * beta_bc_t_t0


def beta_dc_fcm(fcm: float) -> float:
    """Calculate the mean compressive strength factor
        βf_cm(fcm) for drying creep coefficient.

    fib Model Code 2020, eq. (14.6-63)

    Args:
        fcm (float): Mean compressive strength at an age of 28 days fcm in MPa.

    Returns:
        float: Mean compressive strength factor βf_cm(f_cm).

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if fcm <= 0:
        raise ValueError('Mean compressive strength f_cm must be positive.')

    return 412 / fcm**1.4


def beta_RH_c(RH: float, h: float) -> float:
    """Calculate the relative humidity factor βRH(RH)
        for drying creep coefficient.

    fib Model Code 2020, eq. (14.6-64)

    Args:
        RH (float): Relative humidity of the ambient environment in %.
        h (float): Notional size of the member h in mm.

    Returns:
        float: Relative humidity factor βRH(RH).

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if RH < 0 or RH > 100:
        raise ValueError('Relative humidity RH must be between 0 and 100.')
    if h <= 0:
        raise ValueError('Notional size h must be positive.')

    return (1 - RH / 100) / (0.1 * h / 100) ** (1 / 3)


def beta_dc_t0(t0_adj: float) -> float:
    """Calculate the time factor βt0(t0) for drying creep coefficient.

    fib Model Code 2020, eq. (14.6-65)

    Args:
        t0_adj (float): Adjusted age at loading t0,adj in days.

    Returns:
        float: Time factor βt0(t0).

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if t0_adj <= 0:
        raise ValueError('Adjusted age at loading t0_adj must be positive.')

    return 1 / (0.1 + t0_adj**0.2)


def gamma_t0(t0_adj: float) -> float:
    """Calculate the factor γ for drying creep development with time.

    fib Model Code 2020, eq. (14.6-66b)

    Args:
        t0_adj (float): Adjusted age at loading t0,adj in days.

    Returns:
        float: Factor γ.

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if t0_adj <= 0:
        raise ValueError('t0_adj must be positive.')

    return 1 / (2.3 + 3.5 / math.sqrt(t0_adj))


def beta_h(h: float, fcm: float) -> float:
    """Calculate the factor beta_h for drying creep development with time.

    fib Model Code 2020, eqs. (14.6-66c), (14.6-66d)

    Args:
        h (float): Notional size of the member h in mm.
        fcm (float): Mean compressive strength at an age of 28 days fcm in MPa.

    Returns:
        float: Factor beta_h.

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if h <= 0 or fcm <= 0:
        raise ValueError(
            'Notional size h and mean compressive'
            + ' strength fcm must be positive.'
        )

    alpha_fcm = (35 / fcm) ** 0.5
    return min(1.5 * h + 250 * alpha_fcm, 1500 * alpha_fcm)


def beta_dc_t_t0(t: float, t0: float, gamma_t0: float, beta_h: float) -> float:
    """Calculate the development of drying creep with time βdc(t, t0).

    fib Model Code 2020, eq. (14.6-66a)

    Args:
        t (float): Age of concrete at the moment considered in days.
        t0 (float): Age of concrete at loading t0 in days.
        h (float): Notional size of the member h in mm.
        gamma_t0 (float): gamma_t0 coefficient.
        beta_h (float): beta_h coefficient.


    Returns:
        float: Development of drying creep with time βdc(t, t0).

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if t <= t0:
        raise ValueError(
            'Age of concrete t must be greater than age at loading t0.'
        )
    if gamma_t0 <= 0:
        raise ValueError('Notional size h must be positive.')

    return ((t - t0) / (beta_h + (t - t0))) ** gamma_t0


def phi_dc_t_t0(
    beta_dc_fm: float,
    beta_RH: float,
    beta_dc_t0: float,
    beta_dc_t_t0: float,
) -> float:
    """Calculate the drying creep coefficient φdc(t, t0)
        for concrete under sustained loading.

    fib Model Code 2020, eq. (14.6-62)

    Args:
        beta_dc_fm (float): beta_dc_fm coefficient.
        beta_RH (float): beta_RH coefficient.
        beta_dc_t0 (float): beta_dc_t0 coefficient.
        beta_dc_t_t0 (float): beta_dc_t_t0 coefficient.

    Returns:
        float: Drying creep coefficient φdc(t, t0).

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    return beta_dc_fm * beta_RH * beta_dc_t0 * beta_dc_t_t0


def phi_l_t_t0(phi_t_t0: float, rho: float, concrete_grade: str) -> float:
    """Calculate the creep coefficient φl for lightweight aggregate concrete.

    fib Model Code 2020, eq. (14.6-67)

    Args:
        phi_t_t0 (float): Creep coefficient φ(t, t0)
            according to Eq. (14.6-58).
        rho (float): Oven-dry density of
            lightweight aggregate concrete in kg/m3.
        concrete_grade (str): Grade of concrete (e.g., 'LC12', 'LC16', 'LC20').

    Returns:
        float: Creep coefficient φl for lightweight aggregate concrete.

    Raises:
        ValueError: If input values are not within
            valid ranges or invalid grade.
    """
    if rho <= 0:
        raise ValueError('Oven-dry density ρ must be positive.')
    if phi_t_t0 < 0:
        raise ValueError('Creep coefficient φ(t, t0) must not be negative.')

    # Calculate ηE
    eta_E = (rho / 2200) ** 2

    # Calculate the initial creep coefficient for lightweight concrete
    phi_l = eta_E * phi_t_t0

    # Adjust for concrete grades LC12 and LC16
    if concrete_grade in ['LC12', 'LC16']:
        phi_l *= 1.3

    return phi_l


def t0_adj(t0_T: float, strength_class: Literal['CS', 'CN', 'CR']) -> float:
    """Calculate the adjusted age at loading t0,adj for
        concrete under sustained loading.

    fib Model Code 2020, eq. (14.6-68)

    Args:
        t0_T (float): Age of concrete at loading in
            days adjusted according to Eq. (14.6-80).
        strength_class (str): Strength development
            class of concrete ('CS', 'CN', 'CR').

    Returns:
        float: Adjusted age at loading t0,adj in days.

    Raises:
        ValueError: If input values are not
            within valid ranges or invalid class.
    """
    if t0_T <= 0:
        raise ValueError('Age of concrete at loading t0_T must be positive.')

    if strength_class == 'CS':
        alpha = -1
    if strength_class == 'CN':
        alpha = 0
    if strength_class == 'CR':
        alpha = 1

    return max(0.5, t0_T * (9 / (2 + t0_T**1.2) + 1) ** alpha)


def phi_sigma_t_t0(phi_t_t0: float, sigma_c: float, fcm_t0: float) -> float:
    """Calculate the non-linear creep coefficient φσ(t, t0)
        for concrete under high stress.

    fib Model Code 2020, eq. (14.6-69)

    Args:
        phi_t_t0 (float): Linear creep coefficient φ(t, t0)
            according to Eq. (14.6-58).
        sigma_c (float): Applied stress in MPa.
        fcm_t0 (float): Mean compressive strength of concrete
            at the age of loading t0 in MPa.

    Returns:
        float: Non-linear creep coefficient φσ(t, t0).

    Raises:
        ValueError: If input values are not within valid ranges or
            if stress-strength ratio is not within 0.4 < kσ ≤ 0.6.
    """
    if phi_t_t0 < 0:
        raise ValueError(
            'Linear creep coefficient φ(t, t0) must not be negative.'
        )
    if fcm_t0 <= 0:
        raise ValueError(
            'Mean compressive strength f_cm(t0) must be positive.'
        )

    # Calculate the stress-strength ratio
    k_sigma = abs(sigma_c) / fcm_t0

    # Calculate the non-linear creep coefficient
    return phi_t_t0 * math.exp(1.50 * (k_sigma - 0.4))


def eps_cs_t_ts(eps_cbs_t: float, eps_cds_t_ts: float) -> float:
    """Calculate the total shrinkage strain εcs(t, ts).

    fib Model Code 2020, eq. (14.6-70)

    Args:
        eps_cbs_t (float): Basic shrinkage strain εcbs(t).
        eps_cds_t_ts (float): Drying shrinkage strain εcds(t, ts).

    Returns:
        float: Total shrinkage strain εcs(t, ts).

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if eps_cbs_t < 0:
        raise ValueError(
            'Basic shrinkage strain εcbs(t) must not be negative.'
        )
    if eps_cds_t_ts < 0:
        raise ValueError(
            'Drying shrinkage strain εcds(t, ts) must not be negative.'
        )

    return eps_cbs_t + eps_cds_t_ts


def eps_cbs_t(eps_cbs0: float, beta_bs_t: float) -> float:
    """Calculate the basic shrinkage strain εcbs(t).

    fib Model Code 2020, eq. (14.6-71)

    Args:
        eps_cbs0 (float): Basic shrinkage strain coefficient εcbs0.
        beta_bs_t (float): Development function for basic shrinkage βbs(t).

    Returns:
        float: Basic shrinkage strain εcbs(t).

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if eps_cbs0 < 0:
        raise ValueError(
            'Basic shrinkage strain coefficient εcbs0 must not be negative.'
        )
    if beta_bs_t < 0:
        raise ValueError('Development function βbs(t) must not be negative.')

    return eps_cbs0 * beta_bs_t


def eps_cds_t_ts(
    eps_cds0: float, beta_RH: float, beta_ds_t_ts: float
) -> float:
    """Calculate the drying shrinkage strain εcds(t, ts).

    fib Model Code 2020, eq. (14.6-72)

    Args:
        eps_cds0 (float): Drying shrinkage strain coefficient εcds0.
        beta_RH (float): Relative humidity factor βRH(RH).
        beta_ds_t_ts (float): Development function for
            drying shrinkage βds(t - ts).

    Returns:
        float: Drying shrinkage strain εcds(t, ts).

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if eps_cds0 < 0:
        raise ValueError(
            'Drying shrinkage strain coefficient εcds0 must not be negative.'
        )
    if beta_RH < 0 or beta_RH > 1:
        raise ValueError(
            'Relative humidity factor βRH must be between 0 and 1.'
        )
    if beta_ds_t_ts < 0:
        raise ValueError(
            'Development function βds(t - ts) must not be negative.'
        )

    return eps_cds0 * beta_RH * beta_ds_t_ts


def eps_cbs_0(fcm: float, alpha_bs: float) -> float:
    """Calculate the basic notional shrinkage coefficient εcbs0(fcm).

    fib Model Code 2020, eq. (14.6-73)

    Args:
        fcm (float): Mean compressive strength at the age of
            28 days fcm in MPa.
        alpha_bs (float): Coefficient dependent on
            the strength development class of concrete.

    Returns:
        float: Basic notional shrinkage coefficient εcbs0(fcm).

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if fcm <= 0:
        raise ValueError('Mean compressive strength f_cm must be positive.')
    if alpha_bs < 0:
        raise ValueError('Coefficient alpha_bs must not be negative.')

    return -alpha_bs * (((0.1 * fcm) / (6 + 0.1 * fcm)) ** 2.5) * 10e-6


def beta_bs_t(t: float) -> float:
    """Calculate the time function for basic shrinkage βbs(t).

    fib Model Code 2020, eq. (14.6-74)

    Args:
        t (float): Age of concrete in days.

    Returns:
        float: Time function for basic shrinkage βbs(t).

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if t < 0:
        raise ValueError('Age of concrete t must not be negative.')

    return 1 - math.exp(-0.2 * t**0.5)


# Alpha bs, alpha_ds_1 and alpha_ds_2 coefficients
alpha_coefficients = {
    'CS': (800, 3, 0.013),
    'CN': (700, 4, 0.012),
    'CR': (600, 6, 0.012),
}


def alpha_bs(strength_class: Literal['CS', 'CN', 'CR']) -> tuple:
    """Retrieve the coefficient αbs strength development class of concrete.

    Args:
        strength_class (str): The strength development
            class of concrete ('CS', 'CN', 'CR').

    Returns:
        float: the value of αbs.

    Raises:
        ValueError: If the strength development class is not recognized.
    """
    # Convert input to uppercase to avoid case sensitivity issues
    strength_class = strength_class.upper()

    # Retrieve the coefficients based on the provided
    # strength development class
    if strength_class in alpha_coefficients:
        return alpha_coefficients[strength_class][0]
    raise ValueError(
        "Invalid strength development class. Choose 'CS', 'CN', or 'CR'."
    )


def alpha_ds_1(strength_class: Literal['CS', 'CN', 'CR']) -> tuple:
    """Retrieve the coefficient αds1 strength development class of concrete.

    Args:
        strength_class (str): The strength development
            class of concrete ('CS', 'CN', 'CR').

    Returns:
        float: the value of αds1.

    Raises:
        ValueError: If the strength development class is not recognized.
    """
    # Convert input to uppercase to avoid case sensitivity issues
    strength_class = strength_class.upper()

    # Retrieve the coefficients based on the provided
    # strength development class
    if strength_class in alpha_coefficients:
        return alpha_coefficients[strength_class][1]
    raise ValueError(
        "Invalid strength development class. Choose 'CS', 'CN', or 'CR'."
    )


def alpha_ds_2(strength_class: Literal['CS', 'CN', 'CR']) -> tuple:
    """Retrieve the coefficient αds2 strength development class of concrete.

    Args:
        strength_class (str): The strength development
            class of concrete ('CS', 'CN', 'CR').

    Returns:
        float: the value of αds2.

    Raises:
        ValueError: If the strength development class is not recognized.
    """
    # Convert input to uppercase to avoid case sensitivity issues
    strength_class = strength_class.upper()

    # Retrieve the coefficients based on the provided
    # strength development class
    if strength_class in alpha_coefficients:
        return alpha_coefficients[strength_class][2]
    raise ValueError(
        "Invalid strength development class. Choose 'CS', 'CN', or 'CR'."
    )


def eps_cds_0_fcm(fcm: float, alpha_ds1: float, alpha_ds2: float) -> float:
    """Calculate the notional drying shrinkage coefficient εcds0(fcm).

    fib Model Code 2020, eq. (14.6-75)

    Args:
        fcm (float): Mean compressive strength at the age of 28 days in MPa.
        alpha_ds1 (float): Coefficient dependent on the strength
            development class of concrete.
        alpha_ds2 (float): Coefficient dependent on the strength
            development class of concrete.

    Returns:
        float: Notional drying shrinkage coefficient εcds0(fcm).

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if fcm <= 0:
        raise ValueError('Mean compressive strength f_cm must be positive.')

    return ((220 + 110 * alpha_ds1) * math.exp(-alpha_ds2 * fcm)) * 10e-6


def beta_RH_s(RH: float, RH_eq: float) -> float:
    """Calculate the relative humidity factor βRH(RH) for shrinkage.

    fib Model Code 2020, eq. (14.6-76)

    Args:
        RH (float): Relative humidity of the ambient atmosphere in %.
        RH_eq (float): value of equivalent humidity.

    Returns:
        float: Relative humidity factor βRH(RH).

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if RH < 20 or RH > 100:
        raise ValueError('Relative humidity RH must be between 0 and 100.')

    ratio = RH / RH_eq
    if 20 <= RH <= RH_eq:
        return -1.55 * (1 - ratio**3)
    if RH_eq < RH < 100:
        return -1.55 * (1 - ratio**2)
    if RH == 100:
        return -1.55 * (1 - ratio**2) + 0.25

    raise ValueError('Relative humidity RH is out of defined range.')


def beta_ds_t_ts(t: float, ts: float, h: float) -> float:
    """Calculate the time-development function
        for drying shrinkage βds(t - ts).

    fib Model Code 2020, eq. (14.6-77)

    Args:
        t (float): Age of concrete in days.
        ts (float): Age of concrete at the beginning of drying in days.
        h (float): Notional size of the member in mm.

    Returns:
        float: Time-development function for drying shrinkage βds(t - ts).

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if t < ts:
        raise ValueError(
            'Concrete age t must be greater than or '
            + 'equal to the beginning of drying ts.'
        )
    if h <= 0:
        raise ValueError('Notional size h must be positive.')

    t_diff = t - ts
    return (t_diff / (0.0035 * h**2 + t_diff)) ** 0.5


def RH_eq(fcm: float) -> float:
    """Calculate the equivalent relative humidity RH_eq for a given
    fcm.

    Args:
        fcm (float): Mean compressive strength at the age of 28 days in MPa.
    """
    if fcm <= 0:
        raise ValueError('Mean compressive strength f_cm must be positive.')

    return min(99, 99 * (35 / fcm) ** 0.1)


def eps_lcs_t_ts(
    eps_cs_t_ts: float,
    concrete_grade: Literal[
        'LC8',
        'LC12',
        'LC16',
        'LC20',
        'LC25',
        'LC30',
        'LC35',
        'LC40',
        'LC45',
        'LC50',
        'LC55',
        'LC60',
        'LC70',
        'LC80',
    ],
) -> float:
    """Calculate the shrinkage of lightweight aggregate concrete εlcs(t, ts).

    fib Model Code 2020, eq. (14.6-79)

    Args:
        eps_cs_t_ts (float): Shrinkage strain for normal
            weight concrete εcs(t, ts).
        concrete_grade (str): Grade of lightweight aggregate concrete
            (e.g., 'LC8', 'LC12', 'LC16', 'LC20').

    Returns:
        float: Shrinkage strain for lightweight aggregate concrete εlcs(t, ts).

    Raises:
        ValueError: if invalid concrete grade.
    """
    # Define eta values for different concrete grades
    eta_values = {
        'LC8': 1.5,
        'LC12': 1.5,
        'LC16': 1.5,
        'LC20': 1.2,
        'LC25': 1.2,
        'LC30': 1.2,
        'LC35': 1.2,
        'LC40': 1.2,
        'LC45': 1.2,
        'LC50': 1.2,
        'LC55': 1.2,
        'LC60': 1.2,
        'LC70': 1.2,
        'LC80': 1.2,
    }

    # Convert input to uppercase to avoid case sensitivity issues
    concrete_grade = concrete_grade.upper()

    # Retrieve the eta value based on the provided concrete grade
    if concrete_grade in eta_values:
        eta = eta_values[concrete_grade]
    else:
        raise ValueError(
            'Invalid concrete grade. Choose from '
            + "'LC8', 'LC12', 'LC16', 'LC20', and higher."
        )

    # Calculate shrinkage for lightweight aggregate concrete
    return eta * eps_cs_t_ts


def t_T_maturity(
    temperature_intervals: List[Tuple[float, float]],
) -> float:
    """Calculate the temperature-adjusted concrete age tT.

    fib Model Code 2020, eq. (14.6-80)

    Args:
        temperature_intervals (List[Tuple[float, float]]): A list of tuples
            where each tuple contains
            (Δti, T(Δti)), where Δti is the number of days
            where a temperature T prevails, and T(Δti) is the
            mean temperature in °C during the time period Δti.

    Returns:
        float: Temperature-adjusted concrete age tT in days.

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if not temperature_intervals:
        raise ValueError('Temperature intervals list must not be empty.')

    t_T = 0.0

    for delta_t, T in temperature_intervals:
        if delta_t < 0:
            raise ValueError('Number of days Δti must not be negative.')
        if T < -273.15:
            raise ValueError(
                'Temperature T(Δti) must be above absolute zero (-273°C).'
            )

        t_T += delta_t * math.exp(13.65 - 4000 / (273 + T))

    return t_T


def eps_c_T(delta_T: float, alpha_T: str) -> float:
    """Calculate the thermal strain εcT due to thermal expansion of concrete.

    fib Model Code 2020, eq. (14.6-81)

    Args:
        delta_T (float): Change in temperature in Kelvin (K).
        alpha_T (str): thermal expecansion coefficient.

    Returns:
        float: Thermal strain εcT.

    Raises:
        ValueError: If input values are not within
            valid ranges or invalid concrete type.
    """
    if alpha_T < 0:
        raise ValueError('alpha_T must not be negative.')

    # Calculate thermal strain
    return alpha_T * delta_T


def alpha_T(concrete_type: Literal['normal', 'lightweight']) -> float:
    """Return the thermal expansion coefficient as a function
        of the concrete type.

    Args:
        concrete_type (str): concrete type.

    Returns:
        float: the thermal expansion coefficient in K-1
    """
    return 8e-6 if concrete_type == 'lightweight' else 10e-6


def fcm_T(fcm: float, T: float) -> float:
    """Calculate the compressive strength fcm(T) for
        normal weight concrete at a given temperature.

    fib Model Code 2020, eq. (14.6-82)

    Args:
        fcm (float): Compressive strength in MPa at T = 20°C.
        T (float): Temperature in °C.

    Returns:
        float: Compressive strength fcm(T) at temperature T in MPa.

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if fcm <= 0:
        raise ValueError('Compressive strength fcm must be positive.')

    return fcm * (1.06 - 0.003 * T)


def flcm_T(flcm: float, T: float) -> float:
    """Calculate the compressive strength flcm(T)
        for lightweight aggregate concrete at a given temperature.

    fib Model Code 2020, eq. (14.6-83)

    Args:
        flcm (float): Compressive strength in MPa at T = 20°C.
        T (float): Temperature in °C.

    Returns:
        float: Compressive strength flcm(T) at temperature T in MPa.

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if flcm <= 0:
        raise ValueError('Compressive strength flcm must be positive.')

    return flcm * (1.04 - 0.002 * T)


def fcm_hT_T(fcm_href_Tref: float, D_fc: float, S_fc: float) -> float:
    """Calculate the compressive strength fcm(h_T, T) under drying conditions.

    fib Model Code 2020, eq. (14.6-84)

    Args:
        fcm_href_Tref (float): Compressive strength in MPa
            at Tref = 20°C at reference storage conditions.
        D_fc (float): Drying effect parameter.
        S_fc (float): Thermal dilation incompatibility factor.

    Returns:
        float: Compressive strength fcm(h_T, T) in MPa.

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if fcm_href_Tref <= 0:
        raise ValueError('Compressive strength f_cm_ref must be positive.')

    return fcm_href_Tref * (1 + D_fc + S_fc)


def D_fc(h_Tref: float) -> float:
    """Calculate the drying effect parameter D_fc.

    fib Model Code 2020, eq. (14.6-85)

    Args:
        T (float): Temperature in °C.
        h_Tref (float): Relative humidity of the
            concrete pores at temperature T.

    Returns:
        float: Drying effect parameter D_fc.

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if h_Tref < 0 or h_Tref > 1:
        raise ValueError('Relative humidity h_T must be between 0 and 1.')

    return -0.32 * h_Tref**26 + 0.32


def hT_Tref(hT: float, K_hT: float, T: float, T_ref: float = 20.0) -> float:
    """Calculate the relative humidity in the concrete pores h_Tref.

    fib Model Code 2020, eq. (14.6-86)

    Args:
        hT (float): Relative humidity at reference temperature T_ref.
        K_hT (float): Coefficient for computing the relative humidity.
        T (float): Temperature in °C.
        T_ref (float): Reference temperature in °C, default is 20°C.

    Returns:
        float: Relative humidity h_Tref in the concrete pores at temperature T.

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if hT < 0 or hT > 1:
        raise ValueError('Relative humidity h_ref must be between 0 and 1.')

    return hT - K_hT * (T - T_ref)


def K_hT(hT: float) -> float:
    """Calculate the humidity effect coefficient K_hT.

    fib Model Code 2020, eq. (14.6-87)

    Args:
        h_ref (float): Relative humidity at reference temperature T_ref.

    Returns:
        float: Humidity effect coefficient K_hT.

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if hT < 0 or hT > 1:
        raise ValueError('Relative humidity h_ref must be between 0 and 1.')

    return 0.0135 * hT * (1 - hT) / (1.25 - hT)


def S_fc(T: float, alpha_fc: float) -> float:
    """Calculate the thermal dilation incompatibility factor S_fc.

    fib Model Code 2020, eq. (14.6-88)

    Args:
        T (float): Temperature in °C.
        alpha_fc (float): Alpha coefficient.

    Returns:
        float: Thermal dilation incompatibility factor S_fc.

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    ratio = (T - 20) / 100
    return 0.35 * ratio**2 + alpha_fc * ratio


def alpha_fc(h_Tref: float) -> float:
    """Calculate the alpha coefficient alpha_fc.

    fib Model Code 2020, eq. (14.6-89)

    Args:
        h_Tref (float): Relative humidity at reference temperature T_ref.

    Returns:
        float: Alpha coefficient alpha_fc.

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if h_Tref < 0 or h_Tref > 1:
        raise ValueError('Relative humidity h_ref must be between 0 and 1.')

    return 2.3e-9 * math.exp(19.1 * h_Tref) - 0.55 * math.exp(0.64 * h_Tref)


def fctm_T(fctm: float, T: float) -> float:
    """Calculate the uniaxial tensile strength fctm(T)
    for concrete at a given temperature.

    fib Model Code 2020, eq. (14.6-90)

    Args:
        fctm (float): Uniaxial tensile strength in MPa at T = 20°C.
        T (float): Temperature in °C.

    Returns:
        float: Uniaxial tensile strength fctm(T) at temperature T in MPa.

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if fctm <= 0:
        raise ValueError('Uniaxial tensile strength fctm must be positive.')
    return fctm * (1.16 - 0.008 * T)


def GF_T(GF: float, T: float, concrete_type: Literal['dry', 'mass']) -> float:
    """Calculate the fracture energy GF(T) for dry
        concrete at a given temperature.

    fib Model Code 2020, eq. (14.6-93a)

    Args:
        GF (float): Fracture energy in N/m at T = 20°C.
        T (float): Temperature in °C.
        concrete_type (str): Concrete type.

    Returns:
        float: Fracture energy G_F(T) at temperature T in N/m.

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if GF <= 0:
        raise ValueError('Fracture energy GF must be positive.')

    return (
        GF * (1.06 - 0.003 * T)
        if concrete_type == 'dry'
        else GF * (1.12 - 0.006 * T)
    )


def fctm_hT_T(fctm_ref: float, D_fct: float, S_fct: float) -> float:
    """Calculate the uniaxial tensile strength f_ctm(h_T, T) under
        combined effects of temperature and drying.

    fib Model Code 2020, eq. (14.6-94)

    Args:
        fctm_ref (float): Uniaxial tensile strength in MPa at Tref = 20°C.
        D_fct (float): Drying effect parameter.
        S_fct (float): Thermal dilation incompatibility factor.

    Returns:
        float: Uniaxial tensile strength f_ctm(h_T, T) in MPa.

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if fctm_ref <= 0:
        raise ValueError(
            'Uniaxial tensile strength f_ctm_ref must be positive.'
        )

    return fctm_ref * (1 + D_fct + S_fct)


def D_fct(T: float, hT_Tref: float) -> float:
    """Calculate the drying effect parameter D_fct.

    fib Model Code 2020, eq. (14.6-95)

    Args:
        T (float): Temperature in °C.
        hT_Tref (float): Relative humidity of the
            concrete pores at temperature T.

    Returns:
        float: Drying effect parameter D_fct.

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if T < 20 or T > 80:
        raise ValueError('Temperature T must be in the range 20°C ≤ T ≤ 80°C.')
    if hT_Tref < 0 or hT_Tref > 1:
        raise ValueError('Relative humidity h_T must be between 0 and 1.')

    return -0.44 * hT_Tref**4 + 0.44


def S_fct(T: float, alpha_fct: float) -> float:
    """Calculate the thermal dilation incompatibility factor S_fct.

    fib Model Code 2020, eq. (14.6-96)

    Args:
        T (float): Temperature in °C.
        alpha_fct (float): Alpha coefficient.

    Returns:
        float: Thermal dilation incompatibility factor S_fct.

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if T < 20 or T > 80:
        raise ValueError('Temperature T must be in the range 20°C ≤ T ≤ 80°C.')

    ratio = (T - 20) / 100
    return 0.45 * ratio**2 + alpha_fct * ratio


def alpha_fct(hT_Tref: float) -> float:
    """Calculate the alpha coefficientalpha_fct.

    fib Model Code 2020, eq. (14.6-97)

    Args:
        hT_Tref (float): Relative humidity at reference temperature T_ref.

    Returns:
        float: Alpha coefficient alpha_fct.

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if hT_Tref < 0 or hT_Tref > 1:
        raise ValueError('Relative humidity h_ref must be between 0 and 1.')

    return 7.6e-5 * math.exp(8.5 * hT_Tref) - 0.72 * math.exp(0.48 * hT_Tref)


def E_ci_T(E_ci: float, T: float) -> float:
    """Calculate the modulus of elasticity E_ci(T)
        for normal weight concrete at a given temperature.

    fib Model Code 2020, eq. (14.6-98)

    Args:
        E_ci (float): Modulus of elasticity in MPa at T = 20°C.
        T (float): Temperature in °C.

    Returns:
        float: Modulus of elasticity E_ci(T) at temperature T in MPa.

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if E_ci <= 0:
        raise ValueError('Modulus of elasticity E_ci must be positive.')
    if T < 0 or T > 80:
        raise ValueError('Temperature T must be in the range 0°C ≤ T ≤ 80°C.')

    return E_ci * (1.06 - 0.003 * T)


def El_ci_T(E_l_ci: float, T: float) -> float:
    """Calculate the modulus of elasticity E_l_ci(T)
        for lightweight aggregate concrete at a given temperature.

    fib Model Code 2020, eq. (14.6-99)

    Args:
        E_lci (float): Modulus of elasticity in MPa at T = 20°C.
        T (float): Temperature in °C.

    Returns:
        float: Modulus of elasticity E_l_ci(T) at temperature T in MPa.

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if E_l_ci <= 0:
        raise ValueError('Modulus of elasticity E_lci must be positive.')
    if T < 0 or T > 80:
        raise ValueError('Temperature T must be in the range 0°C ≤ T ≤ 80°C.')

    return E_l_ci * (1.04 - 0.002 * T)


def E_cm_hT_T(E_cm_ref: float, D_Ec: float, S_Ec: float) -> float:
    """Calculate the modulus of elasticity E_cm(hT, T)
        under combined effects of temperature and drying.

    fib Model Code 2020, eq. (14.6-100)

    Args:
        E_cm_ref (float): Modulus of elasticity in MPa at
            Tref = 20°C at reference storage conditions.
        D_Ec (float): Drying effect parameter.
        S_Ec (float): Thermal dilation incompatibility factor.

    Returns:
        float: Modulus of elasticity E_cm(h_T, T) in MPa.

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if E_cm_ref <= 0:
        raise ValueError('Modulus of elasticity E_cm_ref must be positive.')

    return E_cm_ref * (1 + D_Ec + S_Ec)


def D_Ec(T: float, hT_Tref: float) -> float:
    """Calculate the drying effect parameter D_Ec.

    fib Model Code 2020, eq. (14.6-101)

    Args:
        T (float): Temperature in °C.
        hT_Tref (float): Relative humidity of the
            concrete pores at temperature T.

    Returns:
        float: Drying effect parameter D_Ec.

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if T < 20 or T > 80:
        raise ValueError('Temperature T must be in the range 20°C ≤ T ≤ 80°C.')
    if hT_Tref < 0 or hT_Tref > 1:
        raise ValueError('Relative humidity h_T must be between 0 and 1.')

    if hT_Tref == 1.0:
        return 0

    return -0.12 * hT_Tref**19 - 0.14


def S_Ec(T: float, alpha_Ec: float) -> float:
    """Calculate the thermal dilation incompatibility factor S_Ec.

    fib Model Code 2020, eq. (14.6-102)

    Args:
        T (float): Temperature in °C.
        alpha_Ec (float): Alpha coefficient for elastic modulus.

    Returns:
        float: Thermal dilation incompatibility factor S_Ec.

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if T < 20 or T > 80:
        raise ValueError('Temperature T must be in the range 20°C ≤ T ≤ 80°C.')

    ratio = (T - 20) / 100
    return alpha_Ec * ratio


def alpha_Ec(hT_Tref: float) -> float:
    """Calculate the alpha coefficient alpha_Ec for the modulus of elasticity.

    fib Model Code 2020, eq. (14.6-103)

    Args:
        h_ref (float): Relative humidity at reference temperature T_ref.

    Returns:
        float: Alpha coefficient alpha_Ec.

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if hT_Tref < 0 or hT_Tref > 1:
        raise ValueError('Relative humidity h_ref must be between 0 and 1.')

    return 2.3e-5 * math.exp(9.0 * hT_Tref) - 0.2 * math.exp(0.7 * hT_Tref)


def beta_h_T(beta_h: float, T: float) -> float:
    """Calculate the temperature-dependent coefficient
        for creep time-development βh,T.

    fib Model Code 2020, eqs. (14.6-104) and (14.6-105)

    Args:
        beta_h (float): Coefficient βh according to Eq. (14.6-66c).
        T (float): Temperature in °C.

    Returns:
        float: Temperature-dependent coefficient βh,T.

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if beta_h <= 0:
        raise ValueError('Coefficient βh must be positive.')
    if T < 0 or T > 80:
        raise ValueError('Temperature T must be in the range 0°C ≤ T ≤ 80°C.')

    beta_T = math.exp(1500 / ((273 + T) - 5.12))
    return beta_h * beta_T


def phi_bc_T(phi_bc: float, T: float) -> float:
    """Calculate the temperature-dependent basic creep coefficient φbc,T.

    fib Model Code 2020, eq. (14.6-106) and (14.6-108)

    Args:
        phi_bc (float): Basic creep coefficient φbc according to Eq. (14.6-59).
        T (float): Temperature in °C.

    Returns:
        float: Temperature-dependent basic creep coefficient φbc,T.

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if phi_bc <= 0:
        raise ValueError('Basic creep coefficient φbc must be positive.')
    if T < 0 or T > 80:
        raise ValueError('Temperature T must be in the range 0°C ≤ T ≤ 80°C.')

    phi_T = math.exp(0.015 * (T - 20))
    return phi_bc * phi_T


def phi_dc_T(phi_dc: float, T: float) -> float:
    """Calculate the temperature-dependent drying creep coefficient φdc,T.

    fib Model Code 2020, eq. (14.6-107) and (14.6-108)

    Args:
        phi_dc (float): Drying creep coefficient φdc
            according to Eq. (14.6-62).
        T (float): Temperature in °C.

    Returns:
        float: Temperature-dependent drying creep coefficient φdc,T.

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if phi_dc <= 0:
        raise ValueError('Drying creep coefficient φdc must be positive.')
    if T < 0 or T > 80:
        raise ValueError('Temperature T must be in the range 0°C ≤ T ≤ 80°C.')

    phi_T = math.exp(0.015 * (T - 20))
    return phi_dc * phi_T**1.2


def phi_t_t0_T(phi_t_t0: float, T: float) -> float:
    """Creep coefficient ΔφT,trans for taking into consideration
        increase in temperature while member is under load.

    fib Model Code 2020, eq. (14.6-109) and (14.6-110)

    Args:
        phi_t_t0 (float): base creep coefficient.
        T (float): Temperature in °C.

    Returns:
        float: Transient thermal creep coefficient ΔφT,trans.

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if T < 0 or T > 80:
        raise ValueError('Temperature T must be in the range 0°C ≤ T ≤ 80°C.')

    delta_phi_T_trans = 0.0004 * (T - 20) ** 2
    return phi_t_t0 + delta_phi_T_trans


def alpha_sT(T: float, h: float) -> float:
    """Calculate the temperature-dependent coefficient
        for drying shrinkage time development αsT(T).

    fib Model Code 2020, eq. (14.6-111)

    Args:
        T (float): Temperature in °C.
        h (float): Notional size parameter.

    Returns:
        float: Temperature-dependent coefficient αsT(T).

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if T < 0 or T > 80:
        raise ValueError('Temperature T must be in the range 0°C ≤ T ≤ 80°C.')
    if h <= 0:
        raise ValueError('Notional size parameter h must be positive.')

    return 0.035 * h**2 * math.exp(-0.06 * (T - 20))


def beta_RH_T(beta_RH: float, beta_sT: float) -> float:
    """Calculate the temperature-dependent coefficient
        for drying shrinkage magnitude βRH,T.

    fib Model Code 2020, eq. (14.6-112)

    Args:
        beta_RH (float): Coefficient βRH.
        beta_sT (float): Coefficient βsT.

    Returns:
        float: Temperature-dependent coefficient βRH,T.

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if beta_RH <= 0 or beta_sT <= 0:
        raise ValueError('Coefficients βRH and βsT must be positive.')

    return beta_RH * beta_sT


def beta_sT(RH: float, T: float) -> float:
    """Calculate the temperature-dependent coefficient replacing βRH.

    fib Model Code 2020, eq. (14.6-113)

    Args:
        RH (float): Relative humidity of the ambient environment in %.
        T (float): Temperature in °C.

    Returns:
        float: Temperature-dependent coeSfficient βRH,T.

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if RH < 0 or RH > 100:
        raise ValueError('Relative humidity RH must be between 0 and 100.')
    if T < 0 or T > 80:
        raise ValueError('Temperature T must be in the range 0°C ≤ T ≤ 80°C.')

    return 1 + (4 / (103 - RH)) * ((T - 20) / 40)


def beta_RH(RH: float, RH_T: float) -> float:
    """Calculate the coefficient βRH,T for relative humidity.

    fib Model Code 2020, eq. (14.6-114)

    Args:
        RH (float): Relative humidity of the ambient environment in %.
        RH_T (float): Relative humidity in % obtained from eq. (14.6-115)

    Returns:
        float: Coefficient βRH.

    Raises:
        ValueError: If input values are not within valid ranges.
    """
    if RH < 0 or RH > 100:
        raise ValueError('Relative humidity RH must be between 0 and 100.')

    if 40 <= RH < RH_T:
        return -1.55 * (1 - (RH / 100) ** 3)
    if RH >= RH_T:
        return 0.25
    return ValueError('RH has not a valid value.')


def RH_T(fcm: float, T: float) -> float:
    """Computes the relative humidity in % adjusted for the temperature T.

    fib Model Code 2020, eq. (14.6-115), (14.6-116) and  (14.6-117)

    Args:
        fcm (float): mean compressive strength of concrete in MPa.
        T (float): Temperature in °C.

    Returns:
        float: Value of RH_T.

    Raise:
        ValueError: if values are not within a valid range.
    """
    if fcm < 0:
        raise ValueError('fcm cannot be negative.')
    if T < 0 or T > 80:
        raise ValueError('Temperature T must be in the range 0°C ≤ T ≤ 80°C.')

    beta_s1 = min((35 / fcm) ** 0.1, 1)
    beta_s1_T = ((T - 20) / 25) ** 3

    return min(99 * beta_s1 - beta_s1_T, 100)
