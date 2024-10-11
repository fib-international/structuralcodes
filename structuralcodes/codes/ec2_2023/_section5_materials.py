"""Functions from Section 5 of EN 1992-1-1:2023."""

import math
import typing as t

import numpy as np
import scipy.interpolate

from structuralcodes.codes import mc2010

VALID_STRENGTH_DEV_CLASSES = ('CS', 'CN', 'CR', 'SLOW', 'NORMAL', 'RAPID')

DUCTILITY_CLASSES = {
    'A': {
        'epsuk': 2.5e-2,
        'k': 1.05,
    },
    'B': {
        'epsuk': 5.0e-2,
        'k': 1.08,
    },
    'C': {
        'epsuk': 7.5e-2,
        'k': 1.15,
    },
}


def fcm(fck: float, delta_f: float = 8.0) -> float:
    """Determines the mean strength of concrete from its characteristic value.

    EN 1992-1-1:2023, Table 5.1.

    Args:
        fck (float): Is the characteristic compressive strength in MPa.

    Keyword Args:
        delta_s (float): The increment in MPa to compute the mean compressive
            strength.

    Returns:
        float: The mean compressive strength in MPa
    """
    return mc2010.fcm(fck, delta_f)


def fctm(fck: float) -> float:
    """Compute the mean concrete tensile strength from the characteristic
    compressive strength.

    EN 1992-1-1:2023, Table 5.1.

    Args:
        fck (float): The characteristic compressive strength in MPa.

    Returns:
        float: The mean tensile strength in MPa.
    """
    if abs(fck) <= 50:
        return 0.3 * math.pow(abs(fck), 2 / 3)
    return 1.1 * math.pow(abs(fck), 1 / 3)


def fctk_5(fctm: float) -> float:
    """Compute the 5% mean concrete tensile strength fractile.

    EN 1992-1-1:2023, Table 5.1.

    Args:
        fctm (float): The mean concrete tensile strength in MPa.

    Returns:
        float: The 5% mean concrete tensile strength fractile in MPa.
    """
    return abs(fctm) * 0.7


def fctk_95(fctm: float) -> float:
    """Compute the 95% mean concrete tensile strength fractile.

    EN 1992-1-1:2023, Table 5.1.

    Args:
        fctm (float): The mean concrete tensile strength in MPa.

    Returns:
        float: The 5% mean concrete tensile strength fractile in MPa.
    """
    return abs(fctm) * 1.3


def Ecm(fcm: float, kE: float = 9500) -> float:
    """Computes the secant modulus between sigma_c=0 and sigma_c=0.4*fcm.

    EN 1992-1-1:2023, Eq. (5.1).

    Args:
        fcm (float): The mean compressive strength in MPa.

    Keyword Args:
        kE (float): Coefficient relating the aggregates used in concrete.
            Default value is 9500, but it can vary from 5000 to 13000.

    Returns:
        float: The secant concrete modulus in MPa.

    Raises:
        ValueError: If fcm is less than 0.
        ValueError: If kE is not between 5000 and 13000.
    """
    if fcm < 0:
        raise ValueError(f'fcm={fcm} cannot be less than 0')
    if kE < 5000:
        raise ValueError(f'kE={kE} cannot be less than 5000')
    if kE > 13000:
        raise ValueError(f'kE={kE} cannot be larger than 13000')

    return kE * math.pow(fcm, 1 / 3)


def hn(Ac: float, u: float) -> float:
    """Computes the notional size of a given concrete cross-section.

    EN 1992-1-1:2023, Table 5.2.

    Args:
        Ac (float): The concrete cross-sectional area in (mm2).
        u (float): The perimeter exposed to drying (mm).

    Returns:
        float: The notional size (mm).

    Raises:
        ValueError: If Ac is less than 0.
        ValueError: If u is less than 0.
    """
    if Ac < 0:
        raise ValueError(f'Ac={Ac} cannot be less than 0')
    if u < 0:
        raise ValueError(f'u={u} cannot be less than 0')

    return 2 * Ac / u


def A_phi_correction_exp(
    hn: float, atm_conditions: t.Literal['dry', 'humid']
) -> float:
    """Computes the correction exponent for the modification for the phi_50y_t0
    with respect the fck value.

    EN 1992-1-1:2023, Table 5.2.

    Args:
        hn (float): The notional size in mm.
        atm_conditions (str): 'dry' or 'humid'.

    Returns:
        float: The correction exponent value.

    Raises:
        ValueError: If hn is less than 100 or greater than 1000.
        ValueError: If atm_conditions is not 'dry' or 'humid'.
    """
    if hn < 100:
        raise ValueError(f'hn={hn} cannot be less than 100')
    if hn > 1000:
        raise ValueError(f'hn={hn} cannot be larger than 1000')

    atm_conditions = atm_conditions.lower().strip()
    x = (100, 200, 500, 1000)
    if atm_conditions == 'dry':
        y = (0.82, 0.79, 0.75, 0.72)
    elif atm_conditions == 'humid':
        y = (0.71, 0.68, 0.66, 0.64)
    else:
        raise ValueError(
            f'atm_conditions={atm_conditions} can only take '
            + '"dry" and "humid" as values'
        )

    interpol = scipy.interpolate.interp1d(x, y)
    return interpol(hn)


def phi_correction_factor(fck: float, A_exponent: float) -> float:
    """Computes the correction factor for the computation of the phi_50y_t0.

    EN 1992-1-1:2023, Table 5.2.

    Args:
        fck (float): Characteristic strength of concrete in MPa.
        A_exponent (float): The A correction exponent value.

    Returns:
        float: The correction factor value.

    Raises:
        ValueError: If fck is not between 12 and 100 MPa.
        ValueError: If A_exponent is not between 0.64 and 0.82.
    """
    if fck < 12:
        raise ValueError(f'fck={fck} cannot be less than 12')
    if fck > 100:
        raise ValueError(f'fck={fck} cannot be larger than 100')
    if A_exponent < 0.64:
        raise ValueError(f'A_exponent={A_exponent} cannot be less than 0.64')
    if A_exponent > 0.82:
        raise ValueError(f'A_exponent={A_exponent} cannot be less than 0.82')

    return math.pow(35 / fck, A_exponent)


def phi_50y_t0(
    t0: float,
    atm_conditions: t.Literal['dry', 'humid'],
    hn: float,
    strength_dev_class: t.Literal['CS', 'CN', 'CR', 'slow', 'normal', 'rapid'],
) -> float:
    """Computes the creep coefficient of plain concrete at 50 years of loading.
    Interpolation is linear between values.

    EN 1992-1-1:2023, Table 5.2.

    Args:
        t0 (float): Age at loading [days].
        atm_conditions (str): 'dry' or 'humid'.
        hn (float): The notional size in mm.
        strength_dev_class (str): Strength development class 'CS', 'CN', 'CR',
            'slow', 'normal' or 'rapid'.

    Returns:
        float: The creep coefficient.

    Raises:
        ValueError: If t0 is less than 1.
        ValueError: If atm_conditions is not 'dry' or 'humid'.
        ValueError: If hn is less than 100 or larger than 1000.
        ValueError: If strength_dev_class is not 'CS', 'CN',
            'CR', 'slow', 'normal' or 'rapid'.
        ValueError: If combination of t0 and hn is out of scope.
    """
    if t0 < 1:
        raise ValueError(f't0={t0} cannot be less than 1')

    atm_conditions = atm_conditions.lower().strip()
    if atm_conditions not in ('dry', 'humid'):
        raise ValueError(
            f'atm_conditions={atm_conditions} must be "dry" or "humid"'
        )

    if hn < 100:
        raise ValueError(f'hn={hn} must be larger or equal than 100')
    if hn > 1000:
        raise ValueError(f'hn={hn} must be less or equal than 1000')

    strength_dev_class = strength_dev_class.upper().strip()
    if strength_dev_class not in VALID_STRENGTH_DEV_CLASSES:
        raise ValueError(
            f'strength_dev_class={strength_dev_class} must be'
            + '"CS", "CN", "CR", "slow", "normal" or "rapid"'
        )

    if strength_dev_class in ('CS', 'SLOW'):
        _t = (3, 10, 32, 91, 365)
    elif strength_dev_class in ('CN', 'NORMAL'):
        _t = (1, 7, 28, 91, 365)
    elif strength_dev_class in ('CR', 'RAPID'):
        _t = (1, 3, 23, 91, 365)

    h_v = (100, 200, 500, 1000)

    if atm_conditions == 'dry':
        values = (
            4.2,
            3.8,
            3.4,
            3.1,
            3.1,
            2.8,
            2.5,
            2.3,
            2.4,
            2.2,
            1.9,
            1.8,
            1.9,
            1.7,
            1.5,
            1.4,
            1.4,
            1.3,
            1.1,
            1.0,
        )
    elif atm_conditions == 'humid':
        values = (
            3.0,
            2.8,
            2.6,
            2.5,
            2.2,
            2.1,
            2.0,
            1.9,
            1.7,
            1.6,
            1.6,
            1.5,
            1.4,
            1.3,
            1.2,
            1.2,
            1.0,
            0.9,
            0.9,
            0.8,
        )

    grid = np.array(np.meshgrid(_t, h_v)).T.reshape(-1, 2)
    p = (t0, hn)

    interp = scipy.interpolate.griddata(grid, values, p, method='linear')
    _phi_50y_t0 = float(interp)

    if math.isnan(_phi_50y_t0) or math.isnan(_phi_50y_t0):
        raise ValueError('Combination of t0, hn out of scope')

    return _phi_50y_t0


def eps_cs_50y(  # noqa: PLR0912
    fck_28: float,
    atm_conditions: t.Literal['dry', 'humid'],
    hn: float,
    strength_dev_class: t.Literal['CS', 'CN', 'CR', 'slow', 'normal', 'rapid'],
) -> float:
    """Computes the nominal total shrinkage in ‰ for concrete after a duration
    of drying of 50 years.

    EN 1992-1-1:2023, Table 5.3.

    Args:
        fck_28 (float): Characteristic strength at 28 days in MPa
        atm_conditions (str): 'dry' or 'humid'.
        hn (float): The notional size in mm.
        strength_dev_class (str): Strength development class 'CS', 'CN', 'CR',
            'slow', 'normal' or 'rapid'.

    Returns:
        float: The nominal shrinkage value in percent.

    Raises:
        ValueError: If fck_28 is less than 20 MPa or larger than 80 MPa.
        ValueError: If atm_conditions is not 'dry' or 'humid'.
        ValueError: If hn is less than 100 or larger than 1000.
        ValueError: If strength_dev_class is not CS', 'CN', 'CR',
            'slow', 'normal' or 'rapid'.
        ValueError: If combination of fck_28 and hn is out of scope.
    """
    if fck_28 < 20:
        raise ValueError(f'fck_28={fck_28} cannot be less than 20')
    if fck_28 > 80:
        raise ValueError(f'fck_28={fck_28} cannot be larger than 80')

    atm_conditions = atm_conditions.lower().strip()
    if atm_conditions not in ('dry', 'humid'):
        raise ValueError(
            f'atm_conditions={atm_conditions} must be "dry" or "humid"'
        )

    if hn < 100:
        raise ValueError(f'hn={hn} must be larger or equal than 100')
    if hn > 1000:
        raise ValueError(f'hn={hn} must be less or equal than 1000')

    strength_dev_class = strength_dev_class.upper().strip()
    if strength_dev_class == 'SLOW':
        strength_dev_class = 'CS'
    elif strength_dev_class == 'NORMAL':
        strength_dev_class = 'CN'
    elif strength_dev_class == 'RAPID':
        strength_dev_class = 'CR'

    if strength_dev_class == 'CS':
        fck_v = (20, 35, 50)
    elif strength_dev_class == 'CN':
        fck_v = (20, 35, 50, 80)
    elif strength_dev_class == 'CR':
        fck_v = (35, 50, 80)
    else:
        raise ValueError(
            f'strength_dev_class={strength_dev_class} '
            + 'must be "CS", "CN", "CR", "slow", "normal" or "rapid"'
        )

    h_v = (100, 200, 500, 1000)

    data = {
        'dry': {
            'CS': (
                0.57,
                0.56,
                0.48,
                0.36,
                0.53,
                0.51,
                0.45,
                0.35,
                0.49,
                0.48,
                0.43,
                0.35,
            ),
            'CN': (
                0.67,
                0.65,
                0.56,
                0.41,
                0.60,
                0.59,
                0.51,
                0.39,
                0.55,
                0.54,
                0.48,
                0.37,
                0.48,
                0.48,
                0.43,
                0.36,
            ),
            'CR': (
                0.76,
                0.74,
                0.65,
                0.48,
                0.67,
                0.66,
                0.58,
                0.44,
                0.55,
                0.54,
                0.49,
                0.39,
            ),
        },
        'humid': {
            'CS': (
                0.33,
                0.32,
                0.28,
                0.21,
                0.31,
                0.31,
                0.27,
                0.22,
                0.30,
                0.29,
                0.27,
                0.23,
            ),
            'CN': (
                0.38,
                0.37,
                0.32,
                0.24,
                0.34,
                0.34,
                0.30,
                0.24,
                0.31,
                0.31,
                0.28,
                0.23,
                0.30,
                0.30,
                0.28,
                0.25,
            ),
            'CR': (
                0.42,
                0.41,
                0.36,
                0.28,
                0.36,
                0.35,
                0.32,
                0.26,
                0.31,
                0.30,
                0.28,
                0.25,
            ),
        },
    }
    values = data.get(atm_conditions).get(strength_dev_class)

    grid = np.array(np.meshgrid(fck_v, h_v)).T.reshape(-1, 2)
    p = (fck_28, hn)

    interp = scipy.interpolate.griddata(grid, values, p, method='linear')
    _eps_cs_50y = float(interp)

    if math.isnan(_eps_cs_50y):
        raise ValueError('Combination of fck_28, hn out of scope')

    return _eps_cs_50y


def eta_cc(fck: float, fck_ref: float = 40) -> float:
    """Computes the factor to measure the difference between the undistributed
    compressive strength of a cylinder and the effective compressive strength
    in a structural member.

    EN 1992-1-1:2023, Eq. (5.4).

    Args:
        fck (float): The characterisitic compressive strength in MPa.

    Keyword Args:
        fck_ref (float, optional): The reference compressive strength MPa.

    Returns:
        float: The value of the factor eta_cc.

    Raises:
        ValueError: If fck is less than 12 MPa.
        ValueError: If fkc_ref is less or equal to 0.
    """
    if fck < 12:
        raise ValueError(f'fck={fck} must be larger or equal than 12 MPa')
    if fck_ref <= 0:
        raise ValueError(f'fck_ref={fck_ref} must be larger than 0')

    return min(math.pow(fck_ref / fck, 1 / 3), 1)


def k_tc(
    t_ref: float,
    t0: float,
    strength_dev_class: t.Literal['CS', 'CN', 'CR', 'slow', 'normal', 'rapid'],
) -> float:
    """Computes the factor for considering the effect of high sustained loads
    and of time of loading on concrete compressive strength.

    EN 1992-1-1:2023, Eq. (5.3).

    Args:
        t_ref (float): The reference time in days.
        t0 (float): Age at loading in days.
        strength_dev_class (str): Strength development class 'CS', 'CN', 'CR',
            'slow', 'normal', 'rapid'.

    Returns:
        float: The factor value.

    Raises:
        ValueError: If t_ref is less than 0.
        ValueError: If t0 is less than 0.
        ValueError: If strength_dev_class is not 'CS', 'CN', 'CR', 'slow',
            'normal', 'rapid'.
    """
    if t_ref < 0:
        raise ValueError(f't_ref={t_ref} must be larger than 0')
    if t0 < 0:
        raise ValueError(f't0={t0} must be larger than 0')

    strength_dev_class = strength_dev_class.upper().strip()

    if strength_dev_class not in VALID_STRENGTH_DEV_CLASSES:
        raise ValueError(
            f'strength_dev_class={strength_dev_class}'
            + f'should can only take {VALID_STRENGTH_DEV_CLASSES}'
        )

    if (
        strength_dev_class in ('CR', 'CN', 'RAPID', 'NORMAL')
        and t_ref <= 28
        and t0 > 90
    ):
        return 1

    if strength_dev_class in ('CS', 'SLOW') and t_ref <= 56 and t0 > 90:
        return 1

    return 0.85


def fcd(fck: float, eta_cc: float, k_tc: float, gamma_c: float) -> float:
    """Computes the value of the design compressive strength of concrete.

    EN 1992-1-1:2023, Eq. (5.3).

    Args:
        fck (float): Characteristic compressive strength in MPa.
        eta_cc (float): Factor for measuring the difference between the
            undistributed compressive strength of a cylinder and the effective
            compressive strength in the real structural member.
        k_tc (float): Factor for taking into consideration high sustained
            loads and of time of loading.
        gamma_c (float): Partial factor of concrete.

    Returns:
        float: The design compressive strength of concrete in MPa.

    Raises:
        ValueError: If fck is less than 12 MPa.
        ValueError: If _etc_cc is not between 0 and 1.
        ValueError: If gamma_c is less than 1.
    """
    if fck < 12:
        raise ValueError(f'fck={fck} must be larger or equal than 12 MPa')
    if eta_cc < 0 or eta_cc > 1:
        raise ValueError(f'eta_cc={eta_cc} must be between 0 and 1')
    if gamma_c < 1:
        raise ValueError(f'gamma_c={gamma_c} must be larger or equal to 1')

    return eta_cc * k_tc * fck / gamma_c


def k_tt(
    t_ref: float,
    strength_dev_class: t.Literal['CS', 'CN', 'CR', 'slow', 'normal', 'rapid'],
) -> float:
    """Computes the factor for considering the effect of high sustained loads
    and of time of loading on concrete tensile strength.

    EN 1992-1-1:2023, Eq. (5.5).

    Args:
        t_ref (float): The reference time in days.
        strength_dev_class (str): Strength development class 'CS', 'CN', 'CR',
            'slow', 'normal' or 'rapid'.

    Returns:
        float: The factor value.

    Raises:
        ValueError: If t_ref is less than 0.
        ValueError: If strength_dev_class is not 'CS', 'CN', 'CR', 'slow',
            'normal' or 'rapid'.
    """
    if t_ref < 0:
        raise ValueError(f't_ref={t_ref} must be larger than 0')

    strength_dev_class = strength_dev_class.upper().strip()

    if strength_dev_class not in VALID_STRENGTH_DEV_CLASSES:
        raise ValueError(
            f'strength_dev_class={strength_dev_class}'
            + 'should can only take "CS", "CN", "CR", "slow"'
            + '"normal", "rapid" as values',
        )

    if strength_dev_class in ('CR', 'CN', 'RAPID', 'NORMAL') and t_ref <= 28:
        return 0.8

    if strength_dev_class in ('CS', 'SLOW') and t_ref <= 56:
        return 0.8

    return 0.7


def fctd(fctk_5: float, k_tt: float, gamma_c: float) -> float:
    """Computes the value of the design tensile strength of concrete.

    EN 1992-1-1:2023, Eq. (5.5).

    Args:
        fctk_5 (float): The 5% mean concrete tensile strength fractile in MPa.
        k_tt (float): The factor for considering the effect of high sustained
            loads and of time of loading on concrete tensile strength.
        gamma_c (float): Partial factor of concrete.

    Returns:
        float: The design tensile strength of concrete in MPa.

    Raises:
        ValueError: If fctk_5 is less than 0.
        ValueError: If gamma_c is less than 1.
    """
    if fctk_5 < 0:
        raise ValueError(f'fctk_5={fctk_5} must be larger or equal to 0')
    if gamma_c < 1:
        raise ValueError(f'gamma_c={gamma_c} must be larger or equal to 1')

    return k_tt * fctk_5 / gamma_c


def eps_c1(fcm: float) -> float:
    """Computes the strain at maximum compressive strength of concrete (fcm)
    for the Sargin constitutive law.

    EN 1992-1-1:2023, Eq. (5.9).

    Args:
        fcm (float): The mean strength of concrete in MPa.

    Returns:
        float: The strain at maximum compressive strength of concrete.

    Raises:
        ValueError: If fcm is less than 12+8MPa.
    """
    if fcm < 20:
        raise ValueError(f'fcm={fcm} must be larger or equal to 12+8MPa')

    return min(0.7 * math.pow(fcm, 1 / 3), 2.8) / 1000


def eps_cu1(fcm: float) -> float:
    """Computes the strain at concrete failure of concrete.

    EN 1992-1-1:2023, Eq. (5.10).

    Args:
        fcm (float): The mean strength of concrete in MPa.

    Returns:
        float: The maximum strength at failure of concrete.

    Raises:
        ValueError: If fcm is less than 12+8MPa.
    """
    if fcm < 20:
        raise ValueError(f'fcm={fcm} must be larger or equal to 12+8MPa')

    return min(2.8 + 14 * (1 - fcm / 108) ** 4, 3.5) / 1000


def k_sargin(
    Ecm: float,
    fcm: float,
    eps_c1: float,
) -> float:
    """Computes the coefficient k for Sargin constitutive law in compression.

    EN 1992-1-1:2003, eq. (5.7)

    Args:
        Ecm (float): the secant modulus between sigma_c=0 and sigma_c=0.4*fcm
            in MPa
        fcm (float): the mean compressive strength of concrete in MPa
        eps_c1 (float): the strain of concrete at stress fcm

    Returns:
        float: the coefficient k for Sargin constitutive law

    Raises:
        ValueError: if Ecm is less or equal to 0
        ValueError: if fcm is less than 12+8MPa
        ValueError: if eps_c1 is less or equal to 0
    """
    if Ecm <= 0:
        raise ValueError(f'Ecm={Ecm} must be larger than 0')
    if fcm < 20:
        raise ValueError(f'fcm={fcm} must be larger or equal to 12+8MPa')
    if eps_c1 <= 0:
        raise ValueError(f'eps_c1={eps_c1} must be larger than 0')
    return 1.05 * Ecm * eps_c1 / fcm


def eps_c2() -> float:
    """The strain at maximum compressive stress of concrete for the
    parabolic-rectangular law.

    EN 1992-1-1:2023, Eq. 8.4

    Returns:
        float: The strain at maximum compressive stress, absolute value, no
        unit.
    """
    return 2.0e-3


def eps_cu2() -> float:
    """The ultimate strain of the parabolic-rectangular law.

    EN 1992-1-1:2023, Eq. 8.4

    Returns:
        float: The ultimate strain, absolute value, no unit.
    """
    return 3.5e-3


def n_parabolic_rectangular() -> float:
    """The exponent in the parabolic-rectangular law.

    EN 1992-1-1:2023, Eq. 8.4
    Returns:
        float: The exponent n, absolute value, no unit.
    """
    return 2.0


def weight_c(concrete_type: t.Literal['nc', 'npc']) -> float:
    """Returns the mean unit weight of concrete in kN/m3.

    EN 1992-1-1:2023, 5.1.6-5.

    Args:
        concrete_type (str): 'nc' for normal concrete, or 'npc' for normal
            plain concrete.

    Returns:
        float: Mean unit weight in kN/m3.

    Raises:
        ValueError: If concrete_type is not 'nc' or 'npc'.
    """
    concrete_type = concrete_type.lower().strip()

    if concrete_type == 'nc':
        return 25
    if concrete_type == 'npc':
        return 24

    raise ValueError(
        f'concrete_type={concrete_type} can only take'
        + '"nc" or "npc" as values'
    )


def alpha_c_th() -> float:
    """Returns the linear coefficient of thermal expansion in 1/Cº for
    concrete.

    EN 1992-1-1:2023, 5.1.6-6.

    Returns:
        float: The linear coefficient of thermal expansion in 1/Cº for
        concrete.
    """
    return 10 * 10e-6


def Es() -> float:
    """Returns the value of the modulus of elasticity for weldable reinforcing
    steel.

    EN 1992-1-1:2023, 5.2.4-3.

    Returns:
        float: Modulus of elasticity in MPa.
    """
    return 200000


def alpha_s_th() -> float:
    """Returns the linear coefficient of thermal expansion in 1/Cº for weldable
    reinforced steel.

    EN 1992-1-1:2023, 5.2.4-5

    Returns:
        float: The linear coefficient of thermal expansion in 1/Cº for weldable
        reinforcement steel.
    """
    return 10 * 10e-6


def weight_s() -> float:
    """Returns the mean unit weight of reinforced steel for the purposes of
    design in kN/m3.

    EN 1992-1-1:2023.2.4-4.

    Returns:
        float: The mean unit weight in kN/m3.
    """
    return 78.5


def fyd(fyk: float, gamma_s: float) -> float:
    """Design value for the yielding stress for welding reinforcing steel.

    EN 1992-1-1:2023, Eq (5.11).

    Args:
        fyk (float): Characteristic yield stress for the steel in MPa.
        gamma_s (float): Safety coefficient.

    Returns:
        float: Design yielding stress for steel in MPa.

    Raises:
        ValueError: If fyk is less than 0.
        ValueError: If gamma_s is less than 1.
    """
    if fyk < 0:
        raise ValueError(f'fyk={fyk} cannot be less than 0')
    if gamma_s < 1:
        raise ValueError(f'gamma_s={gamma_s} must be larger or equal to 1')

    return fyk / gamma_s


def eps_ud(eps_uk: float, gamma_s: float) -> float:
    """Design value for the ultimate limit strain welding reinforcing steel.

    EN 1992-1-1:2023, 5.2.4-2.

    Args:
        eps_uk (float): Characteristic ultimate limit strain.
        gamma_s (float): Safety coefficient.

    Returns:
        float: Design ultimate strain limit.

    Raises:
        ValueError: If eps_uk is less than 0.
        ValueError: If gamma_s is less than 1.
    """
    if eps_uk < 0:
        raise ValueError(f'eps_uk={eps_uk} must be equal or larger to 0')
    if gamma_s < 1:
        raise ValueError(f'gamma_s={gamma_s} must be larger or equal to 1')

    return eps_uk / gamma_s


def sigma_s(
    eps: float, fy: float, k: float, eps_u: float, Es: float = 200000
) -> float:
    """Compute the stress for welded reinforcing steel in MPa for a given
    strain.

    EN 1992-1-1:2023, 5.2.4.

    Args:
        eps (float): The strain value.
        fy (float): The yielding stress in MPa. Use fyd for the design
            strength, and fyk for the characteristic strength.
        k (float): Curve parameter. Ratio between the ultimate stress and the
            yielding stress. k = 1 for horizontal post-elastic branch without
            strain limit.
        eps_u (float): Ultimate strain at failure. Use eps_ud for the design
            ultimate strain, and eps_uk for the characteristic ultimate strain.

    Keyword Args:
        Es (float): The modulus of elasticity for reinforcing steel.

    Returns:
        float: The nominal stress in MPa.

    Raises:
        ValueError: If eps is less than 0 or larger than eps_uk.
        ValueError: If Es is less or equal to 0.
        ValueError: If fy is less or equal to 0.
        ValueError: If k is less than 1.
        ValueError: If eps_u is less or equal to 0 and k > 1.
    """
    if eps < 0:
        raise ValueError(f'eps={eps} must be larger or equal to 0')
    if eps_u <= 0:
        raise ValueError(f'eps_u={eps_u} must be larger than 0')
    if eps > eps_u and k > 1:
        raise ValueError(
            f'eps={eps} must be equal o less than'
            + 'eps_uk={eps_u} when k is larger than 1'
        )
    if Es <= 0:
        raise ValueError(f'Es={Es} must be larger than 0')
    if fy <= 0:
        raise ValueError(f'fyk={fy} must be larger than 0')
    if k < 0:
        raise ValueError(f'k={k} must be larger than 1')

    eps_y = fy / Es

    # If in elastic area
    if eps <= eps_y:
        return eps * Es

    # If in plastic area
    m = fy * (k - 1) / (eps_u - eps_y)
    return fy + m * (eps - eps_y)


def p_steel_stress_params(
    prestress_class: t.Literal[
        'Y1560',
        'Y1670',
        'Y1770',
        'Y1860',
        'Y1770',
        'Y1860',
        'Y1960',
        'Y2060',
        'Y1030',
        'Y1050',
        'Y1100',
        'Y1230',
    ],
    element: t.Literal['W', 'S', 'B'],
) -> t.Tuple[float, float]:
    """Computes the stress-diagram parameters fp01k and fpk.

    EN 1992-1-1:2023, 5.3.3.

    Args:
        prestress_class (str): Possible values: Y1560, Y1670,
            Y1770, Y1860, Y1770, Y1860, Y1960, Y2060,
            Y1030, Y1050, Y1100 and Y1230
        element (str): Element type, 'W' for Wires, 'S' for Strands, and 'B'
            for Bars.

    Returns:
        Tuple(float, float): With the value of fp01k and fpk in MPa.

    Raises:
        ValueError: If combination of prestress_class and element is not a
            possible value from the range.
    """
    prestress_class = prestress_class.upper().strip()
    element = element.upper().strip()
    data = {
        'W': {
            'Y1570': (1380, 1570),
            'Y1670': (1470, 1670),
            'Y1770': (1550, 1770),
            'Y1860': (1650, 1860),
        },
        'S': {
            'Y1770': (1560, 1770),
            'Y1860': (1640, 1860),
            'Y1960': (1740, 1960),
            'Y2060': (1830, 2060),
        },
        'B': {
            'Y1030': (835, 1030),
            'Y1050': (950, 1050),
            'Y1100': (900, 1100),
            'Y1230': (1080, 1230),
        },
    }

    if element not in data:
        raise ValueError(f'element={element} has not a valid value')

    material = data[element]

    if prestress_class not in material:
        raise ValueError(
            f'prestress_class={prestress_class} has not a valid value'
        )

    return material[prestress_class]


def fpd(fp01k: float, gamma_P: float) -> float:
    """Computes the design value for the prestressing steel stress.

    EN 1992-1-1:2023, 5.3.3.

    Args:
        fp01k (float): The 0.1% proof stress in MPa.
        gamma_P (float): The safety coefficient.

    Returns:
        float: The design value for the design prestressing steel stress in
        MPa.

    Raises:
        ValueError: If fp01k is less than 0.
        ValueError: If gamma_P is less than 1.
    """
    if fp01k < 0:
        raise ValueError(f'fp01k={fp01k} must be larger or equal to 0')
    if gamma_P < 1:
        raise ValueError(f'gamma_P={gamma_P} must be larger or equal than 1')

    return fp01k / gamma_P


def sigma_p(
    eps: float,
    fpy: float,
    fpu: float,
    eps_u: float = 0.035,
    Ep: float = 190000,
) -> float:
    """Computes the stress for prestressing steel as a function of the strain.

    EN 1992-1-1:2023, 5.3.3.

    Args:
        eps (float): Strain value.
        fpy (float): Yielding stress of the steel in MPa. Use fd for design
            stress values, and fp01k for nominal stress values.
        fpu (float): The maximum stress at eps_u in MPa. Use fpd for design
            stress values, fpk for nominal stress values, and fpu == fpy for
            horizontal post-elastic branch without strain limit.

    Keyword Args:
        eps_u (float): Ultimate strain. Use eps_uk = 0.035 for nominal ultimate
            strain, and eps_ud for design ultimate strain.
        Ep (float): Modulus of elasticity of prestressing steel in MPa.

    Raises:
        ValueError: If eps is less than 0 or larger than eps_u.
        ValueError: If fpy is less or equal to 0.
        ValueError: If fpu is less than fpy.
        ValueError: If eps_u is lower or equal to 0.
        ValueError: If _Ep is less or equal to 0.
    """
    if eps < 0:
        raise ValueError(f'eps={eps} must be larger or equal to 0')
    if eps_u <= 0:
        raise ValueError(f'eps_u={eps_u} must be larger than 0')
    if eps_u < eps and fpy != fpu:
        raise ValueError(f'eps_u={eps_u} must be larger than eps={eps}')
    if fpy <= 0:
        raise ValueError(f'fpy={fpy} must be larger than 0')
    if fpu < fpy:
        raise ValueError(f'fpu={fpy} must be larger or equal to fpy={fpy}')
    if Ep <= 0:
        raise ValueError(f'Ep={Ep} must be larger than 0')

    eps_y = fpy / Ep

    # If elastic
    if eps <= eps_y:
        return eps * Ep

    # If plastic
    m = (fpu - fpy) / (eps_u - eps_y)
    return fpy + m * (eps - eps_y)


def reinforcement_duct_props(
    fyk: float,
    ductility_class: t.Literal['A', 'B', 'C'],
) -> t.Dict[str, float]:
    """Return a dict with the minimum characteristic ductility properties for
    reinforcement ductility class.

    EUROCODE 2 1992-1-1:2023, Tab. 5.5.

    Args:
        fyk (float): The characteristic yield strength.
        ductility_class (Literal['A', 'B', 'C']): The reinforcement ductility
            class designation.

    Returns:
        Dict[str, float]: A dict with the characteristik strain value at the
        ultimate stress level (epsuk), and the characteristic ultimate stress
        (ftk).
    """
    duct_props = DUCTILITY_CLASSES.get(ductility_class.upper(), None)
    if duct_props is None:
        raise ValueError(
            'The no properties was found for the provided ductility class '
            f'({ductility_class}).'
        )
    return {
        'epsuk': duct_props['epsuk'],
        'ftk': duct_props['k'] * fyk,
    }
