"""Functions from Section 5 of EN 1992-1-1:2023."""

import math
import typing as t
import numpy as np
import scipy.interpolate

from structuralcodes.codes import mc2010


def fcm(fck: float, delta_f: float = 8.0) -> float:
    """Determines the mean strength of concrete from its characteristic value.

    EN 1992-1-1:2023, Table 5.1

    Args:
        fck (float): is the characteristic compressive strength in MPa.

    Keyword Args:
        delta_s (float): the increment in MPa to compute the mean
            compressive strength

    Returns:
        float: the mean compressive strength in MPa
    """
    return mc2010.fcm(fck, delta_f)


def fctm(fck: float) -> float:
    """Compute the mean concrete tensile strength from the characteristic
    compressive strength.

    EN 1992-1-1:2023, Table 5.1

    Args:
        fck (float): the characteristic compressive strength in MPa.

    Returns:
        float: the mean tensile strength in MPa.
    """
    if abs(fck) <= 50:
        return 0.3 * math.pow(abs(fck), 2 / 3)
    return 1.1 * math.pow(abs(fck), 1 / 3)


def fctk_5(_fctm: float) -> float:
    """Compute the 5% mean concrete tensile strength fractile.

    EN 1992-1-1:2023, Table 5.1

    Args:
        _fctm (float): the mean concrete tensile strength in MPa

    Returns:
        float: the 5% mean concrete tensile strength fractile in MPa
    """
    return abs(_fctm) * 0.7


def fctk_95(_fctm: float) -> float:
    """Compute the 95% mean concrete tensile strength fractile.

    EN 1992-1-1:2023, Table 5.1

    Args:
        _fctm (float): the mean concrete tensile strength in MPa

    Returns:
        float: the 5% mean concrete tensile strength fractile in MPa
    """
    return abs(_fctm) * 1.3


def Ecm(_fcm: float, kE: float = 9500) -> float:
    """Computes the secant modulus between sigma_c=0 and sigma_c=0.4*fcm.

    EN 1992-1-1:2023, Eq. (5.1)

    Args:
        _fcm (float): the mean compressibe strength in MPa

    Keyword Args:
        kE (float): coefficient relating the aggregates used
            in concrete. Default value is 9500, but it can vary
            from 5000 to 13000.

    Returns:
        float: the secant concrete modulus in MPa

    Raises:
        ValueError: if _fcm is less than 0
        ValueError: if kE is not between 5000 and 13000
    """
    if _fcm < 0:
        raise ValueError(f'_fcm={_fcm} cannot be less than 0')
    if kE < 5000:
        raise ValueError(f'kE={kE} cannot be less than 5000')
    if kE > 13000:
        raise ValueError(f'kE={kE} cannot be larger than 13000')

    return kE * math.pow(_fcm, 1 / 3)


def hn(Ac: float, u: float) -> float:
    """Computes the notional size of a given concrete cross-section.

    EN 1992-1-1:2023, Table 5.2

    Args:
        Ac (float): the concrete cross-sectional area in (mm2)
        u (float): the perimeter exposed to drying (mm)

    Returns:
        float: the notional size (mm)

    Raises:
        ValueError: if Ac is less than 0
        ValueError: if u is less than 0
    """
    if Ac < 0:
        raise ValueError(f'Ac={Ac} cannot be less than 0')
    if u < 0:
        raise ValueError(f'u={u} cannot be less than 0')

    return 2 * Ac / u


def A_phi_correction_exp(_hn: float, atm_conditions: str) -> float:
    """Computes the correction exponent for the modification for the phi_50y_t0
    with respect the fck value.

    EN 1992-1-1:2023, Table 5.2

    Args:
        _hn (float): the notional size in mm
        atm_conditions (str): 'dry' or 'humid'

    Returns:
        (float): the correction exponent value

    Raises:
        ValueError: if _hn is less than 100 or greater than 1000
        ValueError: if atm_conditions is not 'dry' or 'humid'
    """
    if _hn < 100:
        raise ValueError(f'_hn={_hn} cannot be less than 100')
    if _hn > 1000:
        raise ValueError(f'_hn={_hn} cannot be larger than 1000')

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
    return interpol(_hn)


def phi_correction_factor(fck: float, A_exponent: float) -> float:
    """Computes the correction factor for the computation of the phi_50y_t0.

    EN 1992-1-1:2023able 5.2

    Args:
        fck (float): characteristic strength of concrete in MPa
        A_exponent (float): the A correction exponent value

    Returns:
        float: the correction factor value

    Raises:
        ValueError: if fck is not between 12 and 100 MPa
        ValueError: if A_exponent is not between 0.64 and 0.82
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
    t0: float, atm_conditions: str, _hn: float, strength_dev_class: str
) -> float:
    """Computes the creep coefficient of plain concrete at 50 years of loading.
    Interpolation is linear between values.

    EN 1992-1-1:2023, Table 5.2

    Args:
        t0 (float): age at loading [days]
        atm_conditions (str): 'dry' or 'humid'
        _hn (float): the notional size in mm
        strength_dev_class (str): 'CS', 'CN' or 'CR'

    Returns:
        float: the creep coefficient

    Raises:
        ValueError: if t0 is less than 1
        ValueError: if atm_conditions is not 'dry' or 'humid'
        ValueError: if _hn is less than 100 or larger than 1000
        ValueError: if strength_dev_class is not 'CS', 'CN' or 'CR'
        ValueError: if combination of t0 and _hn is out of scope
    """
    if t0 < 1:
        raise ValueError(f't0={t0} cannot be less than 1')

    atm_conditions = atm_conditions.lower().strip()
    if atm_conditions not in ('dry', 'humid'):
        raise ValueError(
            f'atm_conditions={atm_conditions} must be "dry" or "humid"'
        )

    if _hn < 100:
        raise ValueError(f'_hn={_hn} must be larger or equal than 100')
    if _hn > 1000:
        raise ValueError(f'_hn={_hn} must be less or equal than 1000')

    strength_dev_class = strength_dev_class.upper().strip()
    if strength_dev_class not in ('CS', 'CN', 'CR'):
        raise ValueError(
            f'strength_dev_class={strength_dev_class} must be'
            + '"CS", "CN" or "CR"'
        )

    if strength_dev_class == 'CS':
        _t = (3, 10, 32, 91, 365)
    elif strength_dev_class == 'CN':
        _t = (1, 7, 28, 91, 365)
    elif strength_dev_class == 'CR':
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
    p = (t0, _hn)

    interp = scipy.interpolate.griddata(grid, values, p, method='linear')
    _phi_50y_t0 = float(interp)

    if math.isnan(_phi_50y_t0) or math.isnan(_phi_50y_t0):
        raise ValueError('Combination of t0, _hn out of scope')

    return _phi_50y_t0


def eps_cs_50y(
    fck_28: float, atm_conditions: str, _hn: float, strength_dev_class: str
) -> float:
    """Computes the nominal total shrinkage in ‰ for concrete after a duration
    of drying of 50 years.

    EN 1992-1-1:2023, Table 5.3

    Args:
        fck_28 (float): characteristic strength at 28 days in MPa
        atm_conditions (str): 'dry' or 'humid'
        _hn (float): the notional size in mm
        strength_dev_class (str): 'CS', 'CN' or 'CR'

    Returns:
        float: the nominal shrinkage value in ‰

    Raises:
        ValueError: if fck_28 is less than 20 MPa or larger than 80 MPa
        ValueError: if atm_conditions is not 'dry' or 'humid'
        ValueError: if _hn is less than 100 or larger than 1000
        ValueError: if strength_dev_class is not 'CS', 'CN' or 'CR'
        ValueError: if combination of fck_28 and _hn is out of scope
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

    if _hn < 100:
        raise ValueError(f'_hn={_hn} must be larger or equal than 100')
    if _hn > 1000:
        raise ValueError(f'_hn={_hn} must be less or equal than 1000')

    strength_dev_class = strength_dev_class.upper().strip()
    if strength_dev_class == 'CS':
        fck_v = (20, 35, 50)
    elif strength_dev_class == 'CN':
        fck_v = (20, 35, 50, 80)
    elif strength_dev_class == 'CR':
        fck_v = (35, 50, 80)
    else:
        raise ValueError(
            f'strength_dev_class={strength_dev_class} '
            + 'must be "CS", "CN" or "CR"'
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
    p = (fck_28, _hn)

    interp = scipy.interpolate.griddata(grid, values, p, method='linear')
    _eps_cs_50y = float(interp)

    if math.isnan(_eps_cs_50y):
        raise ValueError('Combination of fck_28, _hn out of scope')

    return _eps_cs_50y


def eta_cc(fck: float, fck_ref: float = 40) -> float:
    """Computes the factor to measure the difference between the undistributed
    compressibe strength of a cylinder and the effective compressive strength
    in a structural member.

    EN 1992-1-1:2023, Eq. (5.4)

    Args:
        fck (float): the characterisitic compressive strength in MPa
        fck_ref (float, optional): the reference compressive strength MPa

    Returns:
        float: the value of the factor eta_cc

    Raises:
        ValueError: if fck is less than 12 MPa
        ValueError: if fkc_ref is less or equal to 0
    """
    if fck < 12:
        raise ValueError(f'fck={fck} must be larger or equal than 12 MPa')
    if fck_ref <= 0:
        raise ValueError(f'fck_ref={fck_ref} must be larger than 0')

    return min(math.pow(fck_ref / fck, 1 / 3), 1)


def k_tc(t_ref: float, t0: float, strength_dev_class: str) -> float:
    """Computes the factor for considering the effect of high sustained loads
    and of time of loading on concrete compressive strength.

    EN 1992-1-1:2023, Eq. (5.3)

    Args:
        t_ref (float): the reference time in days
        t0 (float): age at loading in days
        strength_dev_class (str): 'CS', 'CN', 'CR', 'slow', 'normal', 'rapid'

    Returns:
        float: the factor value

    Raises:
        ValueError: if t_ref is less than 0
        ValueError: if t0 is less than 0
        ValueError if strength_dev_class is not
            'CS', 'CN', 'CR', 'slow', 'normal', 'rapid'
    """
    if t_ref < 0:
        raise ValueError(f't_ref={t_ref} must be larger than 0')
    if t0 < 0:
        raise ValueError(f't0={t0} must be larger than 0')

    strength_dev_class = strength_dev_class.upper().strip()

    valid_dev_classes = ('cs', 'cn', 'cr', 'slow', 'normal', 'rapid')
    if strength_dev_class not in valid_dev_classes:
        raise ValueError(
            f'strength_dev_class={strength_dev_class}'
            + f'should can only take {valid_dev_classes}'
        )

    if strength_dev_class.lower() in ('CR', 'CN', 'rapid', 'normal') and t_ref <= 28 and t0 > 90:
        return 1

    if strength_dev_class == 'CS' and t_ref <= 56 and t0 > 90:
        return 1

    return 0.85


def fcd(fck: float, _eta_cc: float, _k_tc: float, gamma_C: float) -> float:
    """Computes the value of the design compressive strength of concrete.

    EN 1992-1-1:2023, Eq. (5.3)

    Args:
        fck (float): characteristic compressive strength in MPa
        _eta_cc (float): factor for measuring the difference between
            the undistributed compressive strength of a cylinder and
            the effective compressive strength in the real structural
            member
        _k_tc (float): factor for taking into consideration high
            sustained loads and of time of loading
        gamma_C (float): partial factor of concrete

    Returns:
        float: the design compressive strength of concrete in MPa

    Raises:
        ValueError: if fck is less than 12 MPa
        ValueError if _etc_cc is not between 0 and 1
        ValueError: if gamma_C is less or equal to 0
    """
    if fck < 12:
        raise ValueError(f'fck={fck} must be larger or equal than 12 MPa')
    if _eta_cc < 0 or _eta_cc > 1:
        raise ValueError(f'_eta_cc={_eta_cc} must be between 0 and 1')
    if gamma_C <= 0:
        raise ValueError(f'gamma_C={gamma_C} must be larger than 0')

    return _eta_cc * _k_tc * fck / gamma_C


def k_tt(t_ref: float, strength_dev_class: str) -> float:
    """Computes the factor for considering the effect of high sustained loads
    and of time of loading on concrete tensile strength.

    EN 1992-1-1:2023, Eq. (5.5)

    Args:
        t_ref (float): the reference time in days
        strength_dev_class (str): 'CS', 'CN' or 'CR'

    Returns:
        float: the factor value

    Raises:
        ValueError: if t_ref is less than 0
        ValueError if strength_dev_class is not 'CS', 'CN', or 'CR'
    """
    if t_ref < 0:
        raise ValueError(f't_ref={t_ref} must be larger than 0')

    strength_dev_class = strength_dev_class.upper().strip()

    if strength_dev_class not in ('CS', 'CN', 'CR'):
        raise ValueError(
            f'strength_dev_class={strength_dev_class}'
            + 'should can only take "CS", "CN" or "CR" as values'
        )

    if strength_dev_class in ('CR', 'CN') and t_ref <= 28:
        return 0.8

    if strength_dev_class == 'CS' and t_ref <= 56:
        return 0.8

    return 0.7


def fctd(_fctk_5: float, _k_tt: float, gamma_C: float) -> float:
    """Computes the value of the design tensile strength of concrete.

    EN 1992-1-1:2023, Eq. (5.5)

    Args:
        fctk_5 (float): the 5% mean concrete tensile strength fractile in MPa
        _k_tt (float): the factor for considering the effect of high sustained
            loads and of time of loading on concrete tensile strength
        gamma_C (float): partial factor of concrete

    Returns:
        float: the design tensile strength of concrete in MPa

    Raises:
        ValueError: if fctk_5 is less than 0
        ValueError: gamma_C is less than 1
    """
    if _fctk_5 < 0:
        raise ValueError(f'fctk_5={_fctk_5} must be larger or equal to 0')
    if gamma_C < 1:
        raise ValueError(f'gamma_C={gamma_C} must be larger or equal to 1')

    return _k_tt * _fctk_5 / gamma_C


def eps_c1(_fcm: float) -> float:
    """Computes the strain at maximum compressive strength of concrete (fcm)
    for the Sargin constitutive law.

    EN 1992-1-1:2023, Eq. (5.9)

    Args:
        _fcm (float): the mean strength of concrete in MPa

    Returns:
        float: the strain at maximum compressive strength of concrete

    Raises:
        ValueError: if _fcm is less than 12+8MPa
    """
    if _fcm < 20:
        raise ValueError(f'_fcm={_fcm} must be larger or equal to 12+8MPa')

    return min(0.7 * math.pow(_fcm, 1 / 3), 2.8) / 1000


def eps_cu1(_fcm: float) -> float:
    """Computes the strain at concrete failure of concrete.

    EN 1992-1-1:2023, Eq. (5.10)

    Args:
        _fcm (float): the mean strength of concrete in MPa

    Returns:
        float: the maximum strength at failure of concrete

    Raises:
        ValueError: if _fcm is less than 12+8MPa
    """
    if _fcm < 20:
        raise ValueError(f'_fcm={_fcm} must be larger or equal to 12+8MPa')

    return min(2.8 + 14 * (1 - _fcm / 108) ** 4, 3.5) / 1000


def sigma_c(
    _Ecm: float, _fcm: float, _eps_c: float, _eps_c1: float, _eps_cu1: float
) -> float:
    """Computes the compressive stress of concrete given a strain eps_c under
    short term uniaxial compression.

    EN 1992-1-1:2023, Eq. (5.6)

    Args:
        _Ecm (float): the secant modulus between sigma_c=0 and
            sigma_c=0.4*fcm in MPa
        _fcm (float): the mean compressive strength of concrete
            in MPa
        _eps_c (float): the strain of concrete
        _eps_c1 (float): the strain of concrete at stress _fcm
        _eps_cu1 (float): the strain at failure of concrete

    Returns:
        float: the compressive stress of concrete in MPa

    Raises:
        ValueError: if _Ecm is less or equal to 0
        ValueError: if _fcm is less than 12+8MPa
        ValueError: if _eps_c is less than 0
        ValueError: if _eps_c1 is less or equal to 0
        ValueError: if _eps_cu1 is less or equal than 0
        ValueError: if _eps_c is larger than _eps_cu1
    """
    if _Ecm <= 0:
        raise ValueError(f'_Ecm={_Ecm} must be larger than 0')
    if _fcm < 20:
        raise ValueError(f'_fcm={_fcm} must be larger or equal to 12+8MPa')
    if _eps_c < 0:
        raise ValueError(f'_eps_c={_eps_c} must be larger or equal to 0')
    if _eps_c1 <= 0:
        raise ValueError(f'_eps_c1={_eps_c1} must be larger than 0')
    if _eps_cu1 < 0:
        raise ValueError(f'_eps_cu1={_eps_cu1} must be larger or equal to 0')
    if _eps_c > _eps_cu1:
        raise ValueError(
            f'Current strain value _eps_c={_eps_c} is larger than'
            + f'the strain at failure _eps_cu1={_eps_cu1}'
        )

    k = 1.05 * _Ecm * _eps_c1 / _fcm
    eta = _eps_c / _eps_c1

    n = k * eta - eta**2
    d = 1 + (k - 2) * eta

    return n / d * _fcm


def weight_c(concrete_type: str) -> float:
    """Returns the mean unit weight of concrete in kN/m3.

    EN 1992-1-1:2023, 5.1.6-5

    Args:
        concrete_type (str):
            'nc' for normal concrete
            'npc' for normal plain concrete

    Returns:
        float: mean unit weight in kN/m3

    Raises:
        ValueError: if concrete_type is not 'nc' or 'npc'
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

    EN 1992-1-1:2023, 5.1.6-6

    Returns:
        float: the linear coefficient of thermal expansion in 1/Cº
            for concrete
    """
    return 10 * 10e-6


def r_steel_stress_strain_params(
    ductility_class: str,
) -> t.Tuple[float, float]:
    """Returns the properties that define the stress-strain diagram for
    reinforced steel k and eps_uk.

    EN 1992-1-1:2023, Table 5.5

    Args:
        ductility_class (str): 'A', 'B' or 'C'

    Returns:
        Tuple(float, float): with the characteristic value k
            and the characteristic strain at maximum force eps_uk

    Raises:
        ValueError: if ductility class is not 'A', 'B' or 'C'
    """
    ductility_class = ductility_class.upper().strip()

    if ductility_class == 'A':
        return (1.05, 0.025)
    if ductility_class == 'B':
        return (1.08, 0.05)
    if ductility_class == 'C':
        return (1.25, 0.075)
    raise ValueError(
        f'ductility_class={ductility_class} can only '
        + 'take as values "A", "B" or "C'
    )


def Es() -> float:
    """Returns the value of the modulus of elasticity for weldable reinforcing
    steel.

    EN 1992-1-1:2023, 5.2.4-3

    Returns:
        float: modulus of elasticity in MPa
    """
    return 200000


def alpha_s_th() -> float:
    """Returns the linear coefficient of thermal expansion in 1/Cº for weldable
    reinforced steel.

    EN 1992-1-1:2023, 5.2.4-5

    Returns:
        float: the linear coefficient of thermal expansion in 1/Cº
            for weldable reinforce steel
    """
    return 10 * 10e-6


def weight_s() -> float:
    """Returns the mean unit weight of reinforced steel for the purposes of
    design in kN/m3.

    EN 1992-1-1:2023.2.4-4

    Returns:
        float: the mean unit weight in kN/m3
    """
    return 78.5


def fyd(fyk: float, gamma_S: float) -> float:
    """Design value for the yielding stress for welding reinforcing steel.

    EN 1992-1-1:2023, Eq (5.11)

    Args:
        fyk (float): characteristic yield stress for the steel in MPa
        gamma_S (float): safety coefficient

    Returns:
        float: design yielding stress for steel in MPa

    Raises:
        ValueError: if fyk is less than 0
        ValueError: if gamma_S is less or equal to 0
    """
    if fyk < 0:
        raise ValueError(f'fyk={fyk} cannot be less than 0')
    if gamma_S <= 0:
        raise ValueError(f'gamma_S={gamma_S} must be larger than 0')

    return fyk / gamma_S


def eps_ud(eps_uk: float, gamma_S: float) -> float:
    """Design value for the ultimate limit strain welding reinforcing steel.

    EN 1992-1-1:2023, 5.2.4-2

    Args:
        eps_uk (float): characteristic ultimate limit
            strain
        gamma_S (float): safety coefficient

    Returns:
        float: design ultimate strain limit

    Raises:
        ValueError: if eps_uk is less than 0
        ValueError: if gamma_S is less or equal to 0
    """
    if eps_uk < 0:
        raise ValueError(f'eps_uk={eps_uk} must be equal or larger to 0')
    if gamma_S < 0:
        raise ValueError(f'gamma_S={gamma_S} must be larger than 0')

    return eps_uk / gamma_S


def sigma_s(
    eps: float, fy: float, k: float, eps_u: float, _Es: float = 200000
) -> float:
    """Compute the stress for welded reinforcing steel in MPa for a given
    strain.

    EN 1992-1-1:2023, 5.2.4

    Args:
        eps (float): the strain value
        fy (float): the yielding stress in MPa.
            fyd: for the design strength
            fyk: for the characteristic strength
        k (float): curve parameter. Ratio between the ultimate
            stress and the yielding stress.
            k = 1 for horizontal post-elastic branch without
            strain limit
        eps_u (float): ultimate strain at failure
            eps_ud: for the design ultimate strain
            eps_uk: for the characteristic ultimate strain

    Keyword Args:
        _Es (float): the modulus of elasticity for reinforcing steel
    Returns:
        float: stre nominal stress in MPa

    Raises:
        ValueError: if eps is less than 0 or larger than eps_uk
        ValueError: if _Es is less or equal to 0
        ValueError: if fy is less or equal to 0
        ValueError: if k is less than 1
        ValueError: if eps_u is less or equal to 0 and k > 1
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
    if _Es <= 0:
        raise ValueError(f'_Es={_Es} must be larger than 0')
    if fy <= 0:
        raise ValueError(f'fyk={fy} must be larger than 0')
    if k < 0:
        raise ValueError(f'k={k} must be larger than 1')

    eps_y = fy / _Es

    # If in elastic area
    if eps <= eps_y:
        return eps * _Es

    # If in plastic area
    m = fy * (k - 1) / (eps_u - eps_y)
    return fy + m * (eps - eps_y)


def p_steel_stress_params(
    prestress_class: str, element: str
) -> t.Tuple[float, float]:
    """Computes the stress-diagram parameters fp01k and fpk.

    EN 1992-1-1:2023, 5.3.3

    Args:
        prestress_class (str): possible values: Y1560, 1670,
            Y1770, Y1860, Y1770, Y1860, Y1960, Y2060,
            Y1030, Y1050, Y1100 and Y1230
        element (str):
            'W' for Wires
            'S' for Strands
            'B' for Bars

    Returns:
        Tuple(float, float): with the value of fp01k and fpk in MPa

    Raises:
        ValueError: if combination of prestress_class and element
            is not a possible value from the range
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

    EN 1992-1-1:2023, 5.3.3

    Args:
        fp01k (float): the 0.1% proof stress in MPa
        gamma_P (float): the safety coefficient

    Returns:
        float: the design value for the design prestressing
            steel stress in MPa

    Raises:
        ValueError: if fp01k is less than 0
        ValueError: if gamma_P is less than 1
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

    EN 1992-1-1:2023, 5.3.3

    Args:
        eps (float): strain value
        fpy (float): yielding stress of the steel in MPa
            fd for design stress values
            fp01k for nominal stress values
        fpu (float): the maximum stress at eps_u in MPa
            fpd: for design stress values
            fpk: for nominal stress values
            fpu == fpy for horizontal post-elastic branch without
            strain limit.

    Keyword Args:
        eps_u (float): ultimate strain.
            eps_uk = 0.035 for nominal ultimate strain
            eps_ud for design ultimate strain
        Ep (float): modulus of elasticity of prestressing steel
            in MPa

    Raises:
        ValueError: if eps is less than 0 or larger than eps_u
        ValueError: if fpy is less or equal to 0
        ValueError: if fpu is less than fpy
        ValueError: if eps_u is lower or equal to 0
        ValueError: if _Ep is less or equal to 0
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
