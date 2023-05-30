"""Functions from Section 5 of FprEN 1992-1-1:2022"""

import math
import numpy as np
import scipy.interpolate

from structuralcodes.codes import mc2010


# 5.1.3 Strength


def fcm(fck: float, delta_f: float = 8.0) -> float:
    """Determines the mean strength of concrete from its characteristic
    value

    FprEN 1992-1-1, Table 5.1

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
    compressive strength

    FprEN 1992-1-1, Table 5.1

    Args:
        fck (float): the characteristic compressive strength in MPa.

    Returns:
        float: the mean tensile strength in MPa.
    """
    if abs(fck) <= 50:
        return 0.3 * math.pow(abs(fck), 2 / 3)
    return 1.1 * math.pow(abs(fck), 1 / 3)


def fctk_5(_fctm: float) -> float:
    """Compute the 5% mean concrete tensile strength fractile

    FprEN 1992-1-1, Table 5.1

    Args:
        _fctm (float): the mean concrete tensile strength in MPa

    Returns:
        float: the 5% mean concrete tensile strength fractile in MPa
    """
    return abs(_fctm) * 0.7


def fctk_95(_fctm: float) -> float:
    """Compute the 95% mean concrete tensile strength fractile

    FprEN 1992-1-1, Table 5.1

    Args:
        _fctm (float): the mean concrete tensile strength in MPa

    Returns:
        float: the 5% mean concrete tensile strength fractile in MPa
    """
    return abs(_fctm) * 1.3


# 5.1.4 Elastic deformation


def Ecm(_fcm: float, kE: float = 9500) -> float:
    """Computes the secant modulus between sigma_c=0 and
    sigma_c=0.4*fcm

    FprEN 1992-1-1, Eq. (5.1)

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
    """Computes the notional size of a given concrete
    cross-section

    FprEN 1992-1-1, Table 5.2

    Args:
        Ac (float): the concrete cross-sectional area in (mm)
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
    """Computes the correction exponent for the modification for
    the phi_50y_t0 with respect the fck value:

    FprEN 1992-1-1, Table 5.2

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
    val = interpol(_hn)
    return val


def phi_correction_factor(fck: float, A_exponent: float) -> float:
    """Computes the correction factor for the computation of the
    phi_50y_t0

    FprEN 1992-1-1, Table 5.2

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
    t0: float, atm_conditions: str, _hn: float, concrete_class: str
) -> float:
    """Computes the creep coefficient of plain concrete at 50 years
    of loading. Interpolation is lineal between values

    FprEN 1992-1-1, Table 5.2

    Args:
        t0 (float): age at loading [days]
        atm_conditions (str): 'dry' or 'humid'
        _hn (float): the notional size in mm
        concrete_class (str): 'CS', 'CN' or 'CR'

    Returns:
        float: the creep coefficient

    Raises:
        ValueError: if t0 is less than 0
        ValueError: if atm_conditions is not 'dry' or 'humid'
        ValueError: if _hn is less than 100 or larger than 1000
        ValueError: if concrete_class is not 'CS', 'CN' or 'CR'
        ValueError: if combination of t0 and _hn is out of scope
    """
    if t0 < 0:
        raise ValueError(f't0={t0} cannot be less than 0')

    atm_conditions = atm_conditions.lower().strip()
    if atm_conditions not in ('dry', 'humid'):
        raise ValueError(
            f'atm_conditions={atm_conditions} must be "dry" or "humid"'
        )

    if _hn < 100:
        raise ValueError(f'_hn={_hn} must be larger or equal than 100')
    if _hn > 1000:
        raise ValueError(f'_hn={_hn} must be less or equal than 1000')

    concrete_class = concrete_class.upper().strip()
    if concrete_class not in ('CS', 'CN', 'CR'):
        raise ValueError(
            f'concrete_class={concrete_class} must be "CS", "CN" or "CR"'
        )

    if concrete_class == 'CS':
        t = (3, 10, 32, 91, 365)
    elif concrete_class == 'CN':
        t = (1, 7, 28, 91, 365)
    elif concrete_class == 'CR':
        t = (1, 3, 23, 91, 365)

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

    grid = np.array(np.meshgrid(t, h_v)).T.reshape(-1, 2)
    p = (t0, _hn)

    interp = scipy.interpolate.griddata(grid, values, p, method='linear')
    _phi_50y_t0 = float(interp)

    if math.isnan(_phi_50y_t0) or math.isnan(_phi_50y_t0):
        raise ValueError('Combination of t0, _hn out of scope')

    return _phi_50y_t0


def eps_cs_50y(
    fck_28: float, atm_conditions: str, _hn: float, concrete_class: str
) -> float:
    """Computes the nominal total shrinkage in ‰ for concrete after
    a duration of drying of 50 years

    FprEN 1992-1-1, Table 5.3

    Args:
        fck_28 (float): characteristic strngth at 28 days in MPa
        atm_conditions (str): 'dry' or 'humid'
        _hn (float): the notional size in mm
        concrete_class (str): 'CS', 'CN' or 'CR'

    Returns:
        float: the nominal shrinkage value in ‰

    Raises:
        ValueError: if fck_28 is less than 20 MPa or larger than 80 MPa
        ValueError: if atm_conditions is not 'dry' or 'humid'
        ValueError: if _hn is less than 100 or larger than 1000
        ValueError: if concrete_class is not 'CS', 'CN' or 'CR'
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

    concrete_class = concrete_class.upper().strip()
    if concrete_class == 'CS':
        fck_v = (20, 35, 50)
    elif concrete_class == 'CN':
        fck_v = (20, 35, 50, 80)
    elif concrete_class == 'CR':
        fck_v = (35, 50, 80)
    else:
        raise ValueError(
            f'concrete_class={concrete_class} must be "CS", "CN" or "CR"'
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
    values = data.get(atm_conditions).get(concrete_class)

    grid = np.array(np.meshgrid(fck_v, h_v)).T.reshape(-1, 2)
    p = (fck_28, _hn)

    interp = scipy.interpolate.griddata(grid, values, p, method='linear')
    _eps_cs_50y = float(interp)

    if math.isnan(_eps_cs_50y) or math.isnan(_eps_cs_50y):
        raise ValueError('Combination of fck_28, _hn out of scope')

    return _eps_cs_50y


def eta_cc(fck: float, fck_ref: float = 40) -> float:
    """Computes the factor to measure the difference between the undistributed
    compressibe strength of a cylinder and the effective compressive strength
    in a structural member

    FprEN 1992-1-1, Eq. (5.4)

    Args:
        fck (float): the characterisitic compressive strength in MPa
        fck_ref (float, optional): the reference compressive strength MPa

    Returns:
        float: the value of the factor eta_cc

    Raises:
        ValueError: if fck is less or equal to 0
        ValueError: if fkc_ref is less or equal to 0
    """
    if fck <= 0:
        raise ValueError(f'fck={fck} must be larger than 0')
    if fck_ref <= 0:
        raise ValueError(f'fck_ref={fck_ref} must be larger than 0')

    return min(math.pow(fck_ref / fck, 1 / 3), 1)


def k_tc(t_ref: float, t0: float, concrete_class: str) -> float:
    """Computes the factor for considering the effect of high sustained
    loads and of time of loading on concrete compressive strength

    FprEN 1992-1-1, Eq. (5.3)

    Args:
        t_ref (float): the reference time in days
        t0 (float): age at loading in days
        concrete_class (str): 'CS', 'CN' or 'CR'

    Returns:
        float: the factor value

    Raises:
        ValueError: if t_ref is less than 0
        ValueError: if t0 is less than 0
        ValueError if concrete_class is not 'CS', 'CN', or 'CR'
    """
    concrete_class = concrete_class.upper().strip()

    if concrete_class not in ('CS', 'CN', 'CR'):
        raise ValueError(
            f'concrete_class={concrete_class}'
            + 'should can only take "CS", "CN" or "CR" as values'
        )

    if concrete_class in ('CR', 'CN') and t_ref <= 28 and t0 <= 90:
        return 1

    if concrete_class == 'CS' and t_ref <= 56 and t0 <= 90:
        return 1

    return 0.85
