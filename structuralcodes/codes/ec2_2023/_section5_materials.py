"""Functions from Section 5 of FprEN 1992-1-1:2022"""

import math
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
