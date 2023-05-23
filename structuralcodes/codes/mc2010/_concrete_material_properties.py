"""A collection of material properties for concrete"""
import math
import numpy as np

_EC0 = 21.5e3  # MPa


def fcm(fck: float, delta_f: float = 8.0) -> float:
    """Compute the mean concrete compressive strength from the characteristic
    strength.

    fib Model Code 2010, Eq. (5.1-1)

    Args:
        fck (float): The characteristic compressive strength in MPa.

    Keyword Args:
        delta_f (float): The difference between the mean and the
        characteristic strength.

    Returns:
        float: The mean compressive strength in MPa.
    """
    return abs(fck) + abs(delta_f)


def fctm(fck: float) -> float:
    """Compute the mean concrete tensile strength from the characteristic
    compressive strength.

    fib Model Code 2010, Eqs. (5.1-3a) and (5.1-3b)

    Args:
        fck (float): The characteristic compressive strength in MPa.

    Returns:
        float: The mean tensile strength in MPa.
    """
    if abs(fck) <= 50:
        return 0.3 * math.pow(abs(fck), 2 / 3)
    return 2.12 * math.log(1 + 0.1 * fcm(fck))


def fctkmin(_fctm: float) -> float:
    """Compute the lower bound value of the characteristic tensile strength
    from the mean tensile strength.

    fib Model Code 2010, Eq. (5.1-4)

    Args:
        _fctm (float): The mean tensile strength in MPa.

    Returns:
        float: Lower bound of the characteristic tensile strength in MPa.
    """
    return 0.7 * _fctm


def fctkmax(_fctm: float) -> float:
    """Compute the upper bound value of the characteristic tensile strength
    from the mean tensile strength.

    fib Model Code 2010, Eq. (5.1-5)

    Args:
        _fctm (float): The mean tensile strength in MPa.

    Returns:
        float: Upper bound of the characteristic tensile strength in MPa.
    """
    return 1.3 * _fctm


def Gf(fck: float) -> float:
    """Compute tensile fracture energy from characteristic compressive
    strength.

    fib Model Code 2010, Eq. (5.1-9)

    Args:
        fck (float): The characteristic compressive strength in MPa.

    Returns:
        float: The tensile fracture energy in N/m.
    """
    return 73 * fcm(fck) ** 0.18


def _check_agg_types(agg_type: str) -> None:
    """Check if the input for the aggregate type is included in the
        model formulation of the modulus of elasticity (5.1.7.2).

    Defined in fib MC 2010 (2013), table 5.1-6.

    args:
        agg_type (str): The type of aggregate used in the concrete.
            The following values can be chosen:
            ['Basalt', 'Quartzite', 'Limestone', 'Sandstone'].

    returns:
        Raises a ValueError if aggregate type is outside of specified
            values.
    """
    AGG_TYPES = ['basalt', 'quartzite', 'limestone', 'sandstone']
    if agg_type.lower() not in AGG_TYPES:
        raise ValueError(
            f"The specified cement type {agg_type} is not"
            " implemented in the creep and shrinkage laws of the fib "
            "Model Code 2010. It can only be one of these: "
            f"{list(AGG_TYPES)}."
        )


def _check_cem_strength_class(cem_class: str) -> None:
    """Check if a cement strength class is used for which required model
        parameters are available.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-73 and table 5.1-12.

    args:
        cem_class (str): The cement strength class that is used.
            The choices are:
                '32.5 N',
                '32.5 R', '42.5 N',
                '42.5 R', '52.5 N', '52.5 R'.

    returns:
        Raises a ValueError if a non-implemented cement strength class
         is used.
    """
    CEM_CLASSES = ["32.5 N", "32.5 R", "42.5 N", "42.5 R", "52.5 N", "52.5 R"]
    if cem_class.upper() not in CEM_CLASSES:
        raise ValueError(
            "Unknown cem_class used. Please choose one of "
            f"the following {list(CEM_CLASSES)}."
        )


def _get_Ecmod_coeffs(cem_class: str, fcm: float, agg_type: str) -> dict:
    """Get the coefficients to calculate the (time-dependent) ccncrete
        modulus of elasticity belonging to the specified cement strength
        class, fcm and aggregate type.

    Defined in fib Model Code 2010 (2013), tables 5.1-6 and 5.1-9.

    args:
        cem_class (str): The cement strength class that is used.
        fcm (float): The mean compressive strength of the concrete in
            MPa.
        agg_type (str): The type of aggregate used in the concrete.

    returns:
        A dictionary with the model parameters S and ALPHA_E.
    """
    _check_cem_strength_class(cem_class)
    _check_agg_types(agg_type)
    if cem_class.upper() in ["32.5 R", "42.5 N"]:
        if fcm > 60:
            S = 0.20
        else:
            S = 0.25
    elif cem_class.upper() in ["42.5 R", "52.5 N", "52.5 R"]:
        S = 0.20
    elif cem_class.upper() in ["32.5 N"]:
        if fcm > 60:
            S = 0.20
        else:
            S = 0.38
    if agg_type.lower() == 'basalt':
        ALPHA_E = 1.2
    elif agg_type.lower() == 'quartzite':
        ALPHA_E = 1.0
    elif agg_type.lower() == 'limestone':
        ALPHA_E = 0.9
    elif agg_type.lower() == 'sandstone':
        ALPHA_E = 0.7

    return {'S': S, 'alpha_e': ALPHA_E}


def _calc_E28(
    fc: float, TABULAR_VALUES: dict, fc_value_type: str = "char"
) -> float:
    """Calculate the modulus of elasticity for normal weight concrete
        at 28 days.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-21.

    args:
        fc (float): The mean or characteristic value of the compressive
            strength of the concrete in MPa.
        TABULAR_VALUES (dict): Material constants for different cement
            types and aggregate types as given by the fib Model Code
            2010.
        fc_value_type (str): The specified value of fc is "mean" or "char"
            (default).

    returns:
        float: The modulus of elasticity for normal weight concrete
            at 28 days in MPa.
    """
    if fc_value_type == "char":
        fcm = fc + 8
    elif fc_value_type == "mean":
        fcm = fc
    else:
        raise ValueError(
            "The value_type of fc can only be one of the "
            "following: 'mean' or 'char'."
        )

    return _EC0 * TABULAR_VALUES['alpha_e'] * (fcm / 10) ** (1 / 3)


def _calc_beta_cc(time: float, TABULAR_VALUES: dict) -> float:
    """Calculate multiplication factor beta_cc, used to determine the
        modulus of elasticity at an arbitrary time.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-51.

    args:
        time (float): The time in days at which the modulus of
            elasticity is to be determined.
        TABULAR_VALUES (dict): Material constants for different cement
            types and aggregate types as given by the fib Model Code
            2010.

    returns:
        float: Multiplication factor beta_cc.
    """
    return np.exp(TABULAR_VALUES['S'] * (1 - np.sqrt(28 / time)))


def _calc_beta_e(beta_cc: float) -> float:
    """Calculate multiplication factor beta_e, used to determine the
        modulus of elasticity at an arbitrary time.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-57.

    args:
        beta_cc (float): multiplication factor as defined in the fib Model
            Code 2010 (2013), Eq. 5.1-51.

    returns:
        float: Multiplication factor beta_e.
    """
    return np.sqrt(beta_cc)


def _calc_E(
    time: float, fc: float, TABULAR_VALUES: dict, fc_value_type: str = "char"
) -> float:
    """Calculate the modulus of elasticity for normal weight concrete
        at time 'time' (not 28 days).

    Defined in fib Model Code 2010 (2013), Eq. 5.1-56.

    args:
        time (float): The time in days at which the modulus of
            elasticity is to be determined.
        fc (float): The mean or characteristic value of the compressive
            strength of the concrete in MPa.
        TABULAR_VALUES (dict): Material constants for different cement
            types and aggregate types as given by the fib Model Code
            2010.
        fc_value_type (str): The specified value of fc is "mean" or "char"
            (default).

    returns:
        float: The modulus of elasticity for normal weight concrete at
            time 'time' (not 28 days) in MPa.
    """
    beta_cc = _calc_beta_cc(time, TABULAR_VALUES)
    beta_e = _calc_beta_e(beta_cc)
    return beta_e * _calc_E28(fc, TABULAR_VALUES, fc_value_type)
