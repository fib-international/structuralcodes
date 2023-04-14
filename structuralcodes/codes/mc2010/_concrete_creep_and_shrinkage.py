"""Creep and shrinkage models MC2010 for concrete

[1] fib Model Code for Concrete Structures 2010 (2013).

TODO:
 - Add the light-weight aggregate formulas
 - Include the additional creep and shrinkage terms for other temperatures (see 5.1.10.7 of [1])
"""

import numpy as np

from structuralcodes.codes.mc2010 import _concrete_material_properties


def _check_fcm(fcm: float) -> None:
    """Check if the mean compressive strength is in the range of
        applicability of the creep and shrinkage models:

    Defined in fib Model Code 2010 (2013), section 5.1.9.4.2.

    args:
        fcm (float): The mean compressive strength of the concrete in
            MPa.

    returns:
        Raises a ValueError if the mean compressive strength is outside
            of the specified range.
    """
    if fcm < 20 or fcm > 130:
        raise ValueError("The specified mean compressive strength is "
            "outside of the range of applicability for the creep and "
            "shrinkage laws given in the fib Model Code 2010. The mean"
            " compressive strength has to be within the range of "
            "20-130 MPa. Current compressive strength is"
            f": {fcm} MPa.")


def _check_initial_stress(sigma: float, fcm: float) -> None:
    """Check if the initial compressive stress (based on load at t0)
        is within the range of applicability.

    Defined in fib Model Code 2010 (2013), section 5.1.9.4.2.

    args:
        sigma (float): The compressive stress applied to the concrete
            at t0 in MPa.
        fcm (float): The mean compressive strength of the concrete in
            MPa. Note that its value is non-negative.

    returns:
        Prints a warning if the initial compressive stress is
            greater than 0.4*fcm or raises an ValueError if the
            compressive stress is greater than 0.6*fcm.
    """
    if abs(sigma) > 0.6*fcm:
        raise ValueError("The stress level exceeds the range of application."
            "Maximum allowable stress is 0.6*fcm. Current stress level "
            "is {round(abs(sigma)/fcm, 3)}*fcm.")
    elif abs(sigma) > 0.4*fcm:
        print("WARNING: Initial stress is too high to consider the "
            "concrete as an aging linear visco-elastic material: "
            f"sigma = {round(abs(sigma)/fcm,3)}*fcm > 0.4*fcm. Nonlinear"
            " creep calculations are performed according to subclause "
            "5.1.9.4.3 (d) of the fib Model Code 2010 to account for "
            "large compressive stresses.")


def _check_age_at_loading(t0: float) -> None:
    """Check if the age of the concrete is greater than the minimum
        concrete age.

    Defined in fib Model Code 2010 (2013), section 5.1.9.4.2.

    args:
        t0 (float): The age of the concrete in days at which the
            loading is applied.

    returns:
        Raises a ValueError if the age of the concrete is too low.
    """
    if t0 < 1:
        raise ValueError("The load is applied too soon to the concrete"
        " in order to calculate the creep and shrinkage behaviour "
        "according to the fib Model Code 2010. The minimum age of the "
        f"concrete is 1 day, whereas according to the input the load is"
        " applied after {t0} day.")


def _check_RH(rh: float) -> None:
    """Check if the given relative humidity is within the range of
        applicability.

    Defined in fib Model Code 2010 (2013), section 5.1.9.4.2.

    args:
        rh (float): The relative humidity of the environment.

    returns:
        Raises a ValueError if the relative humidity is outside of the
            range of applicability.
    """
    if (rh < 0.4 or rh > 1) and (rh < 40 or rh > 100):
        raise ValueError("The specified relative humidity is outside "
            "of the range of applicability to calculate the creep and "
            "shrinkage according to the fib Model Code 2010. The "
            "relative humidity has to be within the range of 0.4-1.0 or "
            f"40-100%. Currently rh={rh}.")


def _check_env_temp(T: float) -> None:
    """Check if the given environmental temperature is within the range
        of applicability.

    Defined in fib Model Code 2010 (2013), section 5.1.9.4.2.

    args:
        T (float): The environmental temperature in degrees Celcius.

    returns:
        Prints a warning if the applied environmental temperature is
            outside of the range of applicability.
    """
    if T < 5 or T > 30:
        print("WARNING: The given environmental temperature is outside"
        " of the applicable range of 5-30 degrees Celcius for the"
        " creep and shrinkage calculations according to the fib Model "
        f"Code 2010, T={T} degrees Celcius. Creep and shrinkage will"
        " be calculated according to subclause 5.1.10 of the fib Model"
        " Code 2010.")


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

    CEM_CLASSES = ["32.5 N", "32.5 R", "42.5 N", "42.5 R", "52.5 N",
                   "52.5 R"]
    if cem_class.upper() not in CEM_CLASSES:
        raise ValueError("Unknown cem_class used. Please choose one of "
            f"the following {list(CEM_CLASSES)}.")


def _get_creep_shrinkage_coeffs(cem_class: str) -> dict:
    """Get shrinkage and creep model coefficients belonging to the
        specified cement strength class.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-73 and table 5.1-12.

    args:
        cem_class (str): The cement strength class that is used.

    returns:
        A dictionary with the model parameters ALPHA, ALPHA_AS,
         ALPHA_DS1, and ALPHA_DS2.
    """
    _check_cem_strength_class(cem_class)
    if cem_class.upper() in ["32.5 R", "42.5 N"]:
        ALPHA = 0
        ALPHA_AS= 700
        ALPHA_DS1 = 4
        ALPHA_DS2 = 0.012
    elif cem_class.upper() in ["42.5 R", "52.5 N", "52.5 R"]:
        ALPHA = 1
        ALPHA_AS = 600
        ALPHA_DS1 = 6
        ALPHA_DS2 = 0.012
    elif cem_class.upper() in ["32.5 N"]:
        ALPHA = -1
        ALPHA_AS = 800
        ALPHA_DS1 = 3
        ALPHA_DS2 = 0.013

    return {'alpha': ALPHA, 'alpha_as': ALPHA_AS, 'alpha_ds1': ALPHA_DS1,
            'alpha_ds2': ALPHA_DS2}


def _calc_temp_corr_age(t0: float, T_cur: float) -> float:
    """Calculate the temperature corrected concrete age in days at t0.

    Defined in fib Model Code 2010 (2013). Eq. 5.1-85 (only for a single
     time value input, as required in Eq. 5.1-73).

    args:
        t0 (float): The age of the concrete in days at which the
            loading is applied.
        T_cur (float): The temperature of the environment during curing
            in degrees Celcius.

    returns:
        float: The temperature corrected age of the concrete in days
            at loading.
    """
    _check_age_at_loading(t0)
    dt = np.ones(t0)
    T_dt = np.ones(t0) * T_cur
    return np.sum(dt * np.exp(13.65 - (4000 / (273 + T_dt))))


def _calc_modified_t0(
        t0: float, T_cur: float, TABULAR_VALUES: dict) -> float:
    """Calculate the modified age at loading (t0) to account for the
        effect of the type of cement and curing temperature on the
        degree of hydration and - in turn - on creep.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-73.

    args:
        t0 (float): The age of the concrete in days at which the
            loading is applied.
        T_cur (float): The temperature of the environment during curing
            in degrees Celcius.
        TABULAR_VALUES (dict): Material constants for different cement
            types and aggregate types as given by the fib Model Code
            2010.
    returns:
        float: The temperature corrected age of the concrete in days at
            loading, accounting for the effect of the cement type.
            For slow hardening concrete, the creep coefficient is
            increased due to the lower modified age at loading.
    """
    _check_age_at_loading(t0)
    t0T = _calc_temp_corr_age(t0, T_cur)
    ALPHA = TABULAR_VALUES['alpha']
    return max(t0T * ((9/(2 + t0T**1.2)) + 1)**ALPHA, 0.5)


def _calc_notional_drying_shrinkage(
        fcm: float, TABULAR_VALUES:dict) -> float:
    """Calculate the notional drying shrinkage.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-80.

    args:
        fcm (float): The mean compressive strength of the concrete in
            MPa.
        TABULAR_VALUES (dict): Material constants for different cement
            types and aggregate types as given by the fib Model Code
            2010.

    returns:
        float: The notional drying shrinkage in mm/mm.
    """
    return ((220 + 110*TABULAR_VALUES['alpha_ds1'])
            * np.exp(-TABULAR_VALUES['alpha_ds2'] * fcm))*1e-6


def _calc_beta_ds(
        time: np.ndarray, ts: float, notional_size: float
    ) -> np.ndarray:
    """Calculate the multiplication factor beta_ds.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-82.

    args:
        time (numpy.ndarray): The different times in days at which the
            shrinkage strain is determined.
        ts (float): Age of the concrete when exposed to the
            environment.
        notional_size (float): The notional size of the considered
            element in mm, defined as 2A/u.

    returns:
        numpy.ndarray: Multiplication factor used for calculating the
            drying shrinkage as a function of time.
    """
    return np.sqrt((time - ts)
           / (0.035*(notional_size)**2 + (time - ts)))


def _calc_beta_s1(fcm: float) -> float:
    """Calculate the correction factor beta_s1.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-83.

    args:
        fcm (float): The mean compressive strength of the concrete in
            MPa.

    returns:
        float: Multiplication factor used when calculating the drying
            shrinkage.
    """
    return min((35/fcm)**0.1, 1.0)


def _calc_beta_RH(rh: float, beta_s1: float) -> float:
    """Calculate the multiplication factor beta_RH.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-81

    args:
        rh (float): The relative humidity of the environment.
        beta_s1 (float): Multiplication factor as calculated in the
            fib Model Code 2010 (2013), Eq. 5.1-83.

    returns:
        float: Multiplication factor used when calculating the drying
            shrinkage.
    """
    if rh < 1:
        if rh >= 0.99*beta_s1:
            return 0.25
        elif 0.4*beta_s1 <= rh < 0.99*beta_s1:
            return -1.55*(1 - rh**3)
        else:
            raise ValueError("The specified rh*beta_s1 is not in the "
                "range of application.")
    else:
        if rh >= 99*beta_s1:
            return 0.25
        elif 40*beta_s1 <= rh < 99*beta_s1:
            return -1.55*(1 - (rh/100)**3)
        else:
            raise ValueError("The specified rh*beta_s1 ratio is not in"
                " the range of application.")


def calc_drying_shrinkage(
        time: np.ndarray, fcm: float, ts: float, notional_size: float,
        rh: float, cem_class: str, agg_type: str) -> np.ndarray:
    """Calculate the drying shrinkage of the concrete element.

    Defined in fib Model Code 2010 (2013), Eqs. 5.1-77 and 5.1-80 - 5.1-83

    args:
        time (numpy.ndarray): The different times in days at which the
            shrinkage strain is determined.
        fcm (float): The mean compressive strength of the concrete in
            MPa.
        ts (float): Age of the concrete when exposed to the
            environment.
        notional_size (float): The notional size of the considered
            element in mm, defined as 2A/u.
        rh (float): The relative humidity of the environment.
        cem_class (str): The cement strength class that is used.
            The choices are:
                strength classes 32.5R and 42.5N;
                strength classes 42.5R, 52.5N and 52.5R;
                strength class 32.5N.
        agg_type (str): The type of aggregate used in the concrete.
            The following values can be chosen:
            ['Basalt', 'Quartzite', 'Limestone', 'Sandstone']

    returns:
        numpy.ndarray: The shrinkage strains for the given times in
            mm/mm.
    """
    _check_RH(rh)
    _check_fcm(fcm)
    _concrete_material_properties. _check_agg_types(agg_type)
    _check_cem_strength_class(cem_class)
    TABULAR_VALUES = _get_creep_shrinkage_coeffs(cem_class)
    eps_csd0 = _calc_notional_drying_shrinkage(fcm, TABULAR_VALUES)
    beta_ds = _calc_beta_ds(time, ts, notional_size)
    beta_s1 = _calc_beta_s1(fcm)
    beta_rh = _calc_beta_RH(rh, beta_s1)
    return eps_csd0*beta_rh*beta_ds


def _calc_notional_autogenous_shrinkage(
        fcm: float, TABULAR_VALUE: dict) -> float:
    """Calculate the notional autogenous shrinkage.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-78.

    args:
        fcm (float): The mean compressive strength of the concrete in
            MPa.
        TABULAR_VALUES (dict): Material constants for different cement
            types and aggregate types as given by the fib Model Code
            2010.

    returns:
        float: The notional autogenous shrinkage in mm/mm.
    """
    return -TABULAR_VALUE['alpha_as']*1e-6 \
        * ((0.1*fcm)/(6 + 0.1*fcm))**2.5


def _calc_beta_au(time: np.ndarray) -> np.ndarray:
    """Calculate multiplication factor beta_au which is used to
        determine the autogenous shrinkage.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-79.

    args:
        time (numpy.ndarray): The different times in days at which the
            autogenous strain is determined.

    returns:
        numpy.ndarray: Multiplication factor that is used to determine
            the autogenous shrinkage.
    """
    return 1 - np.exp(-0.2 * np.sqrt(time))


def calc_autogenous_shrinkage(
        time: np.ndarray, fcm: float, cem_class: str, agg_type: str
    ) -> np.ndarray:
    """Calculate the autogenous shrinkage.

    Defined in fib Model Code 2010 (2013), Eqs. 5.1-76 and 5.1-78 - 5.1-79

    args:
        time (numpy.ndarray): The different times in days at which the
            shrinkage strain is determined.
        fcm (float): The mean compressive strength of the concrete in
            MPa.
        cem_class (str): The cement type that is used.
            The choices are:
                strength classes 32.5R and 42.5N;
                strength classes 42.5R, 52.5N and 52.5R;
                strength class 32.5N.
        agg_type (str): The type of aggregate used in the concrete.
            The following values can be chosen:
            ['Basalt', 'Quartzite', 'Limestone', 'Sandstone']

    returns:
        numpy.ndarray: The autogenous shrinkage strains for the given
            times.
    """
    _concrete_material_properties._check_agg_types(agg_type)
    _check_fcm(fcm)
    _check_cem_strength_class(cem_class)
    TABULAR_VALUES = _get_creep_shrinkage_coeffs(cem_class)
    eps_au0 = _calc_notional_autogenous_shrinkage(fcm, TABULAR_VALUES)
    beta_au = _calc_beta_au(time)
    return eps_au0*beta_au


def _calc_beta_bc_fcm(fcm: float) -> float:
    """Calculate multiplication factor that accounts for the effect of
        the compressive strength of the concrete to calculate the basic
        creep coefficient.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-65.

    args:
        fcm (float): The mean compressive strength of the concrete in
            MPa.

    returns:
        float: Multiplication factor beta_bc_fcm.
    """
    return 1.8/fcm**0.7


def _calc_beta_bc_t(
        time: np.ndarray, t0: float, t0_adj: float) -> np.ndarray:
    """Calculate multiplication factor that accounts for the effect of
        the age of the of the concrete to calculate the basic
        creep coefficient.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-66.

    args:
        time (numpy.ndarray): The different times in days at which the
            basic creep coefficient is determined.
        t0 (float): The age of the concrete in days at which the
            loading is applied.
        t0_adj (float): The temperature corrected age of the concrete
            when the loading is applied in days, as defined in fib
            Model Code 2010 (2013). Eq. 5.1-85.

    returns:
        numpy.ndarray: Multiplication factors beta_bc_t.
    """
    return np.log(((30/t0_adj + 0.035)**2) * (time - t0) + 1)


def calc_basic_creep_coefficient(
        beta_bc_fcm: float, beta_bc_t: np.ndarray) -> np.ndarray:
    """Calculate the basic creep coefficient.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-64.

    args:
        beta_bc_fcm (float): Multiplication factor that accounts for
            the influence of the concrete strength of the creep
            behaviour, as defined in fib Model Code 2010 (2013),
            Eq. 5.1-65.
        beta_bc_t (numpy.ndarray): Multiplication factor that accounts
            for the influence of the age of the concrete of the creep
            behaviour, as defined in fib Model Code 2010 (2013),
            Eq. 5.1-66.

    returns:
        numpy.ndarray: The basic creep coefficient.
    """
    return beta_bc_fcm * beta_bc_t


def _calc_beta_dc_fcm(fcm: float) -> float:
    """Calculate multiplication factor that accounts for the effect of
        the strength of the concrete to calculate the drying creep
        coefficient.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-68.

    args:
        fcm (float): The mean compressive strength of the concrete in
            MPa.

    returns:
        float: Multiplication factor beta_dc_fcm.
    """
    return 412/fcm**1.4


def _calc_beta_dc_RH(rh: float, notional_size: float) -> float:
    """Calculate multiplication factor that accounts for the effect of
        the relative humidity of the environment to calculate the
        drying creep coefficient.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-69.

    args:
        rh (float): The relative humidity of the environment.
        notional_size (float): The notional size of the considered
            element in mm, defined as 2A/u.

    returns:
        float: Multiplication factor beta_RH.
    """
    _check_RH(rh)
    if rh < 1:
        return (1 - rh) / ((0.1*notional_size/100)**(1/3))
    else:
        return (1 - rh/100) / ((0.1*notional_size/100)**(1/3))


def _calc_beta_dc_t0(t0_adj: float) -> float:
    """Calculate multiplication factor that accounts for the effect of
        the (temperature corrected) age of the concrete when loading is
        applied to calculate the drying creep coefficient.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-70.

    args:
        t0_adj (float): The temperature corrected age of the concrete
            when the loading is applied in days, as defined in fib
            Model Code 2010 (2013). Eq. 5.1-85.

    returns:
        float: Multiplication factor beta_dc_t0.
    """
    return 1 / (0.1 + t0_adj**0.2)


def _calc_alpha_fcm(fcm: float) -> float:
    """Calculate multiplication factor that accounts for the effect of
        the strength of the concrete to calculate the drying creep
        coefficient.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-71d.

    args:
        fcm (float): The mean compressive strength of the concrete in
            MPa.

    returns:
        float: Multiplication factor alpha_fcm.
    """
    return np.sqrt(35/fcm)


def _calc_beta_h(notional_size: float, alpha_fcm: float) -> float:
    """Calculate multiplication factor that accounts for the effect of
        the notional size of the to calculate the drying creep
        coefficient.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-71c

    args:
        notional_size (float): The notional size of the considered
            element in mm, defined as 2A/h.
        alpha_fcm (float): Multiplication factor that accounts for the
            effect of the strength of the concrete on the drying creep
            coefficient.

    returns:
        float: Multiplication factor beta_h.
    """
    return min(1.5*notional_size + 250*alpha_fcm, 1500*alpha_fcm)


def _calc_gamma_t0(t0_adj: float) -> float:
    """Calculate exponent that accounts for the effect of the
        (temperature corrected) age of the concrete when loaded. Used
        to calculate the drying creep coefficient.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-71b.

    args:
        t0_adj (float): The temperature corrected age of the concrete
            when the loading is applied in days, as defined in fib
            Model Code 2010 (2013). Eq. 5.1-85.

    returns:
        float: Exponent gamma_t0.
    """
    return 1 / (2.3 + 3.5/np.sqrt(t0_adj))


def _calc_beta_dc_t(
        time: np.ndarray, t0: float, beta_h: float, gamma_t0: float
    ) -> np.ndarray:
    """Calculate multiplication factor that accounts for the different
        considered values of time. Used to calculate the drying creep
        coefficient.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-71a.

    args:
        time (numpy.ndarray): The different times in days at which the
            drying creep coefficient is determined.
        t0 (float): The age of the concrete concrete
            when the loading is applied in days.
        beta_h (float): Multiplication factor that accounts for the effect of
            the notional size, as calculated by Eq. 5.1-71c
        gamma_t0 (float): Exponent that accounts for the effect of the
            (temperature corrected) age of the concrete when loaded,
            as calculated by Eq. 5.1-71b.

    returns:
        numpy.ndarray: Multiplcation factor beta_dc_t for the
            considered values of time.
    """
    return ((time - t0) / (beta_h + (time - t0))) ** gamma_t0


def calc_drying_creep_coefficient(
        beta_dc_fcm: float, beta_dc_RH: float, beta_dc_t0: float,
        beta_dc_t: np.ndarray
    ) -> np.ndarray:
    """Calculate drying creep coefficient.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-67.

    args:
        beta_dc_fcm (float): Multiplication factor that accounts for
            the effect of the strength of the concrete, as calculated
            by Eq. 5.1-68.
        beta_dc_RH (float): Multiplication factor that accounts for the
            effect of the relative humidity of the environment, as
            calculated by Eq. 5.1-69.
        beta_dc_t0 (float): multiplication factor that accounts for
            the effect of the (temperature corrected) age of the
            concrete when loading is applied, as calculated by Eq.
            5.1-70.
        beta_dc_t (numpy.ndarray): multiplication factor that accounts
            for the different considered values of time, as calculated
            by Eq. 5.1-71a.

    returns:
        numpy.ndarray: Drying creep coeffcient.
    """
    return beta_dc_fcm * beta_dc_RH * beta_dc_t0 * beta_dc_t


def _calc_k_sigma(sigma: float, fcm: float) -> float:
    """Calculate the ratio between the applied stress and the mean
        concrete compressive strength.

    Defined in fib Model Code 2010 (2013), Eq. 5.1-74.

    args:
        sigma (float): The compressive stress applied to the concrete
            at ts in MPa.
        fcm (float): The mean compressive strength of the concrete in
            MPa.

    returns:
        float: Absolute value of the ratio between the stress in the
            concrete and the mean concrete strength.
    """
    return abs(sigma/fcm)


def calc_creep_coefficient(
        time: np.ndarray, t0: float, T_cur: float, fcm: float, rh: float,
        notional_size: float, sigma: float, agg_type: str,
        cem_class: str
    ) -> np.ndarray:
    """Calculate the creep coefficient.

    Defined in fib Model Code 2010, Eqs. 5.1-63 and 5.1-64 - 5.1-71

    args:
        time (numpy.ndarray): The different times in days at which the
            drying creep coefficient is determined.
        t0 (float): The age of the concrete in days at which the
            loading is applied.
        T_cur (float): The temperature of the environment during curing
            in degrees Celcius.
        fcm (float): The mean compressive strength of the concrete in
            MPa.
        rh (float): The relative humidity of the environment.
        notional_size (float): The notional size of the considered
            element in mm, defined as 2A/h.
        sigma (float): The compressive stress applied to the concrete
            at ts in MPa.
        agg_type (str): The type of aggregate used in the concrete.
            The following values can be chosen:
            ['Basalt', 'Quartzite', 'Limestone', 'Sandstone']
        cem_class (str): The cement type that is used.
            The choices are:
                strength classes 32.5R and 42.5N;
                strength classes 42.5R, 52.5N and 52.5R;
                strength class 32.5N.

    returns:
        numpy.ndarray: The creep coefficient.
    """
    _check_age_at_loading(t0)
    _check_fcm(fcm)
    _check_RH(rh)
    _check_initial_stress(sigma, fcm)
    _concrete_material_properties._check_agg_types(agg_type)
    _check_cem_strength_class(cem_class)
    TABULAR_VALUES = _get_creep_shrinkage_coeffs(cem_class)
    t0_adj = _calc_modified_t0(t0, T_cur, TABULAR_VALUES)
    # Calculate the basic creep coefficient (phi_bc) (see Eqs. 5.1-64 - 5.1-66 of [1]):
    beta_bc_fcm = _calc_beta_bc_fcm(fcm)
    beta_bc_t = _calc_beta_bc_t(time, t0, t0_adj)
    phi_bc = calc_basic_creep_coefficient(beta_bc_fcm, beta_bc_t)
    # Calculate the drying creep coefficient (phi_dc) (see Eqs. 5.1-67 - 5.1-71 of [1]):
    beta_dc_fcm = _calc_beta_dc_fcm(fcm)
    beta_dc_RH = _calc_beta_dc_RH(rh, notional_size)
    beta_dc_t0 = _calc_beta_dc_t0(t0_adj)
    alpha_fcm = _calc_alpha_fcm(fcm)
    beta_h = _calc_beta_h(notional_size, alpha_fcm)
    gamma_t0 = _calc_gamma_t0(t0_adj)
    beta_dc_t = _calc_beta_dc_t(time, t0, beta_h, gamma_t0)
    phi_dc = calc_drying_creep_coefficient(
            beta_dc_fcm, beta_dc_RH, beta_dc_t0, beta_dc_t)
    # Calculate the creep coefficient (phi) (see Eq. 5.1-63)
    phi = phi_bc + phi_dc
    # Include the effect of high stress if needed (see Eq. 5.1-74 of [1]):
    k_sigma = _calc_k_sigma(sigma, fcm)
    if 0.4 <= k_sigma <= 0.6:
        return np.exp(1.5*(k_sigma-0.4)) * phi
    else:
        return phi


def calc_creep_compliance(
        t0: float, fcm: float, T_cur: float, phi: np.ndarray,
        agg_type: str, cem_class: str
    ) -> np.ndarray:
    """Calculate the creep compliance function.

    Defined in fib Model Code 2010, Eq. 5.1-61.

    args:
        t0 (float): The age of the concrete in days at which the
            loading is applied.
        fcm (float): The mean compressive strength of the concrete in
            MPa.
        T_cur (float): The temperature of the environment during curing
            in degrees Celcius.
        phi (numpy.ndarray): The creep coefficient as defined by fib
            Model Code 2010 (2013), Eq. 5.1-63
        agg_type (str): The type of aggregate used in the concrete.
            The following values can be chosen:
            ['Basalt', 'Quartzite', 'Limestone', 'Sandstone']
        cem_class (str): The cement type that is used.
            The choices are:
                strength classes 32.5R and 42.5N;
                strength classes 42.5R, 52.5N and 52.5R;
                strength class 32.5N.

    returns:
        numpy.ndarray: The creep compliance function.
    """
    _check_age_at_loading(t0)
    _check_fcm(fcm)
    _concrete_material_properties._check_agg_types(agg_type)
    _check_cem_strength_class(cem_class)
    TABULAR_VALUES_CEM = _get_creep_shrinkage_coeffs(cem_class)
    TABULAR_VALUES_EC = _concrete_material_properties._get_Ecmod_coeffs(cem_class, fcm, agg_type)
    t0_adj = _calc_modified_t0(t0, T_cur, TABULAR_VALUES_CEM)
    Eci_t0 = _concrete_material_properties._calc_E(t0_adj, fcm, TABULAR_VALUES_EC, "mean")
    E28 = _concrete_material_properties._calc_E28(fcm, TABULAR_VALUES_EC, "mean")
    return (1 / Eci_t0) + (phi / E28)