"""Concrete material properties according to Tab. 3.1."""

import math

from structuralcodes.codes import mc2010


def fcm(fck: float, delta_f: float = 8) -> float:
    """The mean compressive strength of concrete.

    EN 1992-1-1:2004, Table 3.1.

    Args:
        fck (float): The characteristic compressive strength of concrete in
            MPa.

    Keyword Args:
        delta_f (float): The difference between the mean and the
            characteristic strength.

    Returns:
        float: The mean compressive strength in MPa.
    """
    return mc2010.fcm(fck=abs(fck), delta_f=abs(delta_f))


def fctm(fck: float) -> float:
    """The mean tensile strength of concrete.

    EN 1992-1-1: 2004, Table 3.1.

    Args:
        fck (float): The characteristic compressive strength of concrete in
            MPa.

    Returns:
        float: The mean tensile strength in MPa.
    """
    return mc2010.fctm(fck=abs(fck))


def fctk_5(fctm: float) -> float:
    """The 5% fractile of the tensile strength of concrete.

    EN 1992-1-1: 2004, Table 3.1.

    Args:
        fctm (float): The mean tensile strength of concrete in MPa.

    Returns:
        float: The 5% fractile of the tensile strength in MPa.
    """
    return mc2010.fctkmin(fctm=abs(fctm))


def fctk_95(fctm: float) -> float:
    """The 95% fractile of the tensile strength of concrete.

    EN 1992-1-1: 2004, Table 3.1.

    Args:
        fctm (float): The mean tensile strength of concrete in MPa.

    Returns:
        float: The 95% fractile of the tensile strength in MPa.
    """
    return mc2010.fctkmax(fctm=abs(fctm))


def Ecm(fcm: float) -> float:
    """The secant modulus of concrete.

    EN 1992-1-1:2004, Table 3.1.

    Args:
        fcm (float): The mean compressive strength of concrete in MPa.

    Returns:
        float: The secant modulus of concrete in MPa.
    """
    return 22000.0 * math.pow(abs(fcm) / 10, 0.3)


def eps_c1(fcm: float) -> float:
    """The strain at maximum compressive stress of concrete (fcm) for the
    Sargin constitutive law.

    EN 1992-1-1:2004, Table 3.1.

    Args:
        fcm (float): The mean compressive strength of concrete in MPa.

    Returns:
        float: The strain at maximum compressive stress, absolute value, no
        unit.
    """
    return min(0.7 * math.pow(abs(fcm), 0.31), 2.8) / 1000


def eps_cu1(fck: float) -> float:
    """The ultimate strain for the Sargin constitutive law.

    EN 1992-1-1:2004, Table 3.1.

    Args:
        fck (float): The characteristic compressive strength of concrete in
            MPa.

    Returns:
        float: The ultimate strain, absolute value, no unit.
    """
    fck = abs(fck)
    return (
        3.5 / 1000
        if fck < 50
        else (2.8 + 27 * ((98 - fcm(fck)) / 100) ** 4) / 1000
    )


def k_sargin(
    Ecm: float,
    fcm: float,
    eps_c1: float,
) -> float:
    """Computation of k parameter for Sargin constitutive Law.

    EN 1992-1-1:2004, Eq. (3.14)

    Args:
        Ecm (float): the mean elastic modulus of concrete in MPa.
        fcm (float): the mean compressive strength in MPa.
        eps_c1 (float): the strain corresponding to peak stress.
    """
    return 1.05 * Ecm * abs(eps_c1) / fcm


def eps_c2(fck: float) -> float:
    """The strain at maximum compressive stress of concrete for the
    parabolic-rectangular law.

    EN 1992-1-1:2004, Table 3.1.

    Args:
        fck (float): The characteristic compressive strength of concrete in
            MPa.

    Returns:
        float: The strain at maximum compressive stress, absolute value, no
        unit.
    """
    fck = abs(fck)
    return (
        2.0 / 1000 if fck <= 50 else (2.0 + 0.085 * (fck - 50) ** 0.53) / 1000
    )


def eps_cu2(fck: float) -> float:
    """The ultimate strain of the parabolic-rectangular law.

    EN 1992-1-1:2004, Table 3.1.

    Args:
        fck (float): The characteristic compressive strength of concrete in
            MPa.

    Returns:
        float: The ultimate strain, absolute value, no unit.
    """
    fck = abs(fck)
    return (
        3.5 / 1000
        if fck <= 50
        else (2.6 + 35 * ((90 - fck) / 100) ** 4) / 1000
    )


def n_parabolic_rectangular(fck: float) -> float:
    """The exponent in the parabolic-rectangular law.

    EN 1992-1-1:2004, Table 3.1.

    Args:
        fck (float): The characteristic compressive strength of concrete in
            MPa.

    Returns:
        float: The exponent n, absolute value, no unit.
    """
    fck = abs(fck)
    return 2.0 if fck <= 50 else (1.4 + 23.4 * ((90 - fck) / 100) ** 4)


def eps_c3(fck: float) -> float:
    """The strain at maximum compressive stress of the bi-linear law.

    EN 1992-1-1:2004, Table 3.1.

    Args:
        fck (float): The characteristic compressive strength of concrete in
            MPa.

    Returns:
        float: The strain at maximum compressive stress, absolute value, no
        unit.
    """
    fck = abs(fck)
    return 1.75 / 1000 if fck <= 50 else (1.75 + 0.55 * (fck - 50) / 40) / 1000


def eps_cu3(fck: float) -> float:
    """The ultimate strain of the bi-linear law.

    EN 1992-1-1:2004, Table 3.1.

    Args:
        fck (float): The characteristic compressive strength of concrete in
            MPa.

    Returns:
        float: The ultimate strain, absolute value, no unit.
    """
    return eps_cu2(fck)


def fcd(fck: float, alpha_cc: float, gamma_c: float) -> float:
    """The design compressive strength of concrete.

    EN 1992-1-1:2004, Eq. (3.15)

    Args:
        fck (float): The characteristic compressive strength in MPa.
        alpha_cc (float): A factor for considering long-term effects on the
            strength, and effects that arise from the way the load is applied.
        gamma_c (float): The partial factor of concrete.

    Returns:
        float: The design compressive strength of concrete in MPa
    """
    return abs(alpha_cc) * abs(fck) / abs(gamma_c)
