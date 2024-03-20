"""Concrete material properties according to Tab. 3.1."""

import math

from structuralcodes.codes import mc2010


def fcm(fck: float, delta_f: float = 8) -> float:
    """The mean compressive strength of concrete.

    EN 1992-1-1:2004, Table 3.1.

    Args:
        fck (float): The caracteristic compressive strength of concrete.
    """
    return mc2010.fcm(fck=fck, delta_f=delta_f)


def fctm(fck: float) -> float:
    """The mean tensile strength of concrete.

    EN 1992-1-1: 2004, Table 3.1.

    Args:
        fck (float): The caracteristic compressive strength of concrete.
    """
    return mc2010.fctm(fck=fck)


def fctk_5(_fctm: float) -> float:
    """The 5% fractile of the tensile strength of concrete.

    EN 1992-1-1: 2004, Table 3.1.

    Args:
        _fctm (float): The mean tensile strength of concrete.
    """
    return mc2010.fctkmin(_fctm=_fctm)


def fctk_95(_fctm: float) -> float:
    """The 95% fractile of the tensile strength of concrete.

    EN 1992-1-1: 2004, Table 3.1.

    Args:
        _fctm (float): The mean tensile strength of concrete.
    """
    return mc2010.fctkmax(_fctm=_fctm)


def Ecm(_fcm: float) -> float:
    """The secant modulus of concrete.

    EN 1992-1-1:2004, Table 3.1.

    Args:
        _fcm (float): The mean compressive strength of concrete.
    """
    return 22000.0 * math.pow(_fcm / 10, 0.3)


def eps_c1(_fcm: float) -> float:
    """The strain at maximum compressive stress of concrete (fcm) for the
    Sargin constitutive law.

    EN 1992-1-1:2004, Table 3.1.

    Args:
        _fcm (float): The mean compressive strength of concrete.
    """
    return min(0.7 * math.pow(_fcm, 0.31), 2.8) / 1000


def eps_cu1(fck: float) -> float:
    """The ultimate strain for the Sargin constitutive law.

    EN 1992-1-1:2004, Table 3.1.

    Args:
        fck (float): The characteristic compressive strength of concrete.
    """
    return (
        3.5 / 1000
        if fck < 50
        else (2.8 + 27 * ((98 - fcm(fck)) / 100) ** 4) / 1000
    )


def eps_c2(fck: float) -> float:
    """The strain at maximum compressive stress of concrete for the
    parabolic-rectangular law.

    EN 1992-1-1:2004, Table 3.1.

    Args:
        fck (float): The characteristic compressive strength of concrete.

    """
    return (
        2.0 / 1000 if fck < 50 else (2.0 + 0.085 * (fck - 50) ** 0.53) / 1000
    )


def eps_cu2(fck: float) -> float:
    """The ultimate strain of the parabolic-rectangular law.

    EN 1992-1-1:2004, Table 3.1.

    Args:
        fck (float): The characteristic compressive strength of concrete.
    """
    return (
        3.5 / 1000 if fck < 50 else (2.6 + 35 * ((90 - fck) / 100) ** 4) / 1000
    )


def n(fck: float) -> float:
    """The exponent in the parabolic-rectangular law.

    EN 1992-1-1:2004, Table 3.1.

    Args:
        fck (float): The characteristic compressive strength of concrete.
    """
    return 2.0 if fck < 50 else (1.4 + 23.4 * ((90 - fck) / 100) ** 4)


def eps_c3(fck: float) -> float:
    """The ultimate strain of the bi-linear law.

    EN 1992-1-1:2004, Table 3.1.

    Args:
        fck (float): The characteristic compressive strength of concrete.
    """
    return 1.75 / 1000 if fck < 50 else (1.75 + 0.55 * (fck - 50) / 40) / 1000


def eps_cu3(fck: float) -> float:
    """The ultimate strain of the bi-linear law.

    EN 1992-1-1:2004, Table 3.1.

    Args:
        fck (float): The characteristic compressive strength of concrete.
    """
    return eps_cu2(fck)
