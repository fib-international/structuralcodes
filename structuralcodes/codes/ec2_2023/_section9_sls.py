"""Functions from Section 9 of FprEN 1992-1-1:2023"""


def Ec_eff(fcm: float, phi: float, kE: float = 9500) -> float:
    """Returns de Effective modulus of elasticity from fcm and phi

    FprEN 1992-1-1:2022, Eq. (9.1)

    Args:
        fcm (float): The mean compressive strength in MPa.
        phi (float): The creep coefficient

    Keyword Args:
        kE (float): Constant to account for the type of aggregate.

    Returns:
        float: The effective modulus of elastiticy in MPa."""
    Ecm = kE * fcm ** (1 / 3)
    return 1.05 * Ecm / (1 + phi)


def As_min_y(
    NEd: float, b: float, h: float, fct_eff: float, fyk: float
) -> float:
    """Returns the minimum reinforcement to avoid yielding of steel. Box or T
       sections are to be divided into rectangles

    FprEN 1992-1-1:2022, Eq. (9.4)
    Eq. (9.2) and (9.3) are particular cases of the general equation
    Eq. (9.2) is valid for pure bending, hence NEd=0
    Eq. (9.3) is valid for pure tension. The general expression has an upper
    limit that equals the values of Eq. (9.3)

    Args:
        NEd (float): SLS axial force applied on the section or rectangle
                     (compressions are negative) in kN
        b (float): the width of the section or rectangle in meters
        h (float): the height of the section or rectange in meters
        fct_eff (float): effective tension strength of concrete (can normally
                         be taken as the mean tensile strength) in MPa
        fyk (float): characteristic yield strength of steel in MPa

    Returns:
        As_min_y[0] (float): The minimum tensile reinforcement to avoid
                             yielding of steel on the most tensioned fibre of
                             the rectangle (As_min_y1) in cm2
        As_min_y[1] (float): The minimum tensile reinforcement to avoid
                             yielding of steel on the most tensioned fibre of
                             the rectangle (As_min_y2) in cm2
    """
    As_min_y1 = (
        min(
            max(
                (0.3 * NEd / 1000 + 0.2 * kh(b, h) * fct_eff * b * h) / fyk, 0
            ),
            0.5 * kh(b, h) * fct_eff * b * h / fyk,
        )
        * 1e4
    )
    return As_min_y1, min(max(NEd / fyk * 10 - As_min_y1, 0), As_min_y1)


def kh(b: float, h: float) -> float:
    """Returns factor kh, which reduces the tensile strength of concrete to
    account for imposed restrained deformations due
     to shrinkage

    FprEN 1992-1-1:2022, Eq. (9.5)

    Args:
        b (float): width of the rectangle in meters
        h (float): height of the rectangle in meters

    Returns:
        Factor kh which applies to the tensile resistance of concrete"""
    return min(max(0.8 - 0.6 * (min(b, h) - 0.3), 0.5), 0.8)
