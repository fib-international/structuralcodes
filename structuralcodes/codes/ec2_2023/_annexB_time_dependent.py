"""Functions from Annex B of EN 1992-1-1:2022."""


def alpha_c(fcm_28: float) -> float:
    """Returns the coefficient to obtain the tangent modulus of elasticity from
    the secant modulus of elasticity.

    EN 1992-1-1, Table 5.1

    Args:
        fcm_28 (float): The characteristic compressive strength at 28 days
            in MPa.

    Returns:
        float: Coefficient alphac so that Ec=alphac*Ecm,28.
    """
    return max(1 / (0.8 + 0.2 * fcm_28 / 88), 1)
