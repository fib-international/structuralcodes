"""A collection of shear formulas for concrete"""




def vrd(vrdc: float, vrds: float) -> float:
    """Compute the shear resistance of a web or slab.

    fib Model Code 2010, Eq. (7.3-11)

    Args:
        vrdc (float): Design shear resistance of concrete.
        vrds (float): Design shear resistance of shear reinforcement.

    Returns:
        float: Design shear resistance
    """

    return abs(vrdc) + abs(vrds)


def vrdc(fck: float, z: float, bw: float, gamma_c: float = 1.5,) -> float:
    """The design shear resistance of a web or a slab without shear reinforcement.

    fib Model Code 2010, Eq. (7.3-17)

    Args:
        vck (float): The characteristic compressive strength in MPa.
        z (float): the effective shear depth.
        gamma_c: Material factor.
        bw: 

    Returns:
        float: Design shear resistance without shear reinforcement
    """
    kv = 180/(1000+1.25*z)
    if fck**0.5 > 8:
        fsqr = 8
    else:
        fsqr = fck**0.5

    return (kv*fsqr*z*bw)/gamma_c
""he""
    
