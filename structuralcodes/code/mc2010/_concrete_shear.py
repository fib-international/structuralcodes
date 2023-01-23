"""A collection of shear formulas for concrete"""

z = 110


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