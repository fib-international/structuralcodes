"""A collection of shear formulas for concrete"""
import typing as t
import math


def vrd(fck: float, z: float, bw: float, gamma_c: float, asw: float, sw: float, fywd: float, theta: float) -> float:
    """Compute the shear resistance of a web or slab.

    fib Model Code 2010, Eq. (7.3-11)

    Args:
        vrdc (float): Design shear resistance of concrete.
        vrds (float): Design shear resistance of shear reinforcement.

    Returns:
        float: Design shear resistance
    """

    return abs(vrdc(fck, z, bw, gamma_c)) + abs(vrds(asw, sw, z, bw, fywd, theta))


def vrdc(fck: float, z: float, bw: float, gamma_c: float = 1.5) -> float:
    """The design shear resistance of a web or a slab without
    shear reinforcement.

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

def vrds(asw: float, sw: float, z: float, fywd: float, theta: float, alpha: t.optional[float] = 0.0) -> float:
    """    fib Model Code 2010, Eq. (7.3-29)

    Args:
        asw (float):
        sw (float):
        gamma_c:
        bw:

    Returns:
        float:
    ion level."""

    if alpha == 0.0:
        return (asw/sw)*z*fywd*math.cot(theta)
    else:
        return (asw/sw)*z*fywd*(cot(theta) + cot(alpha)) * sin(alpha)



def vrdmax(fck: float, bw: float, Approx_lvl: float ,theta: float, z: float, alfa: float=0, gamma_c: float = 1.5) -> float:
    """The maximum allowed shear resistance 
    
    fib Model Code 2010, eq. (7.3-26) and (7.3-24)
    
    Args:
        fck (float): The characteristic compressive strength in MPa.
        bw (float): The width.
        theta (float): The incline of the reinforment relative to the beam axis
        
    Returns:
        float: The maximum allowed shear resisThe design shear resistance providance regarled by ss of 
        approximatirrups"""

    if Approx_lvl==1:
        if alfa ==0:
            nfc=(30/fck)**(1/3)
            if nfc > 1:
                nfc=1
            return 0.55*nfc*(fck/gamma_c)*bw*z*math.sin(theta)*math.cos(theta)