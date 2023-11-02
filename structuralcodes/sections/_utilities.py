"""General Utilities"""
import math

def marin_integration(part_type: str, y: list[float], z: list[float], n:int, m: int) -> float:
    """implements Marin's algorithm for integrating a polynomial 
        over a closed polygon 
    part_type (str): to control the integration for parts which are not
        polygons such as reinforcement bars. Also this parameters is used
        to subtract polygons which represent voids
    y,z (lists of floats): coordinates of polygon vertices or coordinates
        of centers of bars in case of reinforecement
    n (int): the degree of the polynomial in the y direction
    m (int): the degree of the polynomial in the z direction
    The order of the polygon vertices is significant. If the points
    are in counterclockwise order the integral is positive, otherwise
    the integral is negative
    """
    if part_type == 'P' or part_type == 'p':
        print('Steel bars')
    else:
        mom = 0
        for i in range(len(y)-1):
            if abs(math.sqrt((y[i+1]-y[i])**2+(z[i+1]-z[i])**2))>1e-5:
                ssj = 0 # sum in j
                aux = y[i]*z[i+1]-y[i+1]*z[i]
                if abs(aux) > 1e-5:
                    for j in range(m+1):
                        ssk = 0 # sum in k
                        for k in range(n+1):
                            ssk = ssk+math.comb(j+k, j)*math.comb(m+n-j-k, n-k)*z[i]**(n-k)*z[i+1]**k
                        x3 = y[i]**(m-j)*y[i+1]**j
                        ssj = ssj+ssk*x3
                    mom = mom+ssj*aux
    return mom/(math.comb(m+n,n)*(m+n+1)*(m+n+2))

