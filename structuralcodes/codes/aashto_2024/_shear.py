## AASHTO LRFD functions ##
import math


def s_xe(sx: float, ag: float) -> float:
    """Determines the crack spacing parameter that is influenced by aggregate
    size.

    AASHTO LRFD 2024 10th Edition, Eq. (5.7.4.3.2-7)

    Args:
        ag (float) = maximum agreggate size(in)

    Keyword Args:
        sx (float):  crack spacing parameter that is taken as the lesser
        between dv or the maximum distance between layers of
        longitudinal crack reinforcement

    Returns:
        The crack spacing parameter in (in)

    Raises:
        ValueError: If ag is less than 0
        ValueError: If sx is less than 0
    """
    if ag < 0:
        raise ValueError(f'ag={ag} cannot be less than 0')
    if sx < 0:
        raise ValueError(f'sx={sx} cannot be less than 0')

    return sx * (1.38 / (ag + 0.63))


def eps_s(V: float, rho_l: float, bw: float, dv: float) -> float:
    """Determines the longitudinal strain.

    AASHTO LRFD 2024 10th Edition, Eq. (5.7.3.4.2-4)

    Args:
        V (float): Assumed shear force in kips
        rho_l (float): Longitudinal reinforcement ratio
        bw (float): Width of the web in (in)
        dv (float): Effective depth of longitudinal reinforcement in (in)

    Returns:
        The longitudinal strain

    Raises:
        ValueError: If V is less than 0
        ValueError: If rho_l is less than 0
        ValueError: If bw is less than 0
        ValueError: If dv is less than 0
    """
    if V < 0:
        raise ValueError(f'V={V} cannot be less than 0')
    if rho_l < 0:
        raise ValueError(f'rho_l={rho_l} cannot be less than 0')
    if bw < 0:
        raise ValueError(f'bw={bw} cannot be less than 0')
    if dv < 0:
        raise ValueError(f'dv={dv} cannot be less than 0')

    return (3.5 * V) / (29000 * rho_l * (bw * dv))


# Calculates the beta factor
def beta(s_xe: float, strain: float) -> float:
    """Determines the shear resistance factor.

    AASHTO LRFD 2024 10th Edition, Eq. (5.7.3.4.2-2)

    Args:
        s_xe (float): The crack spacing parameter that is influenced by the
        aggregate size (in)
        strain (float): The longitudinal strain

    Returns:
        The shear resistance factor

    Raises:
        ValueError: If s_xe is not between 12 and 80 (in)
    """
    if s_xe < 12:
        raise ValueError(f's_xe={s_xe} cannot be less than 12')
    if s_xe > 80:
        raise ValueError(f's_xe={s_xe} cannot be greater than 80')

    return (4.8 / (1 + 750 * strain)) * (51 / (39 + s_xe))


# Calculates the shear stress resistance in MPa
def Vc(beta: float, fc_prime: float, bw: float, d: float) -> float:
    """Determines the shear resistance in kips.

    AASHTO LRFD 2024 10th Edition, Eq. (5.7.3.3-3)

    Args:
        beta (float): The shear resistance factor
        fc_prime (float): The compressive strength of concrete in ksi
        bw (float): The width of the web in (in)
        d (float): The effective depth in (in)

    Returns:
        The shear stress resistance

    Raises:
        ValueError: If beta is less than 0
        ValueError: If fc_prime is less than 0
    """
    if beta < 0:
        raise ValueError(f'beta={beta} cannot be less than 0')
    if fc_prime < 0:
        raise ValueError(f'fc_prime={fc_prime} cannot be less than 0')

    return 0.0316 * beta * math.sqrt(fc_prime) * bw * d  # ksi


# Iterate for convergence
def _converge(
    V: float,
    bw: float,
    dv: float,
    rho_l: float,
    s_xe: float,
    strain: float,
    beta: float,
    fc_prime: float,
    tau_MPa: float,
) -> float:
    """Iterates the initial guess of shear stress to determine a more accurate
    calculation.

    Args:
        V (float): The initial assumed value of shear force in kips
        bw (float): The width of the web in (mm)
        dv (float): The effective depth of the longitudinal reinforcement
        in (mm)
        rho_l (float): The longitudinal reinforcement ratio
        s_xe (float): The crack spacing parameter influenced by aggregate
        size (in)
        strain (float): The longitudinal strain
        beta (float): The shear resistance factor
        fc_prime (float): The compressive strength of concrete
        tau_MPa (float): The initial shear stress resistance based on the
        assumed shear force

    Raises:
        ValueError: If VkN is less than 0
        ValueError: If bw is less than 0
        ValueError: If dv is less than 0
        ValueError: If rho_l is less than 0
        ValueError: If s_xe is not between 12 and 80 (in)
        ValueError: If beta is less than 0
        ValueError: If fc_prime is lsess than 0
    """
    if V < 0:
        raise ValueError(f'VkN={V} cannot be less than 0')
    if bw < 0:
        raise ValueError(f'bw={bw} cannot be less than 0')
    if dv < 0:
        raise ValueError(f'dv={dv} cannot be less than 0')
    if rho_l < 0:
        raise ValueError(f'rho_l={rho_l} cannot be less than 0')
    if s_xe < 12:
        raise ValueError(f's_xe={s_xe} cannot be less than 12')
    if s_xe > 80:
        raise ValueError(f's_xe={s_xe} cannot be greater than 80')
    if beta < 0:
        raise ValueError(f'beta={beta} cannot be less than 0')
    if fc_prime < 0:
        raise ValueError(f'fc_prime={fc_prime} cannot be less than 0')

    error = 1
    while error > 0.001:
        tau_ref = V / ((bw / 1000) * (dv / 1000) * 1000)
        delta = tau_ref - tau_MPa

        """
        If delta is negative, the next guess of shear should be bigger
        If delta is positive, the next guess of shear should be smaller (This
        idea is based on the effects shear has on the strain and beta
        factor calculations)
        """

        if delta < 0:
            V += 0.5
            strain = eps_s(V, rho_l, bw, dv)
            beta = beta(s_xe, strain)
            tau_MPa = Vc(beta, fc_prime)
            tau_ref = V / ((bw / 1000) * (dv / 1000) * 1000)
            error = abs(tau_ref - tau_MPa) / tau_MPa

        if delta > 0:
            V -= 0.5
            strain = eps_s(V, rho_l, bw, dv)
            beta = beta(s_xe, strain)
            tau_MPa = Vc(beta, fc_prime)
            tau_ref = V / ((bw / 1000) * (dv / 1000) * 1000)
            error = abs(tau_ref - tau_MPa) / tau_MPa

    return tau_MPa


def theta(strain: float) -> float:
    """Determines the angle theta in degrees.

    AASHTO LRFD 2024 10th Edition, Eq. (5.7.3.4.2-2)

    Args:
        Strain (float): The longitudinal strain

    Returns:
        The angle theta in degrees
    """
    return 29 + 3500 * strain


def beta_reinforcement(strain: float) -> float:
    """Determines the shear resistance factor when there is minimum transverse
    reinforcment.

    AASHTO LRFD 2024 10th Edition, Eq. (5.7.3.4.2-1)

    Args:
        Strain (float): The longitudinal strain

    Returns:
        The shear resistance factor
    """
    return 4.8 / (1 + 750 * strain)


def Vs(
    Av: float, fy: float, dv: float, bw: float, cot_theta: float, s: float
) -> float:
    """Determines the shear resistance of the transverese reinforcement
    in kips.

    AASHTO LRFD 2024 10th Edition, Eq (C5.7.3.3-1)

    Args:
        Av (float): The transverse reinforcement area (in^2)
        fy (float): Steel reinforcement yielding strength ksi
        dv (float): Effective depth of steel reinforcement (in)
        bw (float): The width of the web in (in)
        cot_theta (float): the cotangent of the angle theta
        s (float): The spacing of the transverse reinforcement (in)

    Returns:
        The shear resistance of the transverse reinforcement

    Raises:
        ValueError: If Av is less than 0
        ValueError: If fy is less than 0
        ValueError: If dv is less than 0
        ValueError: If bw is less than 0
        ValueError: If s is less than 0
    """
    if Av < 0:
        raise ValueError(f'Av={Av} cannot be less than 0')
    if fy < 0:
        raise ValueError(f'fy={fy} cannot be less than 0')
    if dv < 0:
        raise ValueError(f'dv={dv} cannot be less than 0')
    if bw < 0:
        raise ValueError(f'bw={bw} cannot be less than 0')
    if s < 0:
        raise ValueError(f's={s} cannot be less than 0')

    return (Av * fy * dv * cot_theta) / s


def Vn(Vc: float, Vs: float, Vp: float) -> float:
    """Determines the nominal shear resistance in kips.

    AASHTO LRFD 2024 10th Edition, Eq (5.7.3.3-1)

    Args:
        Vc (float): Compressive shear resistance in kips
        Vs (float): Shear resistance of tranverse reinforcement in
        kips
        Vp (float): Prestressing shear stress resistance in kips

    Returns:
        The nominal shear stress resistance
    """
    return Vc + Vs + Vp
