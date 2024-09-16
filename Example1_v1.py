import math

from shapely import Polygon

from structuralcodes import codes, materials
from structuralcodes.geometry import SurfaceGeometry
from structuralcodes.sections._generic import GenericSection
from structuralcodes.sections._reinforcement import add_reinforcement_line

# from structuralcodes.plots.section_plots import draw_section_response,draw_section,get_stress_point

### CREATE THE SECTION ########
codes.set_design_code(design_code='ec2_2004')
concrete = materials.concrete.create_concrete(fck=25)
reinforcemnet = materials.reinforcement.create_reinforcement(
    fyk=500,
    Es=200000,
    density=7850,
    ftk=500,
    epsuk=0.07,
)
poly = Polygon(((0, 0), (350, 0), (350, 500), (0, 500)))
geo = SurfaceGeometry(poly, concrete)
geo = add_reinforcement_line(
    geo, (50, 450), (300, 450), 12, reinforcemnet, n=3
)
geo = add_reinforcement_line(geo, (50, 50), (300, 50), 20, reinforcemnet, n=6)
sec = GenericSection(
    geo,
)
# sec.geometry = sec.geometry.translate(-sec.gross_properties.cy, -sec.gross_properties.cz)
sec.geometry = sec.geometry.translate(-175, -250)
### CREATE THE SECTION ########


def calculate_strain_profile_uniaxial(section: GenericSection, n_ed=0, m_ed=0):
    """'Get the strain plane for horizontal neutral axe.

    Args:
        n_ed [N]: Axial load
        m_ed [N*m]: Bending moment My

    Returns:
        eps_a: strain at (0,0)
        chi: curvature of the section
    """
    # Rotate the section of angle theta
    sec_geom = section.geometry

    def interpolate(x1, x2, y1, y2, x):
        if x2 - x1 == 0:
            return y1
        y = y1 + (y2 - y1) * (x - x1) / (x2 - x1)
        return y

    if m_ed < 0:
        r = section.section_calculator.calculate_bending_strength(
            theta=0, n=n_ed
        )
        m1, m2 = r.m_y, 0
        chi1, chi2 = r.chi_y, 0
        e1, e2 = r.eps_a, 0
    else:
        r = section.section_calculator.calculate_bending_strength(
            theta=math.pi, n=n_ed
        )
        m1, m2 = 0, r.m_y
        chi1, chi2 = 0, -r.chi_y
        e1, e2 = 0, r.eps_a

    chi = interpolate(m1, m2, chi1, chi2, m_ed)

    # Previous position of strain at (0,0)
    strain = [interpolate(m1, m2, e1, e2, m_ed), 0, 0]

    ITMAX = 100
    tolerance = 1000  # (1e-3 mkN)
    iter = 0
    M_ = 1e11

    # Stop if the ultimate curvature is exceeded
    if chi > chi2:
        return e2, chi2
    elif chi < chi1:
        return e1, chi1

    while abs(m_ed - M_) > tolerance and iter < ITMAX:
        strain = section.section_calculator.find_equilibrium_fixed_curvature(
            sec_geom, n_ed, chi, strain[0]
        )
        (
            _,
            My,
            Mz,
            _,
        ) = section.section_calculator.integrator.integrate_strain_response_on_geometry(
            geo=sec_geom,
            strain=strain,
            tri=section.section_calculator.triangulated_data,
        )

        M_ = (My**2 + Mz**2) ** 0.5

        if m_ed < M_:
            m2 = M_
            chi2 = chi
        else:
            m1 = M_
            chi1 = chi

        chi = interpolate(m1, m2, chi1, chi2, m_ed)
        eps_a = strain[0]
        iter += 1
        print(
            'chi1',
            round(chi1 * 1e6),
            'chi',
            round(chi * 1e6),
            'chi2',
            round(chi2 * 1e6),
        )
    if iter == ITMAX:
        return None, None
    else:
        return eps_a, chi


def calculate_strain_profile(
    section: GenericSection, n_ed=0, my_ed=0, mz_ed=0
):
    """'Get the strain plane for a given axial force and biaxial bending .

    Args:
        n_ed [N]: Axial load
        my_ed [N*m]: Bending moment My
        mz_ed [N*m]: Bending moment My
    Returns:
        eps_a: strain at (0,0)
        chiy,chiz: curvatures of the section
    """

    def interpolate(x1, x2, y1, y2, x):
        if x2 - x1 == 0:
            return y1
        y = y1 + (y2 - y1) * (x - x1) / (x2 - x1)
        return y

    def angle(x, y):
        """Obtain the angle of the vector with respect to (1,0) in counterclockwise direction."""
        angle = math.atan2(y, x)

        # check angle between [0, 2π]
        if angle < 0:
            angle += 2 * math.pi

        # Degrees
        angle_deg = math.degrees(angle)

        return angle  # , angle_deg

    # Check if the section can carry the axial load
    section.section_calculator.check_axial_load(n=n_ed)

    # obtain the boundaries of theta (angle of moments) corresponding to the quadrants of the alpha angle
    res = section.section_calculator.calculate_bending_strength(math.pi, n_ed)
    theta_0 = angle(res.m_y, res.m_z)

    res = section.section_calculator.calculate_bending_strength(
        3 * math.pi / 2, n_ed
    )
    theta_1 = angle(res.m_y, res.m_z)

    res = section.section_calculator.calculate_bending_strength(0, n_ed)
    theta_2 = angle(res.m_y, res.m_z)

    res = section.section_calculator.calculate_bending_strength(
        math.pi / 2, n_ed
    )
    theta_3 = angle(res.m_y, res.m_z)

    # Angle of moments: theta
    theta = angle(my_ed, mz_ed)

    # get first attemp of alfa for iterations and quadrant of aplication
    if (theta_0 <= theta < theta_1) or (0 <= theta < theta_1):
        alfa_1 = 0
        alfa_2 = math.pi / 2
        # alfa = interpolate(theta_0, theta_1, alfa_1, alfa_2, theta)
        if theta_0 <= theta_1:  # theta_0>0
            alfa = interpolate(theta_0, theta_1, alfa_1, alfa_2, theta)
        else:  # theta_0<0
            alfa = interpolate(
                theta_0 - 2 * math.pi, theta_1, alfa_1, alfa_2, theta
            )
    elif theta_1 <= theta < theta_2:
        alfa_1 = math.pi / 2
        alfa_2 = math.pi
        alfa = interpolate(theta_1, theta_2, alfa_1, alfa_2, theta)
    elif theta_2 <= theta < theta_3:
        alfa_1 = math.pi
        alfa_2 = 3 * math.pi / 2
        alfa = interpolate(theta_0, theta_1, alfa_1, alfa_2, theta)
    else:
        alfa_1 = 3 * math.pi / 2
        alfa_2 = 2 * math.pi
        alfa = interpolate(theta_0, theta_1, alfa_1, alfa_2, theta)

    ITMAX = 50
    iter = 0
    My = 1e11
    Mz = 1e11

    while ((alfa_2 - alfa_1) > 1e-4) and iter < ITMAX:
        # rotated section
        rotated_geom = section.geometry.rotate(alfa)
        _sec = GenericSection(rotated_geom)
        _sec.section_calculator.integrator = (
            section.section_calculator.integrator
        )

        # Momento resultante que será el MY que aplique a la sección rotada
        m_ed = (my_ed**2 + mz_ed**2) ** 0.5  # * sign
        # Se saca el plano deformación de la sección rotada para momento en Y rotado
        eps_a, chi = calculate_strain_profile_uniaxial(_sec, n_ed, m_ed)

        # Se descomponen las curvarutas, que pertenecerán a la sección original
        chiy = abs(chi) * math.cos(alfa)
        chiz = abs(chi) * math.sin(alfa)
        # draw_section(section)
        N, My, Mz = section.section_calculator.integrate_strain_profile(
            (eps_a, chiy, chiz)
        )
        print('------------------')
        print(
            'alfa_1',
            round(alfa_1 * 180 / math.pi),
            'alfa',
            round(alfa * 180 / math.pi),
            'alfa_2',
            round(alfa_2 * 180 / math.pi),
        )
        print(
            'chiy',
            round(chiy * 1e6, 2),
            'chiz',
            round(chiz * 1e6, 2),
            'chi',
            round(chi * 1e6, 2),
        )
        print(
            'N',
            round(N / 1e3),
            'My',
            round(My / 1e6),
            'Mz',
            round(Mz / 1e6),
        )

        if angle(My, Mz) > angle(my_ed, mz_ed):
            alfa_2 = alfa
        else:
            alfa_1 = alfa

        alfa = 0.5 * (alfa_1 + alfa_2)

    if iter == ITMAX:
        return None, None, None
    else:
        return eps_a, chiy, chiz


### TEST calculate_strain_profile_biaxial  ########
n = 0 * 1e3
my = -310 * 1e6  # -150
mz = 0 * 1e6  # -100


res = sec.section_calculator.calculate_strain_profile(n, my, mz)
print(
    'eps',
    round(res[0] * 1e3, 2),
    '  chiy',
    round(res[1] * 1e6, 2),
    '  chiz',
    round(res[2] * 1e6, 2),
)
forces = sec.section_calculator.integrate_strain_profile(res)
print(
    'N',
    round(forces[0] / 1e3, 2),
    '  My',
    round(forces[1] / 1e6, 2),
    '  Mz',
    round(forces[2] / 1e6, 2),
)

"""

res = calculate_strain_profile(sec, n, my, mz)
print(
    'eps',
    round(res[0] * 1e3, 2),
    '  chiy',
    round(res[1] * 1e6, 2),
    '  chiz',
    round(res[2] * 1e6, 2),
)
"""
