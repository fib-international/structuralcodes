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
    """'Get the strain plane for uniaxial bending

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
    strain = [0, 0, 0]
    strain[0] = interpolate(m1, m2, e1, e2, m_ed)

    ITMAX = 200
    tolerance = 1000  # (1e-3 mkN)
    iter = 0
    M_ = 1e11

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

    if iter == ITMAX:
        return None, None
    else:
        return eps_a, chi


def calculate_strain_profile_biaxial(
    section: GenericSection, n_ed=0, my_ed=0, mz_ed=0
):
    """'Test values of the neutral fibre angle to get the bending strain plane My Mz using the uniaxial strain plane function.

    Args:
        n_ed [N]: Axial load
        m_ed [N*m]: Bending moment My

    Returns:
        eps_a: strain at (0,0)
        chi: curvature of the section
    """
    # Check if the section can carry the axial load
    section.section_calculator.check_axial_load(n=n_ed)

    if my_ed >= 0 and mz_ed >= 0:
        alfa_1 = 0
        alfa_2 = math.pi / 2
        alfa_0 = 0
    elif my_ed < 0 and mz_ed >= 0:
        alfa_1 = math.pi / 2
        alfa_2 = math.pi
        alfa_0 = math.pi
    elif my_ed < 0 and mz_ed < 0:
        alfa_1 = math.pi
        alfa_2 = 3 * math.pi / 2
        alfa_0 = math.pi
    else:
        alfa_1 = 3 * math.pi / 2
        alfa_2 = 2 * math.pi
        alfa_0 = 2 * math.pi

    # Tantea primer valor de alfa asumiendo que la resultante de momentos sigue el mismo vector que la resultante de curvaturas
    if my_ed == 0:
        alfa = math.pi / 2 if mz_ed >= 0 else 3 * math.pi / 2
        alfa_1 = 1 / 4 * math.pi if mz_ed >= 0 else 5 / 4 * math.pi
        alfa_2 = 3 / 4 * math.pi if mz_ed >= 0 else 7 / 4 * math.pi
        # sign = 1 if mz_ed > 0 else -1
    else:
        alfa = math.atan(mz_ed / my_ed) + alfa_0
        # sign = np.sign(my_ed)

    ITMAX = 50
    tolerance = 1000  # (1e-3 mkN)
    iter = 0
    My = 1e11
    Mz = 1e11

    while ((alfa_2 - alfa_1) > 1e-4) and iter < ITMAX:
        # Giro de la sección un ángulo alfa
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

        if my_ed == 0:
            if Mz / My > 1e11:
                alfa_2 = alfa
            else:
                alfa_1 = alfa
        elif Mz / My > mz_ed / my_ed:
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
my = 20 * 1e6  # -150
mz = 1 * 1e6  # -100


res = calculate_strain_profile_biaxial(sec, n, my, mz)
print(
    'eps',
    round(res[0] * 1e3, 2),
    '  chiy',
    round(res[1] * 1e6, 2),
    '  chiz',
    round(res[2] * 1e6, 2),
)
"""res = calculate_strain_profile_uniaxial(sec, n, my)
print()
print(res[0] * 1e3)
print(res[1] * 1e6)"""
