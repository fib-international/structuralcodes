import math

import numpy as np
from shapely.geometry import LineString, Polygon

from structuralcodes.codes import _CODE, get_design_codes, set_design_code
from structuralcodes.geometry import CompoundGeometry, SurfaceGeometry
from structuralcodes.materials.concrete import create_concrete
from structuralcodes.materials.constitutive_laws import (
    ParabolaRectangle,
    Sargin,
    UserDefined,
)
from structuralcodes.sections._generic import GenericSection


def calculate_strain_profile_batch(sec: GenericSection, n_ed, my_ed):
    """'Get the moment curvature diagram. Then interpolate to get curvature for my_d.

    Args:
        n_ed [N]: Axial load
        my_ed [N*m]: List of bending moment My

    Returns:
        eps_a: strain at (0,0)
        chi: curvature of the section
    """
    my_ed = np.array(my_ed, dtype=float)
    result_neg = sec.section_calculator.calculate_moment_curvature(0, n=n_ed)
    result_pos = sec.section_calculator.calculate_moment_curvature(
        math.pi, n=n_ed
    )
    chi_y = np.vstack((result_neg.chi_y, result_pos.chi_y)).reshape(-1)
    #  results.m_y in (Nmm)
    m_y = (np.vstack((result_neg.m_y, result_pos.m_y))).reshape(-1)
    eps_axial = np.vstack(
        (result_neg.eps_axial, result_pos.eps_axial)
    ).reshape(-1)
    # sorts results M-chi
    index = np.argsort(chi_y)
    chi_y = chi_y[index]
    m_y = m_y[index]
    eps_axial = eps_axial[index]

    chi_y_ = np.interp(my_ed, m_y, chi_y)  #  results.m_y in (Nmm)
    epsa_ = np.interp(my_ed, m_y, eps_axial)
    return epsa_, chi_y_


def calculate_strain_profile(section: GenericSection, n_ed=0, m_ed=0):
    """'Get the stress in a given point (y,z) inside the cross section.

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

    # Check if the section can carry the axial load
    section.section_calculator.check_axial_load(n=n_ed)

    if m_ed < 0:
        r = section.section_calculator.calculate_bending_strength(
            theta=0, n=n_ed
        )
        m1, m2 = r.m_y, 0
        chi1, chi2 = r.chi_y, 0
    else:
        r = section.section_calculator.calculate_bending_strength(
            theta=math.pi, n=n_ed
        )
        m1, m2 = 0, r.m_y
        chi1, chi2 = 0, -r.chi_y

    chi = interpolate(m1, m2, chi1, chi2, m_ed)

    # Previous position of strain at (0,0)
    strain = [0, 0, 0]

    ITMAX = 200
    tolerance = 1000  # (1e-3 mkN)
    iter = 0
    My = 1e11

    while abs(m_ed - My) > tolerance and iter < ITMAX:
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

        if My > m_ed:
            m2 = My
            chi2 = chi
        else:
            m1 = My
            chi1 = chi

        chi = interpolate(m1, m2, chi1, chi2, m_ed)
        eps_a = strain[0]
        iter += 1

    if iter == ITMAX:
        return None, None
    else:
        return eps_a, chi


def calculate_cracking_moment(section: GenericSection, n=0, plot=False):
    """Calculate the cracking moment of a R.C section in N*mm.

    Args:
        n [N]: Axial external load

    Returns:
        m_cracking_pos: positive My_cracking
        m_cracking_neg: negative My_cracking
    """

    def modify_consitutive_law(geom: SurfaceGeometry):
        mat = geom.material
        gamma_c = 1.5  # !!!! Temporary. gamma_c should be readable from ConstitutiveLaw class (TODO)
        if (
            _CODE is None
        ):  # should be the code of the material used to create the section
            set_design_code(get_design_codes()[1])
        if isinstance(mat, ParabolaRectangle):
            strains = np.linspace(mat._eps_u, 0, 20)
            aux_concrete = create_concrete(abs(mat._fc * gamma_c), gamma_c=1)
        elif isinstance(mat, Sargin):
            strains = np.linspace(mat._eps_cu1, 0, 20)
            aux_concrete = create_concrete(abs(mat._fc * gamma_c), gamma_c=1)
        else:  # UserDefined
            strains = np.linspace(-0.0035, 0, 20)
            eps_max, eps_min = mat.get_ultimate_strain()
            s1 = np.linspace(eps_min, eps_max, 100)
            s2 = mat.get_stress(s1)
            aux_concrete = create_concrete(abs(s2.min() * gamma_c), gamma_c=1)

        stress = aux_concrete.constitutive_law.get_stress((strains))
        # stress = mat.get_stress((strains))
        strains = strains.tolist()
        stress = stress.tolist()

        fctm = aux_concrete.fctm
        Ec = aux_concrete.constitutive_law.get_tangent(0)[0]
        tensile_strain_limit = fctm / Ec
        strains.append(tensile_strain_limit)
        stress.append(fctm)
        ec_const = UserDefined(
            strains, stress, 'constitutive_law_with_tension'
        )
        geom.material = ec_const
        if plot:
            from structuralcodes.plots.section_plots import (
                draw_constitutive_law,
            )

            draw_constitutive_law(ec_const)

    import copy

    sec = copy.deepcopy(section)
    m_cracking_pos, m_cracking_neg = 0, 0
    if sec.geometry.reinforced_concrete:
        if isinstance(sec.geometry, SurfaceGeometry):
            modify_consitutive_law(sec.geometry)
        elif isinstance(sec.geometry, CompoundGeometry):
            for i in range(len(sec.geometry.geometries)):
                modify_consitutive_law(sec.geometry.geometries[i])
        else:
            raise ValueError(
                'geometry must be SurfaceGeometry or CompoundGeometry'
            )

        m_cracking_pos = sec.section_calculator.calculate_bending_strength(
            theta=math.pi, n=n
        ).m_y
        m_cracking_neg = sec.section_calculator.calculate_bending_strength(
            theta=0, n=n
        ).m_y

        # PLOT STRESS DIAGRAM
        if plot:
            from structuralcodes.plots.section_plots import (
                draw_section_response,
            )

            eps, chi = calculate_strain_profile(sec, 0, m_cracking_pos)
            draw_section_response(
                sec,
                eps,
                chi,
                title=f'M cracking pos {round(m_cracking_pos/1e6,2)} mkN',
            )

            eps, chi = calculate_strain_profile(sec, 0, m_cracking_neg)
            draw_section_response(
                sec,
                eps,
                chi,
                title=f'M cracking neg {round(m_cracking_neg/1e6,2)} mkN',
            )

    return m_cracking_pos, m_cracking_neg


def calculate_width_at_z(geometry: CompoundGeometry | GenericSection, z):
    """Calculate the width of a section or a section part at a certain level of z.

    Args:
           geometry: the CompoundGeometry (part) or the whole section (GenericSection)
           z: z value to evaluate the width in mm

    Returns:
           Returns the width at z value
    """

    def width_fibre(geometry, z):
        # Create a horizontal line at the specified y-coordinate
        min_y, min_z, max_y, max_z = geometry.bounds
        horizontal_line = LineString([(min_y, z), (max_y, z)])

        # Find the intersection of the horizontal line with the geometry
        intersection = geometry.intersection(horizontal_line)

        # If the intersection is a MultiLineString, get the total length
        if intersection.is_empty:
            return 0
        elif intersection.geom_type == 'LineString':
            return intersection.length
        elif intersection.geom_type == 'MultiLineString':
            return sum(segment.length for segment in intersection)
        else:
            return 0

    width = 0
    if isinstance(geometry, Polygon):
        width = width_fibre(geometry, z)
    elif isinstance(geometry, CompoundGeometry):
        for g in geometry.geometries:
            width += width_fibre(g.polygon, z)
    elif isinstance(geometry, GenericSection):
        for g in geometry.geometry.geometries:
            width += width_fibre(g.polygon, z)

    return width


def effective_depth(section: GenericSection, neg_bending: bool = False):
    """Calculate the width of a section or a section part at a certain level of z.

    Args:
           section (GenericSection).
           neg_bending: If negative=True the function return 'd' for negative bending moments
           z: z value to evaluate the width in mm

    Returns:
           Returns the efective depth 'd' in mm
    """
    if not section.geometry.reinforced_concrete:
        raise TypeError('Not a reinforced concrete section')
    min_y, max_y, min_z, max_z = section.geometry.calculate_extents()
    res = section.section_calculator.calculate_bending_strength()
    neutral_axe = -res.eps_a / res.chi_y
    print('na. ', neutral_axe)
    sum_a = 0
    sum_az = 0
    mid_z = 0.5 * (min_z + max_z)

    for pg in section.geometry.point_geometries:
        pg_z = pg._point.y
        if (neg_bending and pg_z < mid_z) or (
            not neg_bending and pg_z > mid_z
        ):
            sum_a += pg._area
            sum_az += pg._area * pg_z

    if (
        sum_a == 0
    ):  # To handle the case where sum_a is zero to avoid division by zero
        return 0

    if neg_bending:
        return max_z - sum_az / sum_a
    else:
        return sum_az / sum_a
