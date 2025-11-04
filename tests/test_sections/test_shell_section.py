"""Tests for the Shell Section."""

import math

import numpy as np
import pytest

from structuralcodes.geometry import (
    RectangularGeometry,
    add_reinforcement_line,
)
from structuralcodes.geometry._shell_geometry import (
    ShellGeometry,
    ShellReinforcement,
)
from structuralcodes.materials.basic import GenericMaterial
from structuralcodes.materials.concrete import ConcreteEC2_2004
from structuralcodes.materials.constitutive_laws import (
    ConcreteSmearedCracking,
    ConstantPoissonReduction,
    Elastic2D,
    ElasticPlastic,
    GeneralVecchioCollins,
    ParabolaRectangle,
)
from structuralcodes.materials.reinforcement import ReinforcementEC2_2004
from structuralcodes.sections import GenericSection, ShellSection

# Membrane and bending-strain parameter ranges
eps_x = np.linspace(0.0, 1.0e-3, 2)
eps_y = np.linspace(0.0, 1.0e-3, 2)
eps_xy = np.linspace(0.0, 1.0e-3, 2)
chi_x = np.linspace(0.0, 1.0e-6, 2)
chi_y = np.linspace(0.0, 1.0e-6, 2)
chi_xy = np.linspace(0.0, 1.0e-6, 2)


@pytest.mark.parametrize('chi_xy', chi_xy)
@pytest.mark.parametrize('chi_y', chi_y)
@pytest.mark.parametrize('chi_x', chi_x)
@pytest.mark.parametrize('eps_xy', eps_xy)
@pytest.mark.parametrize('eps_y', eps_y)
@pytest.mark.parametrize('eps_x', eps_x)
@pytest.mark.parametrize(
    'Ec, nu, thickness',
    [(30000, 0.20, 200), (30000, 0.20, 600)],
)
def test_integrate_strain_profile(
    Ec,
    nu,
    thickness,
    eps_x,
    eps_y,
    eps_xy,
    chi_x,
    chi_y,
    chi_xy,
):
    """Elastic plate: strains → stress-resultants."""
    constitutive_law = Elastic2D(Ec, nu)
    material = GenericMaterial(density=2500, constitutive_law=constitutive_law)
    shell = ShellSection(ShellGeometry(thickness, material), mesh_size=1e-3)

    # Numerically integrated forces/moments
    R = shell.section_calculator.integrate_strain_profile(
        np.array([eps_x, eps_y, eps_xy, chi_x, chi_y, chi_xy])
    )

    # Analytical forces/moments
    A_membrane = Ec * thickness / (1 - nu**2)  # Coefficient in Hooke's law
    A_bending = Ec * thickness**3 / 12 / (1 - nu**2)  # Plate bending stiffness
    A_membrane_shear = A_membrane * (1 - nu) / 2
    A_bending_shear = A_bending * (1 - nu) / 2

    Nx = A_membrane * (eps_x + nu * eps_y)
    Ny = A_membrane * (eps_y + nu * eps_x)
    Nxy = A_membrane_shear * eps_xy
    Mx = A_bending * (chi_x + nu * chi_y)
    My = A_bending * (chi_y + nu * chi_x)
    Mxy = A_bending_shear * chi_xy

    expected = np.array([Nx, Ny, Nxy, Mx, My, Mxy])

    assert np.allclose(R, expected)


@pytest.mark.parametrize('nu', ((0, 0.2)))
def test_integrate_strain_profile_stiffness_matrix(nu):
    """Elastic plate: stiffness matrix."""
    E = 30000
    t = 200
    constitutive_law = Elastic2D(E, nu)
    material = GenericMaterial(density=2500, constitutive_law=constitutive_law)
    shell = ShellSection(ShellGeometry(t, material=material), mesh_size=1e-3)
    K = shell.section_calculator.integrate_strain_profile(
        np.array([1e-3] * 3 + [1e-6] * 3), integrate='modulus'
    )

    A11 = E * t / (1 - nu**2)
    A12 = nu * A11
    A33 = E * t / (2 * (1 + nu))
    A44 = E * t**3 / 12 / (1 - nu**2)
    A14 = nu * A44
    A55 = E * t**3 / 12 / (2 * (1 + nu))

    A = np.array(
        [
            [A11, A12, 0, 0, 0, 0],
            [A12, A11, 0, 0, 0, 0],
            [0, 0, A33, 0, 0, 0],
            [0, 0, 0, A44, A14, 0],
            [0, 0, 0, A14, A44, 0],
            [0, 0, 0, 0, 0, A55],
        ]
    )

    assert np.allclose(K, A)


def test_wrong_integrator():
    """Unknown keyword must raise ValueError."""
    constitutive_law = Elastic2D(30000, 0.20)
    material = GenericMaterial(density=2500, constitutive_law=constitutive_law)
    shell = ShellSection(ShellGeometry(200, material=material))
    with pytest.raises(ValueError):
        shell.section_calculator.integrate_strain_profile(
            np.zeros(6), integrate='tangent'
        )


# Loads for reverse strain-solution test
nx = np.linspace(-1e5, 1e5, 2)
ny = np.linspace(-1e5, 1e5, 2)
nxy = np.linspace(-1e5, 1e5, 2)
mx = np.linspace(-1e8, 1e8, 2)
my = np.linspace(-1e8, 1e8, 2)
mxy = np.linspace(-1e8, 1e8, 2)


@pytest.mark.parametrize('mxy', mxy)
@pytest.mark.parametrize('my', my)
@pytest.mark.parametrize('mx', mx)
@pytest.mark.parametrize('nxy', nxy)
@pytest.mark.parametrize('ny', ny)
@pytest.mark.parametrize('nx', nx)
@pytest.mark.parametrize('Ec, nu, t', [(30000, 0.20, 200), (30000, 0.20, 600)])
def test_elastic_strain_profile(Ec, nu, t, nx, ny, nxy, mx, my, mxy):
    """Loads → strains for an isotropic plate."""
    constitutive_law = Elastic2D(Ec, nu)
    material = GenericMaterial(density=2500, constitutive_law=constitutive_law)
    shell = ShellSection(ShellGeometry(t, material=material), mesh_size=1e-3)
    eps = shell.section_calculator.calculate_strain_profile(
        nx, ny, nxy, mx, my, mxy
    )

    sig_x = nx / t
    sig_y = ny / t
    tau_xy = nxy / t
    eps_x = (sig_x - nu * sig_y) / Ec
    eps_y = (sig_y - nu * sig_x) / Ec
    eps_xy = 2 * (1 + nu) * tau_xy / Ec
    D = Ec * t**3 / 12
    chi_x = (mx - nu * my) / D
    chi_y = (my - nu * mx) / D
    chi_xy = 2 * (1 + nu) * mxy / D
    expected = np.array([eps_x, eps_y, eps_xy, chi_x, chi_y, chi_xy])

    assert np.allclose(eps, expected)


def test_default_equals_explicit_mesh_size():
    """Default mesh_size (0.01) equals explicit 0.01."""
    t = 200
    constitutive_law = Elastic2D(30000, 0.20)
    material = GenericMaterial(density=2500, constitutive_law=constitutive_law)
    shell_0 = ShellSection(ShellGeometry(t, material=material))
    K0 = shell_0.section_calculator.integrate_strain_profile(
        np.zeros(6), integrate='modulus'
    )
    shell_1 = ShellSection(ShellGeometry(t, material=material), mesh_size=0.01)
    K1 = shell_1.section_calculator.integrate_strain_profile(
        np.zeros(6), integrate='modulus'
    )
    assert np.allclose(K0, K1)


@pytest.mark.parametrize('invalid', [-0.2, 0.0, 1.5])
def test_invalid_mesh_size_raises(invalid):
    """mesh_size outside (0,1] → ValueError."""
    t = 200
    constitutive_law = Elastic2D(30000, 0.20)
    material = GenericMaterial(density=2500, constitutive_law=constitutive_law)
    shell = ShellSection(ShellGeometry(t, material), mesh_size=invalid)
    with pytest.raises(ValueError):
        shell.section_calculator.integrate_strain_profile(
            np.zeros(6), integrate='stress'
        )


@pytest.mark.parametrize(
    'nx,nxy,expected',
    [
        (
            0,
            500,
            np.array([1.43163e-3, 1.919842e-3, 3.466593e-3, 0, 0, 0]),
        ),
        (
            -500,
            500,
            np.array([0.618194e-3, 1.476430e-3, 2.065709e-3, 0, 0, 0]),
        ),
        (
            500,
            500,
            np.array([2.446863e-3, 2.283832e-3, 4.894237e-3, 0, 0, 0]),
        ),
    ],
)
def test_parabola_zero_initial_nu(nx, nxy, expected):
    """Test parabola rectangle with nu = 0."""
    uniaxial_compression = ParabolaRectangle(fc=45)
    strength_reduction = GeneralVecchioCollins(c_1=0.8, c_2=100)
    poisson_reduction = ConstantPoissonReduction(initial_nu=0)
    parabola_rectangle = ConcreteSmearedCracking(
        uniaxial_compression=uniaxial_compression,
        strength_reduction_lateral_cracking=strength_reduction,
        poisson_reduction=poisson_reduction,
    )
    elastic_plastic = ElasticPlastic(200000, 500)
    concrete = GenericMaterial(
        density=2500, constitutive_law=parabola_rectangle
    )
    reinforcement = ReinforcementEC2_2004(
        fyk=500,
        Es=200000,
        ftk=500,
        epsuk=3e-2,
        constitutive_law=elastic_plastic,
    )
    geo = ShellGeometry(350, concrete)
    Asx1 = ShellReinforcement(-132, 1, 200, 16, reinforcement, 0)
    Asx2 = ShellReinforcement(132, 1, 200, 16, reinforcement, 0)
    Asy1 = ShellReinforcement(-118, 1, 200, 12, reinforcement, np.pi / 2)
    Asy2 = ShellReinforcement(118, 1, 200, 12, reinforcement, np.pi / 2)
    geo.add_reinforcement([Asx1, Asx2, Asy1, Asy2])
    section = ShellSection(geo, mesh_size=0.5)
    calculator = section.section_calculator
    strain = calculator.calculate_strain_profile(
        nx, 0, nxy, 0, 0, 0, max_iter=150
    )

    assert np.allclose(strain, expected)


@pytest.mark.parametrize(
    'nx,nxy,expected',
    [
        (
            0,
            500,
            np.array([1.431675e-3, 1.919788e-3, 3.466592e-3, 0, 0, 0]),
        ),
        (
            -500,
            500,
            np.array([0.618204e-3, 1.476423e-3, 2.065719e-3, 0, 0, 0]),
        ),
        (
            500,
            500,
            np.array([2.446899e-3, 2.283762e-3, 4.894200e-3, 0, 0, 0]),
        ),
    ],
)
def test_parabola_initial_nu(nx, nxy, expected):
    """Test parabola rectangle with nu = 0.2."""
    uniaxial_compression = ParabolaRectangle(fc=45)
    strength_reduction = GeneralVecchioCollins(c_1=0.8, c_2=100)
    poisson_reduction = ConstantPoissonReduction(initial_nu=0.2)
    parabola_rectangle = ConcreteSmearedCracking(
        uniaxial_compression=uniaxial_compression,
        strength_reduction_lateral_cracking=strength_reduction,
        poisson_reduction=poisson_reduction,
    )
    elastic_plastic = ElasticPlastic(200000, 500)
    concrete = GenericMaterial(
        density=2500, constitutive_law=parabola_rectangle
    )
    reinforcement = ReinforcementEC2_2004(
        fyk=500,
        Es=200000,
        ftk=500,
        epsuk=3e-2,
        constitutive_law=elastic_plastic,
    )
    geo = ShellGeometry(350, concrete)
    Asx1 = ShellReinforcement(-132, 1, 200, 16, reinforcement, 0)
    Asx2 = ShellReinforcement(132, 1, 200, 16, reinforcement, 0)
    Asy1 = ShellReinforcement(-118, 1, 200, 12, reinforcement, np.pi / 2)
    Asy2 = ShellReinforcement(118, 1, 200, 12, reinforcement, np.pi / 2)
    geo.add_reinforcement([Asx1, Asx2, Asy1, Asy2])
    section = ShellSection(geo, mesh_size=0.5)
    strain = section.section_calculator.calculate_strain_profile(
        nx,
        0,
        nxy,
        0,
        0,
        0,
        max_iter=200,
    )

    assert np.allclose(strain, expected)


def test_exceed_max_iterations():
    """Test that the maximum number of iterations is exceeded."""
    uniaxial_compression = ParabolaRectangle(fc=45)
    strength_reduction = GeneralVecchioCollins(c_1=0.8, c_2=100)
    poisson_reduction = ConstantPoissonReduction(initial_nu=0)
    parabola_rectangle = ConcreteSmearedCracking(
        uniaxial_compression=uniaxial_compression,
        strength_reduction_lateral_cracking=strength_reduction,
        poisson_reduction=poisson_reduction,
    )
    elastic_plastic = ElasticPlastic(200000, 500)
    concrete = GenericMaterial(
        density=2500, constitutive_law=parabola_rectangle
    )
    reinforcement = ReinforcementEC2_2004(
        fyk=500,
        Es=200000,
        ftk=500,
        epsuk=3e-2,
        constitutive_law=elastic_plastic,
    )
    geo = ShellGeometry(350, concrete)
    Asx1 = ShellReinforcement(-132, 1, 200, 16, reinforcement, 0)
    Asx2 = ShellReinforcement(132, 1, 200, 16, reinforcement, 0)
    Asy1 = ShellReinforcement(-118, 1, 200, 12, reinforcement, np.pi / 2)
    Asy2 = ShellReinforcement(118, 1, 200, 12, reinforcement, np.pi / 2)
    geo.add_reinforcement([Asx1, Asx2, Asy1, Asy2])
    section = ShellSection(geo)

    with pytest.raises(
        StopIteration, match='Maximum number of iterations reached'
    ):
        section.section_calculator.calculate_strain_profile(
            -1000,
            0,
            1000,
            0,
            0,
            0,
        )


@pytest.mark.parametrize(
    'axial_force, with_reinforcement',
    (
        (-1e6, False),
        (-1e6, True),
        (1e5, True),
    ),
)
def test_compare_uniaxial_with_generic_section_reinforcement(  # noqa: PLR0915
    axial_force: float, with_reinforcement: bool
):
    """Compare the uniaxial response of the shell section with the generic
    section, with reinforcement.
    """
    # Arrange
    # Geometry
    width = 1000
    height = 350
    cover = 35
    diameter = 16
    spacing = 200

    # Material parameters
    fc = 45  # Concrete compressive strength
    nu = 0.0  # Poisson's ratio
    fyk = 500  # Steel yield strength
    Es = 200000  # Steel modulus of elasticity

    # Create a GenericSection
    reinforcement = ReinforcementEC2_2004(
        fyk=fyk,
        Es=Es,
        epsuk=3e-2,
        ftk=fyk,
        constitutive_law='elasticperfectlyplastic',
    )
    parabola_rectangle = ParabolaRectangle(fc=fc)
    concrete_for_generic = ConcreteEC2_2004(
        fck=fc, constitutive_law=parabola_rectangle
    )
    generic_geo = RectangularGeometry(
        height=height, width=width, material=concrete_for_generic
    )

    if with_reinforcement:
        # Bottom reinforcement
        generic_geo = add_reinforcement_line(
            generic_geo,
            (
                -width / 2 + cover + diameter / 2,
                -height / 2 + cover + diameter / 2,
            ),
            (
                width / 2 - cover - diameter / 2,
                -height / 2 + cover + diameter / 2,
            ),
            diameter=diameter,
            material=reinforcement,
            s=spacing,
        )

        # Top reinforcement
        generic_geo = add_reinforcement_line(
            generic_geo,
            (
                -width / 2 + cover + diameter / 2,
                height / 2 - cover - diameter / 2,
            ),
            (
                width / 2 - cover - diameter / 2,
                height / 2 - cover - diameter / 2,
            ),
            diameter=diameter,
            material=reinforcement,
            s=spacing,
        )

    generic_sec = GenericSection(geometry=generic_geo, integrator='fiber')

    # Create a ShellSection
    uniaxial_compression = ParabolaRectangle(fc=fc)
    strength_reduction = GeneralVecchioCollins(c_1=0.8, c_2=100)
    poisson_reduction = ConstantPoissonReduction(initial_nu=nu)
    smeared_cracking = ConcreteSmearedCracking(
        uniaxial_compression=uniaxial_compression,
        strength_reduction_lateral_cracking=strength_reduction,
        poisson_reduction=poisson_reduction,
    )
    concrete_for_shell = GenericMaterial(
        density=2500, constitutive_law=smeared_cracking
    )
    shell_geo = ShellGeometry(material=concrete_for_shell, thickness=height)

    if with_reinforcement:
        # Create reinforcement for the shell section
        z_btm = -height / 2 + cover + diameter / 2  # -132
        z_top = height / 2 - cover - diameter / 2  # 132

        # Bottom reinforcement
        shell_rein_btm = ShellReinforcement(
            z_btm, 1, spacing, diameter, reinforcement, 0
        )

        # Top reinforcement
        shell_rein_top = ShellReinforcement(
            z_top, 1, spacing, diameter, reinforcement, 0
        )

        # Add reinforcement to the shell geometry
        shell_geo.add_reinforcement([shell_rein_btm, shell_rein_top])

    shell_sec = ShellSection(geometry=shell_geo, mesh_size=0.5)

    # Calculate "exact" solution
    if axial_force > 0:
        exact_longitudinal_strain = axial_force / (
            generic_sec.gross_properties.area_reinforcement * Es
        )
    else:
        exact_longitudinal_strain = 0
        num_iter = 0
        area = width * height
        while True:
            num_iter += 1
            if num_iter >= 40:
                break
            exact_longitudinal_strain = axial_force / (
                area
                * concrete_for_generic.constitutive_law.get_secant(
                    exact_longitudinal_strain
                )
                + generic_sec.gross_properties.area_reinforcement * Es
            )

    # Act
    # Calculate strains
    generic_strain = generic_sec.section_calculator.calculate_strain_profile(
        n=axial_force, my=0.0, mz=0.0
    )
    shell_strain = shell_sec.section_calculator.calculate_strain_profile(
        nx=axial_force / width,
        ny=0.0,
        nxy=0.0,
        mx=0.0,
        my=0.0,
        mxy=0.0,
    )

    # Compare with exact solution
    # Note that the tolerance in isclose is set rather loose because the
    # iterative solver of the shell section converges slower than the generic
    # section.
    assert math.isclose(
        shell_strain[0],
        exact_longitudinal_strain,
        rel_tol=1e-4,
    )

    # Compare GenericSection and ShellSection
    assert math.isclose(
        generic_strain[0],
        shell_strain[0],
        rel_tol=1e-4,
    )


@pytest.mark.parametrize('nu', [0.0, 0.2])
@pytest.mark.parametrize(
    'strain',
    [
        (-2e-3, 0.0, 0.0),
        (-2e-3, -2e-3, 0.0),
        (0.0, 0.0, 0.75e-3),
        (-2e-3, 0.0, 0.75e-3),
        (-2e-3, -2e-3, 0.75e-3),
    ],
)
def test_compare_constitutive_law_and_section(strain, nu):
    """Compare the stress calculated with a constitutive law with the stress
    resultant from the section.
    """
    # Arrange
    fck = 45
    thickness = 450
    uniaxial_compression = ParabolaRectangle(fc=fck)
    strength_reduction = GeneralVecchioCollins(c_1=0.8, c_2=100)
    poisson_reduction = ConstantPoissonReduction(initial_nu=nu)
    constitutive_law = ConcreteSmearedCracking(
        uniaxial_compression=uniaxial_compression,
        strength_reduction_lateral_cracking=strength_reduction,
        poisson_reduction=poisson_reduction,
    )
    concrete = GenericMaterial(density=2500, constitutive_law=constitutive_law)
    shell_geometry = ShellGeometry(thickness=thickness, material=concrete)
    shell_section = ShellSection(geometry=shell_geometry)

    # Act
    stress = constitutive_law.get_stress(eps=strain)
    stress_resultant = (
        shell_section.section_calculator.integrate_strain_profile(
            strain=[*strain, 0.0, 0.0, 0.0],
            integrate='stress',
        )
    )

    # Assert
    assert np.allclose(stress * thickness, stress_resultant[:3])
