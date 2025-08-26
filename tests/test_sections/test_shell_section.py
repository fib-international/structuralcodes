"""Tests for the Shell Section."""

import numpy as np
import pytest

from structuralcodes.geometry._shell_geometry import (
    ShellGeometry,
    ShellReinforcement,
)
from structuralcodes.materials.constitutive_laws import (
    Elastic2D,
    ElasticPlastic,
    ParabolaRectangle2D,
)
from structuralcodes.sections import ShellSection

# Membrane and bending-strain parameter ranges
eps_x = np.linspace(0.0, 1.0e-3, 2)
eps_y = np.linspace(0.0, 1.0e-3, 2)
eps_xy = np.linspace(0.0, 1.0e-3, 2)
chi_x = np.linspace(0.0, 1.0e-6, 2)
chi_y = np.linspace(0.0, 1.0e-6, 2)
chi_xy = np.linspace(0.0, 1.0e-6, 2)


@pytest.mark.parametrize('eps_x', eps_x)
@pytest.mark.parametrize('eps_y', eps_y)
@pytest.mark.parametrize('eps_xy', eps_xy)
@pytest.mark.parametrize('chi_x', chi_x)
@pytest.mark.parametrize('chi_y', chi_y)
@pytest.mark.parametrize('chi_xy', chi_xy)
@pytest.mark.parametrize(
    'Ec, nu, thickness',
    [(30000, 0.20, 200), (30000, 0.20, 600)],
)
def test_integrate_strain_profile(
    Ec, nu, thickness, eps_x, eps_y, eps_xy, chi_x, chi_y, chi_xy
):
    """Elastic plate: strains → stress-resultants."""
    material = Elastic2D(Ec, nu)
    shell = ShellSection(ShellGeometry(thickness, material))

    # Numerically integrated forces/moments
    R = shell.section_calculator.integrate_strain_profile(
        np.array([eps_x, eps_y, eps_xy, chi_x, chi_y, chi_xy])
    )

    # Analytical forces/moments
    A = Ec * thickness / (1 - nu**2)
    B = Ec * thickness**3 / 12 / (1 - nu**2)
    Nx = A * (eps_x + nu * eps_y)
    Ny = A * (eps_y + nu * eps_x)
    Nxy = Ec * thickness / (2 * (1 + nu)) * eps_xy
    Mx = B * (chi_x + nu * chi_y)
    My = B * (chi_y + nu * chi_x)
    Mxy = Ec * thickness**3 / 12 / (2 * (1 + nu)) * chi_xy

    expected = np.array([Nx, Ny, Nxy, Mx, My, Mxy])

    assert np.allclose(R, expected, rtol=1e-2, atol=1e-2)


def test_integrate_strain_profile_tangent(E=30000, nu=0.20, t=200):
    """Elastic plate: tangent stiffness matrix."""
    shell = ShellSection(ShellGeometry(t, Elastic2D(E, nu)))
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

    assert np.allclose(K, A, rtol=1e-2)


def test_wrong_integrator():
    """Unknown keyword must raise ValueError."""
    shell = ShellSection(ShellGeometry(200, Elastic2D(30000, 0.20)))
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


@pytest.mark.parametrize('nx', nx)
@pytest.mark.parametrize('ny', ny)
@pytest.mark.parametrize('nxy', nxy)
@pytest.mark.parametrize('mx', mx)
@pytest.mark.parametrize('my', my)
@pytest.mark.parametrize('mxy', mxy)
@pytest.mark.parametrize('Ec, nu, t', [(30000, 0.20, 200), (30000, 0.20, 600)])
def test_elastic_strain_profile(Ec, nu, t, nx, ny, nxy, mx, my, mxy):
    """Loads → strains for an isotropic plate."""
    shell = ShellSection(ShellGeometry(t, Elastic2D(Ec, nu)))
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

    assert np.allclose(eps, expected, rtol=1e-4, atol=1e-3)


def test_default_equals_explicit_mesh_size(t=200):
    """Default mesh_size (0.01) equals explicit 0.01."""
    shell_0 = ShellSection(ShellGeometry(t, Elastic2D(30_000, 0.20)))
    K0 = shell_0.section_calculator.integrate_strain_profile(
        np.zeros(6), integrate='modulus'
    )
    shell_1 = ShellSection(
        ShellGeometry(t, Elastic2D(30_000, 0.20)), mesh_size=0.01
    )
    K1 = shell_1.section_calculator.integrate_strain_profile(
        np.zeros(6), integrate='modulus'
    )
    assert np.allclose(K0, K1)


@pytest.mark.parametrize('invalid', [-0.2, 0.0, 1.5])
def test_invalid_mesh_size_raises(invalid, t=200):
    """mesh_size outside (0,1] → ValueError."""
    shell = ShellSection(
        ShellGeometry(t, Elastic2D(30000, 0.20)), mesh_size=invalid
    )
    with pytest.raises(ValueError):
        shell.section_calculator.integrate_strain_profile(
            np.zeros(6), integrate='stress'
        )


def test_parabola_section():
    """ParabolaRectangle2D with mesh_size 0.5."""
    concrete = ParabolaRectangle2D(35, nu=0)
    reinforcement = ElasticPlastic(200000, 500)

    Asx1 = ShellReinforcement(-157, 1, 300, 16, reinforcement, 0)

    geo = ShellGeometry(400, concrete)
    geo.add_reinforcement([Asx1])

    section = ShellSection(geo, mesh_size=0.1)
    calculator = section.section_calculator
    strain = calculator.calculate_strain_profile(-1000, 0, 0, 0, 0, 0)

    assert np.allclose(
        strain,
        np.array([-7.20105642e-05, 0, 0, 1.07581081e-08, 0, 0]),
        rtol=1e-6,
        atol=1e-7,
    )


@pytest.mark.parametrize(
    'nx,nxy,expected,tol',
    [
        (
            0,
            1000,
            np.array([3.046432e-05, 3.079979e-05, 1.914840e-04, 0, 0, 0]),
            1e-4,
        ),
        (
            -1000,
            1000,
            np.array([-5.285263e-05, 2.741669e-05, 1.649706e-04, 0, 0, 0]),
            0.001,
        ),
        (
            1000,
            1000,
            np.array([1.330415e-04, 2.785210e-05, 2.186547e-04, 0, 0, 0]),
            0.001,
        ),
    ],
)
def test_parabola_cracked(nx, nxy, expected, tol):
    """Test parabola rectangle with nu = 0."""
    concrete = ParabolaRectangle2D(45, nu=0)
    reinforcement = ElasticPlastic(200000, 500)
    geo = ShellGeometry(350, concrete)
    Asx1 = ShellReinforcement(-132, 1, 200, 16, reinforcement, 0)
    Asx2 = ShellReinforcement(132, 1, 200, 16, reinforcement, 0)
    Asy1 = ShellReinforcement(-118, 1, 200, 12, reinforcement, np.pi / 2)
    Asy2 = ShellReinforcement(118, 1, 200, 12, reinforcement, np.pi / 2)
    geo.add_reinforcement([Asx1, Asx2, Asy1, Asy2])
    section = ShellSection(geo)
    calculator = section.section_calculator
    strain = calculator.calculate_strain_profile(
        nx, 0, nxy, 0, 0, 0, initial=True, tol=tol
    )

    assert np.allclose(
        strain,
        expected,
        rtol=1e-6,
        atol=1e-7,
    )


@pytest.mark.parametrize(
    'nx,nxy,expected,tol',
    [
        (
            0,
            1000,
            np.array([2.923437e-05, 2.961959e-05, 2.150748e-04, 0, 0, 0]),
            1e-4,
        ),
        (
            -1000,
            1000,
            np.array([-5.332211e-05, 3.874431e-05, 1.897691e-04, 0, 0, 0]),
            0.001,
        ),
        (
            1000,
            1000,
            np.array([1.302065e-04, 1.461867e-05, 2.411031e-04, 0, 0, 0]),
            0.001,
        ),
    ],
)
def test_parabola_uncracked(nx, nxy, expected, tol):
    """Test parabola rectangle with nu = 0.2."""
    concrete = ParabolaRectangle2D(45)
    reinforcement = ElasticPlastic(200000, 500)
    geo = ShellGeometry(350, concrete)
    Asx1 = ShellReinforcement(-132, 1, 200, 16, reinforcement, 0)
    Asx2 = ShellReinforcement(132, 1, 200, 16, reinforcement, 0)
    Asy1 = ShellReinforcement(-118, 1, 200, 12, reinforcement, np.pi / 2)
    Asy2 = ShellReinforcement(118, 1, 200, 12, reinforcement, np.pi / 2)
    geo.add_reinforcement([Asx1, Asx2, Asy1, Asy2])
    section = ShellSection(geo)
    calculator = section.section_calculator
    strain = calculator.calculate_strain_profile(
        nx, 0, nxy, 0, 0, 0, initial=True, tol=tol
    )

    assert np.allclose(
        strain,
        expected,
        rtol=1e-6,
        atol=1e-7,
    )


def test_exceed_max_iterations():
    """Test that the maximum number of iterations is exceeded."""
    concrete = ParabolaRectangle2D(45)
    reinforcement = ElasticPlastic(200000, 500)
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
