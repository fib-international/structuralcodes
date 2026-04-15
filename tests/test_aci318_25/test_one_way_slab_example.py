"""End-to-end validation: one-way slab design with both paths.

Problem: 4000 psi concrete, Gr 60 rebar, 20 ft span, one end continuous.
Design per 12 in. strip.
"""

import math
import sys

import pytest

sys.modules['triangle'] = type(sys)('triangle')

from shapely.geometry import Polygon  # noqa: E402

import structuralcodes  # noqa: E402
from structuralcodes.codes import aci318_25  # noqa: E402
from structuralcodes.geometry import (  # noqa: E402
    CompoundGeometry,
    PointGeometry,
    SurfaceGeometry,
)
from structuralcodes.materials.concrete import create_concrete  # noqa: E402
from structuralcodes.materials.concrete._concreteACI318_25 import (  # noqa: E402
    ConcreteACI318_25,
)
from structuralcodes.materials.reinforcement import (  # noqa: E402
    create_reinforcement,
)
from structuralcodes.materials.reinforcement._reinforcementACI318_25 import (  # noqa: E402
    ReinforcementACI318_25,
)
from structuralcodes.sections import BeamSection  # noqa: E402

FC = 27.58
FY = 420.0
SPAN = 6096.0
B = 305.0


@pytest.fixture(autouse=True)
def _reset_design_code():
    yield
    structuralcodes.set_design_code(None)


class TestPathAClosedForm:
    """Path A: closed-form ACI equations."""

    def test_min_thickness(self):
        """Minimum slab thickness per ACI 318-25 Table 7.3.1.1."""
        h = aci318_25.min_thickness(SPAN, 'one_end_continuous')
        assert math.isclose(h, SPAN / 24, rel_tol=1e-6)
        assert h > 200

    def test_flexure_design(self):
        """Flexure design gives phi*Mn >= Mu (tension-controlled section)."""
        h = 254.0  # 10 in.
        d = h - 19 - 16 / 2  # 227 mm

        Mu = 40e6  # N-mm
        phi = 0.9
        As = aci318_25.As_required(Mu, phi, FY, FC, B, d)
        As_min = aci318_25.As_min_slab(FY, B, h)
        As_design = max(As, As_min)

        a = aci318_25.stress_block_depth_sr(As_design, FY, FC, B)
        c = aci318_25.neutral_axis_depth(a, aci318_25.beta1(FC))
        eps_t = aci318_25.eps_t_from_c(c, d)
        assert aci318_25.As_max_check(eps_t, FY)
        assert aci318_25.phi_flexure(eps_t, FY) == 0.9

        Mn = aci318_25.Mn_singly_reinforced(As_design, FY, FC, B, d)
        assert phi * Mn >= Mu

    def test_shear_check(self):
        """Shear check confirms no stirrups required for the given Vu."""
        d = 227.0
        rho_w = 0.009

        Vc = aci318_25.Vc_detailed(FC, B, d, rho_w)
        phi_Vc = aci318_25.phi_shear() * Vc

        Vu = 30000.0  # N
        assert not aci318_25.shear_reinforcement_required(Vu, phi_Vc)


class TestPathBSectionIntegrator:
    """Path B: section integrator with ACI materials."""

    def test_section_analysis(self):
        """Section integrator returns positive gross area for valid section."""
        concrete = ConcreteACI318_25(
            fck=FC, constitutive_law='parabolarectangle'
        )
        steel = ReinforcementACI318_25(
            fyk=FY,
            Es=200000,
            ftk=550,
            epsuk=0.05,
            constitutive_law='elasticperfectlyplastic',
        )

        poly = Polygon([(0, 0), (B, 0), (B, 254), (0, 254)])
        surf = SurfaceGeometry(poly, concrete)
        bar1 = PointGeometry(point=(100, 27), diameter=16, material=steel)
        bar2 = PointGeometry(point=(205, 27), diameter=16, material=steel)
        section_geo = CompoundGeometry([surf, bar1, bar2])

        section = BeamSection(section_geo, integrator='marin')
        props = section.gross_properties
        assert props.area > 0

    def test_factory_round_trip(self):
        """Factory functions return ACI 318-25 types under the ACI code."""
        structuralcodes.set_design_code('aci318_25')
        c = create_concrete(fck=FC)
        assert isinstance(c, ConcreteACI318_25)

        r = create_reinforcement(fyk=FY, Es=200000, ftk=550, epsuk=0.05)
        assert isinstance(r, ReinforcementACI318_25)


class TestCrossCheck:
    """Cross-check between paths."""

    def test_mn_agreement(self):
        """Closed-form Mn and integrator Mn should agree within 5%."""
        As = 400.0  # mm2
        d = 227.0
        h = 254.0

        # Path A
        Mn_closed = aci318_25.Mn_singly_reinforced(As, FY, FC, B, d)

        # Path B
        concrete = ConcreteACI318_25(
            fck=FC, constitutive_law='parabolarectangle'
        )
        steel = ReinforcementACI318_25(
            fyk=FY,
            Es=200000,
            ftk=550,
            epsuk=0.05,
            constitutive_law='elasticperfectlyplastic',
        )

        poly = Polygon([(0, 0), (B, 0), (B, h), (0, h)])
        surf = SurfaceGeometry(poly, concrete)
        bar1 = PointGeometry(point=(100, 27), diameter=16, material=steel)
        bar2 = PointGeometry(point=(205, 27), diameter=16, material=steel)
        section_geo = CompoundGeometry([surf, bar1, bar2])

        section = BeamSection(section_geo, integrator='marin')
        calc = section.section_calculator
        strain = calc.find_equilibrium_fixed_pivot(
            geom=section.geometry,
            n=0,
            yielding=True,
        )
        N, My, Mz, data = (
            calc.integrator.integrate_strain_response_on_geometry(
                geo=section.geometry,
                strain=strain,
                integrate='stress',
            )
        )
        Mn_integrator = abs(My)

        assert math.isclose(Mn_closed, Mn_integrator, rel_tol=0.05)
