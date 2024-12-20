import os
import sys

import numpy as np
import pytest
from shapely import Polygon

from structuralcodes.geometry import CompoundGeometry, SurfaceGeometry
from structuralcodes.materials.concrete import ConcreteEC2_2004
from structuralcodes.materials.reinforcement import ReinforcementEC2_2004
from structuralcodes.sections import GenericSection
from structuralcodes.sections._reinforcement import add_reinforcement_line

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(os.path.dirname(current_dir))
sys.path.append(parent_dir)

from SLS_section_response import (
    calculate_cracking_moment,
    calculate_strain_profile,
    calculate_strain_profile_batch,
    calculate_width_at_z,
    effective_depth,
)


def test_calculate_width_at_z():
    """Test width_at_z in Polygon, CompoundGeometry and GenericSection."""
    # Test with a simple polygon
    polygon = Polygon([(0, 0), (4, 0), (4, 4), (0, 4)])
    assert calculate_width_at_z(polygon, 2) == 4

    # Test with a polygon where z is outside the bounds
    assert calculate_width_at_z(polygon, 5) == 0

    # Test with a compound geometry
    compound_geometry = CompoundGeometry(
        [
            SurfaceGeometry(
                Polygon([(0, 0), (2, 0), (2, 2), (0, 2)]), ConcreteEC2_2004(25)
            ),
            SurfaceGeometry(
                Polygon([(3, 0), (5, 0), (5, 2), (3, 2)]), ConcreteEC2_2004(25)
            ),
        ]
    )
    assert calculate_width_at_z(compound_geometry, 1) == 4

    # Test with a generic section
    generic_section = GenericSection(compound_geometry)
    assert calculate_width_at_z(generic_section, 1) == 4

    # Test with an empty polygon
    empty_polygon = Polygon()
    assert calculate_width_at_z(empty_polygon, 1) == 0

    # Test with an reinforced concrete section
    concrete = ConcreteEC2_2004(25)
    reinforcemnet = ReinforcementEC2_2004(
        fyk=500,
        Es=200000,
        density=7850,
        ftk=500,
        epsuk=0.07,
    )
    poly = Polygon(((0, 0), (0, 500), (350, 500), (350, 0)))
    geo = SurfaceGeometry(poly, concrete)
    geo = add_reinforcement_line(
        geo, (50, 50), (300, 50), 12, reinforcemnet, n=3
    )
    geo = add_reinforcement_line(
        geo, (50, 450), (300, 450), 20, reinforcemnet, n=6
    )
    generic_section = GenericSection(geo)
    assert calculate_width_at_z(generic_section, 1) == 350


def test_effective_depth():
    """Test effective_depth 'd' in GenericSection."""
    concrete = ConcreteEC2_2004(25)
    reinforcemnet = ReinforcementEC2_2004(
        fyk=500,
        Es=200000,
        density=7850,
        ftk=500,
        epsuk=0.07,
    )
    poly = Polygon(((0, 0), (0, 500), (350, 500), (350, 0)))
    geo = SurfaceGeometry(poly, concrete)
    geo = add_reinforcement_line(
        geo, (50, 60), (300, 60), 12, reinforcemnet, n=3
    )
    geo = add_reinforcement_line(
        geo, (50, 450), (300, 450), 20, reinforcemnet, n=6
    )
    section = GenericSection(geo)

    # Test for positive bending (neg_bending=False)
    effective_d = effective_depth(section, neg_bending=False)
    assert effective_d == 450, f'Expected 450, got {effective_d}'

    # Test for negative bending (neg_bending=True)
    effective_d_neg = effective_depth(section, neg_bending=True)
    assert effective_d_neg == 440, f'Expected 440, got {effective_d_neg}'


def test_calculate_cracking_moment():
    """Test cracking moments in reinforced concrete section."""
    concrete = ConcreteEC2_2004(25)
    reinforcemnet = ReinforcementEC2_2004(
        fyk=500,
        Es=200000,
        density=7850,
        ftk=500,
        epsuk=0.07,
    )
    # Create section
    poly = Polygon(((0, 0), (350, 0), (350, 500), (0, 500)))
    geo = SurfaceGeometry(poly, concrete)
    geo = add_reinforcement_line(
        geo, (50, 50), (300, 50), 12, reinforcemnet, n=3
    )
    geo = add_reinforcement_line(
        geo, (50, 450), (300, 450), 20, reinforcemnet, n=6
    )
    sec = GenericSection(geo)
    mcr_pos, mcr_neg = calculate_cracking_moment(sec, plot=False)
    mcr_pos = mcr_pos * 1e-6
    mcr_neg = mcr_neg * 1e-6
    assert mcr_pos == pytest.approx(
        46.286, abs=1
    ), f'Expected approximately 46.286, got {mcr_pos}'
    assert mcr_neg == pytest.approx(
        -41.889, abs=1
    ), f'Expected approximately -41.889, got {mcr_neg}'


def test_calculate_strain_profile():
    """Test strain profile given N and M."""
    concrete = ConcreteEC2_2004(25)
    reinforcemnet = ReinforcementEC2_2004(
        fyk=500,
        Es=200000,
        density=7850,
        ftk=500,
        epsuk=0.07,
    )
    # Create section
    poly = Polygon(((0, 0), (350, 0), (350, 500), (0, 500)))
    geo = SurfaceGeometry(poly, concrete)
    geo = add_reinforcement_line(
        geo, (50, 50), (300, 50), 12, reinforcemnet, n=3
    )
    geo = add_reinforcement_line(
        geo, (50, 450), (300, 450), 20, reinforcemnet, n=6
    )
    sec = GenericSection(geo)
    eps, chi = calculate_strain_profile(sec, 0, 200 * 1e6)  # m = 200 mkN
    eps = eps * 1e3
    chi = chi * 1e6
    assert eps == pytest.approx(
        -1.017, abs=1e-2
    ), f'Expected approximately -1.017 mm/m, got {eps}'
    assert chi == pytest.approx(
        5.322, abs=1e-2
    ), f'Expected approximately 5.322 1/km, got {chi}'


def test_calculate_strain_profile_batch():
    """Test strain profile given N and np.array of M."""
    concrete = ConcreteEC2_2004(25)
    reinforcemnet = ReinforcementEC2_2004(
        fyk=500,
        Es=200000,
        density=7850,
        ftk=500,
        epsuk=0.07,
    )
    # Create section
    poly = Polygon(((0, 0), (350, 0), (350, 500), (0, 500)))
    geo = SurfaceGeometry(poly, concrete)
    geo = add_reinforcement_line(
        geo, (50, 50), (300, 50), 12, reinforcemnet, n=3
    )
    geo = add_reinforcement_line(
        geo, (50, 450), (300, 450), 20, reinforcemnet, n=6
    )
    sec = GenericSection(geo)

    my_ed = np.array([-50, 50, 100]) * 1e6  # mkN
    n_ed = 0 * 1e3  # kN
    eps_batch, chi_batch = calculate_strain_profile_batch(sec, n_ed, my_ed)
    for i, eps1 in enumerate(eps_batch):
        chi1 = chi_batch[i]
        eps2, chi2 = calculate_strain_profile(sec, n_ed, my_ed[i])
        print(eps1, eps2)
        print(chi1, chi2)
        assert eps1 == pytest.approx(
            eps2, rel=1e-2, abs=1e-8
        ), f'Expected approximately {eps2} , got {eps1}'
        assert chi1 == pytest.approx(
            chi2, rel=1e-2, abs=1e-8
        ), f'Expected approximately {chi2} , got {chi1}'
