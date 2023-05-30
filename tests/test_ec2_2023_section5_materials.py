"""Tests for _section5_materials module"""
import math

import pytest

from structuralcodes.codes.ec2_2023 import _section5_materials


@pytest.mark.parametrize(
    'fck, delta, expected',
    [(12, 8, 20), (-16, 8, 24), (25, -10, 35), (-55, -15, 70)],
)
def test_fck(fck, delta, expected):
    """Test the fck function"""
    assert math.isclose(
        _section5_materials.fcm(fck, delta_f=delta), expected, rel_tol=0.005
    )


@pytest.mark.parametrize(
    'test_input1, expected',
    [
        (12, 1.572),
        (16, 1.905),
        (20, 2.21),
        (25, 2.565),
        (30, 2.896),
        (35, 3.21),
        (40, 3.509),
        (45, 3.795),
        (50, 4.072),
        (55, 4.183),
        (60, 4.306),
        (70, 4.533),
        (80, 4.74),
        (90, 4.93),
        (100, 5.106),
    ],
)
def test_fctm(test_input1, expected):
    """Test the fctm function."""
    assert math.isclose(
        _section5_materials.fctm(test_input1), expected, rel_tol=0.05
    )


@pytest.mark.parametrize(
    'fctm, expected',
    [
        (1.6, 1.12),
        (-2.2, 1.54),
        (-4.1, 2.87),
    ],
)
def test_fctk_5(fctm, expected):
    """Test the fctk_5 function"""
    assert math.isclose(
        _section5_materials.fctk_5(fctm), expected, rel_tol=0.005
    )


@pytest.mark.parametrize(
    'fctm, expected',
    [
        (1.6, 2.08),
        (-2.2, 2.86),
        (-4.1, 5.33),
    ],
)
def test_fctm_95(fctm, expected):
    """Test the fctk_95 function"""
    assert math.isclose(
        _section5_materials.fctk_95(fctm), expected, rel_tol=0.005
    )


@pytest.mark.parametrize(
    'fcm, kE, expected',
    [
        (20, 9500, 25786.9673576516),
        (24, 5000, 14422.4957030741),
        (43, 13000, 45544.1747850274),
        (53, 10000, 37562.8575422107),
        (88, 6000, 26687.7610868318),
    ],
)
def test_Ecm(fcm, kE, expected):
    """Test the Ecm function"""
    assert math.isclose(
        _section5_materials.Ecm(fcm, kE=kE), expected, rel_tol=10e-5
    )


@pytest.mark.parametrize(
    'fcm, kE',
    [
        (-20, 9500),
        (24, -5000),
        (43, 50000),
        (53, 3000),
    ],
)
def test_Ecm_raises_errors(fcm, kE):
    """Test the Ecm function raises errors"""
    with pytest.raises(ValueError):
        _section5_materials.Ecm(fcm, kE)


@pytest.mark.parametrize(
    'Ac, u, expected',
    [
        (400, 20, 40),
        (500, 30, 33.3333333333333),
        (200, 15, 26.6666666666667),
        (100, 40, 5),
    ],
)
def test_hn(Ac, u, expected):
    """Test the hn function"""
    assert math.isclose(_section5_materials.hn(Ac, u), expected, rel_tol=10e-5)


@pytest.mark.parametrize('Ac, u', [(-40, 20), (40, -20), (-40, -20)])
def test_hn_raises_errors(Ac, u):
    """Test the hn function raises errors"""
    with pytest.raises(ValueError):
        _section5_materials.hn(Ac, u)


@pytest.mark.parametrize(
    'hn, atm_conditions, expected',
    [
        (100, 'dry', 0.82),
        (150, 'dry', 0.805),
        (800, 'dry', 0.732),
        (200, 'humid', 0.68),
        (400, 'humid', 0.66666667),
        (800, 'humid', 0.648),
    ],
)
def test_A_phi_correction_exp(hn, atm_conditions, expected):
    """Test the A_phi_correction_exp function"""
    assert math.isclose(
        _section5_materials.A_phi_correction_exp(hn, atm_conditions),
        expected,
        rel_tol=10e-5,
    )


@pytest.mark.parametrize(
    'hn, atm_conditions', [(100, '1231'), (80, 'dry'), (1100, 'humid')]
)
def test_A_phi_correction_exp_raises_errors(hn, atm_conditions):
    """Test A_phi_correction_exp raises errors"""
    with pytest.raises(ValueError):
        _section5_materials.A_phi_correction_exp(hn, atm_conditions)


@pytest.mark.parametrize(
    't0, atm_conditions, _hn, concrete_class, expected',
    [
        (10, 'dry', 500, 'CS', 2.5),
        (28, 'humid', 200, 'CN', 1.6),
        (91, 'dry', 750, 'CR', 1.45),
        (60, 'humid', 600, 'CS', 1.41016),
    ],
)
def test_phi_50y_t0(t0, atm_conditions, _hn, concrete_class, expected):
    """Test phi_50y_t0 function"""
    assert math.isclose(
        _section5_materials.phi_50y_t0(
            t0, atm_conditions, _hn, concrete_class
        ),
        expected,
        rel_tol=10e-5,
    )


@pytest.mark.parametrize(
    't0, atm_conditions, _hn, concrete_class',
    [
        (-1, 'dry', 500, 'CS'),
        (50, 'ASDF', 500, 'CS'),
        (50, 'dry', 50, 'CS'),
        (50, 'dry', 1500, 'CS'),
        (50, 'dry', 500, 'ASD'),
        (1, 'dry', 100, 'CS'),
    ],
)
def test_phi_50y_t0_raises_errors(t0, atm_conditions, _hn, concrete_class):
    """Test phi_50y_t0 raises errors"""
    with pytest.raises(ValueError):
        _section5_materials.phi_50y_t0(t0, atm_conditions, _hn, concrete_class)


@pytest.mark.parametrize(
    'fck_28, atm_conditions, _hn, concrete_class, expected',
    [
        (35, 'dry', 500, 'CS', 0.45),
        (50, 'humid', 1000, 'CN', 0.23),
        (80, 'dry', 200, 'CR', 0.54),
        (40, 'humid', 300, 'CS', 0.29),
        (25, 'dry', 800, 'CN', 0.46333),
    ],
)
def test_eps_cs_50y(fck_28, atm_conditions, _hn, concrete_class, expected):
    """Test phi_50y_t0 function"""
    assert math.isclose(
        _section5_materials.eps_cs_50y(
            fck_28, atm_conditions, _hn, concrete_class
        ),
        expected,
        rel_tol=10e-5,
    )


@pytest.mark.parametrize(
    'fck_28, atm_conditions, _hn, concrete_class',
    [
        (15, 'dry', 500, 'CS'),
        (90, 'dry', 500, 'CS'),
        (50, 'ASDF', 500, 'CS'),
        (50, 'dry', 50, 'CS'),
        (50, 'dry', 1500, 'CS'),
        (50, 'dry', 500, 'ASD'),
        (1, 'dry', 100, 'CS'),
    ],
)
def test_eps_cs_50y_raises_errors(fck_28, atm_conditions, _hn, concrete_class):
    """Test eps_cs_50y raises errors"""
    with pytest.raises(ValueError):
        _section5_materials.eps_cs_50y(
            fck_28, atm_conditions, _hn, concrete_class
        )


@pytest.mark.parametrize(
    'fck, fck_ref, expected',
    [(60, 40, 0.873580464736299), (40, 45, 1), (60, 50, 0.941036028881029)],
)
def test_eta_cc(fck, fck_ref, expected):
    """Test eta_cc function"""
    assert math.isclose(
        _section5_materials.eta_cc(fck, fck_ref),
        expected,
        rel_tol=10e-5,
    )


@pytest.mark.parametrize('fck, fck_ref', [(-10, 40), (0, 20), (30, -10)])
def test_eta_cc_raises_errors(fck, fck_ref):
    """Test eta_cc raises errors"""
    with pytest.raises(ValueError):
        _section5_materials.eta_cc(fck, fck_ref)


def test_k_tc_raises_errors(t_ref, t0, concrete_class):
    """Test k_tc taises errors"""
    with pytest.raises(ValueError):
        _section5_materials.k_tc(t_ref, t0, concrete_class)
