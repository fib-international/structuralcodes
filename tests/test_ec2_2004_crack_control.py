"""Tests for EUROCODE 2-1-1:2004 Chapter 7.3 Crack Control"""
import math

import pytest
from structuralcodes.codes.ec2_2004 import _crack_control


@pytest.mark.parametrize(
    'test_exposure_class, test_load_combination, expected',
    [
        ('X0', 'f', 0.2),
        ('x0', 'F', 0.2),
        ('X0', 'qp', 0.4),
        ('x0', 'QP', 0.4),
        ('XC2', 'f', 0.2),
        ('xc2', 'F', 0.2),
        ('XC3', 'f', 0.2),
        ('xc3', 'F', 0.2),
        ('XC4', 'f', 0.2),
        ('xc4', 'F', 0.2),
        ('XC2', 'qp', 0.3),
        ('xc2', 'QP', 0.3),
        ('XC3', 'qp', 0.3),
        ('xc3', 'QP', 0.3),
        ('XC4', 'qp', 0.3),
        ('xc4', 'QP', 0.3),
        ('XD1', 'qp', 0.3),
        ('xd1', 'QP', 0.3),
        ('XD2', 'qp', 0.3),
        ('xd2', 'QP', 0.3),
        ('XS1', 'qp', 0.3),
        ('xs1', 'QP', 0.3),
        ('XS2', 'qp', 0.3),
        ('xs2', 'QP', 0.3),
        ('XS3', 'qp', 0.3),
        ('xs3', 'QP', 0.3),
    ],
)
def test_w_max_returns_expected_values(
    test_exposure_class, test_load_combination, expected
):
    """Test that the w_max function returns expected values"""
    w_max = _crack_control.w_max(test_exposure_class, test_load_combination)
    assert w_max == expected


@pytest.mark.parametrize(
    'test_exposure_class, test_load_combination',
    [('dummy1', 'f'), ('dummy2', 'qp'), ('XD1', 'dummy3'), ('XS1', 'dummy4')],
)
def test_w_max_not_valid_input_raises_valueerror(
    test_exposure_class, test_load_combination
):
    """Test that not valid input returns ValueError"""
    with pytest.raises(ValueError):
        _crack_control.w_max(test_exposure_class, test_load_combination)


@pytest.mark.parametrize(
    'h, expected',
    [
        (200, 1),
        (300, 1),
        (800, 0.65),
        (1000, 0.65),
        (400, 0.93),
        (500, 0.86),
        (600, 0.79),
        (700, 0.72),
    ],
)
def test_k_crack_min_steel_area_returns_expected_values(h, expected):
    """Test the k_crack_min_steel_area function"""
    k = _crack_control.k_crack_min_steel_area(h)
    assert math.isclose(k, expected)


def test_k_crack_min_steel_area_raises_valueerror():
    """Test that not valid input returns ValueError exeption"""
    with pytest.raises(ValueError):
        h = -100
        _crack_control.k_crack_min_steel_area(h)


def test_kc_crack_min_steel_area_pure_tension_returns_expected_values():
    """Test the kc_crack_min_steel_area_pure_tension function"""
    assert 1 == _crack_control.kc_crack_min_steel_area_pure_tension()


@pytest.mark.parametrize(
    'h, b, fct_eff, n_ed, expected',
    [
        (600, 400, 3, 20, 0.3925926),
        (600, 400, 3, -20, 0.4166667),
        (400, 200, 4, 3, 0.397500),
        (200, 50, 5, -80, 1),
        (200, 50, 5, 80, 0),
    ],
)
def test_kc_crack_min_steel_area_rectangular_returns_expected_values(
    h, b, fct_eff, n_ed, expected
):
    """Test the kc_crack_min_steel_area_rectangular"""
    kc = _crack_control.kc_crack_min_steel_area_rectangular(
        h, b, fct_eff, n_ed
    )
    assert math.isclose(kc, expected, rel_tol=0.000001)


def test_kc_crack_min_steel_area_rectangular_rasies_valueerror():
    """Test the kc_crack_min_steel_area_rectangular raises Value
    Error for not correct input values for b and h"""
    with pytest.raises(ValueError):
        _crack_control.kc_crack_min_steel_area_rectangular(
            h=-100, b=100, fct_eff=100, n_ed=10
        )
        _crack_control.kc_crack_min_steel_area_rectangular(
            h=100, b=-100, fct_eff=100, n_ed=10
        )


@pytest.mark.parametrize(
    'f_cr, a_ct, fct_eff, expected',
    [
        (30, 10000, 5, 0.54),
        (20, 5000, 3, 1.2),
        (55, 7500, 4, 1.65),
        (55, 50000, 4, 0.5),
    ],
)
def test_kc_crack_min_steel_area_flanges(f_cr, a_ct, fct_eff, expected):
    """Test the kc_crack_min_steel_area_flanges function"""
    kc = _crack_control.kc_crack_min_steel_area_flanges(f_cr, a_ct, fct_eff)
    assert math.isclose(kc, expected, rel_tol=0.000001)


@pytest.mark.parametrize(
    'a_ct, s_steel, fct_eff, k, kc, expected',
    [
        (10000, 500, 3, 1, 1, 60),
        (80000, 500, 5, 0.65, 0.5, 260),
        (80000, 400, 4, 0.9, 0.75, 540),
    ],
)
def test_crack_min_steel_area_returns_expected_values(
    a_ct, s_steel, fct_eff, k, kc, expected
):
    """Test the crack_min_steel_area returns expected values"""
    as_min = _crack_control.crack_min_steel_area(a_ct, s_steel, fct_eff, k, kc)
    assert math.isclose(as_min, expected, rel_tol=0.000001)
