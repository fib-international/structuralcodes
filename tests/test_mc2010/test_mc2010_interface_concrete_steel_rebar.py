"""Test for the functions of _interface_concrete_steel_rebar."""

import math
import warnings

import pytest

from structuralcodes.codes.mc2010 import _interface_concrete_steel_rebar as bss


@pytest.mark.parametrize(
    'bond, expected',
    [
        ('Good', 1.0),
        ('Other', 0.7),
    ],
)
def test_eta2(bond, expected):
    """Test eta_2 function."""
    assert bss.eta_2(bond) == expected


@pytest.mark.parametrize('invalid_bond', ['Bad', None, 1])
def test_eta_2_invalid_input(invalid_bond):
    """Test eta_2 warnings."""
    with pytest.raises(
        ValueError,
        match="Invalid bond condition. Must be 'good' or 'other'.",
    ):
        bss.eta_2(invalid_bond)


@pytest.mark.parametrize(
    'n_t, A_st, n_b, phi, s_t, expected',
    [
        (2, 314.1593, 1, 20, 200, 0.05),
        (4, 490.8739, 10, 25, 300, 0.02618),
    ],
)
def test_K_tr(n_t, A_st, n_b, phi, s_t, expected):
    """Test K_tr function."""
    assert math.isclose(
        bss.K_tr(n_t, A_st, n_b, phi, s_t), expected, rel_tol=1e-5
    )


@pytest.mark.parametrize(
    'eta_2, f_cm, phi, c_min, c_max, k_m, K_tr, expected',
    [
        (
            1.0,
            15,
            20,
            10,
            30,
            12,
            0.05,
            8.900446,
        ),
        (
            0.7,
            30,
            20,
            30,
            30,
            6,
            0.06,
            7.485077,
        ),
        (
            0.7,
            120,
            20,
            20,
            90,
            0,
            0.04,
            8.185118,
        ),
    ],
)
def test_tau_bu_split(eta_2, f_cm, phi, c_min, c_max, k_m, K_tr, expected):
    """Test tau_bu_split function."""
    assert math.isclose(
        bss.tau_bu_split(eta_2, f_cm, phi, c_min, c_max, k_m, K_tr),
        expected,
        rel_tol=1e-3,
    )


@pytest.mark.parametrize(
    'bond, f_cm, expected',
    [
        ('Good', 30, 13.69306394),
        ('Other', 21.3, 5.76899038),
        ('Good', 55, 18.54049622),
        ('Other', 60, 9.682458366),
    ],
)
def test_tau_bmax(bond, f_cm, expected):
    """Test tau_bmax function."""
    assert math.isclose(bss.tau_bmax(bond, f_cm), expected, rel_tol=1e-2)


@pytest.mark.parametrize('invalid_bond, f_cm', [('Bad', 20), ('poor', 20)])
def test_tau_bmax_invalid_input(invalid_bond, f_cm):
    """Test tau_bmax warnings."""
    with pytest.raises(
        ValueError,
        match="Invalid bond condition. Must be 'good' or 'other'.",
    ):
        bss.tau_bmax(invalid_bond, f_cm)


@pytest.mark.parametrize(
    'bond, expected',
    [
        ('Good', 1.0),
        ('Other', 1.8),
    ],
)
def test_s_1_valid(bond, expected):
    """Test s_1."""
    assert bss.s_1(bond) == expected


@pytest.mark.parametrize(
    'bond',
    [
        'Excellent',
    ],
)
def test_s_1_invalid(bond):
    """Test s_1 warnings."""
    with pytest.raises(
        ValueError, match="Invalid bond condition. Must be 'good' or 'other'."
    ):
        bss.s_1(bond)


@pytest.mark.parametrize(
    'bond, expected',
    [
        ('Good', 2.0),
        ('Other', 3.6),
    ],
)
def test_s_2_valid(bond, expected):
    """Test s_2."""
    assert bss.s_2(bond) == expected


@pytest.mark.parametrize(
    'bond',
    [
        'Excellent',
    ],
)
def test_s_2_invalid(bond):
    """Test s_2 warnings."""
    with pytest.raises(
        ValueError, match="Invalid bond condition. Must be 'good' or 'other'."
    ):
        bss.s_2(bond)


@pytest.mark.parametrize(
    'failmod, bond, confin, c_clear, s_1, expected',
    [
        ('PO', 'Good', 'Unconfined', 10.0, 1.0, 10.0),
        ('PO', 'Other', 'Stirrups', 15.0, 1.8, 15.0),
        ('SP', 'Good', 'Unconfined', 10.0, 1.0, 1.2),
        ('SP', 'Good', 'Stirrups', 10.0, 1.0, 5),
        ('SP', 'Other', 'Unconfined', 15.0, 1.8, 2.16),
        ('SP', 'Other', 'Stirrups', 15.0, 1.8, 7.5),
    ],
)
def test_s_3_valid(failmod, bond, confin, c_clear, s_1, expected):
    """Test s_3."""
    assert math.isclose(
        bss.s_3(failmod, bond, confin, c_clear, s_1), expected, rel_tol=1e-3
    )


@pytest.mark.parametrize(
    'failmod, bond, confin, c_clear, s_1, warning_message',
    [
        (
            'Splitting',
            'Good',
            'Unconfined',
            10,
            1.0,
            "Invalid failmod value. Must be 'PO' or 'SP'.",
        ),
        (
            'PO',
            'Bad',
            'Unconfined',
            10,
            1.0,
            "Invalid bond condition. Must be 'good' or 'other'.",
        ),
        (
            'SP',
            'Good',
            'confin',
            10,
            1.0,
            "Invalid confinement value. Must be 'unconfined' or 'stirrups'.",
        ),
    ],
)
def test_s_3_warnings(failmod, bond, confin, c_clear, s_1, warning_message):
    """Test s_3 warnings."""
    with pytest.raises(ValueError, match=warning_message):
        bss.s_3(failmod, bond, confin, c_clear, s_1)


@pytest.mark.parametrize(
    'tau_bmax, tau_bu_split, alpha, s_1, expected',
    [
        (13.7, 11.0, 0.4, 1.2, 0.6932),
        (20.1, 19, 0.4, 1.0, 0.868748),
    ],
)
def test_s_tau_bu_split(tau_bmax, tau_bu_split, alpha, s_1, expected):
    """Test s_tau_bu_split."""
    assert math.isclose(
        bss.s_tau_bu_split(tau_bmax, tau_bu_split, alpha, s_1),
        expected,
        rel_tol=1e-3,
    )


@pytest.mark.parametrize(
    'f_cm, phi, l_b, c_min, c_max, k_m, K_tr, expected',
    [
        (30, 16, 200, 20, 40, 12, 0.03, 370.1602),
        (23, 12, 400, 20, 40, 0, 0.03, 513.146),
        (40, 32, 500, 20, 40, 6, 0.03, 297.027),
    ],
)
def test_f_stm(f_cm, phi, l_b, c_min, c_max, k_m, K_tr, expected):
    """Test f_stm."""
    assert math.isclose(
        bss.f_stm(f_cm, phi, l_b, c_min, c_max, k_m, K_tr),
        expected,
        rel_tol=1e-3,
    )


@pytest.mark.parametrize(
    'f_cm, phi, l_b, c_min, c_max, k_m, K_tr',
    [
        (
            12,
            16,
            60,
            20,
            30,
            0,
            0.03,
        ),
        (
            30,
            20,
            200,
            5,
            40,
            0,
            0,
        ),
        (
            30,
            20,
            200,
            20,
            100,
            0,
            0,
        ),
        (
            30,
            20,
            200,
            20,
            40,
            1,
            0.1,
        ),
    ],
)
def test_f_stm_warnings(f_cm, phi, l_b, c_min, c_max, k_m, K_tr):
    """Test f_stm warnings."""
    raises = pytest.raises(UserWarning)
    with raises:
        warnings.filterwarnings(action='error', category=UserWarning)
        bss.f_stm(f_cm, phi, l_b, c_min, c_max, k_m, K_tr)


@pytest.mark.parametrize(
    'f_y, l_b, phi, expected',
    [
        (500, 160, 16, 12.5),
        (435, 150, 30, 21.75),
    ],
)
def test_tau_yield(f_y, l_b, phi, expected):
    """Test tau_yield."""
    assert math.isclose(
        bss.tau_yield(f_y, l_b, phi),
        expected,
        rel_tol=1e-3,
    )
