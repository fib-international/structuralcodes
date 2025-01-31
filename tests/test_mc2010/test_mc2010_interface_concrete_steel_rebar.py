import math
import warnings

import pytest
from structuralcodes.codes.mc2010 import _interface_concrete_steel_rebar as BSS


@pytest.mark.parametrize(
    "f_cm, bond, expected_exception",
    [
        (30, 'Good', None),
        (30, 'Other', None),
        (30.0, 'Bad', ValueError),  # Value error
        (30, 'good', ValueError),  # Value error
    ],
)
def test_bond_stress_slip(f_cm, bond, expected_exception):
    """Test BondStressSlip function"""
    if expected_exception:
        with pytest.raises(expected_exception):
            BSS.BondStressSlip(f_cm, bond)
    else:
        bond_stress_slip = BSS.BondStressSlip(f_cm, bond)
        assert bond_stress_slip.f_cm == f_cm
        assert bond_stress_slip.bond == bond


@pytest.mark.parametrize(
    "f_cm, bond, expected_tau_bmax",
    [
        (30, 'Good', 13.69306394),
        (21.3, 'Other', 5.76899038),
        (55, 'Good', 18.54049622),
        (60, 'Other', 9.682458366),
    ],
)
def test_tau_bmax(f_cm, bond, expected_tau_bmax):
    bond_stress_slip = BSS.BondStressSlip(f_cm=f_cm, bond=bond)
    assert math.isclose(
        bond_stress_slip.tau_bmax(), expected_tau_bmax, rel_tol=1e-2
    )


@pytest.mark.parametrize(
    "f_cm, bond, phi, c_min, c_max, k_m, K_tr, expected_tau_bu_split, expected_warnings",
    [
        (
            15,
            "Good",
            20,
            10,
            30,
            12,
            0.05,
            9.203289553,
            [
                "Warning: Eq.(6.1-19) is valid for 15 MPa < f_cm < 110 MPa",
                "Warning: Eq.(6.1-19) is valid for 0.5 < c_min / phi < 3.5",
            ],
        ),
        (
            30,
            "Other",
            20,
            30,
            30,
            6,
            0.06,
            7.303392633,
            [
                "Warning: Eq.(6.1-19) is valid for 1.0 < c_max / c_min < 5.0",
                "Warning: Eq.(6.1-19) is valid for K_tr <= 0.05",
            ],
        ),
        (
            120,
            "Other",
            20,
            20,
            90,
            0,
            0.04,
            8.185118287,
            ["Warning: Eq.(6.1-19) is valid for 15 MPa < f_cm < 110 MPa"],
        ),
        (30, 'Good', 16, 20, 40, 6, 0.02, 9.32212033, []),
        (40, 'Other', 16, 30, 40, 12, 0.02, 8.081104719, []),
        (90, 'Other', 32, 30, 35, 6, 0, 5.961183984, []),
    ],
)
def test_tau_bu_split(
    f_cm,
    bond,
    phi,
    c_min,
    c_max,
    k_m,
    K_tr,
    expected_tau_bu_split,
    expected_warnings,
):
    bss = BSS.BondStressSlip(f_cm, bond)
    with warnings.catch_warnings(record=True) as w:
        assert math.isclose(
            bss.tau_bu_split(phi, c_min, c_max, k_m, K_tr),
            expected_tau_bu_split,
            rel_tol=1e-3,
        )
        assert len(w) == len(expected_warnings)
        for i, warning in enumerate(w):
            assert issubclass(warning.category, UserWarning)
            assert str(warning.message) == expected_warnings[i]


@pytest.mark.parametrize(
    "f_cm, bond, tau_bu_split, expected_s_tau_bu_split",
    [
        (30, 'Good', 10.54002878, 0.51981),
        (85, 'Good', 11.93113642, 0.19279),
        (20, 'Other', 4.432754843, 1.007846),
        (50, 'Other', 5.948520129, 0.66817),
    ],
)
def test_s_tau_bu_split(f_cm, bond, tau_bu_split, expected_s_tau_bu_split):
    assert math.isclose(
        BSS.BondStressSlip(f_cm=f_cm, bond=bond).s_tau_bu_split(tau_bu_split),
        expected_s_tau_bu_split,
        rel_tol=1e-3,
    )


@pytest.mark.parametrize(
    "f_ym, l_b, phi, expected_tau_yield",
    [
        (500, 160, 16, 12.5),
        (435, 150, 30, 21.75),
    ],
)
def test_tau_yield(f_ym, l_b, phi, expected_tau_yield):
    assert math.isclose(
        BSS.BondStressSlip.tau_yield(f_ym, l_b, phi),
        expected_tau_yield,
        rel_tol=1e-3,
    )


def test_tau_yield_default():
    assert math.isclose(
        BSS.BondStressSlip.tau_yield(500), 500 / 20, rel_tol=1e-3
    )


def test_tau_yield_invalid():
    with pytest.raises(ValueError):
        BSS.BondStressSlip.tau_yield(500, l_b=5)
