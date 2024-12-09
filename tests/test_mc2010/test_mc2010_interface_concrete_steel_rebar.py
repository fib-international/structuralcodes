import math
import pytest
from structuralcodes.codes.mc2010 import _interface_concrete_steel_rebar as BSS


def test_BondStressSlip_invalid():
    with pytest.raises(ValueError):
        BSS.BondStressSlip(f_cm=30, bond='poor')
    with pytest.raises(ValueError):
        BSS.BondStressSlip(f_cm=30, bond='moderate')


@pytest.mark.parametrize("f_cm, bond, expected_tau_bmax", [
    (30, 'Good', 13.69306),
    (21.3, 'Other', 5.76899),
    (55, 'Good', 18.5405),
    (60, 'Other', 9.682458),
])
def test_tau_bmax(f_cm, bond, expected_tau_bmax):
    bond_stress_slip = BSS.BondStressSlip(f_cm=f_cm, bond=bond)
    assert math.isclose(bond_stress_slip.tau_bmax(), expected_tau_bmax, rel_tol=1e-3)


@pytest.mark.parametrize("f_cm, bond, phi, c_min, c_max, k_m, K_tr, expected_tau_bu_split", [
    (30, 'Good', 16, 20, 40, 6, 0.02, 9.32212),
    (40, 'Other', 16, 30, 40, 12, 0.02, 8.081104719),
    (23.5, 'Good', 32, 30, 30, 12, 0.05, 9.649550592),
    (90, 'Other', 32, 30, 35, 6, 0, 5.961183984),
])
def test_tau_bu_split(f_cm, bond, phi, c_min, c_max, k_m, K_tr, expected_tau_bu_split):
    BSS1 = BSS.BondStressSlip(f_cm=f_cm, bond=bond)
    assert math.isclose(BSS1.tau_bu_split(phi, c_min, c_max, k_m, K_tr), expected_tau_bu_split, rel_tol=1e-3)


@pytest.mark.parametrize("f_cm, bond, tau_bu_split, expected_s_tau_bu_split", [
    (30, 'Good', 10.54002878, 0.51981),
    (85, 'Good', 11.93113642, 0.19279),
    (20, 'Other', 4.432754843, 1.007846),
    (50, 'Other', 5.948520129, 0.66817),
])
def test_s_tau_bu_split(f_cm, bond, tau_bu_split, expected_s_tau_bu_split):
    assert math.isclose(BSS.BondStressSlip(f_cm=f_cm, bond=bond).s_tau_bu_split(tau_bu_split), expected_s_tau_bu_split,
                        rel_tol=1e-3)


@pytest.mark.parametrize("f_ym, l_b, phi, expected_tau_yield", [
    (500, 160, 16, 12.5),
    (435, 150, 30, 21.75),
])
def test_tau_yield(f_ym, l_b, phi, expected_tau_yield):
    assert math.isclose(BSS.BondStressSlip.tau_yield(f_ym, l_b, phi), expected_tau_yield, rel_tol=1e-3)


def test_tau_yield_default():
    assert math.isclose(BSS.BondStressSlip.tau_yield(500), 500 / 20, rel_tol=1e-3)


def test_tau_yield_invalid():
    with pytest.raises(ValueError):
        BSS.BondStressSlip.tau_yield(500, l_b=5)
