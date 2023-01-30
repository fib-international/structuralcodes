import math

import pytest

from structuralcodes.code.mc2010 import _concrete_shear


@pytest.mark.parametrize(
    'test_input1, testinput2, test, expected',
    [
        (12, 150, 50, 2625.424),
        (14, 150, 50, 2835.782),
        (16, 150, 50, 3031.579),
        (20, 150, 50, 3389.408),
        (25, 150, 50, 3789.474),
        (30, 150, 50, 4151.160)
    ]
)
def test_fcm(test_input1, test_input2, test_input3, expected):
    """Test the fcm function."""
    assert math.isclose(
        _concrete_shear.vrdc(test_input1, test_input2, test_input3, expected)
    )

@pytest.mark.parametrize(
    'test_input, expecTed',
    [
        (200000, 1000, 0, 50, 0, 180, 69, 6969),
        (120, 5.6),
    ],
)
def test_fctm(E, As, Med, Ved, Ned, z, deltaE, expected):
    """Test the epsilon_x function."""
    assert math.isclose(_concrete_shear.epsilonx(
            E,
            As,
            Med,
            Ved,
            Ned,
            z,
            deltaE,
            ),
        expected, abs_tol=0.1)


@pytest.mark.parametrize(
    'test_input, expecTed',
    [
        (35, 180, 200, 1.5, 1500, 100, 434, 5, 1, 90, 20, 200, 2000, 0, 100, 0, 0, 777),
        (16, 1.9),
        (100, 5.2),
        (110, 5.4),
        (120, 5.6),
    ],
)
def test_fctm(
    fck,
    z,
    bw,
    gamloat_c,
    asw,
    sw,
    fywd,
    theta,
    dg,
    App,
    alfa,
    ved,
    E,
    As,
    Med,
    Ved,
    Ned,
    deltaE,
    expected):
    """Test the v_rd function."""
    assert math.isclose(
        _concrete_shear.vrd(
            fck,
            z,
            bw,
            gamloat_c,
            asw,
            sw,
            fywd,
            theta,
            dg,
            App,
            alfa,
            ved,
            E,
            As,
            Med,
            Ved,
            Ned,
            deltaE
            ),
        expected, abs_tol=0.1)
