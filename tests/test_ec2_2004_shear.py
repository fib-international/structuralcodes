"""Test for the ec2.shear module"""
import math
import pytest

from structuralcodes.codes.ec2_2004 import shear
from structuralcodes.codes.ec2_2004.shear import _k, _rho_L, _sigma_cp, _theta


@pytest.mark.parametrize(
    'd, expected',
    [
        (100, 2.0),
        (150, 2.0),
        (200, 2.0),
        (250, 1.89),
        (450, 1.667),
        (1000, 1.45),
    ],
)
def test_k(d, expected):
    """Test the _k function."""
    assert math.isclose(_k(d), expected, rel_tol=0.01)


@pytest.mark.parametrize(
    'Asl, bw, d, expected',
    [
        (100, 100, 250, 0.004),
        (200, 100, 250, 0.008),
        (400, 100, 250, 0.016),
        (1000, 100, 250, 0.02),
        (1500, 100, 250, 0.02),
    ],
)
def test_rho_L(Asl, bw, d, expected):
    """Test the _rho_L function"""
    assert math.isclose(_rho_L(Asl, bw, d), expected, rel_tol=0.01)


@pytest.mark.parametrize(
    'Ned, Ac, fcd, expected',
    [
        (0, 100 * 300, 13.3, 0.0),
        (1e4, 100 * 300, 13.3, 0.3333),
        (5e4, 100 * 300, 13.3, 1.6667),
        (1e5, 100 * 300, 13.3, 2.66),
        (1e5, 500 * 500, 13.3, 0.4),
        (1e8, 500 * 500, 25.0, 5.0),
    ],
)
def test_sigma_cp(Ned, Ac, fcd, expected):
    """Test the sigma_cp function."""
    assert math.isclose(_sigma_cp(Ned, Ac, fcd), expected, rel_tol=0.01)


def test_alpha_l():
    """To be included later."""
    # pass


@pytest.mark.parametrize(
    'fck, d, expected',
    [
        (20.0, 250, 0.407),
        (20.0, 450, 0.337),
        (37.5, 250, 0.557),
        (37.5, 450, 0.461),
    ],
)
def test_vmin(fck, d, expected):
    """Test the vmin function."""
    assert math.isclose(shear.vmin(fck, d), expected, rel_tol=0.01)


# Expand test below.
@pytest.mark.parametrize(
    'fck, d, Asl, bw, Ned, Ac, k1, gamma_c, expected',
    [
        (20, 250, 200, 100, 5e4, 30000, 0.15, 1.5, 20538),
        (37.5, 450, 1000, 500, 1e8, 250000, 0.15, 1.5, 283334),
        (20, 250, 0, 100, 5e4, 30000, 0.15, 1.5, 16425),
        (37.5, 450, 0, 500, 1e8, 250000, 0.15, 1.5, 272475),
    ],
)
def test_VRdc(fck, d, Asl, bw, Ned, Ac, k1, gamma_c, expected):
    """Test the VRdc function."""
    assert math.isclose(
        shear.VRdc(fck, d, Asl, bw, Ned, Ac, k1, gamma_c),
        expected,
        rel_tol=0.01,
    )


def test_Vrdc_prin_stress():
    """To be included later."""
    # pass


@pytest.mark.parametrize(
    'bw, d, fck, expected',
    [
        (100, 250, 20, 91770.0),
        (100, 250, 37.5, 159375.0),
        (500, 450, 37.5, 1434375.0),
    ],
)
def test_VEdmax_unreinf(bw, d, fck, expected):
    """test VEedmax_unreinf function."""
    assert math.isclose(
        shear.VEdmax_unreinf(bw, d, fck), expected, rel_tol=0.01
    )


@pytest.mark.parametrize(
    'fck, expected',
    [
        (20, 0.552),
        (30, 0.528),
        (37.5, 0.51),
        (40, 0.504),
        (50, 0.48),
        (60, 0.456),
        (70, 0.432),
        (90, 0.384),
        (100, 0.36),
        (120, 0.312),
    ],
)
def test_v(fck, expected):
    """Test the v function."""
    assert math.isclose(shear.v(fck), expected, rel_tol=0.01)


@pytest.mark.parametrize(
    'theta, cot_min, cot_max, expected',
    [(20, 1, 2.5, None), (100, 1, 2.5, None)],
)
def test_theta(theta, cot_min, cot_max, expected):
    """Test the theta function."""
    with pytest.raises(ValueError) as exc_info:
        assert _theta(theta, cot_min, cot_max) == expected
    assert str(exc_info.value).startswith("Wrong value for theta is chosen.")


@pytest.mark.parametrize(
    'Asw, s, z, theta, fyk, alpha, gamma_s, expected',
    [
        (2 * 1 / 4 * math.pi * 16**2, 200, 4880, 45, 400, 90, 1.15, 3413e3),
        (2 * 1 / 4 * math.pi * 20**2, 200, 4880, 45, 400, 90, 1.15, 5332e3),
        # Check if ValueError is raised.
        (2 * 1 / 4 * math.pi * 20**2, 200, 1757, 100, 400, 90, 1.15, 1),
        (2 * 1 / 4 * math.pi * 16**2, 200, 4880, 45, 400, 45, 1.15, 4825e3),
        (2 * 1 / 4 * math.pi * 20**2, 200, 4880, 45, 400, 45, 1.15, 7537e3),
        # Check if ValueError is raised.
        (2 * 1 / 4 * math.pi * 20**2, 200, 1757, 100, 400, 45, 1.15, 1),
    ],
)
def test_VRds(Asw, s, z, theta, fyk, alpha, gamma_s, expected):
    """Test the VRds function."""
    try:
        assert math.isclose(
            shear.VRds(Asw, s, z, theta, fyk, alpha, gamma_s),
            expected,
            rel_tol=0.01,
        )
    except ValueError:
        with pytest.raises(ValueError) as exc_info:
            assert math.isclose(
                shear.VRds(Asw, s, z, theta, fyk, alpha, gamma_s),
                expected,
                rel_tol=0.01,
            )
        assert str(exc_info.value).startswith(
            "Wrong value for theta is chosen."
        )


@pytest.mark.parametrize(
    'bw, z, fck, theta, Ned, Ac, gamma_c, alpha, expected',
    [
        (100, 300, 20, 45, 100e3, 100 * 400, 1.5, 90, 131152),
        (100, 300, 20, 21.8, 100e3, 100 * 400, 1.5, 90, 90445),
        (100, 300, 20, 45, 100e3, 100 * 400, 1.5, 45, 262304),
        (100, 300, 20, 21.8, 100e3, 100 * 400, 1.5, 45, 126620),
    ],
)
def test_VRdmax(bw, z, fck, theta, Ned, Ac, gamma_c, alpha, expected):
    """Test the VRdmax function."""
    assert math.isclose(
        shear.VRdmax(bw, z, fck, theta, Ned, Ac, gamma_c, alpha),
        expected,
        rel_tol=0.01,
    )


@pytest.mark.parametrize(
    'Ned, Ac, fcd, expected',
    [
        (-1, 100 * 400, 20 / 1.5, 1),
        (0, 100 * 400, 20 / 1.5, 1),
        (100e3, 100 * 400, 20 / 1.5, 1.1875),
        (250e3, 100 * 400, 20 / 1.5, 1.25),
        (500e3, 100 * 400, 20 / 1.5, 0.15625),
        # Check if ValueError is raised
        (10000e3, 100 * 400, 20 / 1.5, 0.15625),
    ],
)
def test_alpha_cw(Ned, Ac, fcd, expected):
    """Test the alpha_cw function."""
    try:
        assert math.isclose(
            shear.alpha_cw(Ned, Ac, fcd), expected, rel_tol=0.01
        )
    except ValueError:
        with pytest.raises(ValueError) as exc_info:
            assert math.isclose(
                shear.alpha_cw(Ned, Ac, fcd), expected, rel_tol=0.01
            )
        assert str(exc_info.value).startswith("sigma_cp/fcd=")


@pytest.mark.parametrize(
    'fcd, fck, bw, s, fywd, Ned, Ac, alpha, expected',
    [
        (20 / 1.5, 20, 100, 200, 435, 100e3, 100 * 400, 90, 201),
        (20 / 1.5, 20, 100, 200, 435, 100e3, 100 * 400, 45, 284),
    ],
)
def test_Asw_max(fcd, fck, bw, s, fywd, Ned, Ac, alpha, expected):
    """Test the Asw_max function."""
    assert math.isclose(
        shear.Asw_max(fcd, fck, bw, s, fywd, Ned, Ac, alpha),
        expected,
        rel_tol=0.01,
    )
