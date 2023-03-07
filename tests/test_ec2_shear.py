"""Test for the ec2.shear module"""
import math
import pytest

from structuralcodes.codes.ec2 import shear

@pytest.mark.parametrize(
    'd, expected',
    [
        (10, 2.),
        (20, 2.),
        (50, 2.),
        (100, 2.),
        (150, 2.),
        (200, 2.),
        (250, 1.89),
        (500, 1.63),
        (1000, 1.45),
    ],
)

def test_k(d, expected):
    """Test the _k function."""
    assert math.isclose(
        shear._k(d), expected, rel_tol=0.01
    )

@pytest.mark.parametrize(
    'Asl, bw, d, expected',
    [
        (100, 100, 300, 0.00333),
        (150, 100, 300, 0.005),
        (200, 100, 300, 0.00667),
        (250, 100, 300, 0.00833),
        (500, 100, 300, 0.01667),
        (1000, 100, 300, 0.02),
        (1500, 100, 300, 0.02),
    ],
)

def test_rho_L(Asl, bw, d, expected):
    """Test the _rho_L function"""
    assert math.isclose(
        shear._rho_L(Asl, bw, d), expected, rel_tol=0.01
    )

@pytest.mark.parametrize(
    'Ned, Ac, fcd, expected',
    [
        (0, 100*300, 13.3, 0.),
        (1e3, 100*300, 13.3, 0.0333),
        (5e3, 100*300, 13.3, 0.1667),
        (1e4, 100*300, 13.3, 0.3333),
        (5e4, 100*300, 13.3, 1.6667),
        (1e5, 100*300, 13.3, 2.66),
        (1e5, 500*500, 13.3, 0.4),
        (1e5, 500*500, 25., 0.4),
        (1e8, 500*500, 25., 5.0),
    ],
)

def test_sigma_cp(Ned, Ac, fcd, expected):
    """Test the sigma_cp function."""
    assert math.isclose(
        shear._sigma_cp(Ned, Ac, fcd), expected, rel_tol=0.01
    )

def test_alpha_l():
    """To be included later."""
    pass

# Expand test below.
@pytest.mark.parametrize(
    'fck, d, Asl, bw, Ned, Ac, k1, gamma_c, expected',
    [
        (40, 5422, 0.02*563*5422, 563, 7.35*16.292e6, 16.292e6,
            0.15, 1.5, 4323e3),
    ],
)

def test_VRdc(fck, d, Asl, bw, Ned, Ac, k1, gamma_c, expected):
    """Test the VRdc function."""
    assert math.isclose(
        shear.VRdc(fck, d, Asl, bw, Ned, Ac, k1, gamma_c),
        expected, rel_tol=0.01
    )

@pytest.mark.parametrize(
    'fck, d, expected',
    [
        (20., 100, 0.443),
        (20., 300, 0.383),
        (20., 500, 0.326),
        (20., 1000, 0.273),
        (40., 100, 0.626),
        (40., 300, 0.542),
        (40., 500, 0.462),
        (40., 1000, 0.385),
    ],
)

def test_vmin(fck, d, expected):
    """Test the vmin function."""
    assert math.isclose(
        shear.vmin(fck, d), expected, rel_tol=0.01
    )

def test_Vrdc_prin_stress():
    """To be included later."""
    pass

@pytest.mark.parametrize(
    'bw, d, fck, fcd, expected',
    [
        (100, 200, 20, 20/1.5, 73600.),
        (100, 500, 20, 20/1.5, 184000.),
        (300, 700, 20, 20/1.5, 772800.)
    ]
)

def test_VEdmax_unreinf(bw, d, fck, fcd, expected):
    """test VEedmax_unreinf function."""
    assert math.isclose(
        shear.VEdmax_unreinf(bw, d, fck, fcd), expected, rel_tol=.01
    )

@pytest.mark.parametrize(
    'fck, expected',
    [
        (12, .5712),
        (16, .5616),
        (20, .552),
        (25, .54),
        (30, .528),
        (35, .516),
        (40, .504),
        (45, .492),
        (50, .48),
        (55, .468),
        (60, .456),
        (70, .432),
        (80, .408),
        (90, .384),
        (100, .36),
        (110, .336),
        (120, .312),
    ],
)

def test_v(fck, expected):
    """Test the v function."""
    assert math.isclose(
        shear.v(fck), expected, rel_tol=0.01
    )

@pytest.mark.parametrize(
    '_theta, cot_min, cot_max, expected',
    [
        (20, 1, 2.5, None),
        (100, 1, 2.5, None)
    ]
)

def test_theta(_theta, cot_min, cot_max, expected):
    """Test the theta function."""
    with pytest.raises(ValueError) as exc_info:
        math.isclose(
            shear.theta(_theta, cot_min, cot_max), expected,
            rel_tol=0.01)
        assert str(exc_info.value).startswith("Wrong value for theta is chosen.")

@pytest.mark.parametrize(
    'Asw, s, z, _theta, fyk, gamma_s, alpha, expected',
    [
        (2*1/4*math.pi*16**2, 200, 4880, 45, 400, 1.15, 90, 3413e3),
        (2*1/4*math.pi*20**2, 200, 4880, 45, 400, 1.15, 90, 5332e3),
        (2*1/4*math.pi*16**2, 200, 4880, 30, 400, 1.15, 90, 5911e3),
        (2*1/4*math.pi*20**2, 200, 4880, 30, 400, 1.15, 90, 9236e3),
        (2*1/4*math.pi*16**2, 200, 4400, 45, 400, 1.15, 90, 3077e3),
        (2*1/4*math.pi*20**2, 200, 4400, 45, 400, 1.15, 90, 4808e3),
        (2*1/4*math.pi*16**2, 200, 4400, 30, 400, 1.15, 90, 5330e3),
        (2*1/4*math.pi*20**2, 200, 4400, 30, 400, 1.15, 90, 8327e3),
        (2*1/4*math.pi*16**2, 200, 1757, 45, 400, 1.15, 90, 1229e3),
        (2*1/4*math.pi*20**2, 200, 1757, 45, 400, 1.15, 90, 1920e3),
        (2*1/4*math.pi*16**2, 200, 1757, 30, 400, 1.15, 90, 2128e3),
        (2*1/4*math.pi*20**2, 200, 1757, 30, 400, 1.15, 90, 3326e3),
        (2*1/4*math.pi*20**2, 200, 1757, 100, 40, 1.15, 90, 1), # Check if ValueError is raised.
        (2*1/4*math.pi*16**2, 200, 4880, 45, 400, 1.15, 45, 4825e3),
        (2*1/4*math.pi*20**2, 200, 4880, 45, 400, 1.15, 45, 7537e3),
        (2*1/4*math.pi*16**2, 200, 4880, 30, 400, 1.15, 45, 6591e3),
        (2*1/4*math.pi*20**2, 200, 4880, 30, 400, 1.15, 45, 10296e3),
        (2*1/4*math.pi*16**2, 200, 4400, 45, 400, 1.15, 45, 4350e3),
        (2*1/4*math.pi*20**2, 200, 4400, 45, 400, 1.15, 45, 6796e3),
        (2*1/4*math.pi*16**2, 200, 4400, 30, 400, 1.15, 45, 5943e3),
        (2*1/4*math.pi*20**2, 200, 4400, 30, 400, 1.15, 45, 9284e3),
        (2*1/4*math.pi*16**2, 200, 1757, 45, 400, 1.15, 45, 1737e3),
        (2*1/4*math.pi*20**2, 200, 1757, 45, 400, 1.15, 45, 2714e3),
        (2*1/4*math.pi*16**2, 200, 1757, 30, 400, 1.15, 45, 2373e3),
        (2*1/4*math.pi*20**2, 200, 1757, 30, 400, 1.15, 45, 3707e3),
        (2*1/4*math.pi*20**2, 200, 1757, 100, 400, 1.15, 45, 1), # Check if ValueError is raised.
    ],
)

def test_VRds(Asw, s, z, _theta, fyk, gamma_s, alpha, expected):
    """Test the VRds function."""
    try:
        assert math.isclose(
            shear.VRds(Asw, s, z, _theta, fyk, gamma_s, alpha),
            expected, rel_tol=0.01)
    except ValueError:
        with pytest.raises(ValueError) as exc_info:
            math.isclose(
                shear.VRds(Asw, s, z, _theta, fyk, gamma_s, alpha),
                expected, rel_tol=0.01)
            assert str(exc_info.value).startswith("Wrong value for theta is chosen.")

@pytest.mark.parametrize(
    'bw, z, fck, _theta, Ned, Ac, gamma_c, alpha, expected',
    [
        (100, 300, 20, 45, 100e3, 100*400, 1.5, 90, 131152),
        (100, 300, 20, 21.8, 100e3, 100*400, 1.5, 90, 90445),
        (100, 300, 20, 45, 100e3, 100*400, 1.5, 45, 262304),
        (100, 300, 20, 21.8, 100e3, 100*400, 1.5, 45, 126620)
    ]
)

def test_VRdmax(bw, z, fck, _theta, Ned, Ac, gamma_c, alpha, expected):
    """Test the VRdmax function."""
    assert math.isclose(
        shear.VRdmax(bw, z, fck, _theta, Ned, Ac, gamma_c, alpha),
        expected, rel_tol=.01)

@pytest.mark.parametrize(
    'Ned, Ac, fcd, expected',
    [
        (0, 100*400, 20/1.5, 1),
        (100e3, 100*400, 20/1.5, 1.1875),
        (250e3, 100*400, 20/1.5, 1.25),
        (500e3, 100*400, 20/1.5, 0.15625)
    ]
)

def test_alpha_cw(Ned, Ac, fcd, expected):
    """Test the alpha_cw function."""
    assert math.isclose(
        shear.alpha_cw(Ned, Ac, fcd), expected, rel_tol=.01
    )

@pytest.mark.parametrize(
    'fcd, fck, bw, s, fywd, Ned, Ac, alpha, expected',
    [
        (20/1.5, 20, 100, 200, 435, 100e3, 100*400, 90, 201),
        (20/1.5, 20, 100, 200, 435, 100e3, 100*400, 45, 284)
    ]
)

def test_Asw_max(fcd, fck, bw, s, fywd, Ned, Ac, alpha, expected):
    """Test the Asw_max function."""
    assert math.isclose(
        shear.Asw_max(fcd, fck, bw, s, fywd, Ned, Ac, alpha), expected,
        rel_tol=.01
    )
